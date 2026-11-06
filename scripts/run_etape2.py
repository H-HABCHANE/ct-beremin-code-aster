#!/usr/bin/env python3
"""Etape 2 : J et K_J elastoplastiques sur l'eprouvette CT, avec controle independant par l'aire P-delta."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analytic import k_e399, k_from_g_plane_strain, mpa_sqrt_m  # noqa: E402
from aster_run import copy_comm, run, write_export  # noqa: E402
from materiau import LOIS, Loi  # noqa: E402
from msh_to_aster import read_msh22  # noqa: E402
from run_etape1 import W, mesh, parse_table, write_csv  # noqa: E402

A_W = 0.5
A0, B0 = A_W * W, (1.0 - A_W) * W
ARC = math.pi * 0.125 * W
# Couronnes d'integration de J (mm). La derniere englobe la zone plastique de pointe et reste a distance
# du trou (a environ 22 mm de la pointe) et de la face arriere (a 25 mm).
CONTOURS_J = [[0.5, 1.0], [1.0, 2.0], [2.0, 4.0], [4.0, 8.0], [8.0, 16.0]]
INDICE_CONTOUR_LOIN = len(CONTOURS_J) - 1
ETA_E1820 = 2.0 + 0.522 * B0 / W  # facteur eta_pl de ASTM E1820 pour une eprouvette CT
ALPHA = math.sqrt((2 * A0 / B0) ** 2 + 4 * A0 / B0 + 2) - (2 * A0 / B0 + 1)  # Kumar, German, Shih (1981)


def charge_limite(sigma_0: float) -> float:
    """Charge limite P/B en deformations planes : P0 = 1.455 alpha b sigma_0."""
    return 1.455 * ALPHA * B0 * sigma_0


def arc_elements(msh_text: str):
    nodes, elements, groups = read_msh22(msh_text)
    arc_tags = set(groups["TROU_HAUT"])
    return nodes, [conn for tag, _kind, conn in elements if tag in arc_tags]


def deplacement_conjugue(table_rows, nodes, segments):
    """Moyenne de DY le long de l'arc charge (Simpson sur SEG3) : travail = P/B x delta pour un effort uniforme."""
    by_coord = {(round(x, 6), round(y, 6)): tag for tag, (x, y) in nodes.items()}
    par_ordre = defaultdict(dict)
    for row in table_rows:
        tag = by_coord[(round(float(row["COOR_X"]), 6), round(float(row["COOR_Y"]), 6))]
        par_ordre[int(row["NUME_ORDRE"])][tag] = float(row["DY"])
    out = {}
    for ordre, dy in par_ordre.items():
        total, longueur = 0.0, 0.0
        for n1, n2, n3 in segments:
            segment = math.dist(nodes[n1], nodes[n3]) + math.dist(nodes[n3], nodes[n2])
            total += segment * (dy[n1] + 4.0 * dy[n3] + dy[n2]) / 6.0
            longueur += segment
        out[ordre] = total / longueur
    return out


def rayon_zone_plastique(points, pas: float = 0.25) -> float:
    """Rayon de la zone plastique connectee a la pointe.

    `points` contient (distance a la pointe, plastifie ?). On avance par couronnes de largeur `pas` et on
    s'arrete a la premiere couronne maillee sans point plastifie : la plasticite de la face arriere ou du
    trou, separee de la pointe, n'est donc pas comptee.
    """
    couronnes: dict[int, bool] = {}
    for distance, plastifie in points:
        indice = int(distance // pas)
        couronnes[indice] = couronnes.get(indice, False) or plastifie
    rayon = 0.0
    for indice in sorted(couronnes):
        if not couronnes[indice]:
            break
        rayon = (indice + 1) * pas
    return rayon


def run_case(loi: Loi, h_tip: float, instants: list[float], quad: bool = False) -> list[dict]:
    label = f"{loi.nom}_h{h_tip:.3f}_{'q8si' if quad else 't6'}_n{len(instants)}"
    case_dir = ROOT / "build" / "etape2" / label
    case_dir.mkdir(parents=True, exist_ok=True)
    mesh(case_dir, A_W, h_tip, quad=quad)
    ordres_zp = [ordre for ordre in (5, 10, 15, 20, 25, 30) if ordre <= len(instants)] or [len(instants)]
    nodes, segments = arc_elements((case_dir / "ct.msh").read_text())
    params = {
        "E": loi.young, "nu": loi.nu, "courbe": loi.courbe_traction(), "contours": CONTOURS_J,
        "instants": [0.0] + instants, "traction_unitaire": 1.0 / ARC,
        "modelisation": "D_PLAN_SI" if quad else "D_PLAN", "ordres_zp": ordres_zp,
    }
    (case_dir / "params.json").write_text(json.dumps(params))
    comm = copy_comm(ROOT / "aster/ct_elastoplastique.comm", case_dir)
    write_export(case_dir, comm, {20: "ct.mail", 30: "params.json"},
                 {38: "reaction.txt", 39: "g.txt", 41: "arc.txt", 42: "ligament.txt", 43: "zp.txt"}, memory_mb=4000)
    run(case_dir)

    g_par_ordre = defaultdict(list)
    for row in parse_table((case_dir / "g.txt").read_text()):
        g_par_ordre[int(row["NUME_ORDRE"])].append(float(row["G"]))
    delta = deplacement_conjugue(parse_table((case_dir / "arc.txt").read_text()), nodes, segments)
    reaction = {int(r["NUME_ORDRE"]): float(r["DY"]) for r in parse_table((case_dir / "reaction.txt").read_text())}
    ligament = parse_table((case_dir / "ligament.txt").read_text())
    points_zp = defaultdict(list)
    for r in parse_table((case_dir / "zp.txt").read_text()):
        distance = math.hypot(float(r["COOR_X"]) - A0, float(r["COOR_Y"]))
        points_zp[int(r["NUME_ORDRE"])].append((distance, float(r["V1"]) > 1.0e-6))

    ordres = sorted(o for o in g_par_ordre if o > 0)
    premier = ordres[0]
    p1 = float(params["instants"][premier])
    souplesse = delta[premier] / p1  # delta sur le demi-modele
    j_el_unitaire = sum(g_par_ordre[premier]) / len(g_par_ordre[premier]) / p1**2
    sigma_y = float(loi.contrainte(0.002))

    rows, travail, p_prec, d_prec = [], 0.0, 0.0, 0.0
    for ordre in ordres:
        p_sur_b = float(params["instants"][ordre])
        travail += 0.5 * (p_sur_b + p_prec) * (delta[ordre] - d_prec)  # demi-modele
        p_prec, d_prec = p_sur_b, delta[ordre]
        # Modele complet : deux moities, donc deux fois le travail du demi-modele.
        aire_pl = 2.0 * (travail - 0.5 * p_sur_b**2 * souplesse)
        j_el = j_el_unitaire * p_sur_b**2
        j_e1820 = j_el + ETA_E1820 * aire_pl / B0
        j_values = g_par_ordre[ordre]
        j_large = j_values[INDICE_CONTOUR_LOIN]  # couronne de 8 a 16 mm, voir CONTOURS_J
        k_el = k_from_g_plane_strain(j_el, loi.young, loi.nu)
        k_j = k_from_g_plane_strain(j_large, loi.young, loi.nu)
        rayon = rayon_zone_plastique(points_zp[ordre]) if ordre in points_zp else ""
        rows.append({
            "loi": loi.nom, "acier": loi.acier, "h_tip_mm": h_tip, "element": "QUAD8_SI" if quad else "TRIA6",
            "P_sur_B_N_par_mm": p_sur_b, "P_sur_P0": p_sur_b / charge_limite(sigma_y),
            "reaction_N_par_mm": reaction[ordre], "delta_LL_mm": 2.0 * delta[ordre],
            "J_N_par_mm": j_large,
            "dispersion_J_contours_pct": 100.0 * (max(j_values) - min(j_values)) / j_large,
            "J_E1820_N_par_mm": j_e1820, "ecart_J_E1820_pct": 100.0 * (j_e1820 - j_large) / j_large,
            "K_el_MPa_sqrt_m": mpa_sqrt_m(k_el),
            "K_E399_MPa_sqrt_m": mpa_sqrt_m(k_e399(p_sur_b, W, A_W)),
            "K_J_MPa_sqrt_m": mpa_sqrt_m(k_j),
            "rayon_zone_plastique_mm": rayon, "sigma_y_MPa": sigma_y,
            **{f"J_R{r_inf:g}_{r_sup:g}": value for (r_inf, r_sup), value in zip(CONTOURS_J, j_values)},
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        rows = run_case(LOIS["tanguy2005"], 0.2, [100.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0])
        for r in rows:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                   if k in ("P_sur_B_N_par_mm", "P_sur_P0", "reaction_N_par_mm", "J_N_par_mm", "dispersion_J_contours_pct",
                            "ecart_J_E1820_pct", "K_J_MPa_sqrt_m", "K_el_MPa_sqrt_m", "K_E399_MPa_sqrt_m",
                            "rayon_zone_plastique_mm")})
        return
    instants = [100.0 * i for i in range(1, 31)]
    cas = [("tanguy2005", 0.1, False), ("tanguy2005", 0.05, False), ("tanguy2005", 0.025, False),
           ("tanguy2005", 0.05, True), ("tanguy2005_visqueux", 0.05, False), ("bordet2006", 0.05, False)]
    tout = []
    for nom, h_tip, quad in cas:
        rows = run_case(LOIS[nom], h_tip, instants, quad)
        tout += rows
        print(nom, h_tip, quad, f"K_J final = {rows[-1]['K_J_MPa_sqrt_m']:.2f}", flush=True)
    write_csv(ROOT / "results/etape2/resultats.csv", tout)


if __name__ == "__main__":
    main()
