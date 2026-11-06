#!/usr/bin/env python3
"""Etape 4 : calculs pour comparer K_J et Beremin selon la geometrie et la loi d'ecrouissage."""
from __future__ import annotations

import json
import math
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analytic import f_e399, k_e399, k_from_g_plane_strain, mpa_sqrt_m  # noqa: E402
from aster_run import copy_comm, run, write_export  # noqa: E402
from materiau import BORDET, TANGUY, TANGUY_VISQUEUX  # noqa: E402
from run_etape1 import mesh, parse_table, write_csv  # noqa: E402
from run_etape3 import lance  # noqa: E402

MAILLAGE = {"r0": 0.005, "n_theta": 16, "n_r": 36, "q_r": 1.17}

# Calculs a pointe emoussee, grandes deformations : r0 = 5 um, comme le cas de reference de l'etape 3.
BEREMIN = {
    "aW0p3_r0_5um": ({**MAILLAGE, "aW": 0.3}, TANGUY),
    "loi_visqueuse_r0_5um": (MAILLAGE, TANGUY_VISQUEUX),
    "loi_18MND5_r0_5um": (MAILLAGE, BORDET),
}

# K_J de reference pour les nouvelles geometries : fissure vive, petites deformations, CALC_G.
# La couronne la plus eloignee reste a distance du trou : environ 14 mm de la pointe a a/W = 0,3.
K_J_GEOMETRIES = {
    "W500": {"W": 500.0, "aW": 0.5, "h_far": 25.0, "contours": [[1.0, 2.0], [2.0, 4.0], [4.0, 8.0], [8.0, 16.0]]},
    "aW0p3": {"W": 50.0, "aW": 0.3, "h_far": 2.5, "contours": [[1.0, 2.0], [2.0, 4.0], [4.0, 8.0], [6.0, 12.0]]},
}
K_CIBLES = [4.3 * i for i in range(1, 27)]  # K elastiques en MPa.sqrt(m), jusqu'a environ 112


def calcul_k_j(nom: str, conf: dict) -> list[dict]:
    largeur, a_w = conf["W"], conf["aW"]
    dossier = ROOT / "build" / "etape4" / f"KJ_{nom}"
    dossier.mkdir(parents=True, exist_ok=True)
    mesh(dossier, a_w, 0.05, largeur=largeur, h_far=conf["h_far"])
    unite = mpa_sqrt_m(k_e399(1.0, largeur, a_w))
    instants = [0.0] + [round(k / unite, 3) for k in K_CIBLES]
    params = {
        "E": TANGUY.young, "nu": TANGUY.nu, "courbe": TANGUY.courbe_traction(), "contours": conf["contours"],
        "instants": instants, "traction_unitaire": 1.0 / (math.pi * 0.125 * largeur), "modelisation": "D_PLAN",
        "ordres_zp": [len(instants) - 1],
    }
    (dossier / "params.json").write_text(json.dumps(params))
    comm = copy_comm(ROOT / "aster/ct_elastoplastique.comm", dossier)
    write_export(dossier, comm, {20: "ct.mail", 30: "params.json"},
                 {38: "reaction.txt", 39: "g.txt", 41: "arc.txt", 42: "ligament.txt", 43: "zp.txt"},
                 memory_mb=4000, ncpus=4)
    run(dossier, timeout_s=3600)
    par_ordre: dict[int, list[float]] = {}
    for row in parse_table((dossier / "g.txt").read_text()):
        par_ordre.setdefault(int(row["NUME_ORDRE"]), []).append(float(row["G"]))
    lignes = []
    for ordre, valeurs in sorted(par_ordre.items()):
        if ordre == 0:
            continue
        j = valeurs[-1]
        lignes.append({
            "geometrie": nom, "P_sur_B_N_par_mm": instants[ordre], "K_el_E399_MPa_sqrt_m": K_CIBLES[ordre - 1],
            "J_N_par_mm": j, "K_J_MPa_sqrt_m": mpa_sqrt_m(k_from_g_plane_strain(j, TANGUY.young, TANGUY.nu)),
            "ecart_contours_externes_pct": 100.0 * (valeurs[-1] - valeurs[-2]) / j,
        })
    return lignes


def main() -> None:
    print("f(a/W) :", {a: round(f_e399(a), 3) for a in (0.3, 0.5)}, flush=True)
    lignes = []
    for nom, conf in K_J_GEOMETRIES.items():
        lignes += calcul_k_j(nom, conf)
        print(f"K_J {nom} : termine, K_J final = {lignes[-1]['K_J_MPa_sqrt_m']:.2f}", flush=True)
    write_csv(ROOT / "results/etape4/k_j_geometries.csv", lignes)

    def tente(item):
        nom, (geo, loi) = item
        print(f"{nom} : demarre", flush=True)
        try:
            print(lance(nom, geo, loi), flush=True)
        except Exception:
            print(f"{nom} : ECHEC\n{traceback.format_exc(limit=1)}", flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(tente, BEREMIN.items()))


if __name__ == "__main__":
    main()
