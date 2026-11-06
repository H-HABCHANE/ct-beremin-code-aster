#!/usr/bin/env python3
"""Etape 1 : verification de K et de J en elasticite contre la formule ASTM E399."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analytic import g_from_k_plane_strain, k_e399, k_from_g_plane_strain, mpa_sqrt_m  # noqa: E402
from aster_run import copy_comm, run, write_export  # noqa: E402
from msh_to_aster import read_msh22, write_aster  # noqa: E402

W, YOUNG, NU = 50.0, 200000.0, 0.3
P_OVER_B = 400.0  # N/mm d'epaisseur (P = 10 kN, B = 25 mm)
CONTOURS = [[0.5, 1.0], [1.0, 2.0], [2.0, 4.0], [4.0, 8.0], [1.0, 10.0]]
MESH_SIZES = [0.4, 0.2, 0.1, 0.05, 0.025]
RATIOS = [0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]
REF_H = 0.05
FINE_H = 0.025


def mesh(case_dir: Path, a_w: float, h_tip: float, quad: bool = False, largeur: float = W, h_far: float = 2.5) -> int:
    msh = case_dir / "ct.msh"
    subprocess.run(
        ["gmsh", str(ROOT / "mesh/ct_half.geo"), "-2", "-format", "msh2",
         "-setnumber", "aW", str(a_w), "-setnumber", "W", str(largeur), "-setnumber", "h_tip", str(h_tip),
         "-setnumber", "h_far", str(h_far),
         "-setnumber", "quad", "1" if quad else "0",
         "-o", str(msh)],
        check=True, capture_output=True,
    )
    nodes, elements, groups = read_msh22(msh.read_text())
    (case_dir / "ct.mail").write_text(write_aster(nodes, elements, groups))
    return len(nodes)


def parse_table(text: str) -> list[dict[str, str]]:
    """Lit les tables IMPR_TABLE concatenees dans un meme fichier."""
    rows, header = [], None
    for line in text.splitlines():
        cols = line.split()
        if not cols or line.startswith("#"):
            header = None if not cols else header
            continue
        if header is None:
            header = cols
        elif len(cols) == len(header):
            rows.append(dict(zip(header, cols)))
    return rows


def run_case(a_w: float, h_tip: float, chargement: str = "arc") -> dict:
    case_dir = ROOT / "build" / f"aW{a_w:.4f}_h{h_tip:.3f}_{chargement}"
    case_dir.mkdir(parents=True, exist_ok=True)
    n_nodes = mesh(case_dir, a_w, h_tip)
    params = {
        "E": YOUNG, "nu": NU, "contours": CONTOURS, "chargement": chargement, "p_sur_b": P_OVER_B,
        "traction": P_OVER_B / (math.pi * 0.125 * W),  # effort reparti sur le demi-trou superieur
    }
    (case_dir / "params.json").write_text(json.dumps(params))
    comm = copy_comm(ROOT / "aster/ct_elastique.comm", case_dir)
    write_export(case_dir, comm, {20: "ct.mail", 30: "params.json"},
                 {38: "reaction.txt", 39: "g.txt", 40: "energie.txt"})
    run(case_dir)
    g_rows = parse_table((case_dir / "g.txt").read_text())
    if len(g_rows) != len(CONTOURS):
        raise RuntimeError(f"{case_dir}: {len(g_rows)} lignes de G pour {len(CONTOURS)} contours")
    k_ref = k_e399(P_OVER_B, W, a_w)
    g_values = [float(row["G"]) for row in g_rows]
    k1_values = [float(row["K1"]) for row in g_rows]
    g_mean = sum(g_values) / len(g_values)
    k_irwin = k_from_g_plane_strain(g_mean, YOUNG, NU)
    k1_mean = sum(k1_values) / len(k1_values)
    return {
        "a_W": a_w, "h_tip_mm": h_tip, "chargement": chargement, "noeuds": n_nodes,
        "G_moyen_N_par_mm": g_mean,
        "dispersion_G_contours_pct": 100.0 * (max(g_values) - min(g_values)) / g_mean,
        "K_G_MPa_sqrt_m": mpa_sqrt_m(k_irwin),
        "K1_direct_MPa_sqrt_m": mpa_sqrt_m(k1_mean),
        "K_E399_MPa_sqrt_m": mpa_sqrt_m(k_ref),
        "ecart_K_G_E399_pct": 100.0 * (k_irwin - k_ref) / k_ref,
        "ecart_K1_G_pct": 100.0 * (k1_mean - k_irwin) / k_irwin,
        "reaction_N_par_mm": float(parse_table((case_dir / "reaction.txt").read_text())[0]["DY"]),
        "energie_demi_modele_N": float(parse_table((case_dir / "energie.txt").read_text())[0]["TOTALE"]),
        "G_E399_N_par_mm": g_from_k_plane_strain(k_ref, YOUNG, NU),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="un seul cas, pour tester la chaine")
    args = parser.parse_args()
    out_dir = ROOT / "results/etape1"
    if args.smoke:
        print(run_case(0.5, 0.2))
        return
    cases = [(0.5, h, "arc") for h in MESH_SIZES]
    cases += [(r, REF_H, "arc") for r in RATIOS if r != 0.5]
    cases += [(0.5, REF_H, "point")]
    rows = []
    for a_w, h_tip, chargement in cases:
        rows.append(run_case(a_w, h_tip, chargement))
        print(f"a/W={a_w} h={h_tip} {chargement}: ecart K/E399 = {rows[-1]['ecart_K_G_E399_pct']:+.3f} %", flush=True)
    write_csv(out_dir / "verification.csv", rows)

    # Controle independant de CALC_G : G = 2 dU/da a effort impose (U sur le demi-modele).
    # La difference centree a une erreur en da^2 : on verifie sa convergence.
    check = []
    # Chaque longueur de fissure a son propre maillage : le bruit de remaillage augmente quand da diminue.
    for h_tip in (REF_H, FINE_H):
        g_calc = next(r for r in rows if r["h_tip_mm"] == h_tip and r["a_W"] == 0.5)["G_moyen_N_par_mm"]
        for da_ratio in (0.02, 0.01, 0.005, 0.0025):
            low = run_case(0.5 - da_ratio, h_tip)["energie_demi_modele_N"]
            high = run_case(0.5 + da_ratio, h_tip)["energie_demi_modele_N"]
            g_energy = 2.0 * (high - low) / (2.0 * da_ratio * W)
            check.append({"h_tip_mm": h_tip, "da_mm": da_ratio * W, "G_CALC_G": g_calc,
                          "G_energie_dU_da": g_energy, "ecart_pct": 100.0 * (g_energy - g_calc) / g_calc})
            print(h_tip, da_ratio * W, f"{check[-1]['ecart_pct']:+.3f} %", flush=True)
    write_csv(out_dir / "controle_energie.csv", check)


if __name__ == "__main__":
    main()
