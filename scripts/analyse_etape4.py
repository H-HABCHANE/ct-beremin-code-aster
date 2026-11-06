#!/usr/bin/env python3
"""Etape 4 : ce que predisent K_J (approche historique) et Beremin selon la geometrie, la loi et l'epaisseur."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyse_etape3 import interpole_k  # noqa: E402
from beremin import TANGUY, calcule  # noqa: E402
from run_etape1 import write_csv  # noqa: E402
from run_etape3 import COEF_MULT, EPAISSEUR, REPERTOIRE  # noqa: E402

SORTIE = ROOT / "results/etape4"
CIBLES = (0.05, 0.10, 0.50, 0.90)


def k_j_etape2(loi: str):
    rows = [r for r in csv.DictReader((ROOT / "results/etape2/resultats.csv").open())
            if r["loi"] == loi and float(r["h_tip_mm"]) == 0.05 and r["element"] == "TRIA6"]
    return (np.array([0.0] + [float(r["P_sur_B_N_par_mm"]) for r in rows]),
            np.array([0.0] + [float(r["K_J_MPa_sqrt_m"]) for r in rows]))


def k_j_geometrie(nom: str):
    rows = [r for r in csv.DictReader((SORTIE / "k_j_geometries.csv").open()) if r["geometrie"] == nom]
    return (np.array([0.0] + [float(r["P_sur_B_N_par_mm"]) for r in rows]),
            np.array([0.0] + [float(r["K_J_MPa_sqrt_m"]) for r in rows]))


CAS = {
    "CT W=50, a/W=0,5 (référence)": ("r0_5um", k_j_etape2("tanguy2005")),
    "CT W=500, a/W=0,5": ("W500_r0_5um", k_j_geometrie("W500")),
    "CT W=50, a/W=0,3": ("aW0p3_r0_5um", k_j_geometrie("aW0p3")),
}
LOIS = {
    "16MND5, Tanguy et al. (référence)": ("r0_5um", k_j_etape2("tanguy2005")),
    "16MND5 + sur-contrainte de 22 MPa": ("loi_visqueuse_r0_5um", k_j_etape2("tanguy2005_visqueux")),
    "18MND5, Bordet et al.": ("loi_18MND5_r0_5um", k_j_etape2("bordet2006")),
}


def charge(dossier: str, reference, coef: float = COEF_MULT):
    chemin = REPERTOIRE / dossier / "champs.npz"
    if not chemin.exists():
        raise FileNotFoundError(chemin)
    champs = dict(np.load(chemin))
    p_ref, k_ref = reference
    p = champs["instants"][1:]
    sw, pr = calcule(champs, coef, TANGUY)
    return {"champs": champs, "p": p, "k": np.interp(p, p_ref, k_ref), "sw": sw[1:], "pr": pr[1:]}


def k_au_niveau(d, cible):
    return interpole_k(cible, d["pr"], d["k"])


def pr_a_k(d, k):
    if k > d["k"][-1]:
        return math.nan
    return float(np.interp(k, d["k"], d["pr"]))


def pente_weibull(d, k_min=40.0):
    choix = d["k"] >= k_min
    return float(np.polyfit(np.log(d["k"][choix]), np.log(-np.log1p(-d["pr"][choix])), 1)[0])


def main() -> None:
    SORTIE.mkdir(parents=True, exist_ok=True)
    reference = charge(*CAS["CT W=50, a/W=0,5 (référence)"])
    k_c = k_au_niveau(reference, 0.5)

    geometries = []
    for nom, (dossier, table) in CAS.items():
        d = charge(dossier, table)
        ligne = {"cas": nom, "K_J_max": float(d["k"][-1])}
        ligne.update({f"K_J_Pr_{int(c * 100)}pct": k_au_niveau(d, c) for c in CIBLES})
        ligne["Pr_a_K_J_ref_60"] = pr_a_k(d, 60.0)
        ligne["Pr_au_K_Jc_historique"] = pr_a_k(d, k_c)
        ligne["pente_Weibull_Pr_K_J"] = pente_weibull(d)
        geometries.append(ligne)
    write_csv(SORTIE / "geometries.csv", geometries)

    lois = []
    p_ref = reference["p"][np.argmin(np.abs(reference["p"] - 1500.0))]
    k_ref_1500 = float(np.interp(p_ref, reference["p"], reference["k"]))
    k_ref_2500 = float(reference["k"][-1])
    for nom, (dossier, table) in LOIS.items():
        d = charge(dossier, table)
        lois.append({
            "loi": nom,
            "K_J_a_P_1500": float(np.interp(p_ref, d["p"], d["k"])),
            "ecart_K_J_a_P_1500_pct": 100.0 * (float(np.interp(p_ref, d["p"], d["k"])) / k_ref_1500 - 1.0),
            "K_J_a_P_2500": float(d["k"][-1]),
            "ecart_K_J_a_P_2500_pct": 100.0 * (float(d["k"][-1]) / k_ref_2500 - 1.0),
            "sigma_w_a_P_1500": float(np.interp(p_ref, d["p"], d["sw"])),
            "Pr_a_P_1500": float(np.interp(p_ref, d["p"], d["pr"])),
            "K_J_Pr_50pct": k_au_niveau(d, 0.5),
            "ecart_K_J_Pr_50pct": 100.0 * (k_au_niveau(d, 0.5) / k_c - 1.0),
        })
    write_csv(SORTIE / "lois.csv", lois)

    epaisseurs = []
    table = CAS["CT W=50, a/W=0,5 (référence)"][1]
    k_min = 20.0  # MPa.sqrt(m), correction d'epaisseur de la Master Curve (ASTM E1921)
    for b in (6.25, 12.5, 25.0, 50.0, 100.0):
        d = charge("r0_5um", table, COEF_MULT * b / EPAISSEUR)
        k50 = k_au_niveau(d, 0.5)
        epaisseurs.append({
            "B_mm": b, "K_J_Pr_50pct_Beremin": k50,
            "K_J_Pr_50pct_echelle_B_quart": k_c * (EPAISSEUR / b) ** 0.25,
            "K_J_Pr_50pct_Master_Curve": k_min + (k_c - k_min) * (EPAISSEUR / b) ** 0.25,
            "K_J_historique_sans_effet_B": k_c,
        })
    write_csv(SORTIE / "epaisseurs.csv", epaisseurs)

    for fichier, lignes in (("geometries", geometries), ("lois", lois), ("epaisseurs", epaisseurs)):
        print(fichier)
        for ligne in lignes:
            print(" ", {c: (round(v, 3) if isinstance(v, float) else v) for c, v in ligne.items()})


if __name__ == "__main__":
    main()
