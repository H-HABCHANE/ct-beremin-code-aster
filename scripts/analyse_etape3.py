#!/usr/bin/env python3
"""Depouillement de l'etape 3 : sigma_w et Pr en fonction de K_J, controles et sensibilites."""
from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analytic import k_e399, mpa_sqrt_m  # noqa: E402
from beremin import TANGUY, Parametres, calcule  # noqa: E402
from run_etape1 import parse_table, write_csv  # noqa: E402
from run_etape2 import A_W, ETA_E1820, B0, arc_elements, deplacement_conjugue  # noqa: E402
from run_etape3 import CAS, COEF_MULT, REPERTOIRE  # noqa: E402

E, NU = 210000.0, 0.3
SORTIE = ROOT / "results/etape3"
CAS_SENSIBILITE = "r0_5um"  # cas de reference : atteint la charge maximale et converge en maillage


def reference_etape2():
    """K_J(P/B) de l'etape 2 : fissure vive, loi de reference, h = 0,05 mm."""
    rows = [r for r in csv.DictReader((ROOT / "results/etape2/resultats.csv").open())
            if r["loi"] == "tanguy2005" and float(r["h_tip_mm"]) == 0.05 and r["element"] == "TRIA6"]
    return tuple(np.array([0.0] + [float(r[c]) for r in rows])
                 for c in ("P_sur_B_N_par_mm", "K_J_MPa_sqrt_m", "K_el_MPa_sqrt_m"))


def charge(nom: str):
    dossier = REPERTOIRE / nom
    if not (dossier / "champs.npz").exists():
        return None
    champs = dict(np.load(dossier / "champs.npz"))
    nodes, segments = arc_elements((dossier / "ct.msh").read_text())
    delta = deplacement_conjugue(parse_table((dossier / "arc.txt").read_text()), nodes, segments)
    aster = {int(r["NUME_ORDRE"]): float(r["SIGMA_WEIBULL"]) for r in parse_table((dossier / "beremin.txt").read_text())}
    reaction = {int(r["NUME_ORDRE"]): float(r["DY"]) for r in parse_table((dossier / "reaction.txt").read_text())}
    return {"champs": champs, "delta": delta, "aster": aster, "reaction": reaction, "noeuds": len(nodes)}


def j_energetique(p_sur_b, delta, k_el_sur_p):
    """J par l'aire P-delta (meme methode qu'a l'etape 2)."""
    souplesse = delta[1] / p_sur_b[1]
    travail = np.concatenate(([0.0], np.cumsum(0.5 * (p_sur_b[1:] + p_sur_b[:-1]) * np.diff(delta))))
    aire_pl = 2.0 * (travail - 0.5 * p_sur_b**2 * souplesse)
    k_el = k_el_sur_p * p_sur_b * math.sqrt(1000.0)
    return k_el**2 * (1.0 - NU**2) / E + ETA_E1820 * aire_pl / B0


def pente_locale(x, y):
    return np.gradient(np.log(y), np.log(x))


def interpole_k(pr_cible, pr, k):
    """K_J pour une probabilite donnee, par interpolation en ln(-ln(1-Pr)) contre ln K."""
    choix = (pr > 0.0) & (pr < 1.0)
    if not choix.any() or not pr[choix][0] <= pr_cible <= pr[choix][-1]:
        return math.nan
    return float(np.exp(np.interp(np.log(-np.log1p(-pr_cible)), np.log(-np.log1p(-pr[choix])), np.log(k[choix]))))


def test_plasticite_confinee() -> None:
    """Pente de ln(sigma_w^m) contre ln(K), a meme K, pour W = 50 et 500 mm (theorie confinee : 4)."""
    brut = Parametres(TANGUY.m, TANGUY.sigma_u, TANGUY.v0)
    lignes = []
    for nom, largeur in (("r0_5um", 50.0), ("W500_r0_5um", 500.0)):
        dossier = REPERTOIRE / nom
        if not (dossier / "champs.npz").exists():
            print(f"{nom} : pas de resultats")
            continue
        champs = dict(np.load(dossier / "champs.npz"))
        p_sur_b = champs["instants"][1:]
        k = np.array([mpa_sqrt_m(k_e399(p, largeur, A_W)) for p in p_sur_b])
        sigma_w, _pr = calcule(champs, COEF_MULT, brut)
        sw_m = sigma_w[1:] ** TANGUY.m
        pente = pente_locale(k, sw_m)
        for ki, swi, pi in zip(k, sigma_w[1:], pente):
            lignes.append({"W_mm": largeur, "K_MPa_sqrt_m": ki, "sigma_w_MPa": swi, "pente": pi})
    if lignes:
        write_csv(SORTIE / "test_plasticite_confinee.csv", lignes)


def main() -> None:
    test_plasticite_confinee()
    SORTIE.mkdir(parents=True, exist_ok=True)
    p_ref, kj_ref, kel_ref = reference_etape2()
    k_el_sur_p = kel_ref[1] / p_ref[1]
    resultats, synthese = [], []
    donnees = {nom: charge(nom) for nom, geo in CAS.items() if "W" not in geo}
    for nom, d in donnees.items():
        if d is None:
            print(f"{nom} : pas de resultats")
            continue
        c = d["champs"]
        ordres = c["ordres"].astype(int)[1:]
        p_sur_b = c["instants"][1:]
        delta = np.array([0.0] + [d["delta"][int(o)] for o in ordres])
        p_tout = np.concatenate(([0.0], p_sur_b))
        k_j = np.interp(p_sur_b, p_ref, kj_ref)
        j = j_energetique(p_tout, delta, k_el_sur_p)[1:]
        k_j_energie = np.sqrt(j * E / (1.0 - NU**2)) / math.sqrt(1000.0)

        brut = Parametres(TANGUY.m, TANGUY.sigma_u, TANGUY.v0)
        sw_brut, _ = calcule(c, COEF_MULT, brut)
        sw_aster = np.array([d["aster"][int(o)] for o in c["ordres"].astype(int)])
        ecart_aster = np.abs(sw_brut[1:] / sw_aster[1:] - 1.0).max() if sw_aster[-1] else math.nan

        sw, pr = calcule(c, COEF_MULT, TANGUY)
        sw, pr = sw[1:], pr[1:]
        pente = pente_locale(k_j, sw**TANGUY.m)
        for i in range(len(p_sur_b)):
            resultats.append({
                "cas": nom, "P_sur_B_N_par_mm": p_sur_b[i], "K_J_etape2_MPa_sqrt_m": k_j[i],
                "K_J_energie_MPa_sqrt_m": k_j_energie[i], "J_energie_N_par_mm": j[i],
                "reaction_N_par_mm": d["reaction"][int(ordres[i])], "sigma_w_MPa": sw[i], "Pr": pr[i],
                "sigma_w_sans_correction_MPa": sw_brut[1:][i], "sigma_w_POST_BEREMIN_MPa": sw_aster[1:][i],
                "pente_sw_m_K": pente[i],
            })
        synthese.append({
            "cas": nom, "noeuds": d["noeuds"], "r0_mm": CAS[nom]["r0"], "n_pas": len(p_sur_b),
            "P_max_N_par_mm": float(p_sur_b[-1]), "K_J_max": float(k_j[-1]),
            "ecart_max_Python_POST_BEREMIN": ecart_aster,
            "K_J_Pr_5pct": interpole_k(0.05, pr, k_j), "K_J_Pr_50pct": interpole_k(0.5, pr, k_j),
            "K_J_Pr_95pct": interpole_k(0.95, pr, k_j),
            "sigma_w_a_K_J_60": float(np.exp(np.interp(np.log(60.0), np.log(k_j), np.log(sw)))),
            "Pr_a_K_J_60": float(np.interp(60.0, k_j, pr)),
        })
    write_csv(SORTIE / "resultats.csv", resultats)
    write_csv(SORTIE / "synthese.csv", synthese)

    # Sensibilites : un seul calcul mecanique, post-traitements differents.
    d = donnees.get(CAS_SENSIBILITE)
    if d is not None:
        c = d["champs"]
        p_sur_b = c["instants"][1:]
        k_j = np.interp(p_sur_b, p_ref, kj_ref)
        lignes = []
        variantes = {
            "Tanguy 2005 (m=17,8 ; su=2925 ; k=4 ; pc=1e-4)": (TANGUY, COEF_MULT),
            "sans correction exp(-p/k)": (replace(TANGUY, k=None), COEF_MULT),
            "pc = 0": (replace(TANGUY, p_c=0.0), COEF_MULT),
            "m = 16": (replace(TANGUY, m=16.0), COEF_MULT),
            "m = 20": (replace(TANGUY, m=20.0), COEF_MULT),
            "su -5 %": (replace(TANGUY, sigma_u=0.95 * TANGUY.sigma_u), COEF_MULT),
            "su +5 %": (replace(TANGUY, sigma_u=1.05 * TANGUY.sigma_u), COEF_MULT),
            "epaisseur 12,5 mm": (TANGUY, COEF_MULT / 2.0),
            "epaisseur 50 mm": (TANGUY, COEF_MULT * 2.0),
        }
        for libelle, (params, coef) in variantes.items():
            sw, pr = calcule(c, coef, params)
            sw, pr = sw[1:], pr[1:]
            lignes.append({"variante": libelle, "K_J_Pr_5pct": interpole_k(0.05, pr, k_j),
                           "K_J_Pr_50pct": interpole_k(0.5, pr, k_j), "K_J_Pr_95pct": interpole_k(0.95, pr, k_j)})
        write_csv(SORTIE / "sensibilites.csv", lignes)
    print(json.dumps(synthese, indent=1, default=float))


if __name__ == "__main__":
    main()
