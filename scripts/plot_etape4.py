#!/usr/bin/env python3
"""Figure de synthese de l'etape 4."""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyse_etape4 import CAS, charge, k_au_niveau  # noqa: E402

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE, layout="constrained")
for ax in (left, right):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)

for (nom, (dossier, table)), couleur in zip(CAS.items(), (BLUE, ORANGE, AQUA)):
    d = charge(dossier, table)
    left.plot(d["k"], d["pr"], color=couleur, linewidth=2, label=nom)
reference = charge(*CAS["CT W=50, a/W=0,5 (référence)"])
k_c = k_au_niveau(reference, 0.5)
left.axvline(k_c, color=MUTED, linewidth=1, linestyle="--")
left.text(k_c - 1.5, 0.97, f"K_Jc historique = {k_c:.0f}", color=MUTED, ha="right", va="top", fontsize=9)
left.set_xlim(20, 115)
left.set_ylim(0, 1)
left.set_xlabel("K_J (MPa√m)", color=MUTED)
left.set_ylabel("Probabilité de rupture (Beremin)", color=MUTED)
left.set_title("À K_J égal, Beremin distingue les géométries", color=INK, loc="left", fontweight="bold")
left.legend(frameon=False, labelcolor=MUTED, fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 0.9))

# A B = 6,25 mm, le niveau de 50 % n'est pas atteint par le calcul : ce point est exclu.
rows = [r for r in csv.DictReader((ROOT / "results/etape4/epaisseurs.csv").open()) if r["K_J_Pr_50pct_Beremin"] != "nan"]
b = [float(r["B_mm"]) for r in rows]
for colonne, couleur, nom, style in (("K_J_Pr_50pct_Beremin", BLUE, "Beremin", "-"),
                                     ("K_J_Pr_50pct_Master_Curve", ORANGE, "Correction de type Master Curve", "-"),
                                     ("K_J_historique_sans_effet_B", MUTED, "K_J sans effet d'épaisseur", "--")):
    right.plot(b, [float(r[colonne]) for r in rows], color=couleur, linewidth=2 if style == "-" else 1,
               linestyle=style, marker="o" if style == "-" else None, markersize=8,
               markeredgecolor=SURFACE, markeredgewidth=2, label=nom)
right.set_xscale("log")
right.set_xticks(b, ["12,5", "25", "50", "100"])
right.minorticks_off()
right.set_xlabel("Épaisseur B (mm)", color=MUTED)
right.set_ylabel("K_J à Pr = 50 % (MPa√m)", color=MUTED)
right.set_title("Effet d'épaisseur", color=INK, loc="left", fontweight="bold")
right.legend(frameon=False, labelcolor=MUTED, fontsize=9)

fig.savefig(ROOT / "results/etape4/resultats.png", dpi=160)
