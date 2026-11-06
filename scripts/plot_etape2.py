#!/usr/bin/env python3
"""Figure de synthese de l'etape 2."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

rows = list(csv.DictReader((ROOT / "results/etape2/resultats.csv").open()))


def serie(loi, h, element="TRIA6"):
    return [r for r in rows if r["loi"] == loi and float(r["h_tip_mm"]) == h and r["element"] == element]


def col(data, name):
    return [float(r[name]) for r in data]


ref = serie("tanguy2005", 0.05)
fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE, layout="constrained")
for ax in (left, right):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)

x = col(ref, "K_el_MPa_sqrt_m")
left.plot(x, x, color=MUTED, linewidth=1, linestyle="--", label="Élastique (K = K_el)")
for loi, couleur, nom in (("tanguy2005", BLUE, "16MND5, Tanguy et al. (2005)"),
                          ("tanguy2005_visqueux", AQUA, "16MND5 + sur-contrainte visqueuse"),
                          ("bordet2006", ORANGE, "18MND5, Bordet et al. (2006)")):
    d = serie(loi, 0.05)
    left.plot(col(d, "K_el_MPa_sqrt_m"), col(d, "K_J_MPa_sqrt_m"), color=couleur, linewidth=2, label=nom)
left.set_xlabel("K élastique (MPa√m)", color=MUTED)
left.set_ylabel("K_J (MPa√m)", color=MUTED)
left.set_title("La plasticité augmente K_J à forte charge", color=INK, loc="left", fontweight="bold")
left.legend(frameon=False, labelcolor=MUTED, fontsize=9)

right.axhline(0, color=MUTED, linewidth=1)
right.plot(col(ref, "K_J_MPa_sqrt_m"), col(ref, "ecart_J_E1820_pct"), color=BLUE, linewidth=2,
           marker="o", markersize=8, markeredgecolor=SURFACE, markeredgewidth=2, markevery=3)
right.set_xlabel("K_J (MPa√m)", color=MUTED)
right.set_ylabel("Écart J (aire P–δ) / J (CALC_G) (%)", color=MUTED)
right.set_title("Contrôle par l'aire sous la courbe de charge", color=INK, loc="left", fontweight="bold")

fig.savefig(ROOT / "results/etape2/resultats.png", dpi=160)
