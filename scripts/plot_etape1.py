#!/usr/bin/env python3
"""Figure de verification de l'etape 1."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INK, MUTED, GRID, SERIES, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#fcfcfb"

rows = list(csv.DictReader((ROOT / "results/etape1/verification.csv").open()))
arc = [r for r in rows if r["chargement"] == "arc"]
conv = sorted((r for r in arc if float(r["a_W"]) == 0.5), key=lambda r: float(r["h_tip_mm"]))
ratio = sorted((r for r in arc if float(r["h_tip_mm"]) == 0.05), key=lambda r: float(r["a_W"]))

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), facecolor=SURFACE, layout="constrained")
for ax in (left, right):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)

left.plot([float(r["h_tip_mm"]) for r in conv], [float(r["K_G_MPa_sqrt_m"]) for r in conv],
          color=SERIES, linewidth=2, marker="o", markersize=8, markeredgecolor=SURFACE, markeredgewidth=2)
left.set_xscale("log")
left.invert_xaxis()
left.set_xticks([0.4, 0.2, 0.1, 0.05, 0.025], ["0,4", "0,2", "0,1", "0,05", "0,025"])
left.minorticks_off()
left.axhline(float(conv[0]["K_E399_MPa_sqrt_m"]), color=MUTED, linewidth=1, linestyle="--")
left.text(0.4, float(conv[0]["K_E399_MPa_sqrt_m"]) - 0.004, "ASTM E399", color=MUTED, va="top")
left.set_xlabel("Taille de maille en pointe de fissure (mm)", color=MUTED)
left.set_ylabel("K (MPa√m)", color=MUTED)
left.set_title("K converge à a/W = 0,5", color=INK, loc="left", fontweight="bold")

right.axvspan(0.45, 0.55, color=GRID, alpha=0.6, linewidth=0)
right.text(0.5, 0.86, "0,45–0,55 :\ndomaine E399-74", color=MUTED, ha="center", va="top", fontsize=9)
right.axhline(0, color=MUTED, linewidth=1)
right.plot([float(r["a_W"]) for r in ratio], [float(r["ecart_K_G_E399_pct"]) for r in ratio],
           color=SERIES, linewidth=2, marker="o", markersize=8, markeredgecolor=SURFACE, markeredgewidth=2)
right.set_ylim(-1, 1)
right.set_xlabel("a/W", color=MUTED)
right.set_ylabel("Écart à E399 (%)", color=MUTED)
right.set_title("Écart à la formule E399 selon a/W", color=INK, loc="left", fontweight="bold")

fig.savefig(ROOT / "results/etape1/verification.png", dpi=160)
