#!/usr/bin/env python3
"""Figure de synthese de l'etape 3."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"

rows = list(csv.DictReader((ROOT / "results/etape3/resultats.csv").open()))
ssy = list(csv.DictReader((ROOT / "results/etape3/test_plasticite_confinee.csv").open()))


def serie(cas):
    return [r for r in rows if r["cas"] == cas]


def col(data, name):
    return [float(r[name]) for r in data]


fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE, layout="constrained")
for ax in (left, right):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)

for cas, couleur, nom, style in (("r0_2p5um", AQUA, "r0 = 2,5 µm", "-"), ("r0_5um", BLUE, "r0 = 5 µm", "-"),
                                 ("r0_10um", ORANGE, "r0 = 10 µm", "-"),
                                 ("r0_5um_fin", INK, "r0 = 5 µm, maillage fin", ":")):
    d = serie(cas)
    left.plot(col(d, "K_J_etape2_MPa_sqrt_m"), col(d, "Pr"), color=couleur, linewidth=2, linestyle=style, label=nom)
left.set_xlim(0, 115)
left.set_ylim(0, 1)
left.set_xlabel("K_J (MPa√m)", color=MUTED)
left.set_ylabel("Probabilité de rupture", color=MUTED)
left.set_title("Probabilité de rupture prédite à −100 °C", color=INK, loc="left", fontweight="bold")
left.legend(frameon=False, labelcolor=MUTED, fontsize=9)

right.axhline(4, color=MUTED, linewidth=1, linestyle="--")
right.text(112, 4.03, "Théorie confinée : 4", color=MUTED, ha="right", va="bottom", fontsize=9)
for largeur, couleur, nom in (("50.0", BLUE, "W = 50 mm"), ("500.0", ORANGE, "W = 500 mm")):
    d = [r for r in ssy if r["W_mm"] == largeur and float(r["K_MPa_sqrt_m"]) >= 20]
    right.plot(col(d, "K_MPa_sqrt_m"), col(d, "pente"), color=couleur, linewidth=2, label=nom)
right.set_ylim(3.5, 6)
right.set_xlabel("K élastique (MPa√m)", color=MUTED)
right.set_ylabel("Pente d ln(σw^m) / d ln K", color=MUTED)
right.set_title("Test de la plasticité confinée", color=INK, loc="left", fontweight="bold")
right.legend(frameon=False, labelcolor=MUTED, fontsize=9)

fig.savefig(ROOT / "results/etape3/resultats.png", dpi=160)
