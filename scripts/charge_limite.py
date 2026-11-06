#!/usr/bin/env python3
"""Controle numerique de la charge limite : loi parfaitement plastique, effort croissant jusqu'a la perte de convergence."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aster_run import copy_comm, run, write_export  # noqa: E402
from run_etape1 import mesh  # noqa: E402

SIGMA_0, YOUNG = 597.0, 210000.0
A, B0 = 25.0, 25.0  # mm : longueur de fissure et ligament a W = 50 mm, a/W = 0,5


def formule_epri() -> float:
    """Kumar, German et Shih (1981), equations 3-6 et 3-8, deformations planes, par unite d'epaisseur."""
    rapport = 2.0 * A / B0
    eta = math.sqrt(rapport**2 + 2.0 * rapport + 2.0) - (rapport + 1.0)
    return 1.455 * eta * B0 * SIGMA_0


def main() -> None:
    dossier = ROOT / "build/etape4/limite"
    dossier.mkdir(parents=True, exist_ok=True)
    mesh(dossier, 0.5, 0.1, quad=True)
    courbe = [[SIGMA_0 / YOUNG, SIGMA_0], [SIGMA_0 / YOUNG + 10.0, SIGMA_0 * (1.0 + 1.0e-6)]]
    instants = [0.0, 1000.0, 2000.0, 2500.0, 3000.0] + [3000.0 + 25.0 * i for i in range(1, 61)]
    params = {
        "E": YOUNG, "nu": 0.3, "courbe": courbe, "contours": [[8.0, 16.0]], "instants": instants,
        "traction_unitaire": 1.0 / (math.pi * 6.25), "modelisation": "D_PLAN_SI", "ordres_zp": [1],
    }
    (dossier / "params.json").write_text(json.dumps(params))
    comm = copy_comm(ROOT / "aster/ct_elastoplastique.comm", dossier)
    write_export(dossier, comm, {20: "ct.mail", 30: "params.json"},
                 {38: "reaction.txt", 39: "g.txt", 41: "arc.txt", 42: "ligament.txt", 43: "zp.txt"},
                 memory_mb=4000, ncpus=4)
    try:
        run(dossier, timeout_s=900)
    except RuntimeError as erreur:
        print(f"Arret attendu a la perte de convergence : {erreur}")
    journal = (dossier / "docker.log").read_text()
    blocs = re.split(r"Instant de calcul:\s+", journal)[1:]
    converges = [float(b.split()[0]) for b in blocs if "Critère(s) de convergence atteint(s)" in b]
    print(f"Derniere charge convergee : {max(converges):.1f} N/mm")
    print(f"Formule EPRI : {formule_epri():.1f} N/mm")


if __name__ == "__main__":
    main()
