#!/usr/bin/env python3
"""Etape 3 : calculs en grandes deformations a pointe emoussee. Le depouillement est dans analyse_etape3.py."""
from __future__ import annotations

import json
import math
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analytic import f_e399  # noqa: E402
from aster_run import copy_comm, run, write_export  # noqa: E402
from beremin import TANGUY as BEREMIN  # noqa: E402
from materiau import TANGUY  # noqa: E402
from msh_to_aster import read_msh22, write_aster  # noqa: E402

W = 50.0
EPAISSEUR = 25.0  # mm, comme aux etapes 1 et 2
COEF_MULT = 2.0 * EPAISSEUR  # demi-modele : facteur 2 de symetrie, multiplie par l'epaisseur
P_MAX = 2500.0
REPERTOIRE = ROOT / "build" / "etape3"

# Pas proportionnels a la charge : c'est dJ/J qui pilote l'increment de deformation en pointe.
INSTANTS = [0.0]
while INSTANTS[-1] * 1.06 < P_MAX:
    INSTANTS.append(round(max(200.0, INSTANTS[-1] * 1.06), 3))
INSTANTS.append(P_MAX)

# Taille radiale de la premiere maille environ r0/7, pour eviter leur ecrasement a fort emoussement.
CAS = {
    "r0_2p5um": {"r0": 0.0025, "n_theta": 16, "n_r": 36, "q_r": 1.19},
    "r0_5um": {"r0": 0.005, "n_theta": 16, "n_r": 36, "q_r": 1.17},
    "r0_10um": {"r0": 0.010, "n_theta": 16, "n_r": 36, "q_r": 1.15},
    "r0_5um_fin": {"r0": 0.005, "n_theta": 24, "n_r": 48, "q_r": 1.12},
    # Eprouvette 10 fois plus grande, aux memes K : test de la plasticite confinee (sigma_w^m en K^4).
    "W500_r0_5um": {"r0": 0.005, "n_theta": 16, "n_r": 36, "q_r": 1.17, "W": 500.0, "h_far": 25.0},
}


def controle_aire(nodes, elements, W: float = W) -> float:
    """Refuse un maillage incomplet : Gmsh peut ecrire un fichier meme si une surface a echoue."""
    aire = 0.0
    for _tag, kind, conn in elements:
        sommets = {"QUAD8": 4, "TRIA6": 3}.get(kind)
        if sommets is None:
            continue
        pts = [nodes[i] for i in conn[:sommets]]
        signee = 0.5 * sum(pts[i][0] * pts[(i + 1) % sommets][1] - pts[(i + 1) % sommets][0] * pts[i][1]
                           for i in range(sommets))
        if signee <= 0.0:
            raise ValueError("Maille d'aire negative ou nulle")
        aire += signee
    theorique = 1.25 * W * 0.6 * W - math.pi * (0.125 * W) ** 2
    if abs(aire - theorique) / theorique > 0.01:
        raise ValueError(f"Aire du maillage {aire:.1f} mm2 pour {theorique:.1f} mm2 attendus")
    return aire


def lance(nom: str, geo: dict, loi=TANGUY) -> str:
    dossier = REPERTOIRE / nom
    dossier.mkdir(parents=True, exist_ok=True)
    options = [arg for cle, val in geo.items() for arg in ("-setnumber", cle, str(val))]
    subprocess.run(["gmsh", str(ROOT / "mesh/ct_half_emousse.geo"), "-2", "-format", "msh2", *options,
                    "-o", str(dossier / "ct.msh")], check=True, capture_output=True)
    nodes, elements, groupes = read_msh22((dossier / "ct.msh").read_text())
    largeur = geo.get("W", W)
    controle_aire(nodes, elements, largeur)
    # P/B = K sqrt(W) / f(a/W) : memes K elastiques pour toutes les geometries.
    echelle = math.sqrt(largeur / W) * f_e399(0.5) / f_e399(geo.get("aW", 0.5))
    (dossier / "ct.mail").write_text(write_aster(nodes, elements, groupes))
    params = {
        "E": loi.young, "nu": loi.nu, "courbe": loi.courbe_traction(),
        "instants": [round(instant * echelle, 3) for instant in INSTANTS],
        "traction_unitaire": 1.0 / (math.pi * 0.125 * largeur), "modelisation": "D_PLAN_INCO_UPG", "coef_mult": COEF_MULT,
        "m": BEREMIN.m, "v0": BEREMIN.v0, "sigma_u": BEREMIN.sigma_u, "subd_niveau": 4,
    }
    (dossier / "params.json").write_text(json.dumps(params))
    (dossier / "geometrie.json").write_text(json.dumps(geo))
    comm = copy_comm(ROOT / "aster/ct_beremin.comm", dossier)
    write_export(dossier, comm, {20: "ct.mail", 30: "params.json"},
                 {38: "reaction.txt", 41: "arc.txt", 44: "beremin.txt", 50: "champs.npz", 51: "etat.json"},
                 memory_mb=5000, ncpus=4)
    run(dossier, timeout_s=3 * 3600)
    etat = json.loads((dossier / "etat.json").read_text())
    return f"{nom} : termine, {len(nodes)} noeuds, etat {etat}"


def tente(item) -> str:
    nom, geo = item
    print(f"{nom} : demarre", flush=True)
    try:
        message = lance(nom, geo)
    except Exception:
        message = f"{nom} : ECHEC\n{traceback.format_exc(limit=1)}"
    print(message, flush=True)
    return message


def main() -> None:
    choix = sys.argv[1:] or [nom for nom, geo in CAS.items() if "W" not in geo]
    print(f"{len(INSTANTS) - 1} pas de {INSTANTS[1]} a {INSTANTS[-1]} N/mm", flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(tente, [(nom, CAS[nom]) for nom in choix]))


if __name__ == "__main__":
    main()
