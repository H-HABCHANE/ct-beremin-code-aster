import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from materiau import BORDET, TANGUY, TANGUY_VISQUEUX
from run_etape2 import deplacement_conjugue, rayon_zone_plastique


def test_loi_tanguy_valeurs_limites():
    assert TANGUY.contrainte(0.0) == pytest.approx(589.0)
    assert TANGUY.contrainte(1.0e6) == pytest.approx(589.0 + 184.0 + 472.0)


def test_courbes_traction_strictement_croissantes_et_elastiques_au_depart():
    for loi in (TANGUY, TANGUY_VISQUEUX, BORDET):
        courbe = loi.courbe_traction()
        assert courbe[0][0] == pytest.approx(loi.r0 / loi.young)
        assert all(b[0] > a[0] and b[1] > a[1] for a, b in zip(courbe, courbe[1:]))


def test_surcontrainte_visqueuse():
    assert TANGUY_VISQUEUX.r0 - TANGUY.r0 == pytest.approx(84.0 * 1e-4 ** (1 / 7))


def test_deplacement_conjugue_champ_lineaire():
    # Simpson est exact pour un champ lineaire le long d'un segment droit.
    nodes = {1: (0.0, 0.0), 2: (2.0, 0.0), 3: (1.0, 0.0)}
    rows = [{"NUME_ORDRE": "1", "COOR_X": str(x), "COOR_Y": str(y), "DY": str(3.0 * x)} for x, y in nodes.values()]
    assert deplacement_conjugue(rows, nodes, [(1, 2, 3)])[1] == pytest.approx(3.0)


def test_rayon_zone_plastique_ignore_la_face_arriere():
    points = [(0.1, True), (0.3, True), (0.6, True), (0.8, False), (1.1, False), (25.0, True)]
    assert rayon_zone_plastique(points, pas=0.25) == pytest.approx(0.75)


def test_rayon_zone_plastique_vide():
    assert rayon_zone_plastique([(0.1, False), (0.3, True)], pas=0.25) == 0.0
