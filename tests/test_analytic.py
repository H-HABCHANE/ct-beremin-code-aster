import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analytic import f_e399, g_from_k_plane_strain, k_e399, k_from_g_plane_strain


def test_f_e399_valeur_tabulee_a_w_05():
    assert f_e399(0.5) == pytest.approx(9.66, abs=0.01)


def test_k_e399_unites():
    assert k_e399(400.0, 50.0, 0.5) == pytest.approx(400.0 / math.sqrt(50.0) * f_e399(0.5))


def test_irwin_aller_retour():
    g = g_from_k_plane_strain(550.0, 200000.0, 0.3)
    assert k_from_g_plane_strain(g, 200000.0, 0.3) == pytest.approx(550.0)


def test_f_e399_hors_domaine():
    with pytest.raises(ValueError):
        f_e399(0.1)
