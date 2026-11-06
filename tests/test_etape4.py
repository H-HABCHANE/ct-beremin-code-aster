import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyse_etape3 import interpole_k
from analyse_etape4 import pente_weibull, pr_a_k


def distribution(pente=4.0, k0=100.0):
    k = np.linspace(40.0, 120.0, 50)
    return {"k": k, "pr": 1.0 - np.exp(-((k / k0) ** pente))}


def test_pente_weibull_retrouve_la_pente():
    assert pente_weibull(distribution(pente=4.5)) == pytest.approx(4.5)


def test_interpolation_de_k_sur_une_loi_de_weibull():
    d = distribution(pente=4.0, k0=100.0)
    attendu = 100.0 * np.log(2.0) ** 0.25
    assert interpole_k(0.5, d["pr"], d["k"]) == pytest.approx(attendu)


def test_niveau_hors_du_calcul_non_extrapole():
    d = distribution()
    assert np.isnan(interpole_k(0.999999, d["pr"], d["k"]))
    assert np.isnan(pr_a_k(d, 200.0))
