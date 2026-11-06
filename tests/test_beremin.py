import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from beremin import Parametres, contrainte_effective, contrainte_principale_max, weibull

PARAMS = Parametres(m=10.0, sigma_u=2000.0, v0=1.0e-3)


def test_principale_cisaillement_pur():
    assert contrainte_principale_max(np.array(0.0), np.array(0.0), np.array(0.0), np.array(100.0)) == pytest.approx(100.0)


def test_principale_hors_plan_dominante():
    assert contrainte_principale_max(np.array(100.0), np.array(50.0), np.array(300.0), np.array(0.0)) == pytest.approx(300.0)


def test_weibull_contrainte_uniforme_cas_analytique():
    # sigma_w^m = (coef / V0) * aire * s^m pour une contrainte uniforme.
    s, aire, coef = 1500.0, 0.2, 50.0
    sigma_w, pr = weibull(np.array([[s, s]]), np.array([aire / 2, aire / 2]), coef, PARAMS)
    attendu = s * (coef * aire / PARAMS.v0) ** (1.0 / PARAMS.m)
    assert sigma_w[0] == pytest.approx(attendu)
    assert pr[0] == pytest.approx(1.0 - np.exp(-((attendu / PARAMS.sigma_u) ** PARAMS.m)))


def test_weibull_garde_le_maximum_historique():
    sigma = np.array([[1000.0], [1500.0], [800.0]])
    sigma_w, _ = weibull(sigma, np.array([1.0]), 1.0, PARAMS)
    assert sigma_w[2] == pytest.approx(sigma_w[1])
    assert sigma_w[1] > sigma_w[0]


def test_effective_nulle_sans_plasticite_active():
    p = np.array([[0.0], [0.01], [0.01]])
    s = np.array([[500.0], [800.0], [900.0]])
    sortie = contrainte_effective(s, p, Parametres(10.0, 2000.0, 1e-3, k=None, p_c=0.0))
    assert sortie[:, 0].tolist() == [0.0, 800.0, 0.0]


def test_effective_correction_plastique():
    p = np.array([[0.4]])
    sortie = contrainte_effective(np.array([[1000.0]]), p, Parametres(10.0, 2000.0, 1e-3, k=4.0, p_c=1e-4))
    assert sortie[0, 0] == pytest.approx(1000.0 * np.exp(-0.1))


def test_effective_seuil_p_c():
    p = np.array([[5.0e-5]])
    sortie = contrainte_effective(np.array([[1000.0]]), p, Parametres(10.0, 2000.0, 1e-3, k=None, p_c=1e-4))
    assert sortie[0, 0] == 0.0
