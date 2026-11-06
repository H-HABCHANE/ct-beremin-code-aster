"""Contrainte de Weibull et probabilite de rupture (modele de Beremin), a partir des champs aux points de Gauss."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Parametres:
    m: float
    sigma_u: float
    v0: float
    k: float | None = None  # None : pas de correction exp(-p/k)
    p_c: float = 0.0


# Tanguy et al. (2005), tableau 1 : m, sigma_u, V0, k, p_c.
TANGUY = Parametres(m=17.8, sigma_u=2925.0, v0=1.0e-3, k=4.0, p_c=1.0e-4)


def contrainte_principale_max(t_xx, t_yy, t_zz, t_xy):
    """Plus grande valeur propre d'un tenseur plan, en tenant compte de la composante hors plan."""
    centre = 0.5 * (t_xx + t_yy)
    rayon = np.sqrt((0.5 * (t_xx - t_yy)) ** 2 + t_xy**2)
    return np.maximum(centre + rayon, t_zz)


def contrainte_effective(sigma_1, p, params: Parametres):
    """Equation (10) de Tanguy et al. : nulle tant que p <= p_c ou que la plasticite n'evolue pas."""
    p = np.asarray(p, dtype=float)
    sigma_1 = np.asarray(sigma_1, dtype=float)
    dp = np.diff(p, axis=0, prepend=0.0)
    actif = (p > params.p_c) & (dp > 0.0)
    attenuation = 1.0 if params.k is None else np.exp(-p / params.k)
    return np.where(actif, np.maximum(sigma_1, 0.0) * attenuation, 0.0)


def weibull(sigma_eff, poids, coef_mult: float, params: Parametres):
    """Contrainte de Weibull par pas de chargement, avec le maximum historique de l'equation (11).

    sigma_eff a pour forme (pas, points de Gauss). coef_mult = epaisseur x facteur de symetrie.
    Le calcul passe par le logarithme d'un rapport pour eviter tout depassement a grand m.
    """
    maximum = np.maximum.accumulate(sigma_eff, axis=0)
    integrale_reduite = np.sum((maximum / params.sigma_u) ** params.m * poids, axis=1)
    rapport = coef_mult / params.v0 * integrale_reduite  # (sigma_w / sigma_u) ** m
    sigma_w = params.sigma_u * rapport ** (1.0 / params.m)
    probabilite = -np.expm1(-rapport)
    return sigma_w, probabilite


def calcule(champs, coef_mult: float, params: Parametres, facteur_xy: float = 1.0):
    """Applique le modele aux champs exportes par ct_beremin.comm."""
    sigma_1 = contrainte_principale_max(champs["V3"], champs["V4"], champs["V5"], champs["V6"] * facteur_xy)
    sigma_eff = contrainte_effective(sigma_1, champs["V1"], params)
    return weibull(sigma_eff, champs["w"], coef_mult, params)
