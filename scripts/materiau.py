"""Lois d'ecrouissage a -100 C. Chaque valeur renvoie a sa source (voir docs/donnees_materiau.md)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Loi:
    nom: str
    acier: str
    young: float
    nu: float
    r0: float
    termes: tuple[tuple[float, float], ...]  # couples (Q, b) de R = R0 + somme Q (1 - exp(-b p))

    def contrainte(self, p):
        p = np.asarray(p, dtype=float)
        return self.r0 + sum(q * (1.0 - np.exp(-b * p)) for q, b in self.termes)

    def courbe_traction(self, p_max: float = 3.0, n: int = 80) -> list[tuple[float, float]]:
        """Points (deformation totale, contrainte) pour DEFI_FONCTION, de la limite elastique a p_max."""
        p = np.concatenate(([0.0], np.geomspace(1.0e-4, p_max, n)))
        sigma = self.contrainte(p)
        return [(float(s / self.young + pi), float(s)) for pi, s in zip(p, sigma)]


# Tanguy et al. (2005), tableau 1 : E, nu, Q2, b2. R0, Q1, b1 mesures sur leur figure 1 a -100 C.
TANGUY = Loi("tanguy2005", "16MND5", 210000.0, 0.3, 589.0, ((184.0, 16.3), (472.0, 1.7)))

# Meme loi avec la sur-contrainte visqueuse de leur equation (2) a 1e-4 /s : K1 = 84, n1 = 7 (figure 1).
SURCONTRAINTE_1E4 = 84.0 * (1.0e-4) ** (1.0 / 7.0)
TANGUY_VISQUEUX = Loi("tanguy2005_visqueux", "16MND5", 210000.0, 0.3, 589.0 + SURCONTRAINTE_1E4,
                      ((184.0, 16.3), (472.0, 1.7)))

# Bordet et al. (2006), tableau 5, a -100 C. Autre acier (18MND5) : sert a mesurer la sensibilite.
# Les termes cinematiques (C/D, D) sont cumules : equivalent en chargement monotone seulement.
BORDET = Loi("bordet2006", "18MND5", 208900.0, 0.3, 539.6,
             ((324.3, 12.16), (296.4, 0.422), (73.5, 13.5), (2.2, 688.2)))

LOIS = {loi.nom: loi for loi in (TANGUY, TANGUY_VISQUEUX, BORDET)}
