"""Solutions de reference pour l'eprouvette CT (etape 1)."""
from __future__ import annotations

import math


def f_e399(alpha: float) -> float:
    """Fonction de forme ASTM E399, forme etendue de Srawley (1976), 0.2 <= a/W < 1 ici."""
    if not 0.2 <= alpha < 1.0:
        raise ValueError(f"a/W hors du domaine de la formule : {alpha}")
    poly = 0.886 + 4.64 * alpha - 13.32 * alpha**2 + 14.72 * alpha**3 - 5.6 * alpha**4
    return (2.0 + alpha) * poly / (1.0 - alpha) ** 1.5


def k_e399(p_over_b: float, w: float, alpha: float) -> float:
    """K_I en MPa.sqrt(mm) pour P/B en N/mm et W en mm."""
    return p_over_b / math.sqrt(w) * f_e399(alpha)


def g_from_k_plane_strain(k: float, young: float, nu: float) -> float:
    """Irwin en deformations planes : G = K^2 (1 - nu^2) / E."""
    return k**2 * (1.0 - nu**2) / young


def k_from_g_plane_strain(g: float, young: float, nu: float) -> float:
    return math.sqrt(g * young / (1.0 - nu**2))


def mpa_sqrt_m(k_mpa_sqrt_mm: float) -> float:
    return k_mpa_sqrt_mm / math.sqrt(1000.0)
