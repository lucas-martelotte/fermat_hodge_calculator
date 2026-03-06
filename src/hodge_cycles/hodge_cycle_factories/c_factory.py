from ...diagonal_variety import EvenDimensionalDiagonalVariety
from src.utils.sage_imports import PolynomialQuotientRingElement
from ..complete_intersection_factory import CompleteIntersectionFactory


class ExceptionalTypeCFactory(CompleteIntersectionFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        base_permutation: list[int],
        zeta12: PolynomialQuotientRingElement,
        root4of2_root8of3_rr3m1: PolynomialQuotientRingElement,
    ):
        """
        The variable root4of2_root8of3_rr3m1 should be equal to the
        number 2^(1/4) * 3^(1/8) * sqrt(sqrt(3)-1).
        """
        assert variety.dimension == 2  # Only defined for surfaces
        self.base_permutation = base_permutation
        xs, ms = variety.polynomial_variables, variety.exps
        p = base_permutation
        m0, m1, m2, m3 = ms[p[0]], ms[p[1]], ms[p[2]], ms[p[3]]
        assert m0 % 2 == 0 and m1 % 3 == 0 and m2 % 4 == 0 and m3 % 12 == 0
        x0, x1, x2, x3 = (
            xs[p[0]] ** (m0 // 2),
            xs[p[1]] ** (m1 // 3),
            xs[p[2]] ** (m2 // 4),
            xs[p[3]] ** (m3 // 12),
        )
        z3, z4, root3 = zeta12**4, zeta12**3, 2 * zeta12 - zeta12**3
        a, b = root4of2_root8of3_rr3m1, root3 - 1
        f1 = x0 - z4 * (
            x2**2 + (a**3 / 2) * x2 * x3**3 + (a**2 * root3 / 2) * x3**6
        )
        g1 = x0 + z4 * (
            x2**2 + (a**3 / 2) * x2 * x3**3 + (a**2 * root3 / 2) * x3**6
        )
        f2 = x1 - (b * x3**4 + a * x2 * x3)
        g2 = x1 - z3 * (b * x3**4 + a * x2 * x3)
        g2 *= x1 - z3**2 * (b * x3**4 + a * x2 * x3)
        super().__init__(variety, [f1, f2], [g1, g2])
