from ...diagonal_variety import EvenDimensionalDiagonalVariety
from src.utils.sage_imports import PolynomialQuotientRingElement
from ..complete_intersection_factory import CompleteIntersectionFactory


class OrdinaryTypeBBFactory(CompleteIntersectionFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        base_permutation: list[int],
        zeta12: PolynomialQuotientRingElement,
        root3of2: PolynomialQuotientRingElement,
    ):
        assert variety.dimension == 2  # Only defined for surfaces
        self.base_permutation = base_permutation
        xs, ms = variety.polynomial_variables, variety.exps
        p = base_permutation
        m0, m1, m2, m3 = ms[p[0]], ms[p[1]], ms[p[2]], ms[p[3]]
        assert m0 % 3 == 0 and m1 % 3 == 0 and m2 % 6 == 0 and m3 % 6 == 0
        x0, x1, x2, x3 = (
            xs[p[0]] ** (m0 // 3),
            xs[p[1]] ** (m1 // 3),
            xs[p[2]] ** (m2 // 6),
            xs[p[3]] ** (m3 // 6),
        )
        zeta6, zeta4, root3 = zeta12**2, zeta12**3, 2 * zeta12 - zeta12**3
        f1 = x0 + root3of2**2 * x2**2 - zeta6 * (x1 + root3of2**2 * x3**2)
        g1 = x0 + root3of2**2 * x2**2 - zeta6**3 * (x1 + root3of2**2 * x3**2)
        g1 *= x0 + root3of2**2 * x2**2 - zeta6**5 * (x1 + root3of2**2 * x3**2)
        f2 = x2**3 + root3of2 * x0 * x2 - zeta4 * (x3**3 + root3of2 * x1 * x3)
        g2 = (-3) * (
            x2**3
            + root3of2 * x0 * x2
            - zeta4**3 * (x3**3 + root3of2 * x1 * x3)
        )
        zeta_check = (3, 3, 6, 6)
        super().__init__(
            variety, [f1, f2], [g1, g2], zeta_exps_congruence_check=zeta_check
        )
