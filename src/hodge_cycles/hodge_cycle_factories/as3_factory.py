from src.diagonal_variety import EvenDimensionalDiagonalVariety
from ..hodge_cycle_factory import AlgebraicPrimePrimitiveHodgeCycle
from ..complete_intersection_factory import CompleteIntersection
from ..hodge_cycle_factory import HodgeCycleFactory
from src.utils.sage_imports import prod, PolynomialQuotientRingElement
from ..complete_intersection_factory import CompleteIntersectionFactory
from itertools import product as iterprod
from src.utils.auxiliary import get_pairings
from math import gcd


class AokiShiodaType3Factory(CompleteIntersectionFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        base_permutation: list[int],
        rootvof3: PolynomialQuotientRingElement,
        zeta2d: PolynomialQuotientRingElement,
    ):
        """
        We fix v = gcd(m0/3, m1/3, m2/3, m3) > 1.
        This represents the cycle given by equations

        x0^v + x1^v + x3^v = 0
        3^{1/v} * x0 * x1 * x2 - zeta_{2v} x3

        where
        x0 = xs[perm[0]]**(m0 // 3v)
        x1 = xs[perm[1]]**(m1 // 3v)
        x2 = xs[perm[2]]**(m2 // 3v)
        x3 = xs[perm[3]]**(m3 // v)
        """
        assert variety.dimension == 2  # Only defined for surfaces
        self.base_permutation = base_permutation
        xs, ms = variety.polynomial_variables, variety.exps
        p = base_permutation
        m0, m1, m2, m3 = ms[p[0]], ms[p[1]], ms[p[2]], ms[p[3]]
        assert m0 % 3 == 0 and m1 % 3 == 0 and m2 % 3 == 0
        v = gcd(m0//3, m1//3, m2//3, m3)
        assert v > 1
        x0, x1, x2, x3 = (
            xs[p[0]] ** (m0 // (3*v)),
            xs[p[1]] ** (m1 // (3*v)),
            xs[p[2]] ** (m2 // (3*v)),
            xs[p[3]] ** (m3 // v),
        )
        d = variety.degree
        zeta2v = zeta2d ** ((2 * d) // (2 * v))
        zeta3 = zeta2d**((2 * d) // 3)
        f1 = x0**v + x1**v + x2**v
        g1 = (x0**v + zeta3*x1**v + zeta3**2*x2**v)
        g1 *= (x0**v + zeta3**2*x1**v + zeta3*x2**v)
        f2 = rootvof3 * x0 * x1 * x2 - zeta2v * x3
        g2 = prod(rootvof3 * x0 * x1 * x2 - zeta2v**t * x3 for t in range(3, 2*v, 2))
        super().__init__(variety, [f1, f2], [g1, g2])
