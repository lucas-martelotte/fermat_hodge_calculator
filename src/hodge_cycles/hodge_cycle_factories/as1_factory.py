from src.diagonal_variety import EvenDimensionalDiagonalVariety
from ..hodge_cycle_factory import AlgebraicPrimePrimitiveHodgeCycle
from ..complete_intersection_factory import CompleteIntersection
from ..hodge_cycle_factory import HodgeCycleFactory
from src.utils.sage_imports import prod, PolynomialQuotientRingElement
from ..complete_intersection_factory import CompleteIntersectionFactory
from itertools import product as iterprod
from src.utils.auxiliary import get_pairings
from math import gcd


class AokiShiodaType1Factory(CompleteIntersectionFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        base_permutation: list[int],
        rootvof2: PolynomialQuotientRingElement,
        zeta2d: PolynomialQuotientRingElement,
    ):
        """
        We fix u = gcd(m2, m3), v = gcd(m1, u/2).
        So rootvof2 should be 2^(1/v).

        This represents the cycle given by equations

        x1 - (rootvof2 * x2 * x3) = 0
        x0 - zeta4 * (x2**v + x3**v) = 0

        where
        x0 = xs[perm[0]]**(m0 // 2)
        x1 = xs[perm[1]]**(m1 // v)
        x2 = xs[perm[2]]**(m2 // (2 * v))
        x3 = xs[perm[3]]**(m3 // (2 * v))
        """
        assert variety.dimension == 2  # Only defined for surfaces
        self.base_permutation = base_permutation
        xs, ms = variety.polynomial_variables, variety.exps
        p = base_permutation
        m0, m1, m2, m3 = ms[p[0]], ms[p[1]], ms[p[2]], ms[p[3]]
        assert m0 % 2 == 0 and m2 % 2 == 0 and m3 % 2 == 0
        u = gcd(m2, m3)
        assert u % 2 == 0
        v = gcd(m1, u // 2)
        assert v > 1
        x0, x1, x2, x3 = (
            xs[p[0]] ** (m0 // 2),
            xs[p[1]] ** (m1 // v),
            xs[p[2]] ** (m2 // (2 * v)),
            xs[p[3]] ** (m3 // (2 * v)),
        )
        d = variety.degree
        zeta2v = zeta2d ** ((2 * d) // (2 * v))
        zetav, zeta4 = zeta2v**2, zeta2d ** (d // 2)
        f1 = x1 - (rootvof2 * x2 * x3)
        g1 = prod(
            x1 - zetav**t * (rootvof2 * x2 * x3) for t in range(1, v, 1)
        )
        f2 = x0 - zeta4 * (x2**v + x3**v)
        g2 = x0 + zeta4 * (x2**v + x3**v)
        super().__init__(variety, [f1, f2], [g1, g2])
