from src.diagonal_variety import EvenDimensionalDiagonalVariety
from ..hodge_cycle_factory import AlgebraicPrimePrimitiveHodgeCycle
from ..complete_intersection_factory import CompleteIntersection
from ..hodge_cycle_factory import HodgeCycleFactory
from src.utils.sage_imports import prod, PolynomialQuotientRingElement
from ..complete_intersection_factory import CompleteIntersectionFactory
from itertools import product as iterprod
from src.utils.auxiliary import get_pairings
from math import gcd


class AokiShiodaType2BFactory(CompleteIntersectionFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        base_permutation: list[int],
        root2vof2: PolynomialQuotientRingElement,
        zeta2d: PolynomialQuotientRingElement,
    ):
        """
        We fix v = gcd(m0, m1/2, m2/4, m3/4) > 1.
        This represents the cycle given by equations

        x0 - root2vof2^3 * zeta4v * x1 * x2 * x3
        x2^v + x3^v - zeta4 * (x1^v + zeta4 * root2 * x2^v * x3^v)

        where
        x0 = xs[perm[0]]**(m0 // v)
        x1 = xs[perm[1]]**(m1 // 2v)
        x2 = xs[perm[2]]**(m2 // 4v)
        x3 = xs[perm[3]]**(m3 // 4v)
        """
        assert variety.dimension == 2  # Only defined for surfaces
        self.base_permutation = base_permutation
        xs, ms = variety.polynomial_variables, variety.exps
        p = base_permutation
        m0, m1, m2, m3 = ms[p[0]], ms[p[1]], ms[p[2]], ms[p[3]]
        assert m1 % 2 == 0 and m2 % 4 == 0 and m3 % 4 == 0
        v = gcd(m0, m1//2, m2//4, m3//4)
        assert v > 1
        x0, x1, x2, x3 = (
            xs[p[0]] ** (m0 // v),
            xs[p[1]] ** (m1 // (2*v)),
            xs[p[2]] ** (m2 // (4*v)),
            xs[p[3]] ** (m3 // (4*v)),
        )
        d = variety.degree
        zeta4v = zeta2d**((2 * d) // (4 * v))
        zeta4 = zeta2d**((2 * d) // 4)
        zetav = zeta2d**((2 * d) // v)
        root2 = root2vof2**v
        f1 = x0 - root2vof2**3 * zeta4v * x1 * x2 * x3
        g1 = prod(x0 - root2vof2**3 * zetav**t * zeta4v * x1 * x2 * x3 for t in range(1, v, 1))
        f2 = x2**(2*v) + x3**(2*v) - zeta4 * (x1**v + zeta4 * root2 * x2**v * x3**v)
        g2 = x2**(2*v) + x3**(2*v) + zeta4 * (x1**v + zeta4 * root2 * x2**v * x3**v)
        super().__init__(variety, [f1, f2], [g1, g2])
