from src.diagonal_variety import EvenDimensionalDiagonalVariety
from ..hodge_cycle_factory import AlgebraicPrimePrimitiveHodgeCycle
from ..complete_intersection_factory import CompleteIntersection
from ..hodge_cycle_factory import HodgeCycleFactory
from src.utils.sage_imports import prod, PolynomialQuotientRingElement
from itertools import product as iterprod
from src.utils.auxiliary import get_pairings
from math import gcd


class LineFactory(HodgeCycleFactory):
    def __init__(self, variety: EvenDimensionalDiagonalVariety):
        super().__init__(variety)

    def get_hodge_cycles(
        self,
    ) -> dict[str, AlgebraicPrimePrimitiveHodgeCycle]:
        d = self.variety.degree
        xs = self.variety.polynomial_variables
        K_formal = self.variety.formal_base_field
        zeta2d = K_formal.from_str(f"zeta{2*d}")
        ms = self.variety.exps
        line_ids = self._calculate_line_ids()
        id_to_line: dict[str, AlgebraicPrimePrimitiveHodgeCycle] = {}
        for line_id in line_ids:
            fs, gs = [], []
            for i, j, t in line_id:
                u = gcd(ms[i], ms[j])
                fs.append(
                    xs[i] ** (ms[i] // u)
                    - zeta2d ** ((d // u) * t) * xs[j] ** (ms[j] // u)
                )
                gs.append(
                    prod(
                        xs[i] ** (ms[i] // u)
                        - zeta2d ** ((d // u) * (t + 2 * r))
                        * xs[j] ** (ms[j] // u)
                        for r in range(1, u, 1)
                    )
                )
            id_to_line[str(line_id)] = CompleteIntersection(
                self.variety, fs, gs
            )
        return id_to_line

    def _calculate_line_ids(self) -> list[list[tuple[int, int, int]]]:
        """
        Returns all lines inside of the fermat variety of
        dimension n and degree d  in the convention where
        [(i, j, t), (u, v, s), ...] represents the line
        given by the equations
            {x_i^{m_i/u_{ij}} - zeta_{2u_{ij}}^t * x_j^{m_j/u_{ij}} = 0, ...}
        where u_{ij} = gcd(m_i, m_j); if u_{ij} = 1 for some (i,j),
        such a pair is not included as the line will have class = 0
        in the space of primitive hodge cycles.
        """
        n = self.variety.dimension
        ms = self.variety.exps
        u = [[gcd(mi, mj) for mj in ms] for mi in ms]
        pairings = get_pairings(n + 2)
        line_ids: list[list[tuple[int, int, int]]] = []
        for pairing in pairings:
            if any(u[i][j] == 1 for i, j in pairing):
                continue
            roots = iterprod(*(range(1, 2 * u[i][j], 2) for i, j in pairing))
            for root in roots:
                line_ids.append(
                    [(p[0], p[1], r) for p, r in zip(pairing, root)]
                )
        return line_ids
