from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_4_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 4
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(8)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta8 = K_formal.from_str("zeta8")
    zeta4 = zeta8**2
    root2 = zeta8 - zeta8**3
    K_formal.locals["zeta4"] = zeta4
    K_formal.locals["root2"] = root2

    calculator = HodgeCalculator(
        X,
        {
            (1, 1, 3, 3): -1 / K(2),
            (1, 2, 2, 3): -1 / (2 * root2),
            (2, 2, 2, 2): -1 / K(4),
        },
    )
    return X, calculator
