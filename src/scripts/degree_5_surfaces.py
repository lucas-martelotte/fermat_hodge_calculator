from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_5_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 5
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(5)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta5 = K_formal.from_str("zeta5")
    root5 = -(2 * zeta5**3 + 2 * zeta5**2 + 1)

    K_formal.locals["zeta5"] = zeta5
    K_formal.locals["root5"] = root5

    calculator = HodgeCalculator(
        X,
        {
            # Lines
            (1, 2, 3, 4): -1 / root5,
            (2, 2, 3, 3): -2 / (5 + root5),
            (1, 1, 4, 4): -2 / (5 - root5),
        },
    )
    return X, calculator
