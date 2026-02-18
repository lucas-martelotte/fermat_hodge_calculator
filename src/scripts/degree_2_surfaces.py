from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_2_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 2
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(4)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    # zeta3 = formal_field.from_str("zeta3")
    # zeta6 = -(zeta3**2)
    # zeta6 = formal_field.from_str("zeta6")
    # zeta12 = (1 + zeta6) / root3
    # zeta4 = (2 * zeta6 - 1) / root3
    calculator = HodgeCalculator(
        X,
        {
            # Lines ("true lines")
            (1, 1, 1, 1): (-1 / K(4)),
        },
    )
    return X, calculator
