from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_7_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 7
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(7)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta7 = K_formal.from_str("zeta7")
    # sin_a7b7 means sin(a*pi/7) * sin(b*pi/7)
    sin_1717 = (2 - zeta7**6 - zeta7) / 4
    sin_2727 = (2 - zeta7**5 - zeta7**2) / 4
    sin_3737 = (2 - zeta7**4 - zeta7**3) / 4
    sin_1727 = (1 + zeta7**3 - zeta7**2 - zeta7) * zeta7**2 / 4
    sin_1737 = (1 + zeta7**5 - zeta7**4 - zeta7) * zeta7 / 4
    sin_2737 = (1 + zeta7**5 - zeta7**3 - zeta7**2) * zeta7 / 4

    calculator = HodgeCalculator(
        X,
        {
            # Lines
            (1, 1, 6, 6): -1 / (4 * sin_1717),
            (1, 2, 5, 6): -1 / (4 * sin_1727),
            (1, 3, 4, 6): -1 / (4 * sin_1737),
            (2, 2, 5, 5): -1 / (4 * sin_2727),
            (2, 3, 4, 5): -1 / (4 * sin_2737),
            (3, 3, 4, 4): -1 / (4 * sin_3737),
        },
    )
    return X, calculator
