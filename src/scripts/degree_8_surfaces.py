from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    PolynomialRingOverFormalNumberField,
    FormalCyclotomicField,
    FormalNumberField,
)
from src.utils.sage_imports import (
    PolynomialRing,
    QQ,
    gamma,
    prod,
    pi,
    sage_eval,
    I,
    exp,
    sqrt,
)


def get_degree_8_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 8
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = get_degree_8_surface_formal_field()
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta16 = K_formal.from_str("zeta16")
    root4of2 = K_formal.from_str("root4of2")
    zeta8 = zeta16**2
    root2 = zeta8 - zeta8**3
    # sin_a8 means sin(a*pi/8)
    sin_18 = (1 - zeta8) * zeta16**3 / 2
    sin_28 = 1 / root2
    sin_38 = (zeta16 - zeta16**7) / 2

    K_formal.locals["zeta8"] = zeta8

    calculator = HodgeCalculator(
        X,
        {
            (4, 4, 4, 4): -1 / K(4),
            (2, 3, 4, 7): -1 / (root4of2 * 2),
            (1, 4, 5, 6): -root4of2 / K(2),
            (1, 4, 4, 7): -1 / (4 * sin_18),
            (2, 4, 4, 6): -1 / (4 * sin_28),
            (3, 4, 4, 5): -1 / (4 * sin_38),
            (1, 1, 7, 7): -1 / (4 * sin_18**2),
            (2, 2, 6, 6): -1 / (4 * sin_28**2),
            (3, 3, 5, 5): -1 / (4 * sin_38**2),
            (1, 2, 6, 7): -1 / (4 * sin_18 * sin_28),
            (1, 3, 5, 7): -1 / (4 * sin_18 * sin_38),
            (2, 3, 5, 6): -1 / (4 * sin_28 * sin_38),
        },
    )

    # checking the gamma values
    """
    for key, value in calculator.gamma_values.items():
        lhs = (-1 / (4 * pi**2)) * prod([gamma(t / degree) for t in key])
        rhs = sage_eval(
            str(value),
            locals={
                "zeta16": exp(2 * pi * I / 16),
                "root4of2": sqrt(sqrt(2)),
            },
        )
        print(lhs.n() - rhs.n())
    """

    return X, calculator


def get_degree_8_surface_formal_field() -> FormalNumberField:
    Kbase = PolynomialRing(QQ, ["zeta16base", "root4of2base"])
    zeta16base, root4of2base = Kbase.gens()
    Kbase_equations = [
        zeta16base**8 + 1,
        root4of2base**2 + zeta16base**6 - zeta16base**2,
    ]
    return FormalNumberField(Kbase, Kbase_equations)
