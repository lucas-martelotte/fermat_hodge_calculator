from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalNumberField,
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.utils.sage_imports import (
    PolynomialRing,
    QQ,
    sage_eval,
    sqrt,
    prod,
    gamma,
    pi,
    I,
    exp,
)


def get_degree_9_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 9
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = get_degree_9_surface_formal_field()
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta9 = K(K_formal.from_str("zeta9"))
    root3of3 = K(K_formal.from_str("root3of3"))

    sin_1919 = (2 + zeta9**5 + zeta9**2 - zeta9) / 4
    sin_1929 = (1 + zeta9**3 - zeta9 - zeta9**2) * zeta9**3 / 4
    sin_1939 = (1 + zeta9**7 - zeta9 - zeta9**6) * zeta9 / 4
    sin_1949 = (1 + zeta9**2 - zeta9 - zeta9**4) / 4
    sin_2929 = (2 + zeta9 - zeta9**2 + zeta9**4) / 4
    sin_2939 = (1 + zeta9**5 - zeta9**2 - zeta9**3) * zeta9**2 / 4
    sin_2949 = (1 + zeta9 - zeta9**5 - zeta9**2) / 4
    sin_3939 = K(3) / 4
    sin_3949 = (1 + zeta9**7 - zeta9**4 - zeta9**3) * zeta9 / 4
    sin_4949 = (2 - zeta9**4 - zeta9**5) / 4

    calculator = HodgeCalculator(
        X,
        {
            # Lines (on degree 3)
            (3, 3, 6, 6): -1 / (4 * sin_3939),
            # Lines ("true lines")
            (1, 1, 8, 8): -1 / (4 * sin_1919),
            (2, 2, 7, 7): -1 / (4 * sin_2929),
            (4, 4, 5, 5): -1 / (4 * sin_4949),
            (1, 2, 7, 8): -1 / (4 * sin_1929),
            (1, 3, 6, 8): -1 / (4 * sin_1939),
            (1, 4, 5, 8): -1 / (4 * sin_1949),
            (2, 3, 6, 7): -1 / (4 * sin_2939),
            (2, 4, 5, 7): -1 / (4 * sin_2949),
            (3, 4, 5, 6): -1 / (4 * sin_3949),
            # Aoki-shioda
            (1, 4, 6, 7): -1 / root3of3,
            (2, 3, 5, 8): -1 / root3of3**2,
        },
    )

    # checking the gamma values

    """
    for key, value in calculator.gamma_values.items():
        lhs = (-1 / (4 * pi**2)) * prod([gamma(t / degree) for t in key])
        rhs = sage_eval(
            str(value),
            locals={
                "zeta9": exp(2 * pi * I / 9),
                "root3of3": 3 ** (QQ(1) / 3),
            },
        )
        print(lhs.n() - rhs.n())
    """

    return X, calculator


def get_degree_9_surface_formal_field() -> FormalNumberField:
    Kbase = PolynomialRing(QQ, ["zeta9base", "root3of3base"])
    zeta9base, root3of3base = Kbase.gens()
    Kbase_equations = [
        zeta9base**6 + zeta9base**3 + 1,
        root3of3base**3 - 3,
    ]
    return FormalNumberField(Kbase, Kbase_equations)
