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
    prod,
    gamma,
    sage_eval,
    exp,
    pi,
    I,
    sqrt,
)


def get_degree_12_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 12
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = get_degree_12_surface_formal_field()
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta24 = K_formal.from_str("zeta24")
    zeta12 = zeta24**2
    root12of2 = K_formal.from_str("root12of2")
    rr3m1 = K_formal.from_str("rr3m1")  # sqrt(sqrt(3)-1)
    root8of3 = K_formal.from_str("root8of3")
    root2, root3 = root12of2**6, 2 * zeta12 - zeta12**3
    root6, zeta6, zeta4 = root2 * root3, zeta12**2, zeta12**3
    zeta24 = (zeta12 + zeta12**2 - zeta12**3) / root2
    root4of3, root6of2, root4of2 = root8of3**2, root12of2**2, root12of2**3
    root3of2 = root12of2**4
    root3of4 = root12of2**8
    root4of6 = root4of2 * root4of3
    rr3p1 = (zeta24 + zeta24**3 - zeta24**7) * rr3m1  # sqrt(sqrt(3)+1)

    K_formal.locals["zeta12"] = zeta12

    calculator = HodgeCalculator(
        X,
        {
            (1, 1, 11, 11): -1 / (K(2) - root3),
            (1, 2, 10, 11): -root2 / (root3 - 1),
            (1, 3, 9, 11): -1 / (root3 - 1),
            (1, 4, 8, 11): -root2 / (root3 * (root3 - 1)),
            (1, 4, 9, 10): -root12of2 / (root8of3 * rr3m1),
            (1, 5, 7, 11): -K(1),
            (1, 5, 9, 9): -root4of3 / root2,
            (1, 6, 6, 11): -1 / (root2 * (root3 - 1)),
            (1, 6, 7, 10): -1 / root6of2,
            (1, 6, 8, 9): -(root4of2 * rr3p1) / (2 * root8of3),
            (1, 7, 8, 8): -root2 / root3,
            (2, 2, 10, 10): -K(1),
            (2, 3, 8, 11): -1 / (root12of2 * root8of3 * root4of3 * rr3m1),
            (2, 3, 9, 10): -1 / root2,
            (2, 4, 8, 10): -1 / root3,
            (2, 5, 6, 11): -1 / (root3of2**2 * root6of2),
            (2, 5, 7, 10): -root2 / (root3 + 1),
            (2, 5, 8, 9): -(root4of2 * root3of2**2) / (2 * root8of3 * rr3p1),
            (2, 6, 6, 10): -1 / K(2),
            (2, 6, 8, 8): -1 / (root3of2 * root3),
            (3, 3, 7, 11): -1 / (root2 * root4of3),
            (3, 3, 9, 9): -1 / K(2),
            (3, 4, 6, 11): -rr3p1 / (2 * root4of6 * root8of3),
            (3, 4, 7, 10): -(root12of2 * root8of3) / (root3 * rr3p1),
            (3, 4, 8, 9): -1 / root6,
            (3, 5, 7, 9): -1 / (root3 + 1),
            (3, 6, 6, 9): -1 / (2 * root2),
            (3, 6, 7, 8): -1 / (root4of6 * root8of3 * rr3p1),
            (4, 4, 5, 11): -1 / root6,
            (4, 4, 6, 10): -1 / (root3 * root3of4),
            (4, 4, 8, 8): -1 / K(3),
            (4, 5, 6, 9): -root4of2 / (2 * root8of3 * rr3p1),
            (4, 5, 7, 8): -root2 / (root3 * (1 + root3)),
            (4, 6, 6, 8): -1 / (2 * root3),
            (5, 5, 7, 7): -1 / (root3 + 2),
            (5, 6, 6, 7): -1 / (root2 * (root3 + 1)),
            (6, 6, 6, 6): -1 / K(4),
        },
    )

    """
    for key, value in calculator.gamma_values.items():
        lhs = (-1 / (4 * pi**2)) * prod([gamma(t / degree) for t in key])
        rhs = sage_eval(
            str(value),
            locals={
                "zeta24": exp(2 * pi * I / 24),
                "root12of2": 2 ** (QQ(1) / 12),
                "root8of3": 3 ** (QQ(1) / 8),
                "rr3m1": sqrt(sqrt(3) - QQ(1)),
            },
        )
        print(key, lhs.n() - rhs.n())
    """

    return X, calculator


def get_degree_12_surface_formal_field() -> FormalNumberField:
    Kbase = PolynomialRing(
        QQ, ["zeta24base", "root12of2base", "rr3m1base", "root8of3base"]
    )
    zeta24base, root12of2base, rr3m1base, root8of3base = Kbase.gens()
    zeta12base = zeta24base**2
    root3base = 2 * zeta12base - zeta12base**3
    root2base = zeta24base + zeta24base**3 - zeta24base**5
    Kbase_equations = [
        zeta24base**8 - zeta24base**4 + 1,
        root12of2base**6 - root2base,
        rr3m1base**2 - root3base + 1,
        root8of3base**4 - root3base,
    ]
    return FormalNumberField(Kbase, Kbase_equations)
