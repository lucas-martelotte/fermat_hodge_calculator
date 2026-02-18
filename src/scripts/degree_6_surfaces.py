from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator
from src.utils.auxiliary import lcm
from src.utils.sage_imports import PolynomialRing, QQ
from src.formal_algebra import (
    FormalNumberField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_6_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 6
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree

    K_formal = get_degree_6_surface_formal_field()
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)

    zeta12 = K_formal.from_str("zeta12")
    root3 = 2 * zeta12 - zeta12**3
    zeta6 = zeta12**2
    root3of2 = K_formal.from_str("root3of2")

    K_formal.locals["root3"] = root3
    K_formal.locals["zeta6"] = zeta6

    calculator = HodgeCalculator(
        X,
        {
            # Lines (on degree 2)
            (3, 3, 3, 3): -1 / K(4),
            # Lines (on degree 3)
            (2, 2, 4, 4): -1 / K(3),
            # Lines ("true lines")
            (1, 1, 5, 5): -1 / K(1),
            (1, 2, 4, 5): -1 / root3,
            (1, 3, 3, 5): -1 / K(2),
            (2, 3, 3, 4): -1 / (2 * root3),
            # Aoki-shioda
            (1, 3, 4, 4): -1 / (root3 * root3of2),
            (2, 2, 3, 5): -1 / (root3 * root3of2**2),
        },
    )
    return X, calculator


def get_degree_6_surface_formal_field() -> FormalNumberField:
    Kbase = PolynomialRing(QQ, ["zeta12base", "root3of2base"])
    zeta12base, root3of2base = Kbase.gens()
    Kbase_equations = [
        zeta12base**4 - zeta12base**2 + 1,
        root3of2base**3 - 2,
    ]
    return FormalNumberField(Kbase, Kbase_equations)
