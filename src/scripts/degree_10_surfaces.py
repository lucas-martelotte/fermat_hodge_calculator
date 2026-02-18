from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalNumberField,
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.utils.sage_imports import PolynomialRing, QQ


def get_degree_10_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 10
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(20)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)
    # calculator = BasicHodgeCalculator(X)
    # calculator.get_basis_of_primitive_hodge_cycles()

    # exit()

    # J_hodge_alg = calculator.get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma()
    # I = calculator.multi_indexes
    # idxs = set()
    # for j in J_hodge_alg:
    #    idxs.add(tuple(sorted([formi + 1 for formi in I[j]])))
    # for idx in idxs:
    #    print(str(idx) + ": K(0),")
    # exit()

    zeta20 = K_formal.from_str("zeta20")
    zeta10 = zeta20**2
    root5of2 = K_formal.from_str("root5of2")
    root5 = 1 + 2 * zeta10**2 - 2 * zeta10**3

    sin_1_10 = (root5 - 1) / 4
    sin_2_10 = (1 - zeta20**4) * zeta20**3 / 2
    sin_3_10 = (root5 + 1) / 4
    sin_4_10 = (2 * zeta20 - zeta20**3 + zeta20**5 - zeta20**7) / 2
    sin_5_10 = K(1)

    K_formal.locals["zeta10"] = zeta10

    calculator = HodgeCalculator(
        X,
        {
            # Lines
            (1, 1, 9, 9): -1 / (4 * sin_1_10**2),
            (2, 2, 8, 8): -1 / (4 * sin_2_10**2),
            (3, 3, 7, 7): -1 / (4 * sin_3_10**2),
            (4, 4, 6, 6): -1 / (4 * sin_4_10**2),
            (5, 5, 5, 5): -1 / (4 * sin_5_10**2),
            (1, 2, 8, 9): -1 / (4 * sin_1_10 * sin_2_10),
            (1, 3, 7, 9): -1 / (4 * sin_1_10 * sin_3_10),
            (1, 4, 6, 9): -1 / (4 * sin_1_10 * sin_4_10),
            (1, 5, 5, 9): -1 / (4 * sin_1_10 * sin_5_10),
            (2, 3, 7, 8): -1 / (4 * sin_2_10 * sin_3_10),
            (2, 4, 6, 8): -1 / (4 * sin_2_10 * sin_4_10),
            (2, 5, 5, 8): -1 / (4 * sin_2_10 * sin_5_10),
            (3, 4, 6, 7): -1 / (4 * sin_3_10 * sin_4_10),
            (3, 5, 5, 7): -1 / (4 * sin_3_10 * sin_5_10),
            (4, 5, 5, 6): -1 / (4 * sin_4_10 * sin_5_10),
            # Aoki-shiodaS
            (2, 2, 7, 9): -1 / (2 * root5of2 * sin_2_10),
            (1, 5, 6, 8): -1 / (2 * root5of2 * sin_2_10),
            (2, 4, 5, 9): -1 / (2 * root5of2**4 * sin_2_10),
            (2, 5, 6, 7): -1 / (2 * root5of2**2 * sin_4_10),
            (3, 4, 4, 9): -1 / (2 * root5of2**2 * sin_4_10),
            (3, 4, 5, 8): -1 / (2 * root5of2**3 * sin_4_10),
            (1, 6, 6, 7): -1 / (root5of2**3 * sin_4_10),
            (1, 3, 8, 8): -1 / (root5of2**4 * sin_2_10),
        },
    )
    return X, calculator


def get_degree_10_surface_formal_field() -> FormalNumberField:
    Kbase = PolynomialRing(QQ, ["zeta20base", "root5of2base"])
    zeta20base, root5of2base = Kbase.gens()
    Kbase_equations = [
        zeta20base**8 - zeta20base**6 + zeta20base**4 - zeta20base**2 + 1,
        root5of2base**5 - 2,
    ]
    return FormalNumberField(Kbase, Kbase_equations)
