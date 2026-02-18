from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_11_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 11
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(degree)
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

    zeta11 = K_formal.from_str("zeta11")

    sin_1_11_1_11 = (2 - zeta11**1 - zeta11**10) / 4
    sin_2_11_2_11 = (2 - zeta11**2 - zeta11**9) / 4
    sin_3_11_3_11 = (2 - zeta11**3 - zeta11**8) / 4
    sin_4_11_4_11 = (2 - zeta11**4 - zeta11**7) / 4
    sin_5_11_5_11 = (2 - zeta11**5 - zeta11**6) / 4

    # root of 1024 x^5 - 704 x^3 + 88 x - 11 near x = 0.152316
    sin_1_11_2_11 = (1 + zeta11**3 - zeta11**2 - zeta11) * zeta11**4 / 4
    # root of 1024 x^5 - 704 x^3 + 176 x^2 + 44 x - 11 near x = 0.212919
    sin_1_11_3_11 = (1 + zeta11**9 - zeta11**8 - zeta11) * zeta11 / 4
    # root of 1024 x^5 - 704 x^3 - 176 x^2 + 44 x + 11 near x = 0.256273
    sin_1_11_4_11 = (1 + zeta11**5 - zeta11**4 - zeta11) / 4
    # root of 1024 x^5 - 704 x^3 + 88 x - 11 near x = 0.278865
    sin_1_11_5_11 = (1 + zeta11**7 - zeta11**6 - zeta11) * zeta11**2 / 4
    # root of 1024 x^5 - 704 x^3 + 176 x^2 + 44 x - 11 near x = 0.408589
    sin_2_11_3_11 = (1 + zeta11**5 - zeta11**3 - zeta11**2) * zeta11**3 / 4
    # root of 1024 x^5 - 704 x^3 + 88 x + 11 near x = 0.491784
    sin_2_11_4_11 = (1 + zeta11**9 - zeta11**7 - zeta11**2) * zeta11 / 4
    # root of 1024 x^5 - 704 x^3 + 176 x^2 + 44 x - 11 near x = 0.535138
    sin_2_11_5_11 = (1 + zeta11**7 - zeta11**5 - zeta11**2) * zeta11**2 / 4
    # root of 1024 x^5 - 704 x^3 + 88 x + 11 near x = 0.687454
    sin_3_11_4_11 = (1 + zeta11**7 - zeta11**4 - zeta11**3) * zeta11**2 / 4
    # root of 1024 x^5 - 704 x^3 + 88 x - 11 near x = 0.748057
    sin_3_11_5_11 = (1 + zeta11**9 - zeta11**6 - zeta11**3) * zeta11 / 4
    # root of 1024 x^5 - 704 x^3 - 176 x^2 + 44 x + 11 near x = 0.900373
    sin_4_11_5_11 = (1 + zeta11**9 - zeta11 ^ 5 - zeta11 ^ 4) * zeta11 / 4

    calculator = HodgeCalculator(
        X,
        {
            # Lines ("true lines")
            (1, 1, 10, 10): -1 / (4 * sin_1_11_1_11),
            (2, 2, 9, 9): -1 / (4 * sin_2_11_2_11),
            (3, 3, 8, 8): -1 / (4 * sin_3_11_3_11),
            (4, 4, 7, 7): -1 / (4 * sin_4_11_4_11),
            (5, 5, 6, 6): -1 / (4 * sin_5_11_5_11),
            (1, 2, 9, 10): -1 / (4 * sin_1_11_2_11),
            (1, 3, 8, 10): -1 / (4 * sin_1_11_3_11),
            (1, 4, 7, 10): -1 / (4 * sin_1_11_4_11),
            (1, 5, 6, 10): -1 / (4 * sin_1_11_5_11),
            (2, 3, 8, 9): -1 / (4 * sin_2_11_3_11),
            (2, 4, 7, 9): -1 / (4 * sin_2_11_4_11),
            (2, 5, 6, 9): -1 / (4 * sin_2_11_5_11),
            (3, 4, 7, 8): -1 / (4 * sin_3_11_4_11),
            (3, 5, 6, 8): -1 / (4 * sin_3_11_5_11),
            (4, 5, 6, 7): -1 / (4 * sin_4_11_5_11),
        },
    )
    return X, calculator
