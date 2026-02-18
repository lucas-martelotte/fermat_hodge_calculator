from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_9_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 9
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(9)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)
    # calculator = BasicHodgeCalculator(X)
    # J_hodge_alg = calculator.get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma()
    # I = calculator.multi_indexes
    # idxs = set()
    # for j in J_hodge_alg:
    #    idxs.add(tuple(sorted([formi + 1 for formi in I[j]])))
    # print(idxs)

    zeta9 = K_formal.from_str("zeta9")
    root3of3 = K(0)  # TODO

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
            # Lines
            (1, 1, 8, 8): -1 / (4 * sin_1919),
            (2, 2, 7, 7): -1 / (4 * sin_2929),
            (3, 3, 6, 6): -1 / (4 * sin_3939),
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
    return X, calculator
