from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_12_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 12
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(degree)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)
    #calculator = BasicHodgeCalculator(X)
    #calculator.get_basis_of_primitive_hodge_cycles()

    #J_hodge_alg = calculator.get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma()
    #I = calculator.multi_indexes
    #idxs = set()
    #for j in J_hodge_alg:
    #    idxs.add(tuple(sorted([formi + 1 for formi in I[j]])))
    #for idx in idxs:
    #    print(str(idx) + ": K(0),")
    #exit()

    zeta12 = K_formal.from_str("zeta12")
    
    calculator = HodgeCalculator(
        X,
        {
            (5, 5, 7, 7): K(0),
            (2, 5, 6, 11): K(0),
            (3, 6, 7, 8): K(0),
            (2, 5, 7, 10): K(0),
            (4, 4, 5, 11): K(0),
            (4, 4, 6, 10): K(0),
            (1, 6, 8, 9): K(0),
            (3, 5, 7, 9): K(0),
            (3, 4, 7, 10): K(0),
            (1, 6, 7, 10): K(0),
            (2, 6, 8, 8): K(0),
            (2, 2, 10, 10): K(0),
            (4, 5, 6, 9): K(0),
            (4, 5, 7, 8): K(0),
            (1, 3, 9, 11): K(0),
            (3, 4, 6, 11): K(0),
            (3, 6, 6, 9): K(0),
            (1, 6, 6, 11): K(0),
            (3, 4, 8, 9): K(0),
            (6, 6, 6, 6): K(0),
            (2, 3, 8, 11): K(0),
            (2, 3, 9, 10): K(0),
            (3, 3, 9, 9): K(0),
            (1, 7, 8, 8): K(0),
            (1, 5, 9, 9): K(0),
            (1, 1, 11, 11): K(0),
            (2, 6, 6, 10): K(0),
            (1, 4, 8, 11): K(0),
            (1, 4, 9, 10): K(0),
            (1, 2, 10, 11): K(0),
            (3, 3, 7, 11): K(0),
            (4, 6, 6, 8): K(0),
            (2, 5, 8, 9): K(0),
            (5, 6, 6, 7): K(0),
            (4, 4, 8, 8): K(0),
            (2, 4, 8, 10): K(0),
            (1, 5, 7, 11): K(0),
        }
    )
    return X, calculator
