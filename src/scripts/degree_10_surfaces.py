from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.utils.auxiliary import lcm
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)


def get_degree_10_surface(
    exps: tuple[int, ...],
) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
    degree = 10
    assert len(exps) == 4
    assert all(degree % e == 0 for e in exps)
    assert lcm(*exps) == degree
    K_formal = FormalCyclotomicField(9)
    K = K_formal.K
    R_formal = PolynomialRingOverFormalNumberField(
        K_formal, ["x0", "x1", "x2", "x3"]
    )
    X = EvenDimensionalDiagonalVariety(R_formal, exps)
    #calculator = BasicHodgeCalculator(X)
    #calculator.get_basis_of_primitive_hodge_cycles()

    #exit()

    #J_hodge_alg = calculator.get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma()
    #I = calculator.multi_indexes
    #idxs = set()
    #for j in J_hodge_alg:
    #    idxs.add(tuple(sorted([formi + 1 for formi in I[j]])))
    #for idx in idxs:
    #    print(str(idx) + ": K(0),")
    #exit()

    zeta10 = K_formal.from_str("zeta10")
    
    calculator = HodgeCalculator(
        X,
        {
            (2, 3, 7, 8): K(0),
            (3, 3, 7, 7): K(0),
            (2, 2, 7, 9): K(0),
            (1, 2, 8, 9): K(0),
            (1, 1, 9, 9): K(0),
            (2, 5, 6, 7): K(0),
            (4, 4, 6, 6): K(0),
            (1, 5, 5, 9): K(0),
            (1, 5, 6, 8): K(0),
            (2, 4, 6, 8): K(0),
            (1, 4, 6, 9): K(0),
            (5, 5, 5, 5): K(0),
            (2, 4, 5, 9): K(0),
            (3, 4, 6, 7): K(0),
            (3, 4, 5, 8): K(0),
            (2, 5, 5, 8): K(0),
            (1, 6, 6, 7): K(0),
            (3, 5, 5, 7): K(0),
            (3, 4, 4, 9): K(0),
            (1, 3, 8, 8): K(0),
            (4, 5, 5, 6): K(0),
            (2, 2, 8, 8): K(0),
            (1, 3, 7, 9): K(0),
        }
    )
    return X, calculator
