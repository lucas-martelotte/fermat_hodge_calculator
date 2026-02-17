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
    #calculator = BasicHodgeCalculator(X)
    #J_hodge_alg = calculator.get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma()
    #I = calculator.multi_indexes
    #idxs = set()
    #for j in J_hodge_alg:
    #    idxs.add(tuple(sorted([formi + 1 for formi in I[j]])))
    #print(idxs)

    zeta9 = K_formal.from_str("zeta9")
    
    calculator = HodgeCalculator(
        X,
        {
            (2, 2, 7, 7): K(0), 
            (3, 4, 5, 6): K(0), 
            (4, 4, 5, 5): K(0), 
            (1, 2, 7, 8): K(0), 
            (3, 3, 6, 6): K(0), 
            (1, 4, 5, 8): K(0), 
            (1, 4, 6, 7): K(0), 
            (2, 3, 5, 8): K(0), 
            (2, 3, 6, 7): K(0), 
            (2, 4, 5, 7): K(0), 
            (1, 1, 8, 8): K(0), 
            (1, 3, 6, 8): K(0)
        },
    )
    return X, calculator
