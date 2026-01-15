from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    FormalNumberField,
    PolynomialRingOverFormalNumberField,
)
from src.hodge_calculator import BasicHodgeCalculator

d = 7
K_formal = FormalCyclotomicField(2 * d)
R_formal = PolynomialRingOverFormalNumberField(K_formal, ["x0", "x1", "x2", "x3"])
X = EvenDimensionalDiagonalVariety(R_formal, [d] * 4)
calculator = BasicHodgeCalculator(X, override_files=["basis_of_primitive_hodge_cycles"])
calculator.get_data_from_json("basis_of_primitive_hodge_cycles", "all_multi_indexes")
