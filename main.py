from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    FormalNumberField,
    PolynomialRingOverFormalNumberField,
)

K_formal = FormalCyclotomicField(8)
R_formal = PolynomialRingOverFormalNumberField(K_formal, ["x0", "x1", "x2", "x3"])
X = EvenDimensionalDiagonalVariety(R_formal, [4, 4, 4, 4])
