# from sage.libs.pari import pari
# pari.allocatemem(4_000_000_000)
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.hodge_calculator import BasicHodgeCalculator
from time import perf_counter

start_time = perf_counter()

d = 2
K_formal = FormalCyclotomicField(2 * d)


from src.utils.sage_imports import Matrix

zetad = K_formal.K.gen()
A = Matrix(K_formal.K, 5, 1, lambda i, j: zetad**i)


R_formal = PolynomialRingOverFormalNumberField(
    K_formal, ["x0", "x1", "x2", "x3"]
)
X = EvenDimensionalDiagonalVariety(R_formal, tuple([d] * 4))
calculator = BasicHodgeCalculator(
    X, override_files=["basis_of_primitive_hodge_cycles"]
)
rank = calculator.get_rank_of_primitive_hodge_cycles()
print(rank)


end_time = perf_counter()
print(end_time - start_time)


# Still need to fix degree 7 and 8
