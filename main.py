# from sage.libs.pari import pari
# pari.allocatemem(4_000_000_000)
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.scripts import (
    get_degree_2_surface,
    get_degree_3_surface,
    get_degree_4_surface,
    get_degree_5_surface,
    get_degree_6_surface,
)
from src.hodge_calculator import BasicHodgeCalculator
from time import perf_counter

start_time = perf_counter()


# X, calculator = get_degree_2_surface((2, 2, 2, 2))
# X, calculator = get_degree_3_surface((3, 3, 3, 3))
# X, calculator = get_degree_4_surface((4, 4, 4, 4))
# X, calculator = get_degree_5_surface((5, 5, 5, 5))
X, calculator = get_degree_6_surface((6, 6, 6, 6))

R_formal = X.formal_polynomial_ring
K_formal = X.formal_base_field
zetad = K_formal.K.gen()

periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
rank = calculator.get_rank_of_primitive_hodge_cycles()
print(periods)


# end_time = perf_counter()
# print(end_time - start_time)


# Still need to fix degree 7 and 8
