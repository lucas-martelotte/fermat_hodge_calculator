# from sage.libs.pari import pari
# pari.allocatemem(4_000_000_000)
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)

from src.scripts import HodgeCalculatorFactory
from src.hodge_calculator import BasicHodgeCalculator
from time import perf_counter
import json


start_time = perf_counter()


# R_formal = X.formal_polynomial_ring
# K_formal = X.formal_base_field
# zetad = K_formal.K.gen()


# periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
# rank = calculator.get_rank_of_primitive_hodge_cycles()
# print(periods)

hodge_calculator_factory = HodgeCalculatorFactory()


print("Starting degree 10")
X, calculator = hodge_calculator_factory.create((10, 10, 10, 10))
M = calculator.get_intersection_matrix_of_primitive_hodge_cycles()
print("Computed for degree 10")

print("Starting degree 11")
X, calculator = hodge_calculator_factory.create((11, 11, 11, 11))
M = calculator.get_intersection_matrix_of_primitive_hodge_cycles()
print("Computed for degree 11")

print("Starting degree 12")
X, calculator = hodge_calculator_factory.create((12, 12, 12, 12))
M = calculator.get_intersection_matrix_of_primitive_hodge_cycles()
print("Computed for degree 12")


# print(calculator.get_rank_of_primitive_hodge_cycles())


end_time = perf_counter()
print(end_time - start_time)

# Still need to fix degree 7 and 8
