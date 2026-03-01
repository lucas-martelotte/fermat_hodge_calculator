# from sage.libs.pari import pari
# pari.allocatemem(4_000_000_000)
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)

from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.hodge_cycles.hodge_cycle_factories import LineFactory
from src.hodge_calculator import HodgeCalculatorFactory
from src.hodge_calculator import BasicHodgeCalculator
from time import perf_counter
import json

"""
d = 3
hodge_calculator_factory = HodgeCalculatorFactory()
X, calculator = hodge_calculator_factory.create((d, d, d, d))
hodge_period_matrix = calculator.get_period_matrix_of_primitive_hodge_cycles()
line_data = calculator.get_hodge_cycle_factory_data("lines")
line_coords = line_data["coordinates"]
line_periods = line_data["period_matrix"]

print(line_periods)
print()
print(hodge_period_matrix.nrows(), hodge_period_matrix.ncols())
print(line_coords.nrows(), line_coords.ncols())
print((hodge_period_matrix.T * line_coords).T)

print()
print("======================================")
print()
"""




for d in [5]:

    start_time = perf_counter()
    # ======================= #

    print(f"Computing for degree {d}")

    hodge_calculator_factory = HodgeCalculatorFactory()
    X, calculator = hodge_calculator_factory.create((d, d, d, d))
    factory = LineFactory(X)
    # lines = factory.get_hodge_cycles()
    # print(len(lines))
    # for line_id, line in lines.items():
    #    print(line_id)
    calculator.compute_periods_from_hodge_cycle_factory(factory, "lines2")
    # line_data = calculator.get_hodge_cycle_factory_data("lines")
    # print(line_data["coordinates"])

    print(f"Finished degree {d}")
    
    # ======================= #
    end_time = perf_counter()
    print(end_time - start_time)
