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


start_time = perf_counter()
# ======================= #
d = 8

hodge_calculator_factory = HodgeCalculatorFactory()
X, calculator = hodge_calculator_factory.create((d, d, d, d))
factory = LineFactory(X)
# lines = factory.get_hodge_cycles()
# print(len(lines))
# for line_id, line in lines.items():
#    print(line_id)
calculator.compute_periods_from_hodge_cycle_factory(factory, "lines")
# line_data = calculator.get_hodge_cycle_factory_data("lines")
# print(line_data["coordinates"])


# ======================= #
end_time = perf_counter()
print(end_time - start_time)
