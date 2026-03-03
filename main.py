from sage.libs.pari import pari

pari.allocatemem(4_000_000_000)

from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.utils.sage_imports import identity_matrix
from src.utils.auxiliary import sage_matrix_map
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.hodge_cycles.hodge_cycle_factories import LineFactory
from src.hodge_calculator import HodgeCalculatorFactory
from src.hodge_calculator import BasicHodgeCalculator
from src.utils.sage_imports import QQ
from time import perf_counter
import json


start_time = perf_counter()

# degree 8: 1655 seconds (standard methods)
# degree 8: rational methods (11 seconds)
d = 12
hodge_calculator_factory = HodgeCalculatorFactory()
X, calculator = hodge_calculator_factory.create((d, d, d, d))
#periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
intersection = calculator.get_intersection_matrix_of_primitive_hodge_cycles()
print("Intersection calculated!")

exit()

K_formal = X.formal_base_field
print("Importing period matrix.")
P = calculator.get_period_matrix_of_primitive_hodge_cycles()
print("Size of period matrix", P.nrows(), P.ncols())

Prat = K_formal.convert_to_rational_horizontally_blocked_matrix(P)
print("Size of blocked period matrix", Prat.nrows(), Prat.ncols())
print("Computing right-inverse.")
Q = Prat.pseudoinverse()
print("Right-inverse computed. Size: ", Q.nrows(), Q.ncols())
print(Prat * Q == identity_matrix(QQ, P.nrows()))


end_time = perf_counter()
print(end_time - start_time)

exit()

print("Computing right-inverse.")
Q = P.pseudoinverse()
print("Right-inverse computed. Size: ", Q.nrows(), Q.ncols())
print("Exporting...")
data = {"right_inverse": sage_matrix_map(str, Q)}
with open(f"./right_inverse_deg_{d}.json", "w") as f:
    json.dump(data, f)
print("Exported successfully!")


exit()

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
