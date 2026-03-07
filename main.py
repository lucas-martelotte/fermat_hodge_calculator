# from sage.libs.pari import pari
# pari.allocatemem(15_000_000_000)

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
from src.utils.sage_imports import QQ, Matrix
from time import perf_counter
import json
from src.utils.sage_imports import RR, zero_matrix
from src.utils.auxiliary import get_pairings, sage_matrix_map
from src.hodge_cycles.hodge_cycle_factories import (
    AokiShiodaType1Factory,
    AokiShiodaType3Factory,
    AokiShiodaType2BFactory,
    ExceptionalTypeCFactory,
    ExceptionalTypeCBFactory,
)


start_time = perf_counter()

# degree 8: 1655 seconds (standard methods)
# degree 8: rational methods (11 seconds)

# print("Creating mat")
# blocked_matrix = zero_matrix(RR, 600, 900 * 96)
# print("Created mat")

d = 12
hodge_calculator_factory = HodgeCalculatorFactory()
X, calculator = hodge_calculator_factory.create((12, 12, 12, 12))
K_formal = X.formal_base_field

# root3of2 = K_formal.from_str("root3of2")
# zeta12 = K_formal.from_str("zeta12")
# root3 = K_formal.from_str("root3")
root6of2 = K_formal.from_str("root6of2")
root4of3 = K_formal.from_str("root4of3")
zeta24 = K_formal.from_str("zeta24")
alpha = K_formal.from_str("alpha")
zeta12 = zeta24**2
root3of2 = root6of2**2
# factory = AokiShiodaType2BFactory(X, [0, 1, 2, 3], root4of2, zeta16)
factory = ExceptionalTypeCBFactory(X, [0, 1, 2, 3], zeta12, root3of2, alpha)
calculator.compute_periods_from_hodge_cycle_factory(factory, "cb")

end_time = perf_counter()
print(end_time - start_time)

exit()


# line_data = calculator.get_hodge_cycle_factory_data("lines")
# line_coords = line_data["coordinates"]
# ncols, nrows = line_coords.ncols(), line_coords.nrows()
# print({line_coords[i, j] for i in range(nrows) for j in range(ncols)})


exit()

period_matrix = calculator.get_period_matrix_of_primitive_hodge_cycles()
line_data = calculator.get_hodge_cycle_factory_data("lines_wrong")
line_periods = line_data["period_matrix"]

hodge_period_matrix_blocked = (
    K_formal.convert_to_rational_horizontally_blocked_matrix(period_matrix)
)
period_matrix_to_solve_blocked = (
    K_formal.convert_to_rational_horizontally_blocked_matrix(line_periods)
)
period_coords = hodge_period_matrix_blocked.solve_left(
    period_matrix_to_solve_blocked
).T
line_coords = Matrix(QQ, period_coords)  # Solution must lie inside (1/d)*ZZ
line_rank = int(line_coords.rank())

data = {
    "rank": int(line_rank),
    "number_of_cycles": line_data["number_of_cycles"],
    "cycle_ids": line_data["cycle_ids"],
    "period_matrix": sage_matrix_map(str, line_periods),
    "coordinates": sage_matrix_map(str, line_coords),
}

calculator.write_data_on_json("lines", data)

exit()

line_data2 = calculator.get_hodge_cycle_factory_data("lines2")


print(line_periods)

print()

print(line_data2["period_matrix"])
assert line_data == line_data2


"""
line_factory = LineFactory(X)
calculator.compute_periods_from_hodge_cycle_factory(line_factory, "lines2")

hodge_basis = calculator.get_basis_of_primitive_hodge_cycles()
hodge_periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
line_data = calculator.get_hodge_cycle_factory_data("lines2")
line_data_correct = calculator.get_hodge_cycle_factory_data("lines")
assert line_data == line_data_correct

line_coords = line_data["coordinates"]
line_periods = line_data["period_matrix"]

hodge_periods_blocked = K_formal.convert_to_rational_horizontally_blocked_matrix(hodge_periods)
line_periods_blocked = K_formal.convert_to_rational_horizontally_blocked_matrix(line_periods)
print("asserting...")
assert line_coords.T * hodge_periods == line_periods
assert line_coords.T * hodge_periods_blocked == line_periods_blocked
"""

exit()

"""
zeta24base**8 - zeta24base**4 + 1,
root6of2base**3 - root2base,
root4of3base**2 - root3base,
alphabase**2 - root2base * root4of3base * (root3base - 1),
"""

exit()

# period_matrix = calculator.get_period_matrix_of_primitive_hodge_cycles()
# for i in range(period_matrix.nrows()):
#    for j in range(period_matrix.ncols()):
#        input()
#        print(period_matrix[i, j])
# line_data = calculator.get_hodge_cycle_factory_data("lines")
# line_coords = line_data["coordinates"]
# ncols, nrows = line_coords.ncols(), line_coords.nrows()
# print({line_coords[i, j] for i in range(nrows) for j in range(ncols)})
# exit()

line_factory = LineFactory(X)
calculator.compute_periods_from_hodge_cycle_factory(line_factory, "lines2")

hodge_basis = calculator.get_basis_of_primitive_hodge_cycles()
hodge_periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
line_data = calculator.get_hodge_cycle_factory_data("lines2")
line_coords = line_data["coordinates"]
line_periods = line_data["period_matrix"]

hodge_periods_blocked = (
    K_formal.convert_to_rational_horizontally_blocked_matrix(hodge_periods)
)
line_periods_blocked = (
    K_formal.convert_to_rational_horizontally_blocked_matrix(line_periods)
)
print("asserting...")
assert line_coords.T * hodge_periods == line_periods
assert line_coords.T * hodge_periods_blocked == line_periods_blocked

exit()

line_data = calculator.get_hodge_cycle_factory_data("lines")
line_coords = line_data["coordinates"]

ncols, nrows = line_coords.ncols(), line_coords.nrows()
print(nrows, ncols)
print(line_data["rank"])


exit()
print({line_coords[i, j] for i in range(nrows) for j in range(ncols)})


exit()
line_factory = LineFactory(X)
calculator.compute_periods_from_hodge_cycle_factory(line_factory, "lines")

end_time = perf_counter()
print(end_time - start_time)

exit()
hodge_basis = calculator.get_basis_of_primitive_hodge_cycles()
hodge_periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
line_data = calculator.get_hodge_cycle_factory_data("lines")
line_coords = line_data["coordinates"]
line_periods = line_data["period_matrix"]

assert line_coords.T * hodge_periods == line_periods

ncols, nrows = line_coords.ncols(), line_coords.nrows()
print(hodge_basis.nrows(), hodge_basis.ncols())
print(nrows, ncols)
print({line_coords[i, j] for i in range(nrows) for j in range(ncols)})

new_coords = Matrix(
    QQ, nrows, ncols, lambda i, j: line_coords[i, j] - line_coords[i, 0]
)
print({new_coords[i, j] for i in range(nrows) for j in range(ncols)})


# periods = calculator.get_period_matrix_of_primitive_hodge_cycles()
# intersection = calculator.get_intersection_matrix_of_primitive_hodge_cycles()
# print("Intersection calculated!")
exit()

print("Importing period matrix.")
P = calculator.get_period_matrix_of_primitive_hodge_cycles()
print("Size of period matrix", P.nrows(), P.ncols())

Prat = K_formal.convert_to_rational_horizontally_blocked_matrix(P)
print("Size of blocked period matrix", Prat.nrows(), Prat.ncols())
V = Matrix(QQ, 1, Prat.ncols(), lambda i, j: 1)
X = Prat.solve_left(V).T
print("Solved!")


"""
print("Size of blocked period matrix", Prat.nrows(), Prat.ncols())
print("Computing right-inverse.")
Q = Prat.pseudoinverse()
print("Right-inverse computed. Size: ", Q.nrows(), Q.ncols())
print(Prat * Q == identity_matrix(RR, P.nrows()))
I = identity_matrix(RR, P.nrows())
print((Prat * Q - I).norm() < 0.00000001)

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
"""
