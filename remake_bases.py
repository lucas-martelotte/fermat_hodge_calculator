from src.hodge_calculator import HodgeCalculatorFactory
from src.utils.sage_imports import Matrix, ZZ
from src.utils.auxiliary import sage_matrix_map
import json

degree = 3
calculator_factory = HodgeCalculatorFactory()
X, calculator = calculator_factory.create((degree, degree, degree, degree))
basis = calculator.get_basis_of_primitive_hodge_cycles()
D, U, V = basis.smith_form()
# print(D - U * basis * V)  # this should be zero
# What should be done in theory...
new_D = Matrix(
    ZZ, D.nrows(), D.ncols(), lambda i, j: 1 if D[i, j] != 0 else 0
)
new_basis = U.inverse() * D * V.inverse()

base_path = f"data/{degree}_2/{degree}_{degree}_{degree}_{degree}"
with open(f"{base_path}/basis_of_primitive_hodge_cycles.json", "r") as f:
    data = json.load(f)
    data["basis_of_primitive_hodge_cycles"] = sage_matrix_map(int, new_basis)
    with open(
        f"{base_path}/new_basis_of_primitive_hodge_cycles.json", "w"
    ) as g:
        json.dump(data, g)
