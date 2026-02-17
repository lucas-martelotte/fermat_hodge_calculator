from ..diagonal_variety import EvenDimensionalDiagonalVariety
from .basic_hodge_calculator import BasicHodgeCalculator
from ..utils.sage_imports import (
    PolynomialQuotientRingElement,
    PolynomialQuotientRing_generic,
    Matrix_generic_dense,
    NumberFieldElement,
    zero_matrix,
    factorial,
    Matrix,
    prod,
)
from ..utils.auxiliary import sage_matrix_map
from typing import Mapping, Sequence, Any
from tqdm import tqdm
from itertools import product as iterprod
from ..utils.settings import NCORES
from multiprocessing import Pool
from ..formal_algebra import FormalNumberField


class HodgeCalculator(BasicHodgeCalculator):
    def __init__(
        self,
        diagonal_variety: EvenDimensionalDiagonalVariety,
        gamma_values: Mapping[tuple[int, ...], PolynomialQuotientRingElement],
        basepath: str = "data/",
        complete_override: bool = False,
        override_files: list[str] | None = None,
    ):
        self.gamma_values = gamma_values
        super().__init__(
            diagonal_variety,
            basepath=basepath,
            complete_override=complete_override,
            override_files=override_files,
        )
        self.file_to_computation.update(
            {
                "periods_of_primitive_hodge_cycles": self.__compute_periods_of_primitive_hodge_cycles
            }
        )

    def __compute_periods_of_primitive_hodge_cycles(
        self,
    ) -> dict[str, Any]:
        rank_of_hodge_cycles = self.get_rank_of_primitive_hodge_cycles()
        if rank_of_hodge_cycles == 0:
            raise Exception(
                f"Cannot compute periods for variety of exponents {self.variety.exps}, "
                + "space of primitive Hodge cycles has dimension 0!"
            )
        I = self.multi_indexes
        J_hodge = self.get_form_idxs_at_infty_in_middle_hodge_comp()
        hodge_cycle_basis = self.get_basis_of_primitive_hodge_cycles()
        full_period_mat = zero_matrix(
            self.variety.base_field, len(I), len(J_hodge)
        )
        full_period_mat = (
            get_period_matrix_of_prim_hodge_cycles_multi_processing(
                self.variety.formal_base_field,
                self.variety.weights,
                self.variety.degree,
                I,
                J_hodge,
                self.gamma_values,
            )
        )
        hodge_period_mat = (full_period_mat.T * hodge_cycle_basis).T
        return {
            "period_matrix_of_primitive_hodge_cycles": sage_matrix_map(
                str, hodge_period_mat
            ),
        }

    # =============== #
    # === GETTERS === #
    # =============== #

    def get_period_matrix_of_primitive_hodge_cycles(
        self,
    ) -> Matrix_generic_dense:
        """
        Returns the full period matrix whose (i,j)-th coordinate is the
        integration pairing between the i-th element in the basis of primitive
        hodge cycles and the j-th form in 'form_idxs_at_infty_in_middle_hodge_component'.
        """
        eval = self.variety.formal_base_field.from_str
        matrix_raw = self.get_data_from_json(
            "periods_of_primitive_hodge_cycles",
            "period_matrix_of_primitive_hodge_cycles",
        )
        return Matrix(
            self.variety.base_field,
            [[eval(p) for p in row] for row in matrix_raw],
        )


# ============================== #
# === MULTIPROCESSING SCRIPT === #
# ============================== #


def full_pairing(
    zetad: NumberFieldElement,
    weights: list[int],
    form: tuple[int, ...],
    cycle: tuple[int, ...],
    gamma_value: NumberFieldElement | None,
) -> NumberFieldElement:
    """
    This is the full cycle-form pairing (given by integration) including
    the Gamma factors. The cycle is any cycle and the form is written
    in the standard basis, i.e. it is a tuple in the list multi_indexes.
    """
    if gamma_value is None:  # Then it is not algebraic, so period is zero
        return 0
    n = len(weights) - 2
    simple_pair = prod(
        zetad ** (weights[i] * (cycle[i] + 1) * (form[i] + 1))
        - zetad ** (weights[i] * cycle[i] * (form[i] + 1))
        for i in range(n + 2)
    )
    if simple_pair == 0:
        return 0
    return (-1) ** (n + 1) * simple_pair * gamma_value / factorial(n // 2)


def full_pairing_matrix_entry_calculator_worker(
    args: tuple[Any, ...],
) -> tuple[int, int, NumberFieldElement]:
    zetad, weights, form, cycle, i, j, gamma_value = args
    return i, j, full_pairing(zetad, weights, form, cycle, gamma_value)


def get_period_matrix_of_prim_hodge_cycles_multi_processing(
    K_formal: FormalNumberField,
    weights: tuple[int, ...],
    degree: int,
    I: Sequence[tuple[int, ...]],
    J_hodge: Sequence[int],
    gamma_values: Mapping[tuple[int, ...], PolynomialQuotientRingElement],
) -> Matrix_generic_dense:
    if len(I) == 0 or len(J_hodge) == 0:
        raise Exception("Period matrix has zero rows (or columns).")
    period_matrix = zero_matrix(K_formal.K, len(I), len(J_hodge))
    zetad = K_formal.from_str(f"zeta{degree}")
    with Pool(NCORES) as pool:
        gamma_value_key = lambda j: tuple(
            sorted(
                [
                    (formi + 1) * wi
                    for formi, wi in zip(I[J_hodge[j]], weights)
                ]
            )
        )
        tasks = [
            (
                zetad,
                weights,
                I[J_hodge[j]],
                I[i],
                i,
                j,
                gamma_values.get(gamma_value_key(j)),
            )
            for i, j in iterprod(range(len(I)), range(len(J_hodge)))
        ]
        for i, j, val in tqdm(
            pool.imap_unordered(
                full_pairing_matrix_entry_calculator_worker,
                tasks,
                chunksize=max(len(tasks) // (10 * NCORES), 1),
            ),
            desc="Assembling full period matrix",
            total=len(tasks),
        ):
            period_matrix[i, j] = val
    return period_matrix


# gamma_tuple = tuple(sorted([(formi + 1) * wi for formi, wi in zip(form, ws)]))
# gamma_value = self.gamma_values[gamma_tuple]
