from ..diagonal_variety import EvenDimensionalDiagonalVariety
from .basic_hodge_calculator import BasicHodgeCalculator
from ..utils.sage_imports import (
    PolynomialQuotientRingElement,
    PolynomialQuotientRing_generic,
    MPolynomial_element,
    Matrix_rational_dense,
    Matrix_integer_dense,
    Matrix_generic_dense,
    NumberFieldElement,
    identity_matrix,
    zero_matrix,
    factorial,
    binomial,
    Matrix,
    prod,
    Rat,
    QQ,
    ZZ,
)
from ..utils.auxiliary import sage_matrix_map
from typing import Mapping, Sequence, Any
from tqdm import tqdm
from itertools import product as iterprod
from ..utils.settings import NCORES
from multiprocessing import Pool
from ..formal_algebra import FormalNumberField
from ..hodge_cycles import (
    HodgeCycleFactory,
    AlgebraicPrimePrimitiveHodgeCycle,
)


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
                "periods_of_primitive_hodge_cycles": self.__compute_periods_of_primitive_hodge_cycles,
                "intersection_of_primitive_hodge_cycles": self.__compute_intersection_of_primitive_hodge_cycles,
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
        full_period_mat = mp_get_period_matrix_of_prim_hodge_cycles(
            self.variety.formal_base_field,
            self.variety.weights,
            self.variety.degree,
            I,
            J_hodge,
            self.gamma_values,
        )
        hodge_period_mat = (full_period_mat.T * hodge_cycle_basis).T
        return {
            "period_matrix_of_primitive_hodge_cycles": sage_matrix_map(
                str, hodge_period_mat
            ),
        }

    def __compute_intersection_of_primitive_hodge_cycles(self):
        """
        Computes the cup pairing matrix, in which the (i,j)-th entry
        is the cup pairing (given by integrating the cup product) of the
        i-th and j-th n-forms at infinity.

        Also computes the dual period matrix, in which the i-th row contains
        the coefficients of the poincare dual of the i-th element of the
        basis of primitive hodge cycles when written as a linear combination
        of the forms in 'indexes_of_forms_for_computing_primitive_hodge_periods'.
        """
        I = self.multi_indexes
        J_inf = self.get_form_idxs_at_infty()
        J_hodge = self.get_form_idxs_at_infty_in_middle_hodge_comp()
        period_matrix = self.get_period_matrix_of_primitive_hodge_cycles()
        full_cup_pairing_matrix = mp_get_cup_pairing_matrix_of_forms_at_infty(
            I, J_inf, self.variety.exps, self.variety.degree
        )
        N, rel_idxs = len(J_hodge), [J_inf.index(j) for j in J_hodge]
        cup = lambda i, j: full_cup_pairing_matrix[rel_idxs[i]][rel_idxs[j]]
        cup_pairing_matrix_at_middle_hodge_comp = Matrix(QQ, N, N, cup)
        cup_inv = cup_pairing_matrix_at_middle_hodge_comp.inverse()
        dual_period_matrix = period_matrix * cup_inv
        intersection_matrix = (
            dual_period_matrix
            * cup_pairing_matrix_at_middle_hodge_comp
            * dual_period_matrix.T
        )  # (DUAL) PERIOD MATRIX IS SLOWEST BECAUSE NF
        return {
            "cup_pairing_matrix_of_forms_at_infinity": [
                [str(p) for p in row] for row in full_cup_pairing_matrix
            ],
            "cup_pairing_matrix_of_forms_at_infinity_inside_the_middle_hodge_comp": [
                [str(p) for p in row]
                for row in cup_pairing_matrix_at_middle_hodge_comp
            ],
            "dual_period_matrix_of_primitive_hodge_cycles": sage_matrix_map(
                str, dual_period_matrix
            ),
            "intersection_matrix_of_primitive_hodge_cycles": sage_matrix_map(
                str, intersection_matrix
            ),
        }

    # ============================================== #
    # === INTERACTION WITH HODGE CYCLE FACTORIES === #
    # ============================================== #

    def compute_periods_from_hodge_cycle_factory(
        self,
        factory: HodgeCycleFactory,
        filename: str,
        override: bool = False,
    ) -> None:
        """
        Computes the period vector of each cycle produced by
        the input of cycle factory and compiles it all in a
        single matrix and, together with some additional
        info, exports it to memory.
        """
        self.check_folder()
        if self._check_json(filename) and not override:
            return  # File already exists
        nf = self.variety.base_field
        xs = self.variety.polynomial_variables
        id_to_cycle = factory.get_hodge_cycles()
        ids, cycles = list(id_to_cycle.keys()), list(id_to_cycle.values())
        I = self.multi_indexes
        J_hodge = self.get_form_idxs_at_infty_in_middle_hodge_comp()
        ms = self.variety.exps
        # ======================================= #
        period_matrix = mp_get_hodge_cycle_factory_period_matrix(
            nf, I, J_hodge, cycles, ms, xs
        )
        # period_matrix = Matrix(nf, len(cycles), len(J_hodge), 0)
        # for i, j in tqdm(
        #    iterprod(range(len(cycles)), range(len(J_hodge))),
        #    desc="Computing periods",
        #    total=len(cycles) * len(J_hodge),
        # ):
        #    form = I[J_hodge[j]]
        #    form_poly = prod([xi**bi for xi, bi in zip(xs, form)])
        #    period = cycles[i].pairing(form_poly) * prod(ms)
        #    period_matrix[i, j] = period
        # ======================================= #
        cycle_coords = self.solve_from_periods(period_matrix)
        data = {
            "rank": int(cycle_coords.rank()),
            "number_of_cycles": len(ids),
            "cycle_ids": ids,
            "period_matrix": sage_matrix_map(str, period_matrix),
            "coordinates": sage_matrix_map(str, cycle_coords),
        }
        self.write_data_on_json(filename, data)

    def solve_from_periods(
        self, period_matrix_to_solve: Matrix_generic_dense
    ) -> Matrix_rational_dense:
        """
        The input matrix's (i,j)-th coordinate should be the period
        of some to-be-discovered i-th primitive hodge cycle with the
        j-th form from 'form_idx_at_infty_in_middle_hodge_comp'. This
        function uses some linear algebra to discover those cycles,
        and then returns a matrix whose i-th column is the i-th cycle
        writen on the standard Z-basis of primitive hodge cycles.
        """
        hodge_period_matrix = (
            self.get_period_matrix_of_primitive_hodge_cycles()
        )
        K_formal = self.variety.formal_base_field
        period_coords = K_formal.solve_left_in_rationals(
            hodge_period_matrix, period_matrix_to_solve
        ).T
        return Matrix(QQ, period_coords)  # Solution must lie inside (1/d)*ZZ

    def get_hodge_cycle_factory_data(self, filename: str) -> dict[str, Any]:
        """
        Loads the data exported from the function named
        'compute_periods_from_hodge_cycle_factory'.
        """
        data = self.get_json(filename)
        nf = self.variety.base_field
        eval = self.variety.formal_base_field.from_str
        period_matrix = Matrix(
            nf, [[eval(p) for p in row] for row in data["period_matrix"]]
        )
        return {
            "rank": data["rank"],
            "cycle_ids": data["cycle_ids"],
            "period_matrix": period_matrix,
            "coordinates": Matrix(QQ, data["coordinates"]),
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

    def get_cup_pairing_matrix_of_forms_at_infty(
        self,
    ) -> Matrix_rational_dense:
        """
        Returns the full gram matrix of the cup pairing of forms at infinity,
        whose (i,j)-th coordinate corresponds to the integration of the cup
        product of the i-th and j-th forms at infinity.
        """
        matrix_raw = self.get_data_from_json(
            "intersection_of_primitive_hodge_cycles",
            "cup_pairing_matrix_of_forms_at_infinity",
        )
        return Matrix(QQ, [[Rat(p) for p in row] for row in matrix_raw])

    def get_cup_pairing_matrix_of_forms_at_infty_in_middle_hodge_comp(
        self,
    ) -> Matrix_rational_dense:
        """
        Returns the gram matrix of the cup pairing of forms at infinity which
        lie inside the middle hodge component. The (i,j)-th coordinate corresponds
        to the integration of the cup product of the i-th and j-th forms in the
        list 'form_idxs_at_infty_in_middle_hodge_comp'.
        """
        matrix_raw = self.get_data_from_json(
            "intersection_of_primitive_hodge_cycles",
            "cup_pairing_matrix_of_forms_at_infinity_inside_the_middle_hodge_component",
        )
        return Matrix(QQ, [[Rat(p) for p in row] for row in matrix_raw])

    def get_dual_period_matrix_of_primitive_hodge_cycles(
        self,
    ) -> Matrix_generic_dense:
        """
        Returns the dual period matrix of primitive hodge cycles, whose i-th
        row contains the coefficients of the poincare dual of the i-th element
        in the standard basis returned by 'get_basis_of_primitive_hodge_cycles'
        """
        nf = self.variety.base_field
        eval = self.variety.formal_base_field.from_str
        matrix_raw = self.get_data_from_json(
            "intersection_of_primitive_hodge_cycles",
            "dual_period_matrix_of_primitive_hodge_cycles",
        )
        return Matrix(
            nf,
            [[eval(p) for p in row] for row in matrix_raw],
        )

    def get_intersection_matrix_of_primitive_hodge_cycles(
        self,
    ) -> Matrix_rational_dense:
        """
        Returns the intersection pairing matrix whose (i,j)-th entry
        corresponds to the intersection pairing of the i-th and j-th elements
        in the basis returned by 'get_basis_of_primitive_hodge_cycles'.
        """
        matrix_raw = self.get_data_from_json(
            "intersection_of_primitive_hodge_cycles",
            "intersection_matrix_of_primitive_hodge_cycles",
        )
        return Matrix(QQ, matrix_raw)


# ================================================================ #
# === MULTIPROCESSING SCRIPT FOR COMPUTING FULL PAIRING MATRIX === #
# ================================================================ #


def mp_full_pairing(
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


def mp_full_pairing_worker(
    args: tuple[Any, ...],
) -> tuple[int, int, NumberFieldElement]:
    zetad, weights, form, cycle, i, j, gamma_value = args
    return i, j, mp_full_pairing(zetad, weights, form, cycle, gamma_value)


def mp_get_period_matrix_of_prim_hodge_cycles(
    K_formal: FormalNumberField,
    ws: tuple[int, ...],
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
            sorted([(bi + 1) * wi for bi, wi in zip(I[J_hodge[j]], ws)])
        )
        idx_to_val = lambda idx: gamma_values.get(gamma_value_key(idx))
        tasks = [
            (zetad, ws, I[J_hodge[j]], I[i], i, j, idx_to_val(j))
            for i, j in iterprod(range(len(I)), range(len(J_hodge)))
        ]
        for i, j, val in tqdm(
            pool.imap_unordered(
                mp_full_pairing_worker,
                tasks,
                chunksize=max(len(tasks) // (10 * NCORES), 1),
            ),
            desc="Assembling full period matrix",
            total=len(tasks),
        ):
            period_matrix[i, j] = val
    return period_matrix


# =============================================================== #
# === MULTIPROCESSING SCRIPT FOR COMPUTING CUP PAIRING MATRIX === #
# =============================================================== #


def mp_weighted_sum_of_index(
    ms: tuple[int, ...], b: tuple[int, ...], c: int = 1
) -> Rat:
    return sum([Rat((c * (bi + 1)) % mi / mi) for bi, mi in zip(b, ms)])


def mp_cup_pairing(
    form1: tuple[int, ...],
    form2: tuple[int, ...],
    ms: tuple[int, ...],
    d: int,
) -> Rat:
    """
    Computes the cup pairing between two forms n-forms at infinity
    (given by integration of their cup product). It is amazing
    that this gives a rational number!
    """
    n = len(ms) - 2
    A1 = mp_weighted_sum_of_index(ms, form1)
    A2 = mp_weighted_sum_of_index(ms, form2)
    assert A1.is_integral() and A2.is_integral()  # forms are at infinity
    if any(bi1 + bi2 + 2 != mi for bi1, bi2, mi in zip(form1, form2, ms)):
        return 0
    output_num = (-1) ** (binomial(n + 1, 2) + A1) * d * prod(ms)
    output_den = factorial(A1 - 1) * factorial(A2 - 1)
    return output_num / output_den


def mp_cup_pairing_worker(
    args: tuple[Any, ...],
) -> tuple[int, int, NumberFieldElement]:
    i, j, formi, formj, ms, d = args
    return i, j, mp_cup_pairing(formi, formj, ms, d)


def mp_get_cup_pairing_matrix_of_forms_at_infty(
    I: Sequence[tuple[int, ...]],
    J_inf: Sequence[int],
    ms: tuple[int, ...],
    d: int,
) -> Matrix_rational_dense:
    cup_pairing_matrix = zero_matrix(QQ, len(J_inf), len(J_inf))
    with Pool(NCORES) as pool:
        tasks = [
            (i, j, I[J_inf[i]], I[J_inf[j]], ms, d)
            for i, j in iterprod(range(len(J_inf)), range(len(J_inf)))
        ]
        for i, j, val in tqdm(
            pool.imap_unordered(
                mp_cup_pairing_worker,
                tasks,
                chunksize=max(len(tasks) // (10 * NCORES), 1),
            ),
            desc="Assembling cup period matrix",
            total=len(tasks),
        ):
            cup_pairing_matrix[i, j] = val
    return cup_pairing_matrix


# ======================================================================== #
# === MULTIPROCESSING SCRIPT FOR COMPUTING HODGE CYCLE FACTORY PERIODS === #
# ======================================================================== #


def mp_hodge_cycle_factory_period_matrix_worker(
    args: tuple[Any, ...],
) -> tuple[int, int, NumberFieldElement]:
    i, j, cycle, form_poly, ms = args
    return i, j, cycle.pairing(form_poly) * prod(ms)


def mp_get_hodge_cycle_factory_period_matrix(
    K: PolynomialQuotientRing_generic,
    I: Sequence[tuple[int, ...]],
    J_hodge: Sequence[int],
    cycles: list[AlgebraicPrimePrimitiveHodgeCycle],
    ms: tuple[int, ...],
    xs: tuple[MPolynomial_element, ...],
) -> Matrix_generic_dense:
    period_matrix = zero_matrix(K, len(cycles), len(J_hodge))
    with Pool(NCORES) as pool:
        form = lambda j: I[J_hodge[j]]
        form_poly = lambda j: prod([xi**bi for xi, bi in zip(xs, form(j))])
        tasks = [
            (i, j, cycles[i], form_poly(j), ms)
            for i, j in iterprod(range(len(cycles)), range(len(J_hodge)))
        ]
        for i, j, val in tqdm(
            pool.imap_unordered(
                mp_hodge_cycle_factory_period_matrix_worker,
                tasks,
                chunksize=max(len(tasks) // (10 * NCORES), 1),
            ),
            desc="Assembling hodge cycle factory period matrix",
            total=len(tasks),
        ):
            period_matrix[i, j] = val
    return period_matrix


# period_matrix = Matrix(nf, len(cycles), len(J_hodge), 0)
# for i, j in tqdm(
#    iterprod(range(len(cycles)), range(len(J_hodge))),
#    desc="Computing periods",
#    total=len(cycles) * len(J_hodge),
# ):
#    form = I[J_hodge[j]]
#    form_poly = prod([xi**bi for xi, bi in zip(xs, form)])
#    period = cycles[i].pairing(form_poly) * prod(ms)
#    period_matrix[i, j] = period
