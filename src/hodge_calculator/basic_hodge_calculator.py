from ..diagonal_variety import EvenDimensionalDiagonalVariety
from typing import Any, Sequence
from tqdm import tqdm  # type: ignore
from itertools import product as iterprod
from src.utils.sage_imports import (
    CyclotomicField,
    Matrix,
    NumberFieldElement,
    Matrix_integer_dense,
    Rat,
    prod,
    ZZ,
    Vector_integer_dense,
    block_matrix,
    identity_matrix,
)
from src.utils.auxiliary import coprimes
from src.utils.json_manager import JsonManager
from multiprocessing import Pool
from ..utils.settings import NCORES


class BasicHodgeCalculator(JsonManager):
    def __init__(
        self,
        diagonal_variety: EvenDimensionalDiagonalVariety,
        basepath: str = "data/",
        complete_override: bool = False,
        override_files: list[str] | None = None,
    ):
        degree, dimension = (
            diagonal_variety.degree,
            diagonal_variety.dimension,
        )
        filename = "_".join([str(e) for e in diagonal_variety.exps])
        savepath = basepath + f"/{degree}_{dimension}/{filename}"
        super().__init__(savepath, complete_override, override_files)
        self.variety = diagonal_variety
        self.file_to_computation["basis_of_primitive_hodge_cycles"] = (
            self.__compute_basis_of_primitive_hodge_cycles
        )
        self._multi_indexes: Sequence[tuple[int, ...]] | None = None

    @property
    def multi_indexes(self) -> Sequence[tuple[int, ...]]:
        """
        Returns the indexes in lexicographical order, i.e. the
        first elements are (0, ..., 0), (0, ..., 0, 1) and the
        last element is (m_1-2, ..., m_n-2).
        """
        if self._multi_indexes is None:
            ms = self.variety.exps
            self._multi_indexes = list(
                iterprod(*(list(range(mi - 1)) for mi in ms))
            )
        return self._multi_indexes

    def weighted_sum_of_index(self, b: list[int], c: int = 1) -> Rat:
        """
        Compute the value of 'A_{cb}', i.e. the summation given by
        <c(b_0+1)/m_0> + ... + <c(b_{n+1}+1)/m_{m+1}> where <x> stands for
        the fractional part operator (which since b is a tuple of
        positive integers, is just taking b_i modulo m_i).
        """
        ms = self.variety.exps
        return sum([Rat((c * (bi + 1)) % mi / mi) for bi, mi in zip(b, ms)])

    def __compute_basis_of_primitive_hodge_cycles(self) -> dict[str, Any]:
        d, n = self.variety.degree, self.variety.dimension
        k = n // 2 + 1
        H0 = list(range(1, k, 1))
        A = lambda b, c: self.weighted_sum_of_index(b, c)
        I = self.multi_indexes
        at_infty = lambda i: A(I[i], 1).is_integral()
        J_inf = [i for i in range(len(I)) if at_infty(i)]
        J = [i for i in range(len(I)) if A(I[i], 1) in H0 or not at_infty(i)]
        J_hodge = [i for i in range(len(I)) if A(I[i], 1) == k]
        J_hodge_alg = [
            i for i in J_hodge if all(A(I[i], c) == k for c in coprimes(d))
        ]
        basis_of_primitive_hodge_cycles = (
            get_Z_basis_of_kernel_of_pairing_multi_processing(
                self.variety.weights, self.variety.degree, I, J
            )
        )
        dimension_prim_hodge = int(basis_of_primitive_hodge_cycles.rank())
        return {
            "rank_of_space_of_primitive_hodge_cycles": dimension_prim_hodge,
            "form_indexes_at_infinity": J_inf,
            "get_form_indexes_either_not_at_infinity_or_before_the_middle_component": J,
            "form_indexes_at_infinity_inside_the_middle_hodge_component": J_hodge,
            "form_indexes_at_infinity_inside_the_middle_hodge_component_with_algebraic_gamma": J_hodge_alg,
            "basis_of_primitive_hodge_cycles": [
                [int(r) for r in row]
                for row in basis_of_primitive_hodge_cycles
            ],
        }

    # =============== #
    # === GETTERS === #
    # =============== #

    def get_form_idxs_at_infty(self) -> list[int]:
        """
        In the list of all multi indexes, returns all positions
        corresponding to n-forms 'at infinity'.
        """
        return self.get_data_from_json(
            "basis_of_primitive_hodge_cycles", "form_indexes_at_infinity"
        )

    def get_form_idxs_not_at_infty_or_before_the_middle_hodge_comp(
        self,
    ) -> list[int]:
        """
        In the list of all multi indexes, returns all positions
        corresponding to one of the following two:
        1) n-forms which are NOT at infinity
        2) n-forms at infinity lying inside some H^{p,q} with p < n/2
        The hodge cycles can be defined as those whose pairing vanishes
        are those forms.
        """
        return self.get_data_from_json(
            "basis_of_primitive_hodge_cycles",
            "get_form_indexes_either_not_at_infinity_or_before_the_middle_component",
        )

    def get_form_idxs_at_infty_in_middle_hodge_comp(
        self,
    ) -> list[int]:
        """
        In the list of all multi indexes returned by the function
        'get_all_multi_indexes', returns all positions in that list
        which correspond to (n/2)-forms at infinity which lie inside
        the middle component H^{n/2, n/2}. These forms are used to
        compute the periods of hodge cycles.
        """
        return self.get_data_from_json(
            "basis_of_primitive_hodge_cycles",
            "form_indexes_at_infinity_inside_the_middle_hodge_component",
        )

    def get_form_idxs_at_infty_in_middle_hodge_comp_with_alg_gamma(
        self,
    ) -> list[int]:
        """
        In the list of all multi indexes returned by the function
        'get_all_multi_indexes', returns all positions in that list
        which correspond to (n/2)-forms at infinity which lie inside
        the middle component H^{n/2, n/2} and such that the gamma
        that arises when integrating over it is algebraic.
        """
        return self.get_data_from_json(
            "basis_of_primitive_hodge_cycles",
            "form_indexes_at_infinity_inside_the_middle_hodge_component_with_algebraic_gamma",
        )

    def get_basis_of_primitive_hodge_cycles(self) -> Matrix_integer_dense:
        """
        Returns a matrix whose i-th column is the i-th element on a
        basis of primitive hodge cycles of X, represented as an integer
        vector in the standard Z-basis of cycles on Y indexed by
        the list returned by 'get_all_multi_indexes'.
        """
        matrix_raw = self.get_data_from_json(
            "basis_of_primitive_hodge_cycles",
            "basis_of_primitive_hodge_cycles",
        )
        return Matrix(ZZ, matrix_raw)

    def get_rank_of_primitive_hodge_cycles(self) -> int:
        """Returns the rank of the lattice of primitive hodge cycles of X."""
        return self.get_data_from_json(
            "basis_of_primitive_hodge_cycles",
            "rank_of_space_of_primitive_hodge_cycles",
        )


# ============================== #
# === MULTIPROCESSING SCRIPT === #
# ============================== #


def simple_pairing(
    zetad: NumberFieldElement,
    weights: list[int],
    form: tuple[int, ...],
    cycle: tuple[int, ...],
) -> NumberFieldElement:
    """
    Returns a simplified version of the cycle-form pairing (given by
    integration) which only serves the purpose of determining if the
    pairing is equal to zero or not.
    """
    return prod(
        zetad ** (weights[i] * (cycle[i] + 1) * (form[i] + 1))
        - zetad ** (weights[i] * cycle[i] * (form[i] + 1))
        for i in range(len(weights))
    )


def simple_pairing_matrix_entry_calculator_worker(
    args: tuple[Any, ...],
) -> tuple[int, int, Vector_integer_dense]:
    zetad, weights, form, cycle, i, j = args
    vector = simple_pairing(zetad, weights, form, cycle).vector()
    return i, j, Matrix(ZZ, len(vector), 1, vector)


def get_Z_basis_of_kernel_of_pairing_multi_processing(
    weights: tuple[int, ...],
    degree: int,
    I: Sequence[tuple[int, ...]],
    J: Sequence[int],
) -> Matrix_integer_dense:
    if len(J) == 0:  # Then the kernel is everything
        return identity_matrix(ZZ, len(I))
    zetad = CyclotomicField(degree).gen()
    with Pool(NCORES) as pool:
        tasks = [
            (zetad, weights, I[J[i]], I[j], i, j)
            for i, j in iterprod(range(len(J)), range(len(I)))
        ]
        res: dict[tuple[int, int], Any] = {}
        for i, j, val in tqdm(
            pool.imap_unordered(
                simple_pairing_matrix_entry_calculator_worker,
                tasks,
                chunksize=max(len(tasks) // (10 * NCORES), 1),
            ),
            desc="Assembling simple pairing matrix",
            total=len(tasks),
        ):
            res[(i, j)] = val
    M = block_matrix(
        ZZ, [[res[(i, j)] for j in range(len(I))] for i in range(len(J))]
    )
    return Matrix(ZZ, M.right_kernel().basis()).T
