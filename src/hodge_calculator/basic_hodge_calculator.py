from ..diagonal_variety import EvenDimensionalDiagonalVariety
from typing import Any, Callable
from os.path import exists
from os import makedirs
import json
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
    zero_matrix,
    parallel,
    Vector_integer_dense,
    block_matrix,
    identity_matrix
)
from src.utils.auxiliary import coprimes, Z_basis_of_kernel
from src.utils.json_manager import JsonManager
from multiprocessing import Pool
from time import perf_counter


class BasicHodgeCalculator(JsonManager):
    def __init__(
        self,
        diagonal_variety: EvenDimensionalDiagonalVariety,
        basepath: str = "data/",
        complete_override: bool = False,
        override_files: list[str] | None = None,
    ):
        degree, dimension = diagonal_variety.degree, diagonal_variety.dimension
        filename = "_".join([str(e) for e in diagonal_variety.exps])
        savepath = basepath + f"/{degree}_{dimension}/{filename}"
        super().__init__(savepath, complete_override, override_files)
        self.variety = diagonal_variety
        self.file_to_computation["basis_of_primitive_hodge_cycles"] = (
            self.__compute_basis_of_primitive_hodge_cycles
        )

    def __compute_basis_of_primitive_hodge_cycles(self) -> dict[str, Any]:
        d, n = self.variety.degree, self.variety.dimension
        ms, k = self.variety.exps, n // 2 + 1
        H, H0 = list(range(n + 2)), list(range(1, k, 1))
        A = lambda b, c: self.weighted_sum_of_index(b, c)
        I = list(iterprod(*(list(range(mi - 1)) for mi in ms)))
        at_infinity = lambda i: A(I[i], 1).is_integral()
        J_inf = [i for i in range(len(I)) if at_infinity(i)]
        J = [i for i in range(len(I)) if A(I[i], 1) in H0 or not at_infinity(i)]
        J_hodge = [i for i in range(len(I)) if A(I[i], 1) == k]
        J_hodge_alg = [i for i in J_hodge if all(A(I[i], c) == k for c in coprimes(d))]
        K = CyclotomicField(d)
        O, zetad = K.ring_of_integers(), K.gen()
        # M = zero_matrix(O, len(J), len(I))

        basis_of_primitive_hodge_cycles = get_Z_basis_of_kernel_of_pairing(
            self.variety.weights, self.variety.degree, I, J
        )

        # @parallel()
        # def compute_entry(index: tuple[int, int]) -> None:
        #    i, j = index[0], index[1]
        #    M[i, j] = self.simple_pairing(zetad, list(I[J[i]]), list(I[j]))

        # input_list = list(iterprod(range(len(J)), range(len(I))))
        # print(input_list[:10])
        # compute_entry(input_list)  # type: ignore
        # print(M[0, 10])

        # start_time = perf_counter()
        # M = zero_matrix(O, len(J), len(I))
        # for i, j in tqdm(
        #    iterprod(range(len(J)), range(len(I))),
        #    desc="Computing basis of primitive hodge cycles",
        #    total=len(J) * len(I),
        # ):
        #    M[i, j] = self.simple_pairing(zetad, list(I[J[i]]), list(I[j]))
        # basis_of_primitive_hodge_cycles = Z_basis_of_kernel(K, M)
        # end_time = perf_counter()
        # print(end_time - start_time)

        dimension_prim_hodge = int(basis_of_primitive_hodge_cycles.rank())
        return {
            "rank_of_space_of_primitive_hodge_cycles": dimension_prim_hodge,
            "all_multi_indexes": I,
            "form_indexes_at_infinity": J_inf,
            "get_form_indexes_either_not_at_infinity_or_before_the_middle_component": J,
            "form_indexes_at_infinity_inside_the_middle_hodge_component": J_hodge,
            "form_indexes_at_infinity_inside_the_middle_hodge_component_with_algebraic_gamma": J_hodge_alg,
            "basis_of_primitive_hodge_cylces": [
                [int(r) for r in row] for row in basis_of_primitive_hodge_cycles
            ],
        }

    def simple_pairing(
        self, zetad: NumberFieldElement, form: list[int], cycle: list[int]
    ) -> NumberFieldElement:
        """
        Returns a simplified version of the cycle-form pairing (given by
        integration) which only serves the purpose of determining if the
        pairing is equal to zero or not.
        """
        n = self.variety.dimension
        ws = self.variety.weights
        return prod(
            zetad ** (ws[i] * (cycle[i] + 1) * (form[i] + 1))
            - zetad ** (ws[i] * cycle[i] * (form[i] + 1))
            for i in range(n + 2)
        )

    def weighted_sum_of_index(self, b: list[int], c: int = 1) -> Rat:
        """
        Compute the value of 'A_{cb}', i.e. the summation given by
        <c(b_0+1)/m_0> + ... + <c(b_{n+1}+1)/m_{m+1}> where <x> stands for
        the fractional part operator (which since b is a tuple of
        positive integers, is just taking b_i modulo m_i).
        """
        ms = self.variety.exps
        return sum([Rat((c * (bi + 1)) % mi / mi) for bi, mi in zip(b, ms)])

    # =============== #
    # === GETTERS === #
    # =============== #


# ============================== #
# === MULTIPROCESSING SCRIPT === #
# ============================== #

NCORES = 8


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


def get_Z_basis_of_kernel_of_pairing(
    weights: list[int], degree: int, I: list[tuple[int, ...]], J: list[int]
) -> Matrix_integer_dense:
    if len(J) == 0: # Then the kernel is everything
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
    M = block_matrix(ZZ, [[res[(i, j)] for j in range(len(I))] for i in range(len(J))])
    return Matrix(ZZ, M.right_kernel().basis()).T
