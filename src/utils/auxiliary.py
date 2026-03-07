from functools import reduce
from math import gcd
from .sage_imports import (
    Matrix_generic_dense,
    NumberField_generic,
    Matrix_integer_dense,
    block_matrix,
    identity_matrix,
    Matrix,
    ZZ,
)
from typing import Callable
from itertools import product as iterprod, combinations
from math import factorial
from collections import Counter


def lcm(*args):
    return reduce(lambda a, b: a * b // gcd(a, b), args)


def coprimes(n: int):
    """Returns all positive (> 0) integers coprime with n"""
    return [i for i in range(1, n + 1, 1) if gcd(i, n) == 1]


def positive_divisors(n: int) -> list[int]:
    """Returns the positive divisors of an integer"""
    return [k for k in range(1, abs(n) + 1, 1) if n % k == 0]


def Z_basis_of_kernel(
    K: NumberField_generic, M: Matrix_generic_dense
) -> Matrix_integer_dense:
    """
    Computes a Z-basis of the Z-module of vectors with integer
    coefficients inside the kernel of M. The matrix M is assumed
    to have entries inside the ring of integers of K. Returns a
    matrix whose columns are the elements of the basis.
    """
    if M.nrows() == 0:  # Then kernel is everything
        return identity_matrix(K, M.ncols())
    M_components = [
        Matrix(ZZ, M.nrows(), M.ncols(), lambda i, j: M[i, j].vector()[k])
        for k in range(K.degree())
    ]
    M_components = [M for M in M_components if not M.is_zero()]
    M_stacked = block_matrix([[M_component] for M_component in M_components])
    M_stacked_kernel_basis = M_stacked.right_kernel().basis()
    return Matrix(ZZ, M_stacked_kernel_basis).T


def sage_matrix_map(f, m):
    return [[f(m[i, j]) for j in range(m.ncols())] for i in range(m.nrows())]


def return_true(*args) -> bool:
    return True


def get_pairings(n: int) -> list[list[tuple[int, int]]]:
    """
    Returns all pairings of the sequence [0, ..., n-1], i.e.
    all partitions of said sequence in tuples of two elements.
    These are returned ordered in lexicographical order.
    """
    assert n % 2 == 0
    if n == 2:
        return [[(0, 1)]]
    output: list[list[tuple[int, int]]] = []
    for i in range(1, n, 1):
        pairings_rec = get_pairings(n - 2)
        for pairing_rec in pairings_rec:
            adjust = lambda x: x + 2 if x >= i - 1 else x + 1
            adjust_pair = lambda p: (adjust(p[0]), adjust(p[1]))
            output.append([(0, i)] + [adjust_pair(p) for p in pairing_rec])
    return output


def get_all_sums_equal_to_d(es: list[int], d: int) -> list[list[int]]:
    """
    Gets integers e_1, ..., e_n as input and returns all
    tuples (u_1, ..., u_n) of non negative integers u_i >= 0
    such that u_1 * e_1 + ... + u_n * e_n = d.
    """
    if es == []:
        return [[]] if d == 0 else []
    sums: list[list[int]] = []
    for u1 in range(d // es[0] + 1):
        rec = get_all_sums_equal_to_d(es[1:], d - es[0] * u1)
        for rec_sum in rec:
            sums.append([u1] + rec_sum)
    return sums


def signed_tuples(n, l):
    """
    Generator of all n-tuples with exactly l nonzero entries in {±1},
    such that the first nonzero entry is +1.
    """
    indices = range(n)
    for positions in combinations(indices, l):
        # positions of nonzero entries
        for signs in iterprod([-1, 1], repeat=l):
            # check first nonzero = +1 condition
            if signs[0] != 1:
                continue
            t = [0] * n
            for pos, s in zip(positions, signs):
                t[pos] = s
            yield tuple(t)


def count_fixed_permutations(a: tuple) -> int:
    """
    Given a list a = [a1, ..., an], returns the number of
    permutations of a that leaves a unchanged.
    """
    cnt = Counter(a)
    res = 1  # factorial(n)
    for c in cnt.values():
        res *= factorial(c)
    return res


def sage_RREF(
    M: Matrix_generic_dense, inv: Callable = lambda x: x.inverse()
) -> tuple[Matrix_generic_dense, Matrix_generic_dense, list[int]]:
    """
    Returns a pair (R, E) where R in in row echelon form,
    E is invertible and R = EM.
    """
    K = M.base_ring()
    rows, cols = M.nrows(), M.ncols()
    A = Matrix(K, rows, cols, M)
    E_total = identity_matrix(K, rows)
    pivot_row = 0
    pivot_col = 0

    def elementary_swap(n, i, j):
        E = identity_matrix(K, n)
        E[i, i], E[j, j] = 0, 0
        E[i, j], E[j, i] = 1, 1
        return E

    def elementary_scale(n, i, k):
        E = identity_matrix(K, n)
        E[i, i] = k
        return E

    def elementary_add(n, i, j, k):
        # row i += k * row j
        E = identity_matrix(K, n)
        E[i, j] = k
        return E

    while pivot_row < rows and pivot_col < cols:
        # Step 1: Find pivot in current column
        pivot = None
        for r in range(pivot_row, rows):
            if A[r, pivot_col] != 0:
                pivot = r
                break
        # If no pivot in this column → move to next column
        if pivot is None:
            pivot_col += 1
            continue
        # Step 2: Swap the pivot row into place
        if pivot != pivot_row:
            E = elementary_swap(rows, pivot_row, pivot)
            A = E * A
            E_total = E * E_total
        # Step 3: Scale pivot row to make pivot = 1
        pivot_val = A[pivot_row, pivot_col]
        if pivot_val != 1:
            pivot_val_inv = inv(pivot_val)
            E = elementary_scale(rows, pivot_row, pivot_val_inv)
            A = E * A
            E_total = E * E_total
        # Step 4: Eliminate rows below
        for r in range(pivot_row + 1, rows):
            if A[r, pivot_col] != 0:
                factor = -A[r, pivot_col]
                E = elementary_add(rows, r, pivot_row, factor)
                A = E * A
                E_total = E * E_total
        pivot_row += 1
        pivot_col += 1
    pivots: list[int] = []
    curr_pivot_row = 0
    for j in range(cols):
        if A[curr_pivot_row, j] != 0:
            pivots.append(j)
            curr_pivot_row += 1
            if curr_pivot_row >= rows:
                break
    return A, E_total, pivots  # R = E * M
