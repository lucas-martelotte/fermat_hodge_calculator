from itertools import product as iterprod
from typing import Any

from src.utils.sage_imports import (
    PolynomialQuotientRingElement,
    Vector_rational_dense,
    Matrix_rational_dense,
    MPolynomialRing_base,
    Matrix_generic_dense,
    MPolynomial_element,
    NumberFieldElement,
    zero_matrix,
    sage_eval,
    Matrix,
    vector,
    prod,
    Rat,
    RR,
    QQ,
)
from src.utils.auxiliary import return_true
from src.utils.settings import NCORES
from multiprocessing import Pool
from tqdm import tqdm


class FormalNumberField:
    def __init__(
        self,
        Kbase: MPolynomialRing_base,
        Kbase_equations: list[MPolynomial_element],
    ) -> None:
        """
        It is assumed that the input equations are
        [P_1, ..., P_n] where P_i(x_1, ..., x_i) is
        irreducible in (Kbase/(P_1, ..., P_{i-1}))[x_i].

        WARNING! It has to be in this order! The i-th equation
        has to correspond to the i-th variable, otherwise things
        are going to break.
        """
        self.Kbase_ideal = Kbase.ideal(Kbase_equations)
        self.Kbase_equations = Kbase_equations
        self.Kbase = Kbase
        self.Kbase_variables = self.Kbase.gens()

        vars, eqs = self.Kbase_variables, self.Kbase_equations
        degree = lambda eq, var: (
            eq.degree(var) if len(vars) > 1 else eq.degree()
        )
        self.Kbase_degrees = [degree(eq, var) for eq, var in zip(eqs, vars)]
        self.Kbase_subs = [
            var**deg - eq
            for eq, var, deg in zip(eqs, vars, self.Kbase_degrees)
        ]
        self.number_of_components = len(self.Kbase_variables)
        assert all(str(x).endswith("base") for x in self.Kbase_variables)

        varnames = [str(x)[:-4] for x in self.Kbase_variables]
        self.K = self.Kbase.quotient(self.Kbase_ideal, names=varnames)
        self.K.element_class.divides = return_true
        self.K.defining_ideal().is_prime = return_true
        self.K.is_field = return_true
        self.K_variables = self.K.gens()

        self.monomial_indexes = list(
            iterprod(*[range(deg) for deg in self.Kbase_degrees])
        )
        self.Kbase_monomials = [
            prod(xi**mi for xi, mi in zip(self.Kbase_variables, m))
            for m in self.monomial_indexes
        ]
        self.K_monomials = [
            prod(xi**mi for xi, mi in zip(self.K_variables, m))
            for m in self.monomial_indexes
        ]
        self.degree = len(self.K_monomials)
        self.locals = {str(x): x for x in self.K_variables}

    def from_str(self, expr: str) -> PolynomialQuotientRingElement:
        """Parses a string expression as an element of K"""
        local_variables = self.locals
        local_variables["zeta2"] = -1  # exception for degree 2
        return sage_eval(expr, locals=local_variables)

    def vector(
        self, x: PolynomialQuotientRingElement
    ) -> Vector_rational_dense:
        """
        Writes an element x in the quotient as a vector in the
        basis of monomials indexed by self.monomial_indexes.
        """
        x_lift, mons = self.canonical_lift(x), self.Kbase_monomials
        return vector([x_lift.monomial_coefficient(mon) for mon in mons])

    def convert_to_rational_horizontally_blocked_matrix(
        self, M: Matrix_generic_dense
    ) -> Matrix_rational_dense:
        """
        Receives a matrix M with entries in K and converts it to a
        blocked matrix [M_1 ... M_d], where d = self.degree, with
        rational entries. The (i,j)-th entry of the block M_k has
        the k-th coefficient of M[i,j] in the Q-basis of K indexed
        by self.monomial_indexes.
        """
        mons = self.Kbase_monomials
        vars = self.Kbase_variables
        degs = self.Kbase_degrees
        subs = self.Kbase_subs
        return mp_convert_to_rational_horizontally_blocked_matrix(
            M, mons, vars, degs, subs, self.Kbase
        )

    def solve_left_in_rationals(
        self, M: Matrix_generic_dense, V: Matrix_generic_dense
    ) -> Matrix_rational_dense:
        """Solves XM = V, where X is a rational matrix."""
        M_rat = self.convert_to_rational_horizontally_blocked_matrix(M)
        V_rat = self.convert_to_rational_horizontally_blocked_matrix(V)
        return M_rat.solve_left(V_rat)

    def generators_of_right_kernel(
        self, M: Matrix_generic_dense
    ) -> Matrix_generic_dense:
        """Returns a matrix E whose columns generate the right-kernel of M."""
        mons, d = self.K_monomials, self.degree
        M_rat = Matrix(QQ, M.nrows() * d, M.ncols() * d, 0)
        print(M_rat.nrows(), M_rat.ncols())
        for i, j in iterprod(range(M.nrows()), range(M.ncols())):
            for a in range(len(mons)):
                curr_vec = list(self.vector(M[i, j] * mons[a]))
                for b in range(len(curr_vec)):
                    M_rat[i * d + b, j * d + a] = curr_vec[b]
        kernel_rat = Matrix(QQ, M_rat.right_kernel().basis()).transpose()
        kernel = Matrix(self.K, M.ncols(), kernel_rat.ncols(), 0)
        for i, j in iterprod(range(M.ncols()), range(kernel_rat.ncols())):
            kernel[i, j] = sum(
                kernel_rat[i * d + k, j] * mons[k] for k in range(len(mons))
            )
        return kernel

    def canonical_lift(
        self, P: PolynomialQuotientRingElement
    ) -> MPolynomial_element:
        """
        Returns the element P lifted to Kbase in such a way that the degree
        of each variable is strictly less than their respective equations
        (ex.: if the variable zbase has equation zbase**k - F(zbase, ybase, xbase),
        where F has degree strictly less than k in zbase, then the canonical lift
        will ensure that the degree in zbase is also strictly less than k).
        """
        return self._canonical_lift_rec(P.lift())

    def _canonical_lift_rec(
        self, P: MPolynomial_element
    ) -> MPolynomial_element:
        vars, degs, subs = (
            self.Kbase_variables,
            self.Kbase_degrees,
            self.Kbase_subs,
        )
        degree = lambda mon, var: (
            mon.degree(var) if len(vars) > 1 else mon.degree()
        )
        output: MPolynomial_element = self.Kbase(0)
        for mon in P.monomials():
            coeff = P.monomial_coefficient(mon)
            newmon: MPolynomial_element = coeff * mon
            for var, deg, sub in reversed(list(zip(vars, degs, subs))):
                mondeg = degree(mon, var)
                if mondeg < deg:
                    continue
                newmon = prod(
                    vari ** degree(mon, vari) for vari in vars if vari != var
                )
                newmon *= (
                    var ** (mondeg % deg) * sub ** (mondeg // deg) * coeff
                )
                newmon = self._canonical_lift_rec(newmon)
                break
            output += newmon
        return output


# ===================================================================================== #
# === MULTIPROCESSING SCRIPT FOR CONVERTING TO RATIONAL HORIZONTALLY BLOCKED MATRIX === #
# ===================================================================================== #

mp_M: Matrix_generic_dense | None = None  # global matrix
mp_mons: list[MPolynomial_element] | None = None  # global list of monomials
mp_vars: list[MPolynomial_element] | None = None  # global list of variables
mp_degs: list[int] | None = None  # global list of degrees
mp_subs: list[MPolynomial_element] | None = None  # global list of subs
mp_Kbase: MPolynomialRing_base | None = None  # global polynomial ring


def mp_degree(mon: MPolynomial_element, var: MPolynomial_element) -> int:
    assert mp_vars is not None
    return mon.degree(var) if len(mp_vars) > 1 else mon.degree()


def mp_canonical_lift(
    P: PolynomialQuotientRingElement,
) -> MPolynomial_element:
    return mp_canonical_lift_rec(P.lift())


def mp_canonical_lift_rec(P: MPolynomial_element) -> MPolynomial_element:
    global mp_vars, mp_degs, mp_subs, mp_Kbase
    assert mp_vars is not None
    assert mp_degs is not None
    assert mp_subs is not None
    assert mp_Kbase is not None
    output: MPolynomial_element = mp_Kbase(0)
    for mon in P.monomials():
        coeff = P.monomial_coefficient(mon)
        newmon: MPolynomial_element = coeff * mon
        for var, deg, sub in reversed(list(zip(mp_vars, mp_degs, mp_subs))):
            mondeg = mp_degree(mon, var)
            if mondeg < deg:
                continue
            newmon = prod(
                vari ** mp_degree(mon, vari)
                for vari in mp_vars
                if vari != var
            )
            newmon *= var ** (mondeg % deg) * sub ** (mondeg // deg) * coeff
            newmon = mp_canonical_lift_rec(newmon)
            break
        output += newmon
    return output


def mp_convert_entry_to_horizontal_block_worker(
    args: int,
) -> tuple[int, list[list[Rat]]]:
    assert mp_mons is not None
    assert mp_M is not None
    output: list[list[Rat]] = []
    i = args
    for j in range(mp_M.ncols()):
        entry = mp_canonical_lift(mp_M[i, j])
        curr_vec = vector(
            [entry.monomial_coefficient(mon) for mon in mp_mons]
        )
        output.append(curr_vec)
    return i, output


def mp_convert_to_rational_horizontally_blocked_matrix(
    M: Matrix_generic_dense,
    mons: list[MPolynomial_element],  # these are in K_base
    vars: list[MPolynomial_element],
    degs: list[int],
    subs: list[MPolynomial_element],
    Kbase: MPolynomialRing_base,
) -> Matrix_rational_dense:
    degree = len(mons)
    blocked_matrix = zero_matrix(QQ, M.nrows(), M.ncols() * degree)
    global mp_M, mp_mons, mp_vars, mp_degs, mp_subs, mp_Kbase
    mp_M, mp_mons, mp_vars, mp_degs, mp_subs, mp_Kbase = (
        M,
        mons,
        vars,
        degs,
        subs,
        Kbase,
    )
    with Pool(NCORES) as pool:
        tasks = list(range(M.nrows()))
        for i, val in tqdm(
            pool.imap_unordered(
                mp_convert_entry_to_horizontal_block_worker,
                tasks,  # type: ignore #iterprod(range(M.nrows()), range(M.ncols())),
                chunksize=max(M.nrows() // (10 * NCORES), 1),
            ),
            desc="Assembling blocked matrix",
            total=M.nrows(),
        ):
            for j in range(M.ncols()):
                for k in range(degree):
                    blocked_matrix[i, j + M.ncols() * k] = val[j][k]
    return blocked_matrix
