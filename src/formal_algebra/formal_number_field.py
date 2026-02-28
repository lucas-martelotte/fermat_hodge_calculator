from itertools import product as iterprod

from src.utils.sage_imports import (
    PolynomialQuotientRingElement,
    Vector_rational_dense,
    Matrix_rational_dense,
    MPolynomialRing_base,
    Matrix_generic_dense,
    MPolynomial_element,
    zero_matrix,
    sage_eval,
    vector,
    prod,
    QQ,
)
from src.utils.auxiliary import return_true


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
        """
        self.Kbase_ideal = Kbase.ideal(Kbase_equations)
        self.Kbase_equations = Kbase_equations
        self.Kbase = Kbase
        self.Kbase_variables = self.Kbase.gens()
        self.number_of_components = len(self.Kbase_variables)
        assert all(str(x).endswith("base") for x in self.Kbase_variables)

        varnames = [str(x)[:-4] for x in self.Kbase_variables]
        self.K = self.Kbase.quotient(self.Kbase_ideal, names=varnames)
        self.K.element_class.divides = return_true
        self.K.defining_ideal().is_prime = return_true
        self.K.is_field = return_true
        self.K_variables = self.K.gens()

        self.monomial_indexes = list(
            iterprod(*[range(eq.degree()) for eq in Kbase_equations])
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
        x_lift, mons = x.lift(), self.Kbase_monomials
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
        M_rat = zero_matrix(QQ, M.nrows(), M.ncols() * self.degree)
        for i, j in iterprod(range(M.nrows()), range(M.ncols())):
            curr_vec = self.vector(M[i, j])
            for k in range(self.degree):
                M_rat[i, j + M.ncols() * k] = curr_vec[k]
        return M_rat

    def solve_left_in_rationals(
        self, M: Matrix_generic_dense, V: Matrix_generic_dense
    ) -> Matrix_rational_dense:
        """Solves XM = V, where X is a rational matrix."""
        M_rat = self.convert_to_rational_horizontally_blocked_matrix(M)
        V_rat = self.convert_to_rational_horizontally_blocked_matrix(V)
        return M_rat.solve_left(V_rat)
