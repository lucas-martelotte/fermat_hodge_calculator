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
        self.K.element_class.divides = lambda _1, _2: True
        self.K.defining_ideal().is_prime = lambda: True
        self.K.is_field = lambda: True
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
        self._locals = {str(x): x for x in self.K_variables}
