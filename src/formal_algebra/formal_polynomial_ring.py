from src.utils.sage_imports import (
    QQ,
    PolynomialRing,
    MPolynomial_element,
    MPolynomialIdeal,
    sage_eval,
)
from .formal_number_field import FormalNumberField


class PolynomialRingOverFormalNumberField:
    def __init__(self, formal_field: FormalNumberField, varnames: list[str]) -> None:
        self.formal_field = formal_field
        self.R = PolynomialRing(formal_field.K, varnames)
        self.R_variables = self.R.gens()
        self.Rbase = PolynomialRing(
            QQ,
            [str(x) for x in self.formal_field.Kbase_variables]
            + [v + "base" for v in varnames],
        )
        K_vars, Kbase = self.formal_field.K_variables, self.formal_field.Kbase
        self.Rbase_variables, R, Rbase = self.Rbase.gens(), self.R, self.Rbase
        R_vars, Rbase_vars = self.R_variables, self.Rbase_variables
        self.R_to_Rbase_locals = {str(x)[:-4]: x for x in Rbase_vars}
        n_comps = self.formal_field.number_of_components
        self.Kbase_to_Rbase = Kbase.hom(list(Rbase_vars)[:n_comps], Rbase)
        Kbase_ideal_gens = self.formal_field.Kbase_ideal.gens()
        self.Rbase_equations = [self.Kbase_to_Rbase(x) for x in Kbase_ideal_gens]
        self.Rbase_ideal = self.Rbase.ideal(self.Rbase_equations)
        self.Rbase_to_R = self.Rbase.hom(list(K_vars) + list(R_vars), R)

    def lift_polynomial(self, poly: MPolynomial_element) -> MPolynomial_element:
        """
        Lifts a polynomial in R to a polynomial in Rbase. This is of
        course not well-defined: we pick one element in the preimage.
        """
        return sage_eval(str(poly), locals=self.R_to_Rbase_locals)

    def lift_ideal(self, I: MPolynomialIdeal) -> MPolynomialIdeal:
        """
        Lifts the ideal I of R to the corresponding ideal of Rbase
        (which means taking the pullback by the projection map
        induced by Kbase -> K).
        """
        lifted_gens = sorted([self.lift_polynomial(x) for x in I.gens()])
        return self.Rbase.ideal(lifted_gens) + self.Rbase_ideal

    def project_ideal(self, I: MPolynomialIdeal) -> MPolynomialIdeal:
        """
        Projects the ideal I of Rbase to R (which means taking the
        image of I by the projection map).
        """
        projected_gens = sorted([self.Rbase_to_R(x) for x in I.gens()])
        return self.R.ideal([p for p in projected_gens if p != 0])
