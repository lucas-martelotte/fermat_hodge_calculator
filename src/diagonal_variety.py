from .formal_algebra import (
    PolynomialRingOverFormalNumberField,
    FormalNumberField,
)
from .utils.sage_imports import MPolynomialRing_base, prod
from .utils.auxiliary import lcm


class EvenDimensionalDiagonalVariety:
    def __init__(
        self,
        formal_polynomial_ring: PolynomialRingOverFormalNumberField,
        exps: tuple[int, ...],
    ) -> None:
        assert len(exps) % 2 == 0 and len(exps) > 0
        self.exps = exps
        self.degree = lcm(*exps)
        self.weights = tuple([self.degree // e for e in exps])
        self.dimension = len(exps) - 2
        self.formal_polynomial_ring = formal_polynomial_ring
        xs, ms = self.polynomial_variables, exps
        self.defining_polynomial = sum(xi**ei for xi, ei in zip(xs, exps))
        R_formal, R = formal_polynomial_ring, formal_polynomial_ring.R
        det_hess_F = prod(
            mi * (mi - 1) * xi ** (mi - 2) for xi, mi in zip(xs, ms)
        )
        jac_F = R.ideal([xi ** (mi - 1) for xi, mi in zip(xs, ms)])
        self.jacobian_of_defining_polynomial = jac_F
        self.determinant_of_hessian_matrix_of_defining_polynomial = det_hess_F
        quot_jac_F = R_formal.Rbase.quotient(R_formal.lift_ideal(jac_F))
        self.quotient_ring_by_jacobian_of_defining_polynomial = quot_jac_F

    @property
    def formal_base_field(self) -> FormalNumberField:
        return self.formal_polynomial_ring.formal_field

    @property
    def base_field(self) -> MPolynomialRing_base:
        return self.formal_base_field.K

    @property
    def polynomial_ring(self) -> MPolynomialRing_base:
        return self.formal_polynomial_ring.R

    @property
    def polynomial_variables(self) -> tuple:
        return self.polynomial_ring.gens()
