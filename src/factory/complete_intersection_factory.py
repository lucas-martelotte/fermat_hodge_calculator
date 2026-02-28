from .hodge_cycle_factory import (
    AlgebraicPrimePrimitiveHodgeCycle,
    HodgeCycleFactory,
)
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.utils.sage_imports import (
    MPolynomial_element,
    PolynomialQuotientRingElement,
    prod,
    factorial,
    matrix,
)
from itertools import permutations, product as iterprod
from src.utils.auxiliary import count_fixed_permutations
from tqdm import tqdm


class CompleteIntersection(AlgebraicPrimePrimitiveHodgeCycle):
    def __init__(
        self,
        ambient_variety: EvenDimensionalDiagonalVariety,
        fs: list[MPolynomial_element],
        gs: list[MPolynomial_element],
    ):
        """
        Assumes without checking that the input equations
        define a complete intersection.
        """
        super().__init__(ambient_variety, fs)
        F = self.ambient_variety.defining_polynomial
        assert sum(f * g for f, g in zip(fs, gs)) == F
        self.fs, self.gs = fs, gs

    def pairing(
        self, poly: MPolynomial_element
    ) -> PolynomialQuotientRingElement:
        X = self.ambient_variety
        R_formal = X.formal_polynomial_ring
        H = sum([[f, g] for f, g in zip(self.fs, self.gs)], [])
        det_hess_F = X.determinant_of_hessian_matrix_of_defining_polynomial
        xs, n = X.polynomial_variables, X.dimension
        nf = self.ambient_variety.formal_base_field.K
        jac_H = matrix([[fi.derivative(xj) for xj in xs] for fi in H])
        det_jac_H = jac_H.det()
        R_quot = X.quotient_ring_by_jacobian_of_defining_polynomial
        poly_lift = R_formal.lift_polynomial(poly)
        det_jac_H_lift = R_formal.lift_polynomial(det_jac_H)
        form_poly_quot = R_quot(poly_lift)
        det_jac_H_quot = R_quot(det_jac_H_lift)
        prod_quot = form_poly_quot * det_jac_H_quot
        det_hess_F_monomial = det_hess_F.monomials()[0]  # should have only 0
        det_hess_F_coeff = det_hess_F.coefficients()[0]  # should have only 0
        prod_in_R = R_formal.Rbase_to_R(prod_quot.lift())
        prod_coeff = prod_in_R.monomial_coefficient(det_hess_F_monomial)
        c = nf(prod_coeff / det_hess_F_coeff)
        det_hess_F_lift = R_formal.lift_polynomial(det_hess_F)
        det_hess_F_quot = R_quot(det_hess_F_lift)
        assert det_hess_F_quot * R_quot(c.lift()) == prod_quot
        ms = self.ambient_variety.exps
        output_num = (-1) ** (n // 2) * c * prod(mi - 1 for mi in ms)
        output_den = factorial(n // 2)
        return output_num / output_den


class CompleteIntersectionFactory(HodgeCycleFactory):
    def __init__(
        self,
        variety: EvenDimensionalDiagonalVariety,
        fs: list[MPolynomial_element],
        gs: list[MPolynomial_element],
    ):
        super().__init__(variety)
        self.fs, self.gs = fs, gs
        assert len(fs) == len(gs)
        assert (
            sum(fi * gi for fi, gi in zip(fs, gs))
            == variety.defining_polynomial
        )

    def get_hodge_cycles(
        self,
    ) -> dict[str, AlgebraicPrimePrimitiveHodgeCycle]:
        """
        In this case, the identifier for each cycle is the
        element of the automorphism group acting on the cycle C,
        where C is the cycle passed as input, corresponding to
        the key of the identity action.
        """
        X = self.variety
        xs, ms = X.polynomial_variables, X.exps
        id_to_equations: dict[
            str, tuple[list[MPolynomial_element], list[MPolynomial_element]]
        ] = {}
        zeta_exps = iterprod(*[range(ei) for ei in ms])
        perms = [
            p
            for p in permutations(range(len(xs)))
            if all(ms[i] == ms[p[i]] for i in range(len(ms)))
        ]
        for cycle_id in tqdm(
            iterprod(perms, zeta_exps),
            desc="Listing cycles up to automorphism",
            total=prod(ms) * count_fixed_permutations(ms),
        ):
            perm, zeta_exp = cycle_id
            cycle_fs = [self.action(perm, zeta_exp, fi) for fi in self.fs]
            cycle_gs = [self.action(perm, zeta_exp, gi) for gi in self.gs]
            if self.__cycle_belongs_to(id_to_equations, cycle_fs, cycle_gs):
                continue  # check if cycle was already included
            id_to_equations[str(cycle_id)] = (cycle_fs, cycle_gs)
        return {
            key: CompleteIntersection(X, val[0], val[1])
            for key, val in id_to_equations.items()
        }

    def action(
        self,
        perm: tuple[int, ...],
        zeta_exp: tuple[int, ...],
        poly: MPolynomial_element,
    ) -> MPolynomial_element:
        """
        Returns the resulting cycle after action on the complete intersection
        (self.fs, self.gs) by the element of the automorphism group of the
        variety parametrized by a permutation 'perm' of the variables and
        a list of elements (e1, ..., en), each representing the root of
        unity zeta_{d/wi}^ei multiplying the i-th coordinate (after the
        permutation).
        """
        K_formal, d = self.variety.formal_base_field, self.variety.degree
        xs, ws, zd = (
            self.variety.polynomial_variables,
            self.variety.weights,
            K_formal.from_str(f"zeta{d}"),
        )
        locals = {
            xs[i]: xs[perm[i]] * zd ** (ws[i] * zeta_exp[i])
            for i in range(len(xs))
        }
        return poly.subs(locals)

    def __cycle_belongs_to(
        self,
        id_to_equations: dict[
            str, tuple[list[MPolynomial_element], list[MPolynomial_element]]
        ],
        cycle_fs: list[MPolynomial_element],
        cycle_gs: list[MPolynomial_element],
    ) -> bool:
        """
        Returns True if the cycle defined by (cycle_fs, cycle_gs)
        already belongs to the dictionary id_to_equations, in the
        sense that another cycle in this dictionary produces the
        exact same variety. (Obs.: the input cycle is assumed
        to be on the same orbit as the cycles in the
        dictionary.)
        """
        return any(
            self.__cycles_on_the_same_orbit_are_equal(
                c[0], c[1], cycle_fs, cycle_gs
            )
            for c in id_to_equations.values()
        )

    def __cycles_on_the_same_orbit_are_equal(
        self,
        c1_fs: list[MPolynomial_element],
        c1_gs: list[MPolynomial_element],
        c2_fs: list[MPolynomial_element],
        c2_gs: list[MPolynomial_element],
    ) -> bool:
        """
        Returns True if the complete intersections (c1_fs, c1_gs) and
        (c2_fs, c2_gs) (assumed to be on the same orbit under the
        action of the automorphism group of the variety) define
        the same cycle or not.
        """
        if len(c1_fs) == 0:
            return True  # Then all ci_fs and ci_gs are empty
        f1, g1 = c1_fs[0], c1_gs[0]
        for i in range(len(c2_fs)):
            f2, g2 = c2_fs[i], c2_gs[i]
            if set(f1.monomials()) != set(f2.monomials()):
                continue  # Then f1 and f2 cannot coincide
            mon = f1.monomials()[0]
            c1 = f1.monomial_coefficient(mon)
            c2 = f2.monomial_coefficient(mon)
            if f1 != f2 * c1 * c2.inverse():
                continue  # Then they cannot coincide
            # From this point on, f1 == f2.
            if set(g1.monomials()) != set(g2.monomials()):
                continue  # Then they cannot coincide
            mon = g1.monomials()[0]
            c1 = g1.monomial_coefficient(mon)
            c2 = g2.monomial_coefficient(mon)
            if g1 != g2 * c1 * c2.inverse():
                continue  # Then they cannot coincide
            return self.__cycles_on_the_same_orbit_are_equal(
                c1_fs[1:],
                c1_gs[1:],
                c2_fs[:i] + c2_fs[i + 1 :],
                c2_gs[:i] + c2_gs[i + 1 :],
            )
        return False
