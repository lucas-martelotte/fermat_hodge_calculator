from .formal_polynomial_ring import PolynomialRingOverFormalNumberField
from src.utils.auxiliary import (
    sage_RREF,
    get_all_sums_equal_to_d,
)
from src.utils.sage_imports import (
    MPolynomial_element,
    MPolynomialIdeal,
    Matrix,
    QQ,
    vector,
    ZZ,
    prod,
    block_matrix,
)
from itertools import product as iterprod


class FormalIdealDissector:
    """
    This class implements useful functions to do computations with
    weighted homogeneous ideals over formal number fields. It assumes
    without checking that the input ideal is nonzero and weighted
    homogeneous.
    """

    def __init__(
        self,
        formal_polynomial_ring: PolynomialRingOverFormalNumberField,
        weights: tuple[int, ...],
        ideal: MPolynomialIdeal,
    ):
        self.formal_polynomial_ring = formal_polynomial_ring
        self.weights = weights
        self.ideal = ideal
        self._bases: dict[int, list[MPolynomial_element]] = {}
        # dictionary {d:[i1, ..., in]} which tells which elements in
        # the basis of degree d part are new elements
        self._new_elements: dict[int, list[int]] = {}

    def _get_monomials_of_degree(self, deg: int) -> list[MPolynomial_element]:
        ws, xs = self.weights, self.formal_polynomial_ring.R.gens()
        exps = get_all_sums_equal_to_d(list(ws), deg)
        return [prod(xi**ei for xi, ei in zip(xs, exp)) for exp in exps]

    def _degree(self, poly: MPolynomialIdeal) -> int:
        """
        Returns the weighted degree of the input polynomial.
        Assumes the polynomial is nonzero, is over the base
        field of this class and is weighted homogeneous.
        """
        mon = poly.monomials()[0]
        xs, ws = self.formal_polynomial_ring.R.gens(), self.weights
        return sum(mon.degree(xi) * wi for xi, wi in zip(xs, ws))

    def _raw_generators_of_degree(
        self, deg: int
    ) -> list[MPolynomial_element]:
        """
        Returns the generators of the ideal whose
        weighted degree is equal to deg.
        """
        return [p for p in self.ideal.gens() if self._degree(p) == deg]

    def _compute_basis_of_component(self, deg: int) -> None:
        """Returns a basis of the component of degree = deg of the ideal."""
        R_formal = self.formal_polynomial_ring
        old_elements: list[MPolynomial_element] = []
        for i in range(1, deg, 1):
            b1 = self.get_basis_of_new_elements_of_component(i)
            monslst = self._get_monomials_of_degree(deg - i)
            old_elements.extend([p1 * p2 for p1, p2 in iterprod(b1, monslst)])
        maybe_new_elements = self._raw_generators_of_degree(deg)
        all_elements = old_elements + maybe_new_elements
        if len(all_elements) == 0:  # No elements of degree = deg in the ideal
            self._bases[deg] = []
            self._new_elements[deg] = []
            return
        mons: set[MPolynomial_element] = set()
        for raw_element in all_elements:
            mons = mons.union(set(raw_element.monomials()))
        monslst = list(mons)
        raw_matrix = Matrix(
            R_formal.formal_field.K,
            len(monslst),
            len(all_elements),
            lambda i, j: all_elements[j].monomial_coefficient(monslst[i]),  # type: ignore
        )
        _, _, pivots = sage_RREF(raw_matrix)
        basis: list[MPolynomial_element] = []
        for pivot in pivots:
            curr_basis_element = self.formal_polynomial_ring.R(0)
            for i in range(raw_matrix.nrows()):
                curr_basis_element += monslst[i] * raw_matrix[i, pivot]
            basis.append(curr_basis_element)
        self._bases[deg] = basis
        self._new_elements[deg] = [
            i for i in range(len(pivots)) if pivots[i] >= len(old_elements)
        ]

    def get_basis_of_component(self, deg: int) -> list[MPolynomial_element]:
        if deg not in self._bases:
            self._compute_basis_of_component(deg)
        return self._bases[deg]

    def get_basis_of_new_elements_of_component(
        self, deg: int
    ) -> list[MPolynomial_element]:
        if deg not in self._bases:
            self._compute_basis_of_component(deg)
        basis = self._bases[deg]
        return [basis[i] for i in self._new_elements[deg]]

    def get_dimension_of_component(self, deg: int) -> int:
        if deg not in self._bases:
            self._compute_basis_of_component(deg)
        return len(self._bases[deg])

    def get_dimension_of_new_elements_of_component(self, deg: int) -> int:
        if deg not in self._bases:
            self._compute_basis_of_component(deg)
        return len(self._new_elements[deg])
