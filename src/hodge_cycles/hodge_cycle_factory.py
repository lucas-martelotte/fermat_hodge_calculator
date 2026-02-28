from abc import ABC, abstractmethod
from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.utils.sage_imports import (
    MPolynomial_element,
    PolynomialQuotientRingElement,
)


class AlgebraicPrimePrimitiveHodgeCycle(ABC):
    def __init__(
        self,
        ambient_variety: EvenDimensionalDiagonalVariety,
        defining_polynomials: list[MPolynomial_element],
    ):
        """
        Assumes without checking that the cycle is irreducible,
        contained inside the variety and that the dimension is
        equal to half the dimension of the ambient variety.
        """
        super().__init__()
        R = ambient_variety.polynomial_ring
        self.ambient_variety = ambient_variety
        self.defining_polynomials = [R(x) for x in defining_polynomials]

    @abstractmethod
    def pairing(
        self, poly: MPolynomial_element
    ) -> PolynomialQuotientRingElement:
        """
        Computes the integration pairing between itself and the
        form w_{poly} defined in (8.3) of the Hodge Theory book II.
        """
        pass


class HodgeCycleFactory(ABC):
    def __init__(self, variety: EvenDimensionalDiagonalVariety):
        super().__init__()
        self.variety = variety

    @abstractmethod
    def get_hodge_cycles(
        self,
    ) -> dict[str, AlgebraicPrimePrimitiveHodgeCycle]:
        """
        Returns a dictionary where the key is a unique identifier
        of the cycle, which is the value.
        """
        pass
