from src.utils.sage_imports import PolynomialRing, QQ, cyclotomic_polynomial
from .formal_number_field import FormalNumberField


class FormalCyclotomicField(FormalNumberField):
    def __init__(self, order: int) -> None:
        assert order >= 2
        Kbase = PolynomialRing(QQ, [f"zeta{order}base"])
        (zeta,) = Kbase.gens()
        super().__init__(Kbase, [cyclotomic_polynomial(order, zeta)])
