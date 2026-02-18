from ..hodge_calculator import BasicHodgeCalculator, HodgeCalculator
from ..diagonal_variety import EvenDimensionalDiagonalVariety
from ..formal_algebra import (
    PolynomialRingOverFormalNumberField,
    FormalCyclotomicField,
)
from src.utils.auxiliary import lcm
from .gamma_database import GAMMA_DATABASE


class HodgeCalculatorFactory:
    def create_basic(
        self, exps: tuple[int, ...]
    ) -> tuple[EvenDimensionalDiagonalVariety, BasicHodgeCalculator]:
        """
        Creates a BasicHodgeCalculator, together with the variety.
        The advantage is that there is no need to pre-compute gamma
        values to run this routine.
        """
        degree = lcm(*exps)
        varnames = [f"x{i}" for i in range(len(exps))]
        K_formal = FormalCyclotomicField(degree)
        R_formal = PolynomialRingOverFormalNumberField(K_formal, varnames)
        variety = EvenDimensionalDiagonalVariety(R_formal, exps)
        calculator = BasicHodgeCalculator(variety)
        return variety, calculator

    def create(
        self, exps: tuple[int, ...]
    ) -> tuple[EvenDimensionalDiagonalVariety, HodgeCalculator]:
        """
        Creates a full HodgeCalculator, together with the variety.
        The corresponding gamma values and field of definition have
        to be pre-computed.
        """
        dimension, degree = len(exps) - 2, lcm(*exps)
        gamma_structure = GAMMA_DATABASE.get((dimension, degree))
        varnames = [f"x{i}" for i in range(len(exps))]
        K_formal = gamma_structure.K_formal
        R_formal = PolynomialRingOverFormalNumberField(K_formal, varnames)
        variety = EvenDimensionalDiagonalVariety(R_formal, exps)
        gamma_values = gamma_structure.gamma_values
        calculator = HodgeCalculator(variety, gamma_values)
        return variety, calculator
