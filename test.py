# from sage.libs.pari import pari
# pari.allocatemem(15_000_000_000)

from src.diagonal_variety import EvenDimensionalDiagonalVariety
from src.formal_algebra import (
    FormalCyclotomicField,
    PolynomialRingOverFormalNumberField,
)
from src.utils.sage_imports import identity_matrix
from src.utils.auxiliary import sage_matrix_map
from src.hodge_calculator import HodgeCalculator, BasicHodgeCalculator
from src.hodge_cycles.hodge_cycle_factories import LineFactory
from src.hodge_calculator import HodgeCalculatorFactory
from src.hodge_calculator import BasicHodgeCalculator
from src.utils.sage_imports import QQ, Matrix, ZZ, block_matrix
from time import perf_counter
import json
from src.utils.sage_imports import RR, zero_matrix
from src.utils.auxiliary import get_pairings, sage_matrix_map
from src.hodge_cycles.hodge_cycle_factories import AokiShiodaType1Factory, AokiShiodaType3Factory, AokiShiodaType2BFactory

start_time = perf_counter()

hodge_calculator_factory = HodgeCalculatorFactory()

degrees_to_check = [12]

# degree 2 surface
if 2 in degrees_to_check:
    print("Analyzing degree 2 surface...")
    X, calculator = hodge_calculator_factory.create((2, 2, 2, 2))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 3 in degrees_to_check:
    print("Analyzing degree 3 surface...")
    X, calculator = hodge_calculator_factory.create((3, 3, 3, 3))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 4 in degrees_to_check:
    print("Analyzing degree 4 surface...")
    X, calculator = hodge_calculator_factory.create((4, 4, 4, 4))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 5 in degrees_to_check:
    print("Analyzing degree 5 surface...")
    X, calculator = hodge_calculator_factory.create((5, 5, 5, 5))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 6 in degrees_to_check:
    print("Analyzing degree 6 surface...")
    X, calculator = hodge_calculator_factory.create((6, 6, 6, 6))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    as1_data = calculator.get_hodge_cycle_factory_data("as1")
    as1_rank = as1_data["rank"]
    as1_coords = Matrix(ZZ, as1_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 1: {as1_rank}")

    as3_data = calculator.get_hodge_cycle_factory_data("as3")
    as3_rank = as3_data["rank"]
    as3_coords = Matrix(ZZ, as3_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 3: {as3_rank}")

    blocked = block_matrix([[line_coords, as1_coords, as3_coords]])
    blocked_rank = blocked.rank()
    print(f"Total rank of all found cycles: {blocked_rank}")
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 7 in degrees_to_check:
    print("Analyzing degree 7 surface...")
    X, calculator = hodge_calculator_factory.create((7, 7, 7, 7))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()

    print(f"Total degree of primitive hodge cycles: {total_rank}")
    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")

    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 8 in degrees_to_check:
    print("Analyzing degree 8 surface...")
    X, calculator = hodge_calculator_factory.create((8, 8, 8, 8))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")
    
    as1_data = calculator.get_hodge_cycle_factory_data("as1")
    as1_rank = as1_data["rank"]
    as1_coords = Matrix(ZZ, as1_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 1: {as1_rank}")

    as2b_data = calculator.get_hodge_cycle_factory_data("as2b")
    as2b_rank = as2b_data["rank"]
    as2b_coords = Matrix(ZZ, as2b_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 2B: {as2b_rank}")

    blocked = block_matrix([[line_coords, as1_coords, as2b_coords]])
    blocked_rank = blocked.rank()
    print(f"Total rank of all found cycles: {blocked_rank}")
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 9 in degrees_to_check:
    print("Analyzing degree 9 surface...")
    X, calculator = hodge_calculator_factory.create((9, 9, 9, 9))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, 3 * line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")
    print(f"Lines need to be multiplied by 3 to have integer entries.")
    
    as3_data = calculator.get_hodge_cycle_factory_data("as3")
    as3_rank = as3_data["rank"]
    as3_coords = Matrix(ZZ, 3 * as3_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 3: {as3_rank}")
    print(f"Aoki-shioda of type 3 needs to be multiplied by 3 to have integer entries.")
    
    blocked = block_matrix([[line_coords, as3_coords]])
    blocked_rank = blocked.rank()
    print(f"Total rank of all found cycles: {blocked_rank}")
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 10 in degrees_to_check:
    print("Analyzing degree 10 surface...")
    X, calculator = hodge_calculator_factory.create((10, 10, 10, 10))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, 2 * line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")
    print(f"Lines need to be multiplied by 2 to have integer entries.")
    
    as1_data = calculator.get_hodge_cycle_factory_data("as1")
    as1_rank = as1_data["rank"]
    as1_coords = Matrix(ZZ, 2 * as1_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 1: {as1_rank}")
    print(f"Aoki-shioda of type 1 needs to be multiplied by 2 to have integer entries.")
    
    blocked = block_matrix([[line_coords, as1_coords]])
    blocked_rank = blocked.rank()
    print(f"Total rank of all found cycles: {blocked_rank}")
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 11 in degrees_to_check:
    print("Analyzing degree 11 surface...")
    X, calculator = hodge_calculator_factory.create((11, 11, 11, 11))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")
    
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")

if 12 in degrees_to_check:
    print("Analyzing degree 12 surface...")
    X, calculator = hodge_calculator_factory.create((12, 12, 12, 12))
    K_formal = X.formal_base_field
    total_rank = calculator.get_rank_of_primitive_hodge_cycles()
    print(f"Total degree of primitive hodge cycles: {total_rank}")

    line_data = calculator.get_hodge_cycle_factory_data("lines")
    line_rank = line_data["rank"]
    line_coords = Matrix(ZZ, 72 * line_data["coordinates"])
    print(f"Total rank of lines: {line_rank}")
    print(f"Lines need to be multiplied by 72 to have integer entries.")
    #print({line_coords[i, j] for i in range(line_coords.nrows()) for j in range(line_coords.ncols())})
        
    as1_data = calculator.get_hodge_cycle_factory_data("as1")
    as1_rank = as1_data["rank"]
    as1_coords = Matrix(ZZ, 72 * as1_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 1: {as1_rank}")
    print(f"Aoki-shioda of type 1 needs to be multiplied by 72 to have integer entries.")
    #print({as1_coords[i, j] for i in range(as1_coords.nrows()) for j in range(as1_coords.ncols())})
    
    as2b_data = calculator.get_hodge_cycle_factory_data("as2b")
    as2b_rank = as2b_data["rank"]
    as2b_coords = Matrix(ZZ, 72 * as2b_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 2b: {as2b_rank}")
    print(f"Aoki-shioda of type 2B needs to be multiplied by 72 to have integer entries.")
    #print({as2b_coords[i, j] for i in range(as2b_coords.nrows()) for j in range(as2b_coords.ncols())})
    
    as3_data = calculator.get_hodge_cycle_factory_data("as3")
    as3_rank = as3_data["rank"]
    as3_coords = Matrix(ZZ, 72 * as3_data["coordinates"])
    print(f"Total rank of aoki-shioda of type 3: {as3_rank}")
    print(f"Aoki-shioda of type 3 needs to be multiplied by 72 to have integer entries.")
    #print({as3_coords[i, j] for i in range(as3_coords.nrows()) for j in range(as3_coords.ncols())})

    c_data = calculator.get_hodge_cycle_factory_data("c")
    c_rank = c_data["rank"]
    c_coords = Matrix(ZZ, 72 * c_data["coordinates"])
    print(f"Total rank of exceptional cycle of type C: {c_rank}")
    print(f"Exceptional cycle of type C needs to be multiplied by 72 to have integer entries.")
    #print({c_coords[i, j] for i in range(c_coords.nrows()) for j in range(c_coords.ncols())})

    cb_data = calculator.get_hodge_cycle_factory_data("cb")
    cb_rank = cb_data["rank"]
    cb_coords = Matrix(ZZ, 72 * cb_data["coordinates"])
    print(f"Total rank of exceptional cycle of type CB: {cb_rank}")
    print(f"Exceptional cycle of type CB needs to be multiplied by 72 to have integer entries.")
    #print({cb_coords[i, j] for i in range(cb_coords.nrows()) for j in range(cb_coords.ncols())})

    nb_data = calculator.get_hodge_cycle_factory_data("nb")
    nb_rank = nb_data["rank"]
    nb_coords = Matrix(ZZ, 72 * nb_data["coordinates"])
    print(f"Total rank of exceptional cycle of type NB: {nb_rank}")
    print(f"Exceptional cycle of type NB needs to be multiplied by 72 to have integer entries.")
    #print({nb_coords[i, j] for i in range(nb_coords.nrows()) for j in range(nb_coords.ncols())})

    blocked = block_matrix([[line_coords, as1_coords, as2b_coords, as3_coords, c_coords, cb_coords, nb_coords]])
    blocked_rank = blocked.rank()
    print(f"Total rank of all found cycles: {blocked_rank}")
    diag, _, _ = line_coords.smith_form()
    diag_entries = {diag[i, i] for i in range(min(diag.nrows(), diag.ncols()))}
    print(f"Entries in the diagonal of scaled smith normal form: {diag_entries}")
    print("=======================")



end_time = perf_counter()
print(end_time - start_time)