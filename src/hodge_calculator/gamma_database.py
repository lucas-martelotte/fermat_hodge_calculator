# ====================== #
# === WELCOME TO THE === #
# === GAMMA DATABASE === #
# ====================== #

# ========================================================= #
# Here, values of Gamma products correspondings to periods  #
# of Fermat varieties of different dimensions and degrees   #
# are stored. These values are necessary in order to        #
# instantiate the class HodgeCalculator.                    #
# ========================================================= #

from typing import Callable
from dataclasses import dataclass
from src.utils.singleton_metaclass import SingletonMetaclass
from itertools import product as iterprod
from src.utils.auxiliary import coprimes
from src.utils.sage_imports import (
    PolynomialQuotientRingElement,
    PolynomialRing,
    sage_eval,
    gamma,
    prod,
    QQ,
    Rat,
    sqrt,
    exp,
    pi,
    I,
)
from ..formal_algebra import (
    FormalCyclotomicField,
    FormalNumberField,
    PolynomialRingOverFormalNumberField,
)


@dataclass
class GammaStructure:
    K_formal: FormalNumberField
    gamma_values: dict[tuple[int, ...], PolynomialQuotientRingElement]
    embedding: Callable[[str], complex]


class GammaDatabase:  # TODO Fix issue with singleton metaclass
    def __init__(self, check_gamma_values: bool = False) -> None:
        super().__init__()
        self._database: dict[tuple[int, int], GammaStructure] = {}
        self._check_gamma_values = check_gamma_values

    def gamma_product(
        self, gamma_tuple: tuple[int, ...], denominator: int
    ) -> complex:
        output = 1 / (2 * pi * I) ** (QQ(len(gamma_tuple)) / 2)
        output *= prod(gamma(ti / denominator) for ti in gamma_tuple)
        return complex(output.n())

    def add(
        self, key: tuple[int, int], gamma_structure: GammaStructure
    ) -> None:
        """
        A = lambda b, c: self.weighted_sum_of_index(degree, b, c)
        k = dimension // 2 + 1
        multi_idxs = list(
            iterprod(*(list(range(degree - 1)) for _ in range(dimension + 2)))
        )
        J_hodge = [
            i for i in range(len(multi_idxs)) if A(multi_idxs[i], 1) == k
        ]
        J_hodge_alg = [
            i
            for i in J_hodge
            if all(A(multi_idxs[i], c) == k for c in coprimes(degree))
        ]
        # Asserting if all gamma values are being considered
        gamma_values = gamma_structure.gamma_values
        embedding = gamma_structure.embedding
        gamma_values_keys_assert = {
            tuple(sorted([fi + 1 for fi in multi_idxs[j]]))
            for j in J_hodge_alg
        }
        print(gamma_values_keys_assert)
        assert gamma_values_keys_assert == set(gamma_values.keys())
        """
        # Asserting if they all have the correct value
        dimension, degree = key
        if self._check_gamma_values:
            for gamma_tuple, gamma_value in gamma_values.items():
                gamma_value_1 = embedding(str(gamma_value))
                gamma_value_2 = self.gamma_product(gamma_tuple, degree)
                assert len(gamma_tuple) == dimension + 2
                assert abs(gamma_value_2 - gamma_value_1) < 0.00001
        self._database[key] = gamma_structure

    def get(self, key: tuple[int, int]) -> GammaStructure:
        return self._database[key]


GAMMA_DATABASE = GammaDatabase(check_gamma_values=False)

# === degree 2 surface === #
K_formal = FormalCyclotomicField(4)
gamma_values = {(1, 1, 1, 1): (-1 / K_formal.K(4))}
embedding_locals = {"zeta4": exp(2 * pi * I / 4)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 2), GammaStructure(K_formal, gamma_values, embedding))

# === degree 3 surface === #
K_formal = FormalCyclotomicField(3)
zeta3 = K_formal.K.gen()
K_formal.locals["zeta6"] = -(zeta3**2)
gamma_values = {(1, 1, 2, 2): (-1 / K_formal.K(3))}
embedding_locals = {"zeta3": exp(2 * pi * I / 3)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 3), GammaStructure(K_formal, gamma_values, embedding))

# === degree 4 surface === #
K_formal = FormalCyclotomicField(8)
zeta8 = K_formal.K(K_formal.from_str("zeta8"))
zeta4 = zeta8**2
root2 = zeta8 - zeta8**3
K_formal.locals["zeta4"] = zeta4
K_formal.locals["root2"] = root2
gamma_values = {
    # Lines (on degree 2)
    (2, 2, 2, 2): -1 / K_formal.K(4),
    # Lines ("true lines")
    (1, 1, 3, 3): -1 / K_formal.K(2),
    (1, 2, 2, 3): -1 / (2 * root2),
}
embedding_locals = {"zeta8": exp(2 * pi * I / 8)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 4), GammaStructure(K_formal, gamma_values, embedding))

# === degree 5 surface === #
K_formal = FormalCyclotomicField(5)
zeta5 = K_formal.K.gen()
root5 = -(2 * zeta5**3 + 2 * zeta5**2 + 1)
K_formal.locals["root5"] = root5
K_formal.locals["zeta10"] = -(zeta5**3)
gamma_values = {
    # Lines ("true lines")
    (1, 2, 3, 4): -1 / root5,
    (2, 2, 3, 3): -2 / (5 + root5),
    (1, 1, 4, 4): -2 / (5 - root5),
}
embedding_locals = {"zeta5": exp(2 * pi * I / 5)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 5), GammaStructure(K_formal, gamma_values, embedding))

# === degree 6 surface === #
Kbase = PolynomialRing(QQ, ["zeta12base", "root3of2base"])
zeta12base, root3of2base = Kbase.gens()
Kbase_equations = [
    zeta12base**4 - zeta12base**2 + 1,
    root3of2base**3 - 2,
]
K_formal = FormalNumberField(Kbase, Kbase_equations)
zeta12 = K_formal.K(K_formal.from_str("zeta12"))
root3of2 = K_formal.K(K_formal.from_str("root3of2"))
root3 = 2 * zeta12 - zeta12**3
zeta6 = zeta12**2
K_formal.locals["root3"] = root3
K_formal.locals["zeta6"] = zeta6
gamma_values = {
    # Lines (on degree 2)
    (3, 3, 3, 3): -1 / K_formal.K(4),
    # Lines (on degree 3)
    (2, 2, 4, 4): -1 / K_formal.K(3),
    # Lines ("true lines")
    (1, 1, 5, 5): -1 / K_formal.K(1),
    (1, 2, 4, 5): -1 / root3,
    (1, 3, 3, 5): -1 / K_formal.K(2),
    (2, 3, 3, 4): -1 / (2 * root3),
    # Aoki-shioda (24 after permutation)
    (1, 3, 4, 4): -1 / (root3 * root3of2),
    (2, 2, 3, 5): -1 / (root3 * root3of2**2),
}
embedding_locals = {
    "zeta12": exp(2 * pi * I / 12),
    "root3of2": 2 ** (QQ(1) / 3),
}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 6), GammaStructure(K_formal, gamma_values, embedding))

# === degree 7 surface === #
K_formal = FormalCyclotomicField(7)
zeta7 = K_formal.K(K_formal.from_str("zeta7"))
K_formal.locals["zeta14"] = -(zeta7**4)
# sin_a7b7 means sin(a*pi/7) * sin(b*pi/7)
sin_1717 = (2 - zeta7**6 - zeta7) / 4
sin_2727 = (2 - zeta7**5 - zeta7**2) / 4
sin_3737 = (2 - zeta7**4 - zeta7**3) / 4
sin_1727 = (1 + zeta7**3 - zeta7**2 - zeta7) * zeta7**2 / 4
sin_1737 = (1 + zeta7**5 - zeta7**4 - zeta7) * zeta7 / 4
sin_2737 = (1 + zeta7**5 - zeta7**3 - zeta7**2) * zeta7 / 4
gamma_values = {
    # Lines ("true lines")
    (1, 1, 6, 6): -1 / (4 * sin_1717),
    (1, 2, 5, 6): -1 / (4 * sin_1727),
    (1, 3, 4, 6): -1 / (4 * sin_1737),
    (2, 2, 5, 5): -1 / (4 * sin_2727),
    (2, 3, 4, 5): -1 / (4 * sin_2737),
    (3, 3, 4, 4): -1 / (4 * sin_3737),
}
embedding_locals = {"zeta7": exp(2 * pi * I / 7)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 7), GammaStructure(K_formal, gamma_values, embedding))

# === degree 8 surface === #
Kbase = PolynomialRing(QQ, ["zeta16base", "root4of2base"])
zeta16base, root4of2base = Kbase.gens()
root2base = zeta16base**2 - zeta16base**6
Kbase_equations = [
    zeta16base**8 + 1,
    root4of2base**2 - root2base,
]
K_formal = FormalNumberField(Kbase, Kbase_equations)
zeta16 = K_formal.K(K_formal.from_str("zeta16"))
root4of2 = K_formal.K(K_formal.from_str("root4of2"))
zeta8 = zeta16**2
root2 = zeta8 - zeta8**3
# sin_a8 means sin(a*pi/8)
sin_18 = (1 - zeta8) * zeta16**3 / 2
sin_28 = 1 / root2
sin_38 = (zeta16 - zeta16**7) / 2
sin_48 = K_formal.K(1)
K_formal.locals["zeta8"] = zeta8
K_formal.locals["root2"] = root2
gamma_values = {
    # Lines (on degree 2)
    (4, 4, 4, 4): -1 / (4 * sin_48**2),
    # Lines (on degree 4)
    (2, 2, 6, 6): -1 / (4 * sin_28**2),
    (2, 4, 4, 6): -1 / (4 * sin_28 * sin_48),
    # Lines ("true lines")
    (1, 1, 7, 7): -1 / (4 * sin_18**2),
    (3, 3, 5, 5): -1 / (4 * sin_38**2),
    (1, 2, 6, 7): -1 / (4 * sin_18 * sin_28),
    (1, 3, 5, 7): -1 / (4 * sin_18 * sin_38),
    (1, 4, 4, 7): -1 / (4 * sin_18 * sin_48),
    (2, 3, 5, 6): -1 / (4 * sin_28 * sin_38),
    (3, 4, 4, 5): -1 / (4 * sin_38 * sin_48),
    # Aoki-shioda (48 after permutation)
    (2, 3, 4, 7): -1 / (root4of2 * 2),
    (1, 4, 5, 6): -root4of2 / K_formal.K(2),
}
embedding_locals = {
    "zeta16": exp(2 * pi * I / 16),
    "root4of2": 2 ** (QQ(1) / 4),
}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 8), GammaStructure(K_formal, gamma_values, embedding))

# === degree 9 surface === #
Kbase = PolynomialRing(QQ, ["zeta9base", "root3of3base"])
zeta9base, root3of3base = Kbase.gens()
Kbase_equations = [
    zeta9base**6 + zeta9base**3 + 1,
    root3of3base**3 - 3,
]
K_formal = FormalNumberField(Kbase, Kbase_equations)
zeta9 = K_formal.K(K_formal.from_str("zeta9"))
K_formal.locals["zeta18"] = -(zeta9**5)
root3of3 = K_formal.K(K_formal.from_str("root3of3"))
sin_1919 = (2 + zeta9**5 + zeta9**2 - zeta9) / 4
sin_1929 = (1 + zeta9**3 - zeta9 - zeta9**2) * zeta9**3 / 4
sin_1939 = (1 + zeta9**7 - zeta9 - zeta9**6) * zeta9 / 4
sin_1949 = (1 + zeta9**2 - zeta9 - zeta9**4) / 4
sin_2929 = (2 + zeta9 - zeta9**2 + zeta9**4) / 4
sin_2939 = (1 + zeta9**5 - zeta9**2 - zeta9**3) * zeta9**2 / 4
sin_2949 = (1 + zeta9 - zeta9**5 - zeta9**2) / 4
sin_3939 = K_formal.K(3) / 4
sin_3949 = (1 + zeta9**7 - zeta9**4 - zeta9**3) * zeta9 / 4
sin_4949 = (2 - zeta9**4 - zeta9**5) / 4
gamma_values = {
    # Lines (on degree 3)
    (3, 3, 6, 6): -1 / (4 * sin_3939),
    # Lines ("true lines")
    (1, 1, 8, 8): -1 / (4 * sin_1919),
    (2, 2, 7, 7): -1 / (4 * sin_2929),
    (4, 4, 5, 5): -1 / (4 * sin_4949),
    (1, 2, 7, 8): -1 / (4 * sin_1929),
    (1, 3, 6, 8): -1 / (4 * sin_1939),
    (1, 4, 5, 8): -1 / (4 * sin_1949),
    (2, 3, 6, 7): -1 / (4 * sin_2939),
    (2, 4, 5, 7): -1 / (4 * sin_2949),
    (3, 4, 5, 6): -1 / (4 * sin_3949),
    # Aoki-shioda (48 after permutation)
    (1, 4, 6, 7): -1 / root3of3,
    (2, 3, 5, 8): -1 / root3of3**2,
}
embedding_locals = {
    "zeta9": exp(2 * pi * I / 9),
    "root3of3": 3 ** (QQ(1) / 3),
}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 9), GammaStructure(K_formal, gamma_values, embedding))

# === degree 10 surface === #
Kbase = PolynomialRing(QQ, ["zeta20base", "root5of2base"])
zeta20base, root5of2base = Kbase.gens()
Kbase_equations = [
    zeta20base**8 - zeta20base**6 + zeta20base**4 - zeta20base**2 + 1,
    root5of2base**5 - 2,
]
K_formal = FormalNumberField(Kbase, Kbase_equations)
zeta20 = K_formal.from_str("zeta20")
zeta10 = zeta20**2
root5of2 = K_formal.from_str("root5of2")
root5 = 1 + 2 * zeta10**2 - 2 * zeta10**3
sin_1_10 = (root5 - 1) / 4
sin_2_10 = (1 - zeta20**4) * zeta20**3 / 2
sin_3_10 = (root5 + 1) / 4
sin_4_10 = (2 * zeta20 - zeta20**3 + zeta20**5 - zeta20**7) / 2
sin_5_10 = K_formal.K(1)
K_formal.locals["zeta10"] = zeta10
K_formal.locals["root5"] = root5
gamma_values = {
    # Lines (on degree 2)
    (5, 5, 5, 5): -1 / (4 * sin_5_10**2),
    # Lines (on degree 5)
    (2, 2, 8, 8): -1 / (4 * sin_2_10**2),
    (4, 4, 6, 6): -1 / (4 * sin_4_10**2),
    (2, 4, 6, 8): -1 / (4 * sin_2_10 * sin_4_10),
    # Lines ("true lines")
    (1, 1, 9, 9): -1 / (4 * sin_1_10**2),
    (3, 3, 7, 7): -1 / (4 * sin_3_10**2),
    (1, 2, 8, 9): -1 / (4 * sin_1_10 * sin_2_10),
    (1, 3, 7, 9): -1 / (4 * sin_1_10 * sin_3_10),
    (1, 4, 6, 9): -1 / (4 * sin_1_10 * sin_4_10),
    (1, 5, 5, 9): -1 / (4 * sin_1_10 * sin_5_10),
    (2, 3, 7, 8): -1 / (4 * sin_2_10 * sin_3_10),
    (2, 5, 5, 8): -1 / (4 * sin_2_10 * sin_5_10),
    (3, 4, 6, 7): -1 / (4 * sin_3_10 * sin_4_10),
    (3, 5, 5, 7): -1 / (4 * sin_3_10 * sin_5_10),
    (4, 5, 5, 6): -1 / (4 * sin_4_10 * sin_5_10),
    # Aoki-shioda (144 after permutation)
    (2, 2, 7, 9): -1 / (2 * root5of2 * sin_2_10),
    (1, 5, 6, 8): -1 / (2 * root5of2 * sin_2_10),
    (2, 4, 5, 9): -1 / (2 * root5of2**4 * sin_2_10),
    (2, 5, 6, 7): -1 / (2 * root5of2**2 * sin_4_10),
    (3, 4, 4, 9): -1 / (2 * root5of2**2 * sin_4_10),
    (3, 4, 5, 8): -1 / (2 * root5of2**3 * sin_4_10),
    (1, 6, 6, 7): -1 / (root5of2**3 * sin_4_10),
    (1, 3, 8, 8): -1 / (root5of2**4 * sin_2_10),
}
embedding_locals = {
    "zeta20": exp(2 * pi * I / 20),
    "root5of2": 2 ** (QQ(1) / 5),
}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 10), GammaStructure(K_formal, gamma_values, embedding))

# === degree 11 surface === #
K_formal = FormalCyclotomicField(11)
zeta11 = K_formal.K.gen()
K_formal.locals["zeta22"] = -(zeta11**6)
# sin_a_11_b_11 means sin(a*pi/11)*sin(b*pi/11)
sin_1_11_1_11 = (2 - zeta11**1 - zeta11**10) / 4
sin_2_11_2_11 = (2 - zeta11**2 - zeta11**9) / 4
sin_3_11_3_11 = (2 - zeta11**3 - zeta11**8) / 4
sin_4_11_4_11 = (2 - zeta11**4 - zeta11**7) / 4
sin_5_11_5_11 = (2 - zeta11**5 - zeta11**6) / 4
sin_1_11_2_11 = (1 + zeta11**3 - zeta11**2 - zeta11) * zeta11**4 / 4
sin_1_11_3_11 = (1 + zeta11**9 - zeta11**8 - zeta11) * zeta11 / 4
sin_1_11_4_11 = (1 + zeta11**5 - zeta11**4 - zeta11) * zeta11**3 / 4
sin_1_11_5_11 = (1 + zeta11**7 - zeta11**6 - zeta11) * zeta11**2 / 4
sin_2_11_3_11 = (1 + zeta11**5 - zeta11**3 - zeta11**2) * zeta11**3 / 4
sin_2_11_4_11 = (1 + zeta11**9 - zeta11**7 - zeta11**2) * zeta11 / 4
sin_2_11_5_11 = (1 + zeta11**7 - zeta11**5 - zeta11**2) * zeta11**2 / 4
sin_3_11_4_11 = (1 + zeta11**7 - zeta11**4 - zeta11**3) * zeta11**2 / 4
sin_3_11_5_11 = (1 + zeta11**9 - zeta11**6 - zeta11**3) * zeta11 / 4
sin_4_11_5_11 = (1 + zeta11**9 - zeta11**5 - zeta11**4) * zeta11 / 4
gamma_values = {
    # Lines ("true lines")
    (1, 1, 10, 10): -1 / (4 * sin_1_11_1_11),
    (2, 2, 9, 9): -1 / (4 * sin_2_11_2_11),
    (3, 3, 8, 8): -1 / (4 * sin_3_11_3_11),
    (4, 4, 7, 7): -1 / (4 * sin_4_11_4_11),
    (5, 5, 6, 6): -1 / (4 * sin_5_11_5_11),
    (1, 2, 9, 10): -1 / (4 * sin_1_11_2_11),
    (1, 3, 8, 10): -1 / (4 * sin_1_11_3_11),
    (1, 4, 7, 10): -1 / (4 * sin_1_11_4_11),
    (1, 5, 6, 10): -1 / (4 * sin_1_11_5_11),
    (2, 3, 8, 9): -1 / (4 * sin_2_11_3_11),
    (2, 4, 7, 9): -1 / (4 * sin_2_11_4_11),
    (2, 5, 6, 9): -1 / (4 * sin_2_11_5_11),
    (3, 4, 7, 8): -1 / (4 * sin_3_11_4_11),
    (3, 5, 6, 8): -1 / (4 * sin_3_11_5_11),
    (4, 5, 6, 7): -1 / (4 * sin_4_11_5_11),
}
embedding_locals = {"zeta11": exp(2 * pi * I / 11)}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 11), GammaStructure(K_formal, gamma_values, embedding))

# === degree 12 surface === #
Kbase = PolynomialRing(
    QQ, ["zeta24base", "root12of2base", "rr3m1base", "root8of3base"]
)
zeta24base, root12of2base, rr3m1base, root8of3base = Kbase.gens()
zeta12base = zeta24base**2
root3base = 2 * zeta12base - zeta12base**3
root2base = zeta24base + zeta24base**3 - zeta24base**5
Kbase_equations = [
    zeta24base**8 - zeta24base**4 + 1,
    root12of2base**6 - root2base,
    rr3m1base**2 - root3base + 1,
    root8of3base**4 - root3base,
]
K_formal = FormalNumberField(Kbase, Kbase_equations)
zeta24 = K_formal.from_str("zeta24")
root12of2 = K_formal.from_str("root12of2")
rr3m1 = K_formal.from_str("rr3m1")  # sqrt(sqrt(3)-1)
root8of3 = K_formal.from_str("root8of3")
zeta12 = zeta24**2
root2, root3 = root12of2**6, 2 * zeta12 - zeta12**3
root6, zeta6, zeta4 = root2 * root3, zeta12**2, zeta12**3
zeta24 = (zeta12 + zeta12**2 - zeta12**3) / root2
root4of3, root6of2, root4of2 = root8of3**2, root12of2**2, root12of2**3
root3of2 = root12of2**4
root3of4 = root12of2**8
root4of6 = root4of2 * root4of3
rr3p1 = (zeta24 + zeta24**3 - zeta24**7) * rr3m1  # sqrt(sqrt(3)+1)
K_formal.locals["zeta12"] = zeta12
K_formal.locals["rr3p1"] = rr3p1
K_formal.locals["root2"] = root2
K_formal.locals["root3"] = root3
gamma_values = {
    # Lines (on degree 2)
    (6, 6, 6, 6): -1 / K_formal.K(4),
    # Lines (on degree 3)
    (4, 4, 8, 8): -1 / K_formal.K(3),
    # Lines (on degree 4)
    (3, 3, 9, 9): -1 / K_formal.K(2),
    (3, 6, 6, 9): -1 / (2 * root2),
    # Lines (on degree 6)
    (2, 2, 10, 10): -K_formal.K(1),
    (2, 4, 8, 10): -1 / root3,
    (2, 6, 6, 10): -1 / K_formal.K(2),
    (4, 6, 6, 8): -1 / (2 * root3),
    # Lines ("true lines")
    (1, 1, 11, 11): -1 / (K_formal.K(2) - root3),
    (1, 2, 10, 11): -root2 / (root3 - 1),
    (1, 3, 9, 11): -1 / (root3 - 1),
    (1, 4, 8, 11): -root2 / (root3 * (root3 - 1)),
    (1, 5, 7, 11): -K_formal.K(1),
    (1, 6, 6, 11): -1 / (root2 * (root3 - 1)),
    (2, 3, 9, 10): -1 / root2,
    (2, 5, 7, 10): -root2 / (root3 + 1),
    (3, 4, 8, 9): -1 / root6,
    (3, 5, 7, 9): -1 / (root3 + 1),
    (4, 5, 7, 8): -root2 / (root3 * (1 + root3)),
    (5, 5, 7, 7): -1 / (root3 + 2),
    (5, 6, 6, 7): -1 / (root2 * (root3 + 1)),
    # Aoki-shioda (on degree 6)
    (2, 6, 8, 8): -1 / (root3 * root3of2),
    (4, 4, 6, 10): -1 / (root3 * root3of2**2),
    # Aoki-shioda ("true aoki-shioda")
    (1, 4, 9, 10): -root12of2 / (root8of3 * rr3m1),
    (1, 5, 9, 9): -root4of3 / root2,
    (1, 6, 7, 10): -1 / root6of2,
    (1, 6, 8, 9): -(root4of2 * rr3p1) / (2 * root8of3),
    (1, 7, 8, 8): -root2 / root3,
    (2, 3, 8, 11): -1 / (root12of2 * root8of3 * root4of3 * rr3m1),
    (2, 5, 6, 11): -1 / (root3of2**2 * root6of2),
    (2, 5, 8, 9): -(root4of2 * root3of2**2) / (2 * root8of3 * rr3p1),
    (3, 3, 7, 11): -1 / (root2 * root4of3),
    (3, 4, 6, 11): -rr3p1 / (2 * root4of6 * root8of3),
    (3, 4, 7, 10): -(root12of2 * root8of3) / (root3 * rr3p1),
    (3, 6, 7, 8): -1 / (root4of6 * root8of3 * rr3p1),
    (4, 4, 5, 11): -1 / root6,
    (4, 5, 6, 9): -root4of2 / (2 * root8of3 * rr3p1),
}
embedding_locals = {
    "zeta24": exp(2 * pi * I / 24),
    "root12of2": 2 ** (QQ(1) / 12),
    "rr3m1": sqrt(sqrt(3) - QQ(1)),
    "root8of3": 3 ** (QQ(1) / 8),
}
embedding: Callable[[str], complex] = lambda expr: complex(
    sage_eval(expr, locals=embedding_locals).n()
)
GAMMA_DATABASE.add((2, 12), GammaStructure(K_formal, gamma_values, embedding))
