from src.utils.sage_imports import (
    NumberField_absolute,
    NumberFieldHomomorphism_im_gens,
    Matrix,
    Integer,
    Rat,
    QQ,
    prod,
    vector,
    Matrix_rational_dense,
    identity_matrix,
    polygen,
    ZZ,
    sage_eval,
    NumberField,
    PolynomialRing,
    MPolynomial_element,
)
from itertools import product as iterprod
from os.path import exists
from os import makedirs
import json
from .singleton_metaclass import Singleton
from ast import literal_eval


class CompositumManager:
    def __init__(
        self,
        compositum: NumberField_absolute,
        nfs: list[NumberField_absolute],
        embeddings: list[NumberFieldHomomorphism_im_gens],
        coefficients: list[int],
        monomial_indexes: list[tuple],
        pivots: list[int],
        compositum_to_components_matrix: Matrix_rational_dense,
    ) -> None:
        self.compositum = compositum
        self.nfs = nfs
        self.embeddings = embeddings
        self.coefficients = coefficients
        self.monomial_indexes = monomial_indexes
        self.pivots = pivots
        # This is a matrix M that converts an element x in the compositum
        # to a sum of monomial in the components, in such a way that the
        # i-th entry of M(x.vector()) is the coefficient of the i-th
        # monomial in pivots.
        self.compositum_to_components_matrix = compositum_to_components_matrix
        self.formal_ring = PolynomialRing(QQ, [str(nf.gen()) for nf in nfs])
        self.formal_variables = self.formal_ring.gens()

    def decouple(self, x: NumberField_absolute) -> dict[tuple, Rat]:
        """
        Receives an element of the compositum as input and returns a
        dictionary {key: valua} such that the key contains the exponents
        of a monomial where the i-th variable is the generator of the i-th
        number field in the compositum and the value is the coefficient
        of x with respect to that monomial.
        """
        x_coeffs = Matrix(QQ, self.compositum.degree(), 1, x.vector())
        rel_mon_coeffs = vector(
            self.compositum_to_components_matrix * x_coeffs
        )
        output: dict[tuple, Rat] = {}
        for i in range(len(self.pivots)):
            if rel_mon_coeffs[i] != 0:
                output[tuple(self.monomial_indexes[self.pivots[i]])] = (
                    rel_mon_coeffs[i]
                )
        return output

    def to_str(self, x: NumberField_absolute) -> str:
        """
        Writes x as a polynomial in the generators of each
        component in the compositum.
        """
        decoupled = self.decouple(x)
        monomials: list[MPolynomial_element] = []
        FR, xs = self.formal_ring, self.formal_variables
        for exp, coeff in decoupled.items():
            vars = prod([FR(xi**ei) for xi, ei in zip(xs, exp)])
            monomials.append(coeff * vars)
        return str(sum(monomials))

    def from_str(self, expr: str) -> NumberField_absolute:
        return self.compositum(
            sage_eval(
                expr,
                locals={
                    str(nf.gen()): e(nf.gen())
                    for nf, e in zip(self.nfs, self.embeddings)
                },
            )
        )

    def save(self, savepath: str, filename: str) -> None:
        if not exists(savepath):
            makedirs(savepath)
        nfs, ks, fs = self.nfs, self.coefficients, self.embeddings
        comp_gen_str = str(self.compositum.gen())
        data = {
            "compositum": str(self.compositum.polynomial()),
            "compositum_generator_name": comp_gen_str,
            "components": [str(nf.polynomial()) for nf in nfs],
            "components_generator_names": [str(nf.gen()) for nf in nfs],
            "embeddings": [str(f(nf.gen())) for nf, f in zip(nfs, fs)],
            "coefficients": [int(k) for k in ks],
            "monomial_indexes": [list(m) for m in self.monomial_indexes],
            "monomial_pivots": self.pivots,
            "compositum_to_components_matrix": [
                [str(q) for q in row]
                for row in self.compositum_to_components_matrix
            ],
        }
        with open(f"{savepath}/{filename}.json", "w") as f:
            json.dump(data, f)


class CompositumManagerFactory(metaclass=Singleton):

    def from_json(self, filepath: str) -> CompositumManager:
        assert exists(f"{filepath}.json")
        with open(f"{filepath}.json", "r") as f:
            data = json.load(f)
        x = polygen(ZZ, "x")
        comp_gen_name = data["compositum_generator_name"]
        comp = NumberField(
            sage_eval(data["compositum"], locals={"x": x}),
            name=comp_gen_name,
        )
        nfs: list[NumberField_absolute] = []
        for nf, gen_name in zip(
            data["components"], data["components_generator_names"]
        ):
            nfs.append(
                NumberField(sage_eval(nf, locals={"x": x}), name=gen_name)
            )
        embeddings: list[NumberFieldHomomorphism_im_gens] = []
        for nf, e in zip(nfs, data["embeddings"]):
            embeddings.append(
                nf.hom([sage_eval(e, locals={comp_gen_name: comp.gen()})])
            )
        coeffs = data["coefficients"]
        exps = [tuple(idx) for idx in data["monomial_indexes"]]
        pivots = data["monomial_pivots"]
        comp_to_components_matrix = Matrix(
            QQ,
            [
                [sage_eval(x) for x in row]
                for row in data["compositum_to_components_matrix"]
            ],
        )
        return CompositumManager(
            comp,
            nfs,
            embeddings,
            coeffs,
            exps,
            pivots,
            comp_to_components_matrix,
        )

    def from_number_fields(
        self, nfs: list[NumberField_absolute]
    ) -> CompositumManager:
        fs, ks, comp = self.compute_compositum(nfs)
        exps = list(iterprod(*[range(nf.degree()) for nf in nfs]))
        mons = [
            prod(f(nf.gen()) ** e for nf, f, e in zip(nfs, fs, e))
            for e in exps
        ]
        coeff_matrix = Matrix(QQ, [list(mon.vector()) for mon in mons]).T
        pivots = coeff_matrix.pivots()
        compositum_to_components_matrix = coeff_matrix[:, pivots].inverse()
        return CompositumManager(
            comp, nfs, fs, ks, exps, pivots, compositum_to_components_matrix
        )

    def compute_compositum(
        self, nfs: list[NumberField_absolute], varname: str = "h"
    ) -> tuple[
        list[NumberFieldHomomorphism_im_gens],
        list[Integer],
        NumberField_absolute,
    ]:
        # print(f"Starting {varname}")
        if len(nfs) == 1:
            K = NumberField(nfs[0].polynomial(), names=(varname,))
            return [nfs[0].hom([K.gen()], K)], [1], K
        k = (len(nfs) + 1) // 2
        fs1, ks1, comp1 = self.compute_compositum(nfs[:k], varname + "1")
        fs2, ks2, comp2 = self.compute_compositum(nfs[k:], varname + "2")
        comps = comp1.composite_fields(comp2, both_maps=True, names=varname)
        assert len(comps) == 1
        comp, F1, F2, k = comps[0]
        assert isinstance(k, Integer)
        fs = [F1 * f1 for f1 in fs1] + [F2 * f2 for f2 in fs2]
        ks = [k * k1 for k1 in ks1] + ks2
        # print(f"Finished {varname}")
        return fs, ks, comp
