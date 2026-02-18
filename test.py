# type: ignore
from sage.all import *

K = CyclotomicField(11)
R = PolynomialRing(K, ["x"])
(x,) = R.gens()
z11 = exp(2 * pi * I / 11).n()

P = 1024 * x**5 - 704 * x**3 - 176 * x**2 + 44 * x + 11  # some polynomial
roots = P.roots()
for root, mult in roots:
    print()
    print(root)
    print(sage_eval(str(root), locals={"zeta11": z11}))
    print()
