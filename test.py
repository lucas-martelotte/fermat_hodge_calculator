from src.utils.sage_imports import parallel
from time import sleep


@parallel(ncpus=2)
def f(n: int) -> int:
    for i in range(10000000):
        pass
    return n**2


def g(n: int) -> int:
    for i in range(10000000):
        pass
    return n**2


N = 1000000

print("a")
y = f([(i,) for i in range(N)])  # type: ignore
# print(list(y))

print("b")
y = [g(i) for i in range(N)]
# print(y)

print("c")
