# %%

import numpy as np


# def J(x):
#     return x[0] ** 3 * x[1]

# def g(x):
#     return x[0] ** 3 - x[1] - 1

def dJ(x):
    return np.array([3*x[0]**2 * x[1], x[0]**3])

def dg(x):
    return np.array([3 * x[0]**2, -1])


def fdf(f, x, h):
    return (f(x + h) - f(x - h)) / 2


def DL(X):
    return np.array(
        [3 * X[0] ** 2 * (X[1] + X[2]), X[0] ** 3 - X[2], X[0] ** 3 - X[1] - 1]
    )

gradient_descent()


from scipy.optimize import root

res = root(DL, [1, 1, 1])
print(res)
x = res.x
# print(J(x) g(x))
# %%
