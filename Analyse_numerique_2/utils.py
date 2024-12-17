import numpy as np

def newton(a, da, u0, du0, t, x, tol=1e-6, max_iter=100):
    X = x
    for _ in range(max_iter):
        g = a(u0(X)) * t + X - x
        gp = da(u0(X)) * du0(X) * t + 1
        if np.abs(gp) < 1e-4:
            # print("Derivative is too small")
            return None
        Xnew = X - g / gp
        if np.abs(Xnew - X) < tol:
            return Xnew
        X = Xnew
    return None


def build_solution(a, da, u0, du0):
    def u(t, x):
        X = newton(a, da, u0, du0, t, x)
        if X is None:
            return None
        return u0(X)

    return u