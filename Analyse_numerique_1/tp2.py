# %%
import numpy as np
import matplotlib.pyplot as plt


def poisson1D(f, cond1, cond2, nint=20):
    a, alpha = cond1
    b, beta = cond2
    h = (b - a) / nint
    x = np.linspace(a, b, nint + 1)
    Ah = np.eye(nint) * -2 + np.eye(nint, k=-1) + np.eye(nint, k=1)
    Ah[0, 0] = -3 / 2
    Ah[0, 1] = 2
    Ah[0, 2] = -1 / 2
    Ah *= 1 / h**2
    F = f(x[:-1])
    F[0] = 0
    Bc = np.zeros(nint)
    Bc[0] = alpha / h
    Bc[-1] = -beta / h**2

    Uh = np.linalg.solve(Ah, F + Bc)
    last = Uh
    return x, np.append(Uh, beta)


def poisson1D_alt(f, cond1, cond2, nint=20):
    a, alpha = cond1
    b, beta = cond2
    h = (b - a) / nint

    # Adjust grid points to include only nint + 1 points (both boundaries)
    x = np.linspace(a, b, nint + 1)

    # Adjust matrix Ah size to match (nint + 1, nint + 1)
    Ah = np.eye(nint + 1) * -2 + np.eye(nint + 1, k=-1) + np.eye(nint + 1, k=1)

    # Neumann condition at the first row (x = a)
    Ah[0, 0] = -3 / 2
    Ah[0, 1] = 2
    Ah[0, 2] = -1 / 2

    # Dirichlet condition at the last row (x = b)
    Ah[-1, -1] = 1
    Ah[-1, -2] = 0
    Ah /= h**2

    # Define F vector and apply the source term and boundary conditions
    F = f(x)
    F[0] = 0  # Neumann condition at x = a
    F[-1] = 0  # This line is optional and can be set to zero

    # Boundary condition adjustments
    Bc = np.zeros(nint + 1)
    Bc[0] = alpha / h  # Neumann condition at the start
    Bc[-1] = beta / h**2  # Neumann condition at the start
    # Note: No need to set Bc[-1] because Dirichlet is enforced in Ah directly

    # Solve the linear system
    Uh = np.linalg.solve(Ah, F + Bc)
    return x, Uh


# print(x,U)
a, alpha = [0, -5]
b, beta = [3, 3]

sol = lambda x: np.exp(x) + (x - b) * (alpha - np.exp(a)) + beta - np.exp(b)

E = lambda X, Yh: np.max(np.abs(sol(X) - Yh))

nint = 34
x, U = poisson1D(np.exp, [0, -5], [3, 3], nint=nint)
plt.figure(figsize=(5, 15))
ax = plt.subplot(411)
plt.scatter(x, U, label="poisson1D")
X = np.linspace(0, 3, 100)
plt.plot(X, sol(X), label="sol")
x, Ualt = poisson1D_alt(np.exp, [0, -5], [3, 3], nint=nint)
Ualt -= np.abs(Ualt[-1] - 3)
plt.scatter(x, Ualt, label="poisson1D_alt")
plt.legend()

nints = np.logspace(np.log10(20), 3.5, num=10).astype(int)

ax = plt.subplot(412)
ax.semilogy(x, np.abs(sol(x) - U), label="sol - U")
plt.legend()
ax = plt.subplot(413)
ax.semilogy(x, np.abs(sol(x) - Ualt), label="sol - Ualt")
plt.legend()
ax = plt.subplot(414)

hcalc = lambda nint: b - a / nint

ax.loglog(
    [(b - a) / nint for nint in nints],
    [E(*poisson1D(np.exp, [0, -5], [3, 3], nint=nint)) for nint in nints],
    label="sol - U",
)
ax.loglog(
    [(b - a) / nint for nint in nints],
    [E(*poisson1D_alt(np.exp, [0, -5], [3, 3], nint=nint)) for nint in nints],
    label="sol - Ualt",
)
plt.tight_layout()
plt.legend()
