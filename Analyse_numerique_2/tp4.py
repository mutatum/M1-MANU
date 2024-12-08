# %%
from typing import Callable
import numpy as np
from time import sleep
from matplotlib import pyplot as plt


def finite_volume_method(
    m: int,
    ab: tuple[float, float],
    u0: Callable,
    T: float,
    flux_scheme: str = "Lax-Friedrichs",
    CFL: float = 0.5,
) -> np.ndarray:

    def f(u):  # Burgers' equation
        return 0.5 * u * u

    def fprime(u):
        return u

    def F(ul, ur, gamma):
        return 0.5 * (f(ul) + f(ur)) - 0.5 * gamma * (ur - ul)

    def MurmanRoe(ul, ur):
        delta_u = ur - ul
        delta_f = f(ur) - f(ul)
        epsilon = 1e-8
        return np.where(np.abs(delta_u) < epsilon, fprime(ul), delta_f / delta_u)

    a, b = ab
    dx = (b - a) / m
    x = np.linspace(a, b, m + 1)
    C = (x[1:] + x[:-1]) / 2
    u = u0(C)
    t = 0.0

    while t < T:
        dt = CFL * dx / np.max(np.abs(fprime(u)))

        if t + dt > T:
            dt = T - t

        u_extended = np.concatenate(([u[-1]], u, [u[0]]))
        u_L = u_extended[:-1]
        u_R = u_extended[1:]

        if flux_scheme == "Murman-Roe":
            gamma = np.abs(MurmanRoe(u_L, u_R))
        else:
            if flux_scheme == "Lax-Friedrichs":
                gamma_value = np.max(np.abs(fprime(u)))
            elif flux_scheme == "Lax-Wendroff":
                gamma_value = dx / dt
            else:
                raise ValueError(f"Invalid scheme: {flux_scheme}")
            gamma = np.full_like(u_L, gamma_value)

        F_i = F(u_L, u_R, gamma)

        flux_divergence = F_i[1:] - F_i[:-1]

        u -= dt / dx * flux_divergence

        if np.any(np.isnan(u)) or np.any(np.abs(u) > 1e6):
            raise RuntimeError("Solution diverged")

        t += dt
    return C, u


def u0_cours(x):
    res = np.where(x < 0.3, 0.0, np.where(x < 0.7, -1.0, 0.5))
    return res


def solution_exacte_cours(x, t):
    if t > 1:
        raise ValueError("t must be less than 1")
    return np.where(x < t, 1, np.where(x <= 1, (1 - x) / (1 - t), 0))


def solution_td(x, t=2):
    t_etoile = 0.8
    if t < t_etoile:
        return np.where(
            x < 0.3 - t / 2,
            0,
            np.where(
                x < 0.7 - t,
                -1,
                np.where(
                    (-1 <= (x - 0.7) / t) & ((x - 0.7) / t <= 0.5), (x - 0.7) / t, 0.5
                ),
            ),
        )
    return np.where(
        x < -np.sqrt(0.8 * t) + 0.7, 0, np.where(x <= 0.7 + t / 2, (x - 0.7) / t, 0.5)
    )


def newton(a, da, u0, du0, t, x, tol=1e-6, max_iter=100):
    X = x

    for _ in range(max_iter):
        g = a(u0(X)) * t + X - x
        gp = da(u0(X)) * du0(X) * t + 1
        if np.abs(gp) < 1e-6:
            print("Derivative is too small")
            return x
        Xnew = X - g / gp
        if np.abs(Xnew - X) < tol:
            return Xnew
        X = Xnew
    # print("Not converged", Xnew)


def build_solution(a, da, u0, du0):
    def u(t, x):
        X = newton(a, da, u0, du0, t, x)
        if X is None:
            return None
        return u0(X)

    return u


u0_tp = lambda x: np.sin(np.pi * 2 * x)

npoints = 800
interval = [0, 1]
T = 0.16
CFL = 0.8

x_centers, sol_LF = finite_volume_method(
    npoints, interval, u0_tp, T, flux_scheme="Lax-Friedrichs", CFL=CFL
)
_, sol_LW = finite_volume_method(
    npoints, interval, u0_tp, T, flux_scheme="Lax-Wendroff", CFL=CFL
)
_, sol_MR = finite_volume_method(
    npoints, interval, u0_tp, T, flux_scheme="Murman-Roe", CFL=CFL
)
print(x_centers.shape, sol_LF.shape, sol_MR.shape)

x_space = np.linspace(interval[0], interval[1], npoints)
# sol_exact = solution_exacte_cours(x_exact, T)
sol_exact = [
    build_solution(
        lambda u: u, lambda u: 1, lambda x: np.sin(2*np.pi*x), lambda x: 2 * np.pi * np.cos(np.pi * 2 * x)
    )(T, x)
    for x in x_space
]
print(len(sol_exact))
plt.plot(sol_exact)

plt.figure(figsize=(10, 6))
plt.plot(x_centers, sol_LF, label="Lax-Friedrichs")
plt.plot(x_centers, sol_LW, label="Lax-Wendroff")
plt.plot(x_centers, sol_MR, label="Murman-Roe")
plt.plot(x_space, sol_exact, label="Exact Solution", linestyle="--")
plt.xlabel("x")
plt.ylabel("u")
plt.title(
    f"Finite Volume Method Solutions with CFL={CFL} and T={T} for {npoints} points"
)
plt.legend()
plt.grid(True)
plt.show()
# %%
