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
    if t>1:
        raise ValueError("t must be less than 1")
    return np.where(x<t, 1, np.where(x<=1, (1-x)/(1-t), 0))

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

def newton(f, df, x0, tol=1e-10, max_iter=1000):
    x = x0
    for i in range(max_iter):
        x_new = x - f(x) / df(x)
        if np.abs(x_new - x) < tol:
            return x_new
        x = x_new
    raise RuntimeError("Newton's method did not converge")


u0_tp = lambda x: np.sin(np.pi * 2 * x)

npoints = 1000
interval = [-3, 2.5]
T = .6
CFL = 0.5

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

x_exact = np.linspace(interval[0], interval[1], npoints)
sol_exact = solution_exacte_cours(x_exact, T)

plt.figure(figsize=(10, 6))
plt.plot(x_centers, sol_LF, label="Lax-Friedrichs")
plt.plot(x_centers, sol_LW, label="Lax-Wendroff")
plt.plot(x_centers, sol_MR, label="Murman-Roe")
plt.plot(x_exact, sol_exact, label="Exact Solution", linestyle="--")
plt.xlabel("x")
plt.ylabel("u")
plt.title(
    f"Finite Volume Method Solutions with CFL={CFL} and T={T} for {npoints} points"
)
plt.legend()
plt.grid(True)
plt.show()
# %%
