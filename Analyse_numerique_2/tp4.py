# %%
from typing import Callable
import numpy as np
from matplotlib import pyplot as plt


def MurmanRoe(ul, ur, f, df, epsilon=1e-6) -> np.ndarray:
    delta_u = ur - ul
    delta_f = f(ur) - f(ul)

    if isinstance(ul, int):
        return np.abs(delta_f / delta_u)

    mask = np.abs(delta_u) < epsilon
    gamma = np.zeros_like(delta_u)
    gamma[mask] = df(ul)[mask]
    gamma[~mask] = delta_f[~mask] / delta_u[~mask]
    return np.abs(gamma)


def LaxFriedrichs_global(ul, ur, df) -> float:
    return np.max(np.abs(df(np.union1d(ul, ur))))


def LaxFriedrichs(dx, dt) -> float:
    return dx / dt


def LaxWendroff(ul, ur, f, df, dx, dt) -> float:
    return np.abs((dt / dx) * df((ul + ur) / 2) * MurmanRoe(ul, ur, f, df))


def Rusanov(ul, ur, df) -> np.ndarray:
    return np.maximum(np.abs(df(ul)), np.abs(df(ur)))


def finite_volume_method(
    f: Callable,
    df: Callable,
    u0: Callable,
    T: float,
    interval: tuple[float, float],
    m: int,
    periodic: bool = True,
    solution: Callable = None,
    flux_scheme: str = "Lax-Friedrichs_global",
    CFL: float = 0.5,
) -> np.ndarray:

    def F(ul, ur, gamma):
        return 0.5 * (f(ul) + f(ur)) - 0.5 * gamma * (ur - ul)

    a, b = interval
    dx = (b - a) / m
    x = np.linspace(a, b, m + 1)
    C = (x[1:] + x[:-1]) / 2
    u = u0(C)
    t = 0.0

    while t < T:
        g = np.max(np.abs(df(u)))
        if g < 1e-8:
            dt = CFL * dx / 2.34
            "here"
        else:
            dt = CFL * dx / g

        if t + dt > T:
            dt = T - t

        if periodic:
            u_extended = np.concatenate(([u[-1]], u, [u[0]]))
        else:
            u_extended = np.concatenate(([solution(t, a)], u, [solution(t, b)]))
        u_L = u_extended[:-1]
        u_R = u_extended[1:]

        # array gammas
        if flux_scheme == "Murman-Roe":
            gamma = np.abs(MurmanRoe(u_L, u_R, f, df))  # array
        elif flux_scheme == "Rusanov" or flux_scheme == "Lax-Friedrichs_local":
            gamma = Rusanov(u_L, u_R, df)
        elif flux_scheme == "Lax-Wendroff":
            gamma = LaxWendroff(u_L, u_L, f, df, dx, dt)  # single value
        else: # scalar gammas
            if flux_scheme == "Lax-Friedrichs_global":
                gamma_value = LaxFriedrichs_global(u_L, u_R, df)  # single value
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
    # print("Not converged", Xnew)
    return None


def build_solution(a, da, u0, du0):
    def u(t, x):
        X = newton(a, da, u0, du0, t, x)
        if X is None:
            return None
        return u0(X)

    return u

# %% plot_comparison

def plot_comparison(
    f,
    df,
    u0,
    T,
    interval,
    npoints=100,
    CFL=0.2,
    ddf: Callable = None,
    du0: Callable = None,
    periodic: bool = True,
    solution: Callable | tuple = None,
    flux: str = None,
):

    x_centers, sol_MR = finite_volume_method(
        f=f,
        df=df,
        u0=u0,
        T=T,
        interval=interval,
        m=npoints,
        flux_scheme="Murman-Roe",
        CFL=CFL,
        periodic=periodic,
        solution=solution,
    )
    _, sol_LW = finite_volume_method(
        f=f,
        df=df,
        u0=u0,
        T=T,
        interval=interval,
        m=npoints,
        flux_scheme="Lax-Wendroff",
        CFL=CFL,
        periodic=periodic,
        solution=solution,
    )
    _, sol_R = finite_volume_method(
        f=f,
        df=df,
        u0=u0,
        T=T,
        interval=interval,
        m=npoints,
        flux_scheme="Rusanov",
        CFL=CFL,
        periodic=periodic,
        solution=solution,
    )
    _, sol_LF = finite_volume_method(
        f,
        df,
        u0,
        T=T,
        interval=interval,
        m=npoints,
        flux_scheme="Lax-Friedrichs_global",
        CFL=CFL,
        periodic=periodic,
        solution=solution,
    )
    x_solution = np.linspace(interval[0], interval[1], 200)
    plt.figure(figsize=(15, 10))
    if du0 != None and ddf != None:
        u_solution = [build_solution(df, ddf, u0, du0)(T, x) for x in x_solution]
    elif solution != None:
        if isinstance(solution, tuple):
            x_solution, u_solution = solution
        elif callable(solution):
            u_solution = [solution(T, x) for x in x_solution]
        else:
            raise ValueError("Invalid solution")
        plt.plot(x_solution, u_solution, label="Exact Solution", linestyle="--")

    if flux == None:
        flux = f.__name__

    plt.plot(x_centers, sol_LF, "-o", label="Global Lax-Friedrichs")
    plt.plot(x_centers, sol_LW, "-o", label="Lax-Wendroff")
    plt.plot(x_centers, sol_MR, "-o", label="Murman-Roe")
    plt.plot(x_centers, sol_R, "-o", label="Local Lax-Friedrichs / Rusanov")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(f"FV {flux} with CFL={CFL} and T={T} for {npoints} points")
    plt.legend()
    plt.grid(True)
    plt.show()


# %%

def burgers(u):  # Burgers' equation
    return 0.5 * u * u


def dburgers(u):
    return u


def buckley(u):
    return 4 * u**2 / (4 * u**2 + (1 - u) ** 2)


def dbuckley(u):
    return 8 * u * (1 - u) / (5 * u**2 - 2 * u + 1) ** 2


def u0_cours(x):
    return np.where(x < 0.3, 0.0, np.where(x < 0.7, -1.0, 0.5))


def solution_exacte_cours(t, x):
    if t > 1:
        raise ValueError("t must be less than 1")
    return np.where(x < t, 1, np.where(x <= 1, (1 - x) / (1 - t), 0))


def solution_td(t, x):
    t_etoile = 0.8
    if t == 0:
        return u0_cours(x)
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


def u0_buckley(x):
    x = np.array(x)  # Ensure x is a NumPy array
    return np.where((-0.5 < x) & (x < 0), 1.0, 0.0)


u0_tp = lambda x: np.sin(np.pi * 2 * x)
du0_tp = lambda x: 2 * np.pi * np.cos(2 * np.pi * x)


# %% Burgers' equation

plot_comparison(
    burgers,
    dburgers,
    u0_tp,
    T=0.3,
    interval=[0, 1],
    npoints=100,
    ddf=lambda _: 1,
    du0=du0_tp,
    CFL=0.5,
    flux="Burgers",
    periodic=True,
)

# %%

plot_comparison(
    burgers,
    dburgers,
    u0_cours,
    T=0.8,
    interval=[-2, 2],
    npoints=100,
    solution=solution_td,
    flux="Buckley",
    periodic=False,
)

# %%
solution_data = np.loadtxt("burgers_t=4tc.dat")
x_solution = solution_data[:, 0]
u_solution = solution_data[:, 1]

plot_comparison(
    burgers,
    dburgers,
    u0_tp,
    T=0.65,
    interval=[x_solution[0], x_solution[-1]],
    # interval=[-1, 1],
    npoints=200,
    solution=(x_solution, u_solution),
    flux="Buckley",
)


# %% Buckley 40 points
solution_data = np.loadtxt("buckley.dat")
x_solution = solution_data[:, 0]
u_solution = solution_data[:, 1]


plot_comparison(
    buckley,
    dbuckley,
    u0_buckley,
    T=0.4,
    interval=[x_solution[0], x_solution[-1]],
    npoints=40,
    solution=(x_solution, u_solution),
    flux="Buckley",
)
# %% Buckley 80 points
solution_data = np.loadtxt("buckley.dat")
x_solution = solution_data[:, 0]
u_solution = solution_data[:, 1]


plot_comparison(
    buckley,
    dbuckley,
    u0_buckley,
    T=0.4,
    interval=[x_solution[0], x_solution[-1]],
    npoints=80,
    solution=(x_solution, u_solution),
    flux="Buckley",
)

# %% Buckley 160 points

plot_comparison(
    buckley,
    dbuckley,
    u0_buckley,
    T=0.4,
    interval=[x_solution[0], x_solution[-1]],
    npoints=160,
    solution=(x_solution, u_solution),
    # ddf=ddbuckley,
    # u0=u0_buckley,
)

# %% Buckley 320 points
solution_data = np.loadtxt("buckley.dat")
x_solution = solution_data[:, 0]
u_solution = solution_data[:, 1]


plot_comparison(
    buckley,
    dbuckley,
    u0_buckley,
    T=0.4,
    CFL=0.1,
    interval=[x_solution[0], x_solution[-1]],
    npoints=320,
    solution=(x_solution, u_solution),
    flux="Buckley",
)
# %%

c = 1.5
plot_comparison(
    f=lambda u: u*c,
    df=lambda u: np.ones_like(u)*c,
    u0=lambda x: np.where(x < 0.5, np.where(x>0, 1.0, 0.0), 0.0),
    T=103,
    interval=[-.5, 1],
    npoints=60,
    CFL=0.9,
    flux="Advection",
    periodic=True,
)

