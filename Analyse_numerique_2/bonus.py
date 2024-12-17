# %%
from typing import Callable
import numpy as np
from matplotlib import pyplot as plt
import utils

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
    flux_schemes_to_plot: list = None,
):
    if flux_schemes_to_plot is None:
        flux_schemes_to_plot = [
            "Lax-Friedrichs_global",
            # "Lax-Wendroff",
            # "Murman-Roe",
            # "Rusanov",
        ]
    solutions_order_1 = {}
    solutions_order_2_RK2 = {}
    solutions_order_2_PC = {}
    for scheme in flux_schemes_to_plot:
        x_centers, sol = finite_volume_method(
            f=f,
            df=df,
            u0=u0,
            T=T,
            interval=interval,
            m=npoints,
            order=1,
            flux_scheme=scheme,
            CFL=CFL,
            periodic=periodic,
            solution=solution,
        )
        solutions_order_1[scheme] = sol
        x_centers, sol = finite_volume_method(
            f=f,
            df=df,
            u0=u0,
            T=T,
            interval=interval,
            m=npoints,
            order=2,
            flux_scheme=scheme,
            scheme="RK2",
            CFL=CFL,
            periodic=periodic,
            solution=solution,
        )
        solutions_order_2_RK2[scheme] = sol
        x_centers, sol = finite_volume_method(
            f=f,
            df=df,
            u0=u0,
            T=T,
            interval=interval,
            m=npoints,
            order=2,
            flux_scheme=scheme,
            scheme='PC',
            CFL=CFL,
            periodic=periodic,
            solution=solution,
        )
        solutions_order_2_PC[scheme] = sol

    x_solution = np.linspace(interval[0], interval[1], 200)
    plt.figure(figsize=(25, 15))

    if du0 is not None and ddf is not None:
        u_solution = [utils.build_solution(df, ddf, u0, du0)(T, x) for x in x_solution]
    elif solution is not None:
        if isinstance(solution, tuple):
            x_solution, u_solution = solution
        elif callable(solution):
            u_solution = [solution(T, x) for x in x_solution]
        else:
            raise ValueError("Invalid solution")
        plt.plot(x_solution, u_solution, label="Exact Solution", linestyle="--")

    if flux is None:
        flux = f.__name__

    for scheme in flux_schemes_to_plot:
        label = scheme.replace("_", " ").title()
        plt.plot(x_centers, solutions_order_1[scheme], "-o", label=label + ' Order 1')
        plt.plot(x_centers, solutions_order_2_RK2[scheme], "-o", label=label + ' Order 2 RK2')
        # plt.plot(x_centers, solutions_order_2_PC[scheme], "-o", label=label + ' Order 2 PC')

    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(f"FV {flux} with CFL={CFL} and T={T} for {npoints} points")
    plt.legend()
    plt.grid(True)
    plt.show()



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


u0_tp  = lambda x: np.sin(np.pi * 2 * x)
du0_tp = lambda x: 2 * np.pi * np.cos(2 * np.pi * x)


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


def compute_gamma_murman_roe(u_L, u_R, f, df, dx, dt):
    return np.abs(MurmanRoe(u_L, u_R, f, df))

def compute_gamma_rusanov(u_L, u_R, f, df, dx, dt):
    return Rusanov(u_L, u_R, df)

def compute_gamma_lax_wendroff(u_L, u_R, f, df, dx, dt):
    return LaxWendroff(u_L, u_R, f, df, dx, dt)

def compute_gamma_lax_friedrichs_global(u_L, u_R, f, df, dx, dt):
    gamma_value = LaxFriedrichs_global(u_L, u_R, df)
    return np.full_like(u_L, gamma_value)

def minmod(x,y,z):
    return np.minimum(0, np.maximum(x,np.maximum(y,z))) + np.maximum(0, np.minimum(x,np.minimum(y,z)))

def slope_limiter(u_padded, dx, alpha):
    x_slope = (u_padded[2:]-u_padded[:-2])/(2*dx)
    y_slope = (u_padded[1:-1]-u_padded[:-2])/dx
    z_slope = (u_padded[2:]-u_padded[1:-1])/dx
    return minmod(x_slope, 2*alpha*y_slope, 2*alpha*z_slope)

flux_schemes = {
    "Murman-Roe": compute_gamma_murman_roe,
    "Rusanov": compute_gamma_rusanov,
    "Lax-Wendroff": compute_gamma_lax_wendroff,
    "Lax-Friedrichs_global": compute_gamma_lax_friedrichs_global,
}

def FV_step(u, f, df, dx, dt, flux_scheme, order=1):

    def F(ul, ur, gamma):
        return 0.5 * (f(ul) + f(ur)) - 0.5 * gamma * (ur - ul)

    u_padded = np.pad(u, order, mode="wrap")

    if order == 2:
        slopes = slope_limiter(u_padded, dx, .5)
        u_L = u_padded[1:-2] + dx * slopes[:-1] / 2
        u_R = u_padded[2:-1] - dx * slopes[1:] / 2
    elif order == 1:
        u_L = u_padded[:-1]
        u_R = u_padded[1:]

    if flux_scheme in flux_schemes:
        gamma = flux_schemes[flux_scheme](u_L, u_R, f, df, dx, dt)
    else:
        raise ValueError(f"Invalid scheme: {flux_scheme}")

    F_i = F(u_L, u_R, gamma)

    flux_divergence = (F_i[1:] - F_i[:-1]) / dx
    return dt * flux_divergence

def finite_volume_method(
    f: Callable,
    df: Callable,
    u0: Callable,
    T: float,
    interval: tuple[float, float],
    m: int,
    order: int,
    periodic: bool = True,
    solution: Callable = None,
    flux_scheme: str = "Lax-Friedrichs_global",
    scheme: str = "RK2",
    CFL: float = 0.5,
) -> np.ndarray:


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
        else:
            dt = CFL * dx / g

        if t + dt > T:
            dt = T - t
        
        if order == 1:
            u = u - FV_step(u, f, df, dx, dt, flux_scheme, order)
        elif order == 2:
            match scheme:
                case "RK2":
                    u1 = u - FV_step(u, f, df, dx, dt, flux_scheme, order)
                    u2 = u1 - FV_step(u1, f, df, dx, dt, flux_scheme, order)
                    u = 0.5 * (u + u2)
                case "PC": 
                    uhalf = u - .5 * FV_step(u, f, df, dx, dt/2, flux_scheme, order)
                    u = u - FV_step(uhalf, f, df, dx, dt, flux_scheme, order)

        if np.any(np.isnan(u)) or np.any(np.abs(u) > 1e6):
            raise RuntimeError("Solution diverged")

        t += dt
    return C, u

# %% Burgers' equation

plot_comparison(
    burgers,
    dburgers,
    u0_tp,
    T=0.6,
    interval=[0, 1],
    npoints=40,
    ddf=lambda _: 1,
    du0=du0_tp,
    CFL=0.2,
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
    flux="Burgers",
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
    T=2/np.pi,
    interval=[x_solution[0], x_solution[-1]],
    # interval=[-1, 1],
    npoints=20,
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
    npoints=2000,
    solution=(x_solution, u_solution),
    flux="Buckley",
)
# %%

c = 1
a,b = 0, 1
L = (b - a)*c
u0 = lambda x: np.where((x % L + a < 0.625) & (x % L + a >= 0.375), 1.0, 0.0)
T = 15
CFL = 0.5
points = 100
# u0 = lambda x: np.where(x < 0.5, np.where(x>0, 1.0, 0.0), 0.0)
x_solution = np.linspace(a, b, 200)
u_solution = u0(x_solution- (T%c)*c)
plot_comparison(
    f=lambda u: u*c,
    df=lambda u: np.ones_like(u)*c,
    u0=u0,
    T=T,
    interval=[a, b],
    npoints=points,
    CFL=CFL,
    solution = (x_solution, u_solution),
    flux="Advection",
    periodic=True,
)


# %%
