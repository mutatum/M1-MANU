# %%
import numpy as np
import visualization as vis
import math
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye, csc_array, kron
from scipy.sparse.linalg import splu, svds


def CoefDF(k, xbar, x):
    """
    Calculate the coefficients of the finite difference approximation for the k-th derivative at a given point xbar.
    Parameters:
    k (int): The order of the derivative.
    xbar (float): The point at which the derivative is to be approximated.
    x (array-like): The array of points used for the finite difference approximation.
    Returns:
    numpy.ndarray: The coefficients of the finite difference approximation for the k-th derivative at xbar.
    Notes:
    - If k is out of bounds (less than 0 or greater than or equal to the length of x), the function will print an error message and return None.
    - The function uses the smallest spacing between points in x and the distance from xbar to the nearest point in x to determine the step size h.
    - The function constructs a matrix A based on the points in x and solves the linear system A * coef = B to find the coefficients.
    """
    x = np.array(x)
    n = len(x)
    A = np.vander(x - xbar, n, increasing=True).T
    b = np.zeros(n)
    b[k] = math.factorial(k)
    h = np.min(x[1:n] - x[0 : n - 1])
    h2 = np.min(np.abs(x - xbar))
    if k < 0 or k >= n:
        raise ValueError(f"k must be between 0 and {n-1}")
    if h2 > 0:
        h = min(h, h2)
    return np.linalg.solve(A, b) * h**k


def build_A_matrix(order, nx, ny, dt, D):
    if order % 2:
        raise ValueError("Order must be even for number of points to be odd")

    n_virtuals = order // 2
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    # Compute normalized coefficients
    coefs = CoefDF(2, 0, np.linspace(-n_virtuals, n_virtuals, 1 + order))
    coefsx = D * (dt / dx**2) * coefs
    coefsy = D * (dt / dy**2) * coefs

    # Construct 1D operators
    Dx = diags(
        coefsx,
        offsets=np.arange(-n_virtuals, n_virtuals + 1),
        shape=(nx + order, nx + order),
    )
    Dy = diags(
        coefsy,
        offsets=np.arange(-n_virtuals, n_virtuals + 1),
        shape=(ny + order, ny + order),
    )

    # Assemble 2D sparse Laplacian using Kronecker product
    Ix = eye(nx + order, format="csc")
    Iy = eye(ny + order, format="csc")
    A = kron(Dx, Iy, format="csc") + kron(Ix, Dy, format="csc")

    # Add identity for implicit time-stepping
    I = eye((nx + order) * (ny + order), format="csc")
    # plt.spy(A, markersize=1)
    # plt.title("Sparsity Pattern of A")
    # plt.show()
    return csc_array(I - A)


def radiator_def(nx, ny, num=1):
    X, Y = np.meshgrid(np.linspace(0, 1, ny), np.linspace(0, 1, nx))
    if num == 1:  # right middle radiator
        radiator = np.where(
            (0.45 <= X) & (X <= 0.55) & (0.9 <= Y) & (Y <= 1), True, False
        )
    elif num == 2:  # middle radiator
        radiator = np.where(
            (0.40 <= X) & (X <= 0.60) & (0.4 <= Y) & (Y <= 0.6), True, False
        )
    elif num == 3:  # left middle radiator
        radiator = np.where(
            (0.0 <= X) & (X <= 0.1) & (0.4 <= Y) & (Y <= 0.6), True, False
        )
    else:  # bottom middle radiator in front of window
        radiator = np.where(
            (0.45 <= X) & (X <= 0.55) & (0 <= Y) & (Y <= 0.1), True, False
        )
    return radiator


def window_def(nx, ny):
    X, Y = np.meshgrid(np.linspace(0, 1, ny), np.linspace(0, 1, nx))
    window = np.where((0.4 <= X) & (X <= 0.6) & (0 <= Y) & (Y <= 1 / ny), True, False)
    return window


def heat_equation_2d(
    nx,
    ny,
    T,
    CFL=0.9,
    D=1,
    T_ext=5.0,
    T_int=20.0,
    T_rad=50.0,
    space_order=4,
    time_points=20,
):
    # Grid setup
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)
    dt_diffusion = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))  # Condition CFL
    dt_time = T / time_points
    dt = min(dt_diffusion, dt_time)
    print(f"dt: {dt}, T: {T}, steps: {math.ceil(T/dt)}")

    # Initialize temperature field
    U = T_int * np.ones((nx, ny), dtype=float)

    radiator = radiator_def(nx, ny, 4)

    window = window_def(nx, ny)
    U[window] = T_ext

    save_interval = T / time_points
    next_save = save_interval
    history = [(0, U.copy())]
    t = 0
    nt = 0

    total_steps = int(T / dt)

    print("About to build A")
    A = build_A_matrix(order=space_order, nx=nx, ny=ny, dt=dt, D=D)

    lu = splu(A)
    phi = np.zeros((nx, ny), dtype=float)
    vis.print_progress_bar(
        0, total_steps, prefix="Progress:", suffix="Complete", length=50
    )
    # Boundary conditions
    BC_coefs = CoefDF(1, 0, np.linspace(0, 1, space_order + 1))
    BC_coefs = -BC_coefs[1:] / BC_coefs[0]
    U_prev = U.copy()
    while t < T:
        # Adjust dt for final step
        if T - t < dt:
            dt = T - t
            A = build_A_matrix(order=space_order, nx=nx, ny=ny, dt=dt,D=D)
            lu = splu(A)
        # Radiator heating
        phi[radiator] = (T_rad - U_prev[radiator]) ** 2 * (T_rad - U[radiator])

        n_virtuals = space_order // 2
        Uflataugmented = np.pad(U, n_virtuals, mode="reflect").flatten()
        phiflataugmented = np.pad(phi, n_virtuals, mode="constant").flatten()

        U = lu.solve(Uflataugmented
                    # + D * phiflataugmented
                ).reshape(
            (nx + space_order, ny + space_order)
        )[n_virtuals:-n_virtuals, n_virtuals:-n_virtuals]
        # 
        # Handle points one away from edges using 4th order FD for Laplacian
        # Calculate coefficients for second derivative

        U[0, 1:-1] = U[1 : space_order + 1, 1:-1].T @ BC_coefs  # bottom
        U[-1, 1:-1] = U[-space_order - 1 : -1, 1:-1].T @ BC_coefs[::-1]  # top
        U[1:-1, 0] = U[1:-1, 1 : space_order + 1] @ BC_coefs[::-1]
        U[1:-1, -1] = U[1:-1, -space_order - 1 : -1] @ BC_coefs  #

        # U[0,1:-1] = U[1,1:-1]
        # U[-1,1:-1] = U[-2,1:-1]
        # U[1:-1,0] = U[1:-1,1]
        # U[1:-1,-1] = U[1:-1,-2]

        corner_coefs = CoefDF(1, 0, np.linspace(0, 1, space_order + 1))
        corner_coefs = -corner_coefs[1:] / corner_coefs[0]
        U[0,0] = .5 * (U[1:space_order+1,0] @ corner_coefs + U[0,1:space_order+1] @ corner_coefs)
        U[0,-1] = .5 * (U[1:space_order+1,-1] @ corner_coefs + U[0,-2:-space_order-2:-1] @ corner_coefs)
        U[-1,0] = .5 * (U[-space_order-1:-1,0] @ corner_coefs[::-1] + U[-1,1:space_order+1] @ corner_coefs)
        U[-1,-1] = .5 * (U[-space_order-1:-1,-1] @ corner_coefs[::-1] + U[-1,-2:-space_order-2:-1] @ corner_coefs)

        
        # U[0, 0] = 0.5 * (U[1, 0] + U[0, 1])
        # U[0, -1] = 0.5 * (U[0, -2] + U[1, -1])
        # U[-1, 0] = 0.5 * (U[-1, 1] + U[-2, 0])
        # U[-1, -1] = 0.5 * (U[-1, -2] + U[-2, -1])
        U[window] = T_ext

        t += dt
        nt += 1
        vis.print_progress_bar(
            nt, total_steps, prefix="Progress:", suffix="Complete", length=50
        )
        if t >= next_save:
            next_save += save_interval
            history.append((t, U.copy()))
        U_prev = U.copy()
    if history[-1][0] != T:
        history.append((T, U.copy()))

    return U, history


# %%

if __name__ == "__main__":

    nx = 40
    ny = 40
    T = 400000.0
    CFL = 0.9
    D = 2.2e-5 # D is the thermal diffusivity in m^2/s
    T_ext = 5.0
    T_int = 20.0
    T_rad = 50.0
    time_points = 50
    solution, history = heat_equation_2d(
        nx=nx,
        ny=ny,
        T=T,
        CFL=CFL,
        D=D,
        T_ext=T_ext,
        T_int=T_int,
        T_rad=T_rad,
        time_points=time_points,
    )
    # print("len: ", len(history))
    # Statistical analysis
    print("Statistical analysis:")
    print("Max temperature: ", np.max(solution))
    print("Min temperature: ", np.min(solution))
    print("Mean temperature: ", np.mean(solution))
    print("Standard deviation: ", np.std(solution))
    vis.plot_interactive(history, np.linspace(0, 1, nx), np.linspace(0, 1, ny))
