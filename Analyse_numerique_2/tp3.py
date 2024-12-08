# %%
import numpy as np
import visualization as vis
import math
import matplotlib.pyplot as plt
from scipy.sparse import eye, lil_array, csc_array
from scipy.sparse.linalg import splu


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
    if nx < order + 2 or ny < order + 2:
        raise ValueError("Grid too small for order of FD scheme")
    n_side = order//2

    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    factor_x = D * dt / dx**2
    factor_y = D * dt / dy**2

    coeffs_centered = CoefDF(2, 0.5, np.linspace(0, 1, order + 1))
    # print(coeffs_centered)

    # bcs = boundaries coefficients
    # 2 + order because decentered FD scheme needs an additional point to reach order
    bcs = [CoefDF(2, i, np.arange(0, order + 2)) for i in range(n_side)]
    # print("bcs", len(bcs), len(bcs[0]))

    A = lil_array((nx * ny, nx * ny), dtype=float)

    print("nx, ny", nx, ny, nx*ny)
    for i in range(nx*ny):
        ix = i // ny
        iy = i % ny
        # print("i, ix, iy", i, ix, iy)
        if iy < n_side:
            # print("condition start X")
            coeffs_y_start = CoefDF(2, iy, np.arange(0, order + 2)) * factor_y
            # print("coeffs_y_start", coeffs_y_start.size)
            for j, coeff in enumerate(coeffs_y_start):
                # print('\t access', i, ix*ny + j, "value", coeff)
                A[i, ix*ny + j] += coeff
        elif iy >= ny - n_side:
            # print("condition end X")
            coeffs_y_end = CoefDF(2, order+2+iy-ny, np.arange(0, order + 2)) * factor_y
            # print(order+2+iy-ny, np.arange(0, order + 2))
            # print("coeffs_y_end", len(coeffs_y_end))
            for j, coeff in enumerate(coeffs_y_end):
                A[i, (ix+1)*ny + j - len(coeffs_y_end)] += coeff
                # print('\t access', i, (ix+1)*ny + j - len(coeffs_y_end), "value", coeff * factor_y)
        else:
            # print("condition center X")
            # print("coeffs_centered", coeffs_centered*factor_y)
            for j, coeff in enumerate(coeffs_centered):
                # print('\t access', i, ix*ny+iy +j -n_side, "value", coeff * factor_y)
                A[i, ix*ny+iy +j -n_side] += coeff * factor_y
        if ix < n_side:
            # print("condition start Y")
            coeffs_x_start = CoefDF(2, ix, np.arange(0, order + 2)) * factor_x
            for j, coeff in enumerate(coeffs_x_start):
                # print("access", i, (ix+j)*ny+iy, (ix+j), "value", coeff * factor_x)
                A[i, j*ny + iy] += coeff
        elif ix >= nx - n_side:
            # print(nx*ny-n_side*ny)
            coeffs_x_end = CoefDF(2, order+2+ix-nx, np.arange(0, order + 2)) * factor_x
            # print(order+2+ix-nx, np.arange(0, order + 2))
            # print("condition end Y")
            for j, coeff in enumerate(coeffs_x_end):
                # print("access", i, i-order-1+j*ny,(ix+j-(len(coeffs_x_end)-1))*ny)
                A[i, (j-len(coeffs_x_end))*ny + iy] += coeff
        else:
            # print("condition center Y")
            for j, coeff in enumerate(coeffs_centered):
                # print("\t access", i, (ix + j - n_side)*ny + iy, "value", coeff * factor_x)
                A[i, (ix + j - n_side)*ny + iy] += coeff * factor_x

    # spy(A) visualisation de la matrice
    I = eye(nx * ny, format="lil")
    return csc_array((I - A))

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
            (0.40 <= X) & (X <= 0.60) & (0 <= Y) & (Y <= 0.2), True, False
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

    radiator = radiator_def(nx, ny, 2)

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
        # phi[radiator] = (T_rad - U_prev[radiator]) ** 2 * (T_rad - U[radiator])
        phi[radiator] = (T_rad - U[radiator])**3

        Uflat= U.flatten()
        phiflat = phi.flatten()

        U= lu.solve(Uflat
                    + D * phiflat
                ).reshape(
            (nx, ny)
        )

        # U[0, 1:-1] = U[1 : space_order + 1, 1:-1].T @ BC_coefs  # bottom
        # U[-1, 1:-1] = U[-space_order - 1 : -1, 1:-1].T @ BC_coefs[::-1]  # top
        # U[1:-1, 0] = U[1:-1, 1 : space_order + 1] @ BC_coefs[::-1]
        # U[1:-1, -1] = U[1:-1, -space_order - 1 : -1] @ BC_coefs  #

        U[0,1:-1] = U[1,1:-1]
        U[-1,1:-1] = U[-2,1:-1]
        U[1:-1,0] = U[1:-1,1]
        U[1:-1,-1] = U[1:-1,-2]

        U[0,0] = U[1,0]
        U[0,-1] = U[1,-1]
        U[-1,0] = U[-2,0]
        U[-1,-1] = U[-2,-1]

        
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
    CFL = 1
    D = 2.2e-5 # D is the thermal diffusivity in m^2/s
    T_ext = 5.0
    T_int = 20.0
    T_rad = 50.0
    time_points = 300
    solution, history = heat_equation_2d(
        nx=nx,
        ny=ny,
        T=T,
        CFL=CFL,
        D=D,
        T_ext=T_ext,
        T_int=T_int,
        T_rad=T_rad,
        space_order = 4,
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

# %%
