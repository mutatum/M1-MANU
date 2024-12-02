# %%
import numpy as np
import visualization as vis
import math
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye, csc_array, kron
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


def build_A_matrix(order, nx, ny, dt):
    """Build system matrix using appropriate stencils for each point"""
    N = nx * ny
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)
    rows, cols, data = [], [], []

    if order >= nx or order >= ny:
        raise ValueError("Order of the stencil is too large for the grid size")
        n_blocks = 7  # Total number of blocks (rows and columns)
block_size = n  # Size of each block

# Predefine blocks
D1 = csc_matrix(np.diag(np.random.rand(block_size)))
D2 = csc_matrix(np.diag(np.random.rand(block_size)))
D3 = csc_matrix(np.diag(np.random.rand(block_size)))
B_blocks = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10]
Y_blocks = [Y1, Y2]
Zero = csc_matrix((block_size, block_size))

blocks = []
for i in range(n_blocks):
    row_blocks = []
    for j in range(n_blocks):
        if i == j:
            # Diagonal blocks
            if i in [0, n_blocks - 1]:
                row_blocks.append(D1)
            elif i in [1, n_blocks - 2]:
                row_blocks.append(D2)
            else:
                row_blocks.append(D3)
        elif i == 0 and 1 <= j <= 5:
            # First row, B blocks
            row_blocks.append(B_blocks[j - 1])
        elif i == 1 and 0 <= j <= 5:
            # Second row, B blocks
            row_blocks.append(B_blocks[5 + j])
        elif 2 <= i <= n_blocks - 3 and abs(i - j) <= 2:
            # Middle rows, Y blocks
            if abs(i - j) == 1:
                row_blocks.append(Y1)
            elif abs(i - j) == 2:
                row_blocks.append(Y2)
            else:
                row_blocks.append(D3)
        elif i == n_blocks - 2 and 1 <= j <= 6:
            # Second last row, B blocks
            row_blocks.append(B_blocks[5 + (n_blocks - j - 1)])
        elif i == n_blocks - 1 and 0 <= j <= 5:
            # Last row, B blocks
            row_blocks.append(B_blocks[j])
        else:
            # Zero blocks
            row_blocks.append(Zero)
    blocks.append(row_blocks)

    A = bmat(blocks, format='csc')
    return A

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
    lx = D * dt / dx**2
    ly = D * dt / dy**2

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
    A = build_A_matrix(order=space_order, nx=nx, ny=ny, dt=dt)
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
            A = build_A_matrix(order=space_order, nx=nx, ny=ny, dt=dt)
            lu = splu(A)
        # Radiator heating
        phi[radiator] = (T_rad - U_prev[radiator]) ** 2 * (T_rad - U[radiator])

        Uflat= U.flatten()
        phiflat = phi.flatten()

        U = lu.solve(Uflat
                    #+ dt * phiflat
                    ).reshape((nx, ny)
        )

        U[0, 1:-1] = U[1 : space_order + 1, 1:-1].T @ BC_coefs  # bottom
        U[-1, 1:-1] = U[-space_order - 1 : -1, 1:-1].T @ BC_coefs[::-1]  # top
        U[1:-1, 0] = U[1:-1, 1 : space_order + 1] @ BC_coefs[::-1]
        U[1:-1, -1] = U[1:-1, -space_order - 1 : -1] @ BC_coefs  #


        U[0, 0] = 0.5 * (U[0, 1] + U[1, 0])
        U[0, -1] = 0.5 * (U[0, -2] + U[1, -1])
        U[-1, 0] = 0.5 * (U[-1, 1] + U[-2, 0])
        U[-1, -1] = 0.5 * (U[-1, -2] + U[-2, -1])
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

    nx = 15
    ny = 15
    T = 0.1
    CFL = .01
    D = 1
    T_ext = 5.0
    T_int = 20.0
    T_rad = 50.0
    time_points = 200
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
    vis.plot_interactive(history, np.linspace(0, 1, nx), np.linspace(0, 1, ny))

# %%
