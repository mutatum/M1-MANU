# %%
import numpy as np
import visualization as vis
from scipy.sparse import csc_array, kron, eye, diags


def add_neumann_bc(L1D, dx):
    # Ghost points
    matrix = L1D.tolil()
    matrix[0, 1] += 4 / 3 / dx**2
    matrix[0, 2] += -1 / 12 / dx**2
    matrix[1, 3] += -1 / 12 / dx**2
    matrix[-1, -2] += 4 / 3 / dx**2
    matrix[-1, -3] += -1 / 12 / dx**2
    matrix[-2, -4] += -1 / 12 / dx**2
    return matrix


def build_2d_laplacian_4th_order_with_neumann_bc(Nx, Ny):

    stencil = [-1 / 12, 4 / 3, -5 / 2, 4 / 3, -1 / 12]
    dx = 1 / (Nx - 1)
    dy = 1 / (Ny - 1)

    Lx = diags(stencil, [-2, -1, 0, 1, 2], shape=(Nx, Nx)) / dx**2
    Ly = diags(stencil, [-2, -1, 0, 1, 2], shape=(Ny, Ny)) / dy**2

    Lx = add_neumann_bc(Lx, dx)
    Ly = add_neumann_bc(Ly, dy)

    L = kron(eye(Ny), Lx) + kron(Ly, eye(Nx))

    return L


def build_step_matrix(Nx, Ny, dt, D):

    L = build_2d_laplacian_4th_order_with_neumann_bc(Nx, Ny)
    A = eye(Nx * Ny) + dt * D * L

    # window def
    X, Y = np.meshgrid(np.linspace(0, 1, Ny), np.linspace(0, 1, Nx))
    window = np.where((0.4 <= X) & (X <= 0.6) & (0.0 == Y), True, False)

    # Carve out window
    A = A.tolil()
    b = np.zeros(Nx * Ny)
    for i in range(Nx * Ny):
        ix, iy = divmod(i, Ny)
        if window[ix, iy]:
            A[i, :] = 0
            b[i] = 5.0

    return csc_array(A), b

def heat(u0, Nx, Ny, T, D, rad_T):

    dx = 1 / (Nx - 1)
    dy = 1 / (Ny - 1)
    dt = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))
    dt = min(dt, T / 10)

    X, Y = np.meshgrid(np.linspace(0, 1, Ny), np.linspace(0, 1, Nx))
    U = u0(X, Y)

    radiator = np.where((0.4 <= X) & (X <= 0.6) & (0.0 <= Y) & (Y <= 0.1), True, False)

    phi = np.zeros_like(U)
    history_manager = vis.HistoryManager(save_frequency=T/30)
    history_manager.add(0, U, force=True)

    t = 0.0
    A, b = build_step_matrix(Nx, Ny, dt, D)
    percent = 0.1
    while t < T:

        phi[radiator] = (rad_T - U[radiator])**3

        if t + dt > T:
            dt = T - t
        U = (A @ U.flatten() + b + D*phi.flatten()).reshape((Nx, Ny))
        t += dt
        
        history_manager.add(t, U)

    history_manager.add(T, U, force=True)
    return U, history_manager

# %%
D = 2.2e-5
Nx, Ny = 60, 60
T = 20000
CFL = 0.55
u0 = lambda x, y: np.ones_like(x) * 20.0


U, history_manager = heat(u0, Nx, Ny, T, D, 50.0)
fig = vis.visualize_heat_equation(history_manager, np.linspace(0,1,Nx), np.linspace(0,1,Ny) , mode='heightmap', dynamic_scaling=False)
print("Statistics:")
print(f"Mean: {U.mean()}, Var: {U.std()}", U.min(), U.max())


# %%
