# %%
import numpy as np
import visualization as vis
from scipy.sparse import csc_array, kron, eye, diags
from scipy.sparse.linalg import gmres, spilu, LinearOperator


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


def build_step_matrix(Nx, Ny, dt, D, implicit=False):
    global window

    L = build_2d_laplacian_4th_order_with_neumann_bc(Nx, Ny)
    if implicit:
        A = eye(Nx * Ny) - D * dt * L
    else:
        print("Explicit")
        A = eye(Nx * Ny) + D * dt * L

    # Carve out window
    A = A.tolil()
    for i in range(Nx * Ny):
        ix, iy = divmod(i, Ny)
        if window[ix, iy]:
            A[i, :] = 0
            A[i, i] = 1

    return csc_array(A)

def explicit_heat(u0, Nx, Ny, T, D, rad_T):
    global window, radiator

    dx = 1 / (Nx - 1)
    dy = 1 / (Ny - 1)
    dt = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))
    dt = min(dt, T / 10)

    U = u0(X, Y)
    U[window] = 5.0
    phi = np.zeros_like(U)

    history_manager = vis.HistoryManager(save_frequency=T/30)

    t = 0.0
    A = build_step_matrix(Nx, Ny, dt, D)
    while t < T:
        if int(t / T * 100) % 5 == 0: print(f"Progress: {t/T*100-1:.1f}% Complete", end="\r")

        phi[radiator] = (rad_T - U[radiator])**3

        if t + dt > T:
            dt = T - t
            A = build_step_matrix(Nx, Ny, dt, D)
        U = (A @ U.flatten() + dt * D * phi.flatten()).reshape((Nx, Ny))
        t += dt
        
        history_manager.add(t, U)

    history_manager.add(T, U, force=True)
    return U, history_manager

def implicit_heat(u0, Nx, Ny, T, D, rad_T):
    global window, radiator

    dx = 1 / (Nx - 1)
    dy = 1 / (Ny - 1)
    dt = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))
    dt = min(dt, T / 10)

    U = u0(X, Y)
    phi = np.zeros_like(U)

    history_manager = vis.HistoryManager(save_frequency=T/30)

    t = 0.0
    A = build_step_matrix(Nx, Ny, dt, D, implicit=True)
    ilu = spilu(A)
    M = LinearOperator(A.shape, ilu.solve)
    while t < T:
        if int(t / T * 100) % 5 == 0: print(f"Progress: {int((t*100)//T)}% Complete", end="\r")

        if t + dt > T:
            dt = T - t
            A = build_step_matrix(Nx, Ny, dt, D, implicit=True)
            ilu = spilu(A)
            M = LinearOperator(A.shape, ilu.solve)

        U[window] = 0
        phi[radiator] = (rad_T - U[radiator])**3
        rhs = U.flatten() + dt* D * phi.flatten()

        Unext, exit_code = gmres(A, rhs, M=M, rtol=1e-6)

        if exit_code != 0:
            print(f"GMRES failed to converge at t = {t:.2f}")
            break

        U = Unext.reshape((Nx, Ny))

        t += dt
        history_manager.add(t, U)

    history_manager.add(T, U, force=True)
    return U, history_manager

# %%
D = 2.2e-5
Nx, Ny = 60, 60
T = 20000.0
CFL = 30.0
u0 = lambda x, y: np.ones_like(x) * 20.0
X, Y = np.meshgrid(np.linspace(0, 1, Ny), np.linspace(0, 1, Nx))
window = np.where((0.4 <= X) & (X <= 0.6) & (0.0 == Y), True, False)
radiator = np.where((0.4 <= X) & (X <= 0.6) & (0.9 <= Y) & (Y <= 1.1), True, False)


U, history_manager = implicit_heat(u0, Nx, Ny, T, D, 50.0)
fig = vis.visualize_heat_equation(history_manager, np.linspace(0,1,Nx), np.linspace(0,1,Ny) , mode='heightmap', dynamic_scaling=False)
print("Statistics:")
print(f"Mean: {U.mean()}, Var: {U.std()}", U.min(), U.max())


# %%
D = 2.2e-5
Nx, Ny = 60, 60
T = 2000.0
CFL = 0.45
u0 = lambda x, y: np.ones_like(x) * 20.0
X, Y = np.meshgrid(np.linspace(0, 1, Ny), np.linspace(0, 1, Nx))
window = np.where((0.4 <= X) & (X <= 0.6) & (0.0 == Y), True, False)
radiator = np.where((0.4 <= X) & (X <= 0.6) & (0.9 <= Y) & (Y <= 1.1), True, False)


U, history_manager = explicit_heat(u0, Nx, Ny, T, D, 50.0)
fig = vis.visualize_heat_equation(history_manager, np.linspace(0,1,Nx), np.linspace(0,1,Ny) , mode='heightmap', dynamic_scaling=False)
print("Statistics:")
print(f"Mean: {U.mean()}, Var: {U.std()}", U.min(), U.max())

# %%
