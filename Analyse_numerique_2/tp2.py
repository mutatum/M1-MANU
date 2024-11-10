# %%
import matplotlib.pyplot as plt
import numpy as np


def schema_chaleur2D(u0, mesh_dimensions, T, D=2.2e-5, CFL=1):
    m, p = mesh_dimensions
    if p < 1 or m < 1:
        print("Incorrect mesh dimensions")
        return None
    dx = 1 / (m - 1)
    dy = 1 / (p - 1)
    dt = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))  # Condition CFL
    lx = D * dt / dx**2
    ly = D * dt / dy**2

    def point_from_index(i, j):
        return i * dx, j * dy

    def apply_u0(u0):
        return lambda i, j: u0(*point_from_index(i, j))

    x = np.linspace(0, 1.0, m)
    y = np.linspace(0, 1.0, p)
    X, Y = np.meshgrid(x, y, indexing="ij")
    Uh = u0(X, Y)
    vmin = np.min(Uh)
    vmax = np.max(Uh)
    t = 0
    it = 0

    plt.figure(figsize=(8, 6))
    cp = plt.contourf(X, Y, Uh, levels=50, cmap='hot', vmin=vmin, vmax=vmax)  # Fixed scale
    plt.colorbar(cp)
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title('Temperature Distribution Over Time')
    
    ## Graphing 3D
    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection="3d")
    # surf = ax.plot_surface(X, Y, Uh, cmap="hot", vmin=vmin, vmax=vmax)
    # ax.set_zlim(vmin, vmax)
    # fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    # ax.set_xlabel("X (m)")
    # ax.set_ylabel("Y (m)")
    # ax.set_zlabel("Temperature")
    # ax.set_title("3D Temperature Distribution Over Time")
    while t < T:
        it += 1
        if T - t < dt:
            dt = T - t
            lx = D * dt / dx**2
            ly = D * dt / dy**2
        t += dt
        Uh[1:-1, 1:-1] = (
            Uh[1:-1, 1:-1]
            + lx * (Uh[2:, 1:-1] + Uh[:-2, 1:-1] - 2 * Uh[1:-1, 1:-1])
            + ly * (Uh[1:-1, 2:] + Uh[1:-1, :-2] - 2 * Uh[1:-1, 1:-1])
        )
        Uh[0, :] = Uh[-1, :] = Uh[:, 0] = Uh[:, -1] = 0.0  # Apply boundary conditions

        if it % 1 == 0:
            plt.clf()  # Clear the current plot
            cp = plt.contourf(X, Y, Uh, levels=50, cmap='hot', vmin=vmin, vmax=vmax)
            plt.colorbar(cp)
            plt.title(f'Temperature Distribution at t = {t:.2f} s')
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            plt.axis('equal')
            plt.pause(0.1)  # Pause for 0.1 seconds
            
            # Graphing 3D
            # ax.cla()  # Clear the current plot
            # surf = ax.plot_surface(X, Y, Uh, cmap="hot", vmin=vmin, vmax=vmax)
            # ax.set_zlim(vmin, vmax)
            # ax.set_xlabel("X (m)")
            # ax.set_ylabel("Y (m)")
            # ax.set_zlabel("Temperature")
            # ax.set_title(f"Temperature Distribution at t = {t:.2f} s")
            # plt.pause(0.01)  # Pause for 0.1 seconds to create animation effect

    return Uh


def square(x, y):
    # Initial condition: set to 1 if within the region [0.4, 0.6] x [0.4, 0.6], else 0
    return np.where((0.4 <= x) & (x <= 0.6) & (0.4 <= y) & (y <= 0.6), 2.0, 0.0)


def circle(x, y):
    # Circular hot spot centered in the domain with radius 0.2
    r = np.sqrt((x-0.5) ** 2 + (y-0.5) ** 2)
    return np.where(r <= 0.6, 2.0, 0.0)


def diagonal(x, y):
    # Diagonal gradient from (0,0) to (L, L)
    return (x + y) / 2


def checkerboard(x, y):
    # Checkerboard pattern with square cells of size 0.02m
    return np.where(((np.floor((x - 0.5) / 0.1) % 2) == (np.floor((y - 0.5) / 0.1) % 2)), 1.0, 0.0)


def random(x, y):
    # Random hot spots scattered throughout the domain
    # np.random.seed(1)  # Seed for reproducibility
    random_spots = (np.random.rand(*x.shape) < 0.05).astype(
        float
    )  # 5% chance for each cell to be hot
    return random_spots * 20


# Run the simulation
T = 100
D = 2.2e-5  # Thermal diffusivity of air in m^2/s
mesh_dimensions = (40, 40)
CFL = 0.45
final_temperature_distribution = schema_chaleur2D(
    random, mesh_dimensions, T=T, D=D, CFL=CFL
)

# Generate the grid for plotting
m, p = mesh_dimensions
x = np.linspace(0, 1, m)
y = np.linspace(0, 1, p)
X, Y = np.meshgrid(x, y)

# Create the 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, final_temperature_distribution, cmap='viridis')

# Customize the plot appearance
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Temperature')
ax.set_title(f"Schéma numérique explicite à T = {T} pour une CFL de {CFL}")

plt.show()

# %%
