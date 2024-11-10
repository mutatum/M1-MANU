# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button
import threading
import time as pytime


def schema_chaleur2D(u0, mesh_dimensions, T, D=2.2e-5, CFL=1, every=1):
    m, p = mesh_dimensions
    if p < 1 or m < 1:
        print("Incorrect mesh dimensions")
        return None
    dx = 1 / (m - 1)
    dy = 1 / (p - 1)
    dt = CFL * (dx**2 * dy**2) / (2 * D * (dx**2 + dy**2))  # Condition CFL
    lx = D * dt / dx**2
    ly = D * dt / dy**2

    x = np.linspace(0, 1.0, m)
    y = np.linspace(0, 1.0, p)
    X, Y = np.meshgrid(x, y, indexing="ij")
    Uh = u0(X, Y)
    vmin = np.min(Uh)
    vmax = np.max(Uh)
    t = 0
    it = 0
    Uh_history = [(t, Uh.copy())]

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

        # Store every iteration
        if it % every == 0:
            Uh_history.append((t, Uh.copy()))

    return Uh, Uh_history


def plot_progression_with_controls(Uh_history, X, Y, plot_type='3d'):
    fig = plt.figure(figsize=(16, 8))
    vmax = np.max(Uh_history[0][1])
    zlim = vmax  # Set z-axis limit to a fixed maximum value
    
    ax = None
    cbar = None
    if plot_type == '3d':
        ax = fig.add_subplot(1, 1, 1, projection='3d')
    elif plot_type == 'height':
        ax = fig.add_subplot(1, 1, 1)

    # Set initial plot
    current_index = 0  # Integer to track the current iteration
    playing = [False]  # Mutable flag to control play/stop
    time, Uh = Uh_history[current_index]
    if plot_type == '3d':
        surf = ax.plot_surface(X, Y, Uh, cmap='hot', vmin=0, vmax=vmax)
        ax.set_zlim(0, zlim)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Temperature')
        ax.set_title(f'3D Temperature Distribution at t = {time:.2f} s')
        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    elif plot_type == 'height':
        cp = ax.contourf(X, Y, Uh, levels=50, cmap='hot', vmin=0, vmax=vmax)
        cbar = fig.colorbar(cp, ax=ax)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'2D Temperature Distribution at t = {time:.2f} s')

    # Function to update the plot
    def update_plot():
        time, Uh = Uh_history[current_index]
        ax.cla()  # Clear the previous plot
        if plot_type == '3d':
            surf = ax.plot_surface(X, Y, Uh, cmap='hot', vmin=0, vmax=vmax)
            ax.set_zlim(0, zlim)  # Keep z-axis limit fixed
            ax.set_xlabel('X (m)')
            ax.set_ylabel('Y (m)')
            ax.set_zlabel('Temperature')
            ax.set_title(f'3D Temperature Distribution at t = {time:.2f} s')
            cbar.update_normal(surf)
        elif plot_type == 'height':
            cp = ax.contourf(X, Y, Uh, levels=50, cmap='hot', vmin=0, vmax=vmax)
            cbar.update_normal(cp)
            ax.set_xlabel('X (m)')
            ax.set_ylabel('Y (m)')
            ax.set_title(f'2D Temperature Distribution at t = {time:.2f} s')
        plt.draw()

    # Button to go to the next iteration
    def next(event):
        nonlocal current_index
        if current_index < len(Uh_history) - 1:
            current_index += 1
            update_plot()

    # Button to go to the previous iteration
    def prev(event):
        nonlocal current_index
        if current_index > 0:
            current_index -= 1
            update_plot()

    # Function to play through the iterations
    def play():
        nonlocal current_index
        while playing[0]:
            if current_index < len(Uh_history) - 1:
                current_index += 1
            else:
                current_index = 0  # Loop back to start
            update_plot()
            pytime.sleep(0.2)  # Pause between frames

    # Button to play/stop the animation
    def play_stop(event):
        if not playing[0]:
            playing[0] = True
            thread = threading.Thread(target=play)
            thread.start()
        else:
            playing[0] = False

    # Adding buttons to the figure
    axprev = plt.axes([0.6, 0.01, 0.1, 0.075])
    axnext = plt.axes([0.71, 0.01, 0.1, 0.075])
    axplay = plt.axes([0.82, 0.01, 0.1, 0.075])
    bnext = Button(axnext, 'Next')
    bnext.on_clicked(next)
    bprev = Button(axprev, 'Previous')
    bprev.on_clicked(prev)
    bplay = Button(axplay, 'Play/Stop')
    bplay.on_clicked(play_stop)

    plt.show()


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
T = 300
D = 2.2e-5  # Thermal diffusivity of air in m^2/s
mesh_dimensions = (100, 100)
CFL = 0.45
final_temperature_distribution, Uh_history = schema_chaleur2D(
    checkerboard, mesh_dimensions, T=T, D=D, CFL=CFL, every=10
)

# Generate the grid for plotting
m, p = mesh_dimensions
x = np.linspace(0, 1, m)
y = np.linspace(0, 1, p)
X, Y = np.meshgrid(x, y)

# Plot the progression of temperature distributions with interactive controls
plot_progression_with_controls(Uh_history, X, Y, plot_type='height')

# %%
