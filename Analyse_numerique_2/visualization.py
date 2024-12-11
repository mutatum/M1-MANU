from time import sleep
import numpy as np
from typing import Optional
import plotly.graph_objects as go
import webbrowser
from pathlib import Path
import tempfile

class HistoryManager:
    def __init__(self, save_frequency=None, max_snapshots=None):
        """
        Initialize history manager with controls for saving frequency.
        
        Args:
            save_frequency (float): Time interval between saves (None = save all steps)
            max_snapshots (int): Maximum number of snapshots to store (None = unlimited)
        """
        self.save_frequency = save_frequency
        self.max_snapshots = max_snapshots
        self.history = []
        self.last_save_time = 0.0

    def should_save(self, t):
        if self.save_frequency is None:
            return True
        return (t - self.last_save_time) >= self.save_frequency

    def add(self, t, U, force=True):
        if not self.should_save(t) or not force:
            return
        
        self.history.append((t, U.copy()))
        self.last_save_time = t
        
        if self.max_snapshots and len(self.history) > self.max_snapshots:
            self.history.pop(0)

    def get_history(self):
        return self.history


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        sleep(0.1)
        print()

import numpy as np
from typing import Union, Callable, Optional, List
import plotly.graph_objects as go

def _format_time(t: float) -> str:
    """
    Formats time `t` into the two most pertinent units (e.g., years, months).
    
    Args:
        t: Time value in seconds.
    
    Returns:
        Formatted string with the most relevant time units.
    """
    units = [
        ("years", 365*24*3600),
        ("months", 30*24*3600),
        ("days", 24*3600),
        ("hours", 3600),
        ("minutes", 60),
        ("seconds", 1),
        ("milliseconds", 1e-3),
        ("microseconds", 1e-6),
        ("nanoseconds", 1e-9)
    ]
    
    remaining = t
    result = []
    for unit, factor in units:
        value = int(remaining // factor)
        if value > 0 or len(result) > 0:  # Skip leading zeros
            result.append(f"{value} {unit}")
            remaining -= value * factor
        if len(result) == 2:  # Limit to two units
            break
    
    return " ".join(result)

def visualize_pde_evolution(history_manager, x_range: np.ndarray, y_range: Optional[np.ndarray] = None,
                            title: str = "PDE Evolution",
                            colorscale: str = 'Viridis',
                            mode: str = 'heightmap',
                            dynamic_scaling: bool = True) -> go.Figure:
    """
    Visualize PDE evolution in 1D or 2D with time control (play, pause, sliders).
    Supports multiple visualization modes: heightmap (3D) or heatmap (2D).
    
    Args:
        history_manager: Instance of HistoryManager containing saved states.
        x_range: Spatial coordinates for x.
        y_range: Spatial coordinates for y (None for 1D problems).
        title: Plot title.
        colorscale: Colormap for the visualization.
        mode: Visualization mode ('heightmap' or 'heatmap').
        dynamic_scaling: Adjust color scale dynamically per frame.
    
    Returns:
        A Plotly Figure object.
    """
    history = history_manager.get_history()
    t_range = np.array([t for t, _ in history])
    is_2d = y_range is not None

    initial_data = history[0][1]

    if dynamic_scaling:
        z_min, z_max = np.min(initial_data), np.max(initial_data)
    else:
        z_min, z_max = np.inf, -np.inf
        for _, data in history:
            z_min = min(z_min, np.min(data))
            z_max = max(z_max, np.max(data))

    if mode == 'heightmap' and is_2d:
        fig = go.Figure(
            data=[go.Surface(z=initial_data, x=x_range, y=y_range, colorscale=colorscale, cmin=z_min, cmax=z_max)]
        )
    elif is_2d:
        fig = go.Figure(
            data=[go.Heatmap(x=x_range, y=y_range, z=initial_data, colorscale=colorscale, zmin=z_min, zmax=z_max)]
        )
    else:
        fig = go.Figure(
            data=[go.Scatter(x=x_range, y=initial_data, mode='lines')]
        )

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y" if is_2d else "Value",
            zaxis_title="Value" if mode == 'heightmap' else None
        )
    )

    frames = []
    for t, data in history:
        formatted_time = _format_time(t)
        if dynamic_scaling:
            z_min, z_max = np.min(data), np.max(data)

        if mode == 'heightmap' and is_2d:
            frames.append(go.Frame(
                data=[go.Surface(z=data, x=x_range, y=y_range, colorscale=colorscale, cmin=z_min, cmax=z_max)],
                name=f"t={formatted_time}"
            ))
        elif is_2d:
            frames.append(go.Frame(
                data=[go.Heatmap(x=x_range, y=y_range, z=data, colorscale=colorscale, zmin=z_min, zmax=z_max)],
                name=f"t={formatted_time}"
            ))
        else:
            frames.append(go.Frame(
                data=[go.Scatter(x=x_range, y=data, mode='lines')],
                name=f"t={formatted_time}"
            ))

    fig.frames = frames

    steps = []
    for i, (t, _) in enumerate(history):
        step = dict(
            method="animate",
            args=[[f"t={_format_time(t)}"], 
                  dict(mode="immediate", frame=dict(duration=0, redraw=True), transition=dict(duration=0))],
            label=f"{_format_time(t)}"
        )
        steps.append(step)

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Time: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                x=0,
                y=1.15,
                xanchor="left",
                yanchor="top",
                buttons=[
                    dict(label="Play", method="animate", 
                         args=[None, {"frame": {"duration": 100, "redraw": True},
                                      "fromcurrent": True, "transition": {"duration": 0}}]),
                    dict(label="Pause", method="animate", 
                         args=[[None], {"frame": {"duration": 0, "redraw": False},
                                        "mode": "immediate"}])
                ]
            )
        ],
        sliders=sliders
    )

    temp_dir = Path(tempfile.gettempdir())
    temp_html = temp_dir / "heat_equation_animation.html"
    
    # Save the figure to a temporary HTML file
    fig.write_html(temp_html)
    
    # Open the default web browser
    webbrowser.open(f'file://{temp_html}')
    
    return fig

def visualize_heat_equation(history_manager, x_range: np.ndarray, y_range: Optional[np.ndarray] = None,
                            mode: str = 'heatmap', dynamic_scaling: bool = True) -> go.Figure:
    """
    Specialized visualization for heat equation solutions (default: heatmap).
    
    Args:
        history_manager: Instance of HistoryManager containing saved states.
        x_range: Spatial coordinates for x.
        y_range: Spatial coordinates for y (None for 1D problems).
        mode: Visualization mode ('heatmap' or 'heightmap').
        dynamic_scaling: Adjust color scale dynamically per frame.
    
    Returns:
        A Plotly Figure object.
    """
    return visualize_pde_evolution(history_manager, x_range, y_range, 
                                   title="Heat Equation Evolution", 
                                   colorscale='RdBu_r', mode=mode, dynamic_scaling=dynamic_scaling)
