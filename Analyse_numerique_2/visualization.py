import plotly.graph_objects as go
import webbrowser
import tempfile
from pathlib import Path
import numpy as np
from time import sleep

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        sleep(0.1)
        print()

def create_heat_visualization(history, X, Y, title="Heat Distribution"):
    # Calculate temperature ranges
    temp_min = min(np.min(state) for _, state in history)
    temp_max = max(np.max(state) for _, state in history)
    
    # Create figure with both surface and heatmap traces
    fig = go.Figure()
    
    surface_data = go.Surface(
        x=X, y=Y, z=history[0][1],
        colorscale='RdBu_r',
        cmin=temp_min,
        cmax=temp_max,
        showscale=True
    )
    
    heatmap_data = go.Heatmap(
        x=X[0,:] if X.ndim == 2 else X,
        y=Y[:,0] if Y.ndim == 2 else Y,
        z=history[0][1],
        colorscale='RdBu_r',
        zmin=temp_min,
        zmax=temp_max,
        showscale=True
    )
    
    fig.add_trace(surface_data)
    fig.add_trace(heatmap_data)
    fig.data[1].visible = False  # Hide heatmap initially
    
    # Create frames for both views
    frames = []
    for t, U in history:
        frame = go.Frame(
            data=[
                go.Surface(x=X, y=Y, z=U, colorscale='RdBu_r', cmin=temp_min, cmax=temp_max),
                go.Heatmap(x=X[0,:] if X.ndim == 2 else X, y=Y[:,0] if Y.ndim == 2 else Y, z=U, colorscale='RdBu_r', zmin=temp_min, zmax=temp_max)
            ],
            name=f't={t:.2f}s'
        )
        frames.append(frame)
    fig.frames = frames
    
    # Update layout with play/pause and view toggle buttons
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X (m)',
            yaxis_title='Y (m)',
            zaxis_title='Temperature (°C)',
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        updatemenus=[
            dict(
                type='buttons',
                showactive=False,
                buttons=[
                    dict(label='Play',
                         method='animate',
                         args=[None, {'frame': {'duration': 50, 'redraw': True},
                                      'fromcurrent': True,
                                      'mode': 'immediate'}]),
                    dict(label='Pause',
                         method='animate',
                         args=[[None], {'frame': {'duration': 0, 'redraw': False},
                                        'mode': 'immediate'}])
                ],
                x=0.1,
                y=0,
                direction='left'
            ),
            dict(
                type='buttons',
                buttons=[
                    dict(label='3D Surface',
                         method='update',
                         args=[{'visible': [True, False]},
                               {'scene': {'visible': True}}]),
                    dict(label='2D Heatmap',
                         method='update',
                         args=[{'visible': [False, True]},
                               {'scene': {'visible': False}}])
                ],
                direction='down',
                x=0,
                y=1
            )
        ],
        sliders=[{
            'active': 0,
            'currentvalue': {'prefix': 'Time: ', 'suffix': 's'},
            'steps': [
                {
                    'method': 'animate',
                    'label': f'{t:.1f}s',
                    'args': [[f't={t:.2f}s'], {
                        'frame': {'duration': 0, 'redraw': True},
                        'mode': 'immediate'
                    }]
                }
                for t, _ in history
            ]
        }]
    )
    
    return fig

def plot_interactive(history, X, Y):
    try:
        # Create visualization
        fig = create_heat_visualization(history, X, Y)
        
        # Create temp file
        temp_dir = Path(tempfile.gettempdir())
        html_path = temp_dir / "heat_visualization.html"
        
        # Save and open
        fig.write_html(str(html_path), auto_open=False)
        webbrowser.open(f'file://{html_path}')
        
    except Exception as e:
        print(f"Error displaying visualization: {e}")
        import traceback
        traceback.print_exc()  # Print full error trace
