"""
3D visualization module for solar wave analysis
Creates 3D ellipsoid models and theta angle distributions
"""

import numpy as np
import plotly.graph_objects as go
import plotly.offline as pyo
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.constants import R_sun


class Visualization3D:
    """Class for creating 3D visualizations"""
    
    @staticmethod
    def create_3d_ellipsoid(ellipse_params, sunpy_map, show_shell=True, 
                          show_normals=False, opacity=0.7):
        """
        Create 3D ellipsoid visualization
        
        Parameters:
        -----------
        ellipse_params : dict
            Ellipse parameters from fitting
        sunpy_map : sunpy.map.Map
            The solar map for coordinate transformation
        show_shell : bool
            Whether to show the ellipsoid shell
        show_normals : bool
            Whether to show normal vectors
        opacity : float
            Opacity of the ellipsoid surface
            
        Returns:
        --------
        plotly.graph_objects.Figure : 3D visualization figure
        """
        
        fig = go.Figure()
        
        # Create the Sun (sphere)
        sun_data = Visualization3D._create_sun_sphere()
        fig.add_trace(go.Surface(
            x=sun_data['x'], y=sun_data['y'], z=sun_data['z'],
            colorscale=[[0, 'orange'], [1, 'yellow']],
            showscale=False,
            name='Sun',
            hoverinfo='skip',
            opacity=0.8
        ))
        
        # Create ellipsoid shell
        if show_shell:
            ellipsoid_data = Visualization3D._create_ellipsoid_shell(
                ellipse_params, sunpy_map
            )
            
            if ellipsoid_data is not None:
                fig.add_trace(go.Surface(
                    x=ellipsoid_data['x'], 
                    y=ellipsoid_data['y'], 
                    z=ellipsoid_data['z'],
                    colorscale='Viridis',
                    opacity=opacity,
                    name='Wave Shell',
                    hovertemplate='<b>Wave Shell</b><br>' +
                                'X: %{x:.2f}<br>' +
                                'Y: %{y:.2f}<br>' +
                                'Z: %{z:.2f}<extra></extra>'
                ))
        
        # Add coordinate system arrows
        Visualization3D._add_coordinate_arrows(fig)
        
        # Layout settings
        fig.update_layout(
            scene=dict(
                xaxis=dict(
                    title='X (Solar Radii)',
                    range=[-2.5, 2.5],
                    showgrid=True,
                    gridcolor='lightgray'
                ),
                yaxis=dict(
                    title='Y (Solar Radii)',
                    range=[-2.5, 2.5],
                    showgrid=True,
                    gridcolor='lightgray'
                ),
                zaxis=dict(
                    title='Z (Solar Radii)',
                    range=[-2.5, 2.5],
                    showgrid=True,
                    gridcolor='lightgray'
                ),
                aspectmode='cube',
                camera=dict(
                    eye=dict(x=1.8, y=1.8, z=1.8),
                    center=dict(x=0, y=0, z=0),
                    up=dict(x=0, y=0, z=1)
                ),
                bgcolor='black'
            ),
            width=1024,
            height=768,
            title=dict(
                text=f'3D Solar Wave Model - {sunpy_map.date.strftime("%Y-%m-%d %H:%M:%S")}',
                x=0.5,
                font=dict(size=16)
            ),
            showlegend=True,
            paper_bgcolor='white',
            plot_bgcolor='white'
        )
        
        return fig
    
    @staticmethod
    def _create_sun_sphere(resolution=50):
        """Create a sphere representing the Sun"""
        theta, phi = np.mgrid[0:np.pi:resolution*1j, 0:2*np.pi:resolution*1j]
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        
        return {'x': x, 'y': y, 'z': z}
    
    @staticmethod
    def _create_ellipsoid_shell(ellipse_params, sunpy_map, resolution=50):
        """Create ellipsoid shell data"""
        try:
            # Extract parameters (convert from arcsec to solar radii)
            x0 = ellipse_params['x0']  # arcsec
            y0 = ellipse_params['y0']  # arcsec
            a = ellipse_params['a']    # arcsec (semi-major axis)
            b = ellipse_params['b']    # arcsec (semi-minor axis)
            angle = ellipse_params['angle']  # degrees
            
            # Convert arcsec to solar radii
            # 1 solar radius ≈ 960 arcsec (approximate)
            arcsec_per_rsun = 960.0
            
            x0_rsun = x0 / arcsec_per_rsun
            y0_rsun = y0 / arcsec_per_rsun
            a_rsun = a / arcsec_per_rsun
            b_rsun = b / arcsec_per_rsun
            
            # Assume the ellipsoid extends from the surface (1 R_sun) outward
            # The height in the z-direction is estimated based on the ellipse parameters
            c_rsun = np.sqrt(a_rsun**2 + b_rsun**2) / 2  # Height estimate
            
            # Create parametric surface for ellipsoid
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            u, v = np.meshgrid(u, v)
            
            # Convert angle to radians
            angle_rad = np.radians(angle)
            
            # Parametric ellipsoid equations
            x_base = a_rsun * np.sin(v) * np.cos(u)
            y_base = b_rsun * np.sin(v) * np.sin(u)
            z_base = c_rsun * np.cos(v)
            
            # Apply rotation around z-axis
            x_rot = x_base * np.cos(angle_rad) - y_base * np.sin(angle_rad)
            y_rot = x_base * np.sin(angle_rad) + y_base * np.cos(angle_rad)
            z_rot = z_base
            
            # Translate to center position and offset from Sun surface
            x_final = x_rot + x0_rsun
            y_final = y_rot + y0_rsun
            z_final = z_rot + 1.0 + c_rsun  # Start from Sun surface (1 R_sun) + height
            
            return {'x': x_final, 'y': y_final, 'z': z_final}
            
        except Exception as e:
            print(f"Error creating ellipsoid shell: {e}")
            return None
    
    @staticmethod
    def _add_coordinate_arrows(fig):
        """Add coordinate system arrows to the figure"""
        # X-axis arrow (red)
        fig.add_trace(go.Scatter3d(
            x=[0, 2], y=[0, 0], z=[0, 0],
            mode='lines+markers',
            line=dict(color='red', width=6),
            marker=dict(size=8, color='red'),
            name='X-axis',
            showlegend=False
        ))
        
        # Y-axis arrow (green)
        fig.add_trace(go.Scatter3d(
            x=[0, 0], y=[0, 2], z=[0, 0],
            mode='lines+markers',
            line=dict(color='green', width=6),
            marker=dict(size=8, color='green'),
            name='Y-axis',
            showlegend=False
        ))
        
        # Z-axis arrow (blue)
        fig.add_trace(go.Scatter3d(
            x=[0, 0], y=[0, 0], z=[0, 2],
            mode='lines+markers',
            line=dict(color='blue', width=6),
            marker=dict(size=8, color='blue'),
            name='Z-axis',
            showlegend=False
        ))
    
    @staticmethod
    def create_theta_distribution(ellipse_params_list, title="Wave Propagation Analysis"):
        """
        Create theta angle distribution plot
        
        Parameters:
        -----------
        ellipse_params_list : list
            List of ellipse parameters from multiple time steps
        title : str
            Title for the plot
            
        Returns:
        --------
        plotly.graph_objects.Figure : Theta distribution figure
        """
        if not ellipse_params_list:
            return None
            
        # Extract timestamps and angles
        timestamps = []
        angles = []
        semi_major = []
        semi_minor = []
        
        for params in ellipse_params_list:
            if 'timestamp' in params and params['timestamp'] is not None:
                timestamps.append(params['timestamp'])
                angles.append(params['angle'])
                semi_major.append(params['a'])
                semi_minor.append(params['b'])
        
        if not timestamps:
            return None
            
        fig = go.Figure()
        
        # Add angle evolution
        fig.add_trace(go.Scatter(
            x=list(range(len(timestamps))),
            y=angles,
            mode='lines+markers',
            name='Propagation Angle',
            line=dict(color='blue', width=2),
            marker=dict(size=8),
            hovertemplate='<b>Time Step:</b> %{x}<br>' +
                         '<b>Angle:</b> %{y:.1f}°<br>' +
                         '<extra></extra>'
        ))
        
        # Add semi-major axis evolution (secondary y-axis)
        fig.add_trace(go.Scatter(
            x=list(range(len(timestamps))),
            y=semi_major,
            mode='lines+markers',
            name='Semi-major Axis',
            line=dict(color='red', width=2),
            marker=dict(size=6),
            yaxis='y2',
            hovertemplate='<b>Time Step:</b> %{x}<br>' +
                         '<b>Semi-major:</b> %{y:.1f} arcsec<br>' +
                         '<extra></extra>'
        ))
        
        # Layout
        fig.update_layout(
            title=dict(
                text=title,
                x=0.5,
                font=dict(size=16)
            ),
            xaxis=dict(
                title='Time Step',
                showgrid=True,
                gridcolor='lightgray'
            ),
            yaxis=dict(
                title='Angle (degrees)',
                side='left',
                showgrid=True,
                gridcolor='lightgray'
            ),
            yaxis2=dict(
                title='Semi-major Axis (arcsec)',
                side='right',
                overlaying='y',
                showgrid=False
            ),
            width=1024,
            height=500,
            hovermode='x unified',
            showlegend=True,
            paper_bgcolor='white',
            plot_bgcolor='white'
        )
        
        return fig
    
    @staticmethod
    def create_wave_evolution_plot(ellipse_params_list, title="Wave Evolution"):
        """
        Create a comprehensive wave evolution plot
        
        Parameters:
        -----------
        ellipse_params_list : list
            List of ellipse parameters from multiple time steps
        title : str
            Title for the plot
            
        Returns:
        --------
        plotly.graph_objects.Figure : Wave evolution figure
        """
        if not ellipse_params_list:
            return None
            
        # Extract data
        timestamps = []
        angles = []
        semi_major = []
        semi_minor = []
        areas = []
        
        for params in ellipse_params_list:
            if 'timestamp' in params and params['timestamp'] is not None:
                timestamps.append(params['timestamp'])
                angles.append(params['angle'])
                a = params['a']
                b = params['b']
                semi_major.append(a)
                semi_minor.append(b)
                # Calculate ellipse area
                areas.append(np.pi * a * b)
        
        if not timestamps:
            return None
            
        time_steps = list(range(len(timestamps)))
        
        # Create subplots
        from plotly.subplots import make_subplots
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Propagation Angle', 'Semi-major Axis', 
                          'Semi-minor Axis', 'Wave Area'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Angle plot
        fig.add_trace(
            go.Scatter(x=time_steps, y=angles, mode='lines+markers',
                      name='Angle', line=dict(color='blue')),
            row=1, col=1
        )
        
        # Semi-major axis plot
        fig.add_trace(
            go.Scatter(x=time_steps, y=semi_major, mode='lines+markers',
                      name='Semi-major', line=dict(color='red')),
            row=1, col=2
        )
        
        # Semi-minor axis plot
        fig.add_trace(
            go.Scatter(x=time_steps, y=semi_minor, mode='lines+markers',
                      name='Semi-minor', line=dict(color='green')),
            row=2, col=1
        )
        
        # Area plot
        fig.add_trace(
            go.Scatter(x=time_steps, y=areas, mode='lines+markers',
                      name='Area', line=dict(color='orange')),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            title_text=title,
            showlegend=False,
            width=1200,
            height=800
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Time Step", row=2, col=1)
        fig.update_xaxes(title_text="Time Step", row=2, col=2)
        fig.update_yaxes(title_text="Degrees", row=1, col=1)
        fig.update_yaxes(title_text="Arcsec", row=1, col=2)
        fig.update_yaxes(title_text="Arcsec", row=2, col=1)
        fig.update_yaxes(title_text="Arcsec²", row=2, col=2)
        
        return fig
    
    @staticmethod
    def save_figure(fig, filename, format='html'):
        """
        Save plotly figure to file
        
        Parameters:
        -----------
        fig : plotly.graph_objects.Figure
            The figure to save
        filename : str
            Output filename
        format : str
            Output format ('html', 'png', 'pdf', etc.)
        """
        try:
            if format.lower() == 'html':
                pyo.plot(fig, filename=filename, auto_open=False)
            else:
                fig.write_image(filename, format=format)
            print(f"Figure saved to {filename}")
        except Exception as e:
            print(f"Error saving figure: {e}")