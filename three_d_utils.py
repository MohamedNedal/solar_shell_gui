import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import tempfile
import webbrowser
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import colors

# Assume the viewer object provides .maps, .current_index, .ellipse_center, .show_shell_checkbox, etc.

def create_3d_ellipsoid(viewer, ellipse_params, sunpy_map, show_shell=True, show_normals=False, show_radials=False, n_lat_radial=10, n_lon_radial=30, show_field_lines=False):
    """
    Create a 3D interactive ellipsoid plot using Plotly.
    Returns a Plotly Figure object.
    """
    theta, phi = np.mgrid[0:np.pi:50j, 0:2*np.pi:50j]
    x_sun = np.sin(theta) * np.cos(phi)
    y_sun = np.sin(theta) * np.sin(phi)
    z_sun = np.cos(theta)

    fig = go.Figure()
    fig.add_trace(go.Surface(x=x_sun, y=y_sun, z=z_sun, colorscale=[[0,'orange'],[1,'orange']], showscale=False, name='Sun', hoverinfo='skip'))

    x0, y0 = ellipse_params['x0'], ellipse_params['y0']
    a_radius, b_radius = ellipse_params['a'], ellipse_params['b']

    hpc_coord = SkyCoord(x0*u.arcsec, y0*u.arcsec,
                         frame='helioprojective',
                         observer=sunpy_map.observer_coordinate,
                         obstime=sunpy_map.date)
    hgs_coord = hpc_coord.transform_to('heliographic_stonyhurst')
    shell_lon_rad = np.deg2rad(hgs_coord.lon.value)
    shell_lat_rad = np.deg2rad(hgs_coord.lat.value)
    xshift = np.cos(shell_lat_rad) * np.cos(shell_lon_rad)
    yshift = np.cos(shell_lat_rad) * np.sin(shell_lon_rad)
    zshift = np.sin(shell_lat_rad)

    b_major = a_radius / sunpy_map.rsun_obs.value
    b_minor = b_radius / sunpy_map.rsun_obs.value
    shell_mesh_res = 50j
    theta_src, phi_src = np.mgrid[0:np.pi:shell_mesh_res, 0:2*np.pi:shell_mesh_res]
    x_src = xshift + b_minor * np.sin(theta_src) * np.cos(phi_src)
    y_src = yshift + b_minor * np.sin(theta_src) * np.sin(phi_src)
    z_src = zshift + b_major * np.cos(theta_src)
    r_shell = np.sqrt(x_src**2 + y_src**2 + z_src**2)
    mask = r_shell > 1
    x_display = np.copy(x_src)
    y_display = np.copy(y_src)
    z_display = np.copy(z_src)
    x_display[~mask] = np.nan
    y_display[~mask] = np.nan
    z_display[~mask] = np.nan

    if show_shell:
        fig.add_trace(go.Surface(x=x_display, y=y_display, z=z_display,
                                 colorscale='Viridis', cmin=0, cmax=90, opacity=1, showscale=True,
                                 name='Shell'))
    fig.update_layout(scene=dict(xaxis=dict(title='X (Solar Radii)'),
                                 yaxis=dict(title='Y (Solar Radii)'),
                                 zaxis=dict(title='Z (Solar Radii)'),
                                 aspectmode='cube'))
    return fig


def show_3d_ellipsoid_plot(viewer):
    """
    Generate and show the 3D ellipsoid in a browser.
    """
    if viewer.ellipse_center is None or not viewer.maps:
        from PyQt5.QtWidgets import QMessageBox
        QMessageBox.warning(viewer, 'Warning', 'Please draw an ellipse and load a map first.')
        return

    x0, y0 = viewer.ellipse_center
    a_radius = viewer.a_slider.value()
    b_radius = viewer.b_slider.value()
    ellipse_params = {'x0': x0, 'y0': y0, 'a': a_radius, 'b': b_radius}

    current_map = viewer.maps[viewer.current_index]

    fig = create_3d_ellipsoid(
        viewer, ellipse_params, current_map,
        show_shell=viewer.show_shell_checkbox.isChecked(),
        show_normals=viewer.show_normals_checkbox.isChecked(),
        show_radials=viewer.show_radials_checkbox.isChecked(),
        n_lat_radial=viewer.n_lat_radial_slider.value(),
        n_lon_radial=viewer.n_lon_radial_slider.value(),
        show_field_lines=viewer.show_field_lines_checkbox.isChecked()
    )

    # Save to temporary HTML and open
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as f:
        temp_filepath = f.name
        pio.write_html(fig, file=f, auto_open=False, include_plotlyjs='cdn')

    webbrowser.open(f'file://{temp_filepath}')


def export_3d_png_silent(viewer):
    """
    Export the current 3D ellipsoid to a PNG silently using Kaleido.
    """
    if viewer.ellipse_center is None or not viewer.maps:
        viewer.label.setText('Draw an ellipse and load a map first.')
        return

    try:
        x0, y0 = viewer.ellipse_center
        a_radius = viewer.a_slider.value()
        b_radius = viewer.b_slider.value()
        ellipse_params = {'x0': x0, 'y0': y0, 'a': a_radius, 'b': b_radius}

        current_map = viewer.maps[viewer.current_index]

        fig = create_3d_ellipsoid(
            viewer, ellipse_params, current_map,
            show_shell=viewer.show_shell_checkbox.isChecked(),
            show_normals=viewer.show_normals_checkbox.isChecked(),
            show_radials=viewer.show_radials_checkbox.isChecked(),
            n_lat_radial=viewer.n_lat_radial_slider.value(),
            n_lon_radial=viewer.n_lon_radial_slider.value(),
            show_field_lines=viewer.show_field_lines_checkbox.isChecked()
        )

        # Export folder
        script_dir = os.path.dirname(os.path.abspath(__file__))
        out_png = os.path.join(script_dir, f'frame_{viewer.current_index+1:03d}.png')

        # Requires kaleido
        import plotly.io as pio
        pio.write_image(fig, out_png, width=1600, height=1200, scale=2)

        viewer.label.setText(f'3D PNG saved silently to {out_png}')
    except Exception as e:
        viewer.label.setText(f'PNG export failed: {e}')