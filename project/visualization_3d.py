"""
visualization_3d.py
-------------------
Standalone function for building an interactive Plotly 3-D figure of the
solar CME ellipsoid shell with optional magnetic field lines.  No Qt dependency.
"""

import numpy as np
import plotly.graph_objects as go
from astropy.coordinates import SkyCoord
from sunpy.sun import constants as const
import astropy.units as u


def create_3d_ellipsoid(
    ellipse_params,
    sunpy_map,
    field_lines=None,
    show_shell=True,
    show_normals=False,
    show_radials=False,
    n_lat_radial=10,
    n_lon_radial=30,
    show_field_lines=False,
):
    """
    Build an interactive Plotly 3-D figure of the solar CME ellipsoid shell.

    The 2-D ellipse (in arcsec) is projected into heliographic coordinates
    and extruded into a 3-D ellipsoid.  The shell surface is coloured by the
    acute angle (theta) between each outward surface normal and the nearest
    magnetic field-line segment.

    Parameters
    ----------
    ellipse_params : dict
        Must contain ``'x0'``, ``'y0'`` (centre in arcsec) and
        ``'a'``, ``'b'`` (semi-axes in arcsec).
    sunpy_map : sunpy.map.GenericMap
        Reference AIA map supplying observer coordinates and ``rsun_obs``.
    field_lines : pfsspy.fieldline.FieldLines or None
        Traced magnetic field lines from a PFSS model.  Required when
        ``show_field_lines=True``.
    show_shell : bool
        Render the coloured ellipsoid surface mesh.
    show_normals : bool
        Draw a sub-sampled set of outward surface normal vectors.
    show_radials : bool
        Draw radial lines from the solar origin out to 2 R_sun.
    n_lat_radial : int
        Number of latitude divisions for the radial-line grid.
    n_lon_radial : int
        Number of longitude divisions for the radial-line grid.
    show_field_lines : bool
        Overlay magnetic field lines colour-coded by polarity
        (red = positive, blue = negative, black = open/neutral).

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Configured 3-D interactive Plotly figure.
    theta_angles : numpy.ndarray
        1-D array of theta angles (degrees) for each visible shell point,
        representing the angle between the local outward normal and the
        nearest field-line segment.
    """
    # --- Sun sphere ---
    t_sph, p_sph = np.mgrid[0:np.pi:50j, 0:2 * np.pi:50j]
    x_sun = np.sin(t_sph) * np.cos(p_sph)
    y_sun = np.sin(t_sph) * np.sin(p_sph)
    z_sun = np.cos(t_sph)

    fig = go.Figure()
    fig.add_trace(go.Surface(
        x=x_sun, y=y_sun, z=z_sun,
        colorscale=[[0, 'orange'], [1, 'orange']],
        showscale=False, name='Sun', hoverinfo='skip',
    ))

    # --- Project ellipse centre to 3-D heliographic frame ---
    x0, y0 = ellipse_params['x0'], ellipse_params['y0']
    a_radius, b_radius = ellipse_params['a'], ellipse_params['b']

    hpc_coord = SkyCoord(
        x0 * u.arcsec, y0 * u.arcsec,
        frame='helioprojective',
        observer=sunpy_map.observer_coordinate,
        obstime=sunpy_map.date,
    )
    hgs_coord = hpc_coord.transform_to('heliographic_stonyhurst')
    shell_lon_rad = np.deg2rad(hgs_coord.lon.value)
    shell_lat_rad = np.deg2rad(hgs_coord.lat.value)

    xshift = np.cos(shell_lat_rad) * np.cos(shell_lon_rad)
    yshift = np.cos(shell_lat_rad) * np.sin(shell_lon_rad)
    zshift = np.sin(shell_lat_rad)

    b_major = a_radius / sunpy_map.rsun_obs.value   # along z-axis (line of sight)
    b_minor = b_radius / sunpy_map.rsun_obs.value   # transverse axes

    # --- Ellipsoid mesh ---
    t_src, p_src = np.mgrid[0:np.pi:50j, 0:2 * np.pi:50j]
    x_src = xshift + b_minor * np.sin(t_src) * np.cos(p_src)
    y_src = yshift + b_minor * np.sin(t_src) * np.sin(p_src)
    z_src = zshift + b_major * np.cos(t_src)

    r_shell = np.sqrt(x_src ** 2 + y_src ** 2 + z_src ** 2)
    mask = r_shell > 1  # only keep points above the photosphere

    x_display = np.where(mask, x_src, np.nan)
    y_display = np.where(mask, y_src, np.nan)
    z_display = np.where(mask, z_src, np.nan)

    # --- Outward surface normals (analytic gradient of the ellipsoid) ---
    nx_all = 2 * (x_src - xshift) / (b_minor ** 2)
    ny_all = 2 * (y_src - yshift) / (b_minor ** 2)
    nz_all = 2 * (z_src - zshift) / (b_major ** 2)
    norm_mag = np.sqrt(nx_all ** 2 + ny_all ** 2 + nz_all ** 2)
    nx_all /= norm_mag
    ny_all /= norm_mag
    nz_all /= norm_mag

    x_outer = x_src[mask].flatten()
    y_outer = y_src[mask].flatten()
    z_outer = z_src[mask].flatten()
    nx_flat = nx_all[mask].flatten()
    ny_flat = ny_all[mask].flatten()
    nz_flat = nz_all[mask].flatten()

    # --- Seed points and directions used for theta computation ---
    all_points_for_theta = []
    all_dirs_for_theta = []

    # --- Optional radial lines ---
    if show_radials:
        radial_length = 2.0
        t_rad, p_rad = np.mgrid[
            0:np.pi:complex(n_lat_radial),
            0:2 * np.pi:complex(n_lon_radial),
        ]
        rdx = (np.sin(t_rad) * np.cos(p_rad)).flatten()
        rdy = (np.sin(t_rad) * np.sin(p_rad)).flatten()
        rdz = np.cos(t_rad).flatten()

        for i in range(len(rdx)):
            r_vec = np.array([rdx[i], rdy[i], rdz[i]])
            r_pts = np.linspace(0, radial_length, 20)
            fig.add_trace(go.Scatter3d(
                x=r_pts * r_vec[0], y=r_pts * r_vec[1], z=r_pts * r_vec[2],
                mode='lines', line=dict(color='gray', width=1),
                showlegend=False, hoverinfo='skip', opacity=0.3,
            ))
            if not show_field_lines or field_lines is None:
                all_points_for_theta.append(r_vec * radial_length)
                all_dirs_for_theta.append(r_vec)

    # --- Optional magnetic field lines ---
    if show_field_lines and field_lines:
        all_points_for_theta = []
        all_dirs_for_theta = []
        for field_line in field_lines:
            color = {0: 'black', -1: 'blue', 1: 'red'}.get(field_line.polarity, 'black')
            coords = field_line.coords
            coords.representation_type = 'cartesian'
            x_f = coords.x / const.radius
            y_f = coords.y / const.radius
            z_f = coords.z / const.radius

            fig.add_trace(go.Scatter3d(
                x=x_f, y=y_f, z=z_f,
                mode='lines', line=dict(color=color, width=2),
                showlegend=False, opacity=0.4,
            ))
            for i in range(len(x_f) - 1):
                pt = np.array([x_f[i].value, y_f[i].value, z_f[i].value])
                dv = np.array([
                    x_f[i + 1].value - x_f[i].value,
                    y_f[i + 1].value - y_f[i].value,
                    z_f[i + 1].value - z_f[i].value,
                ])
                norm = np.linalg.norm(dv)
                if norm > 0:
                    dv /= norm
                all_points_for_theta.append(pt)
                all_dirs_for_theta.append(dv)

    # --- Compute theta angle at each visible shell point ---
    all_points_for_theta = np.array(all_points_for_theta)
    all_dirs_for_theta = np.array(all_dirs_for_theta)
    theta_angles = np.full(len(x_outer), np.nan)

    if len(all_points_for_theta) > 0 and len(x_outer) > 0:
        for i in range(len(x_outer)):
            surf_pt = np.array([x_outer[i], y_outer[i], z_outer[i]])
            normal_vec = np.array([nx_flat[i], ny_flat[i], nz_flat[i]])
            dists = np.linalg.norm(all_points_for_theta - surf_pt, axis=1)
            idx = np.argmin(dists)
            cos_t = np.clip(np.dot(normal_vec, all_dirs_for_theta[idx]), -1.0, 1.0)
            theta_angles[i] = np.degrees(np.arccos(np.abs(cos_t)))

    theta_surface = np.full_like(x_src, np.nan)
    theta_surface[mask] = theta_angles

    # --- Shell surface ---
    if show_shell:
        fig.add_trace(go.Scatter3d(
            x=[xshift], y=[yshift], z=[zshift],
            mode='markers', marker=dict(size=8, color='black'),
            showlegend=False, name='Shell Center',
        ))
        fig.add_trace(go.Surface(
            x=x_display, y=y_display, z=z_display,
            surfacecolor=theta_surface,
            colorscale='Viridis', cmin=0, cmax=90, opacity=1,
            showscale=True,
            colorbar=dict(
                title=dict(text='Theta Angle (degrees)', side='right'),
                len=0.6,
            ),
            name='Shell',
        ))

    # --- Normal vector arrows (sub-sampled) ---
    if show_normals:
        step = max(1, len(x_outer) // 100)
        x_ns = x_outer[::step]
        y_ns = y_outer[::step]
        z_ns = z_outer[::step]
        nx_s = nx_all[mask][::step] * 0.2
        ny_s = ny_all[mask][::step] * 0.2
        nz_s = nz_all[mask][::step] * 0.2
        x_lines, y_lines, z_lines = [], [], []
        for i in range(len(x_ns)):
            x_lines.extend([x_ns[i], x_ns[i] + nx_s[i], None])
            y_lines.extend([y_ns[i], y_ns[i] + ny_s[i], None])
            z_lines.extend([z_ns[i], z_ns[i] + nz_s[i], None])
        fig.add_trace(go.Scatter3d(
            x=x_lines, y=y_lines, z=z_lines,
            mode='lines', line=dict(color='black', width=3),
            showlegend=False, name='Normal Vectors',
        ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-2, 2], title='X (Solar Radii)'),
            yaxis=dict(range=[-2, 2], title='Y (Solar Radii)'),
            zaxis=dict(range=[-2, 2], title='Z (Solar Radii)'),
            aspectmode='cube',
            camera=dict(eye=dict(x=0, y=1.5, z=0)),
        ),
        width=1024, height=768,
        title='Interactive 3D Solar Shell Visualization',
        showlegend=True,
    )
    return fig, theta_angles

