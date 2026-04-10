"""
pfss_model.py
-------------
Standalone function for computing a PFSS (Potential Field Source Surface)
magnetic field model from a GONG synoptic magnetogram.  No Qt dependency.
"""

import numpy as np
import sunpy.map
import pfsspy
import pfsspy.tracing as tracing
from astropy.coordinates import SkyCoord
from sunpy.sun import constants as const
import astropy.units as u


def calculate_pfss_model(gong_filepath, progress_callback=None):
    """
    Load a GONG magnetogram, solve the PFSS model, and trace field lines.

    The source surface is placed at 3 solar radii with 50 radial grid points.
    Field lines are seeded on a 40 × 60 uniform lat/lon grid at 1.05 R_sun.

    Parameters
    ----------
    gong_filepath : str
        Path to a GONG FITS synoptic magnetogram file.
    progress_callback : callable, optional
        Called with a plain status string at key stages, e.g.
        ``progress_callback('Tracing field lines...')``.

    Returns
    -------
    gong_map : sunpy.map.GenericMap
        The loaded (and metadata-patched) GONG map.
    pfss_out : pfsspy.Output
        The computed PFSS solution object.
    field_lines : pfsspy.fieldline.FieldLines
        Traced magnetic field lines from the uniform seed grid.

    Raises
    ------
    Exception
        Re-raises any exception encountered during loading, solving, or tracing.
    """
    def _notify(msg):
        if progress_callback:
            progress_callback(msg)
        print(msg)

    _notify('Loading GONG map...')
    gong_map = sunpy.map.Map(gong_filepath)
    if 'cunit1' not in gong_map.meta:
        gong_map.meta['cunit1'] = u.deg

    nrho = 50
    rss = 3
    pfss_in = pfsspy.Input(gong_map, nrho, rss)

    _notify('Calculating PFSS model...')
    pfss_out = pfsspy.pfss(pfss_in)

    num_footpoints_lat = 40
    num_footpoints_lon = 60
    r_trace = 1.05 * const.radius

    lat = np.linspace(np.radians(-90), np.radians(90), num_footpoints_lat, endpoint=False)
    lon = np.linspace(np.radians(-180), np.radians(180), num_footpoints_lon, endpoint=False)
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')
    lat_flat = lat_grid.ravel() * u.rad
    lon_flat = lon_grid.ravel() * u.rad

    seeds = SkyCoord(lon_flat, lat_flat, r_trace, frame=pfss_out.coordinate_frame)
    tracer = tracing.FortranTracer()

    _notify('Tracing magnetic field lines...')
    field_lines = tracer.trace(seeds, pfss_out)

    return gong_map, pfss_out, field_lines

#J-Plot Code as Viewer.py has too much code already
#def J_plot(viewer):
