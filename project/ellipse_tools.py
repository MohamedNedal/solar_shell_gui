"""
ellipse_tools.py
----------------
Standalone functions for ellipse parameter extraction, unit conversion,
and CSV/text file export.  No Qt dependency.
"""

import os
import csv
from datetime import datetime

import numpy as np
import astropy.units as u


def init_fit_results_csv(filepath):
    """
    Create the CSV log file with a header row if it does not already exist.

    Parameters
    ----------
    filepath : str
        Full path to the CSV file to initialise.
    """
    if not os.path.exists(filepath):
        with open(filepath, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow([
                'timestamp', 'filename', 'x0_arcsec', 'y0_arcsec',
                'a_arcsec', 'b_arcsec', 'a_rsun', 'b_rsun',
                'arcsec_per_pixel', 'notes',
            ])


def extract_ellipse_params(x0, y0, a, b, current_map=None, current_filename='N/A'):
    """
    Compute ellipse parameters in both arcsec and solar-radii units.

    Attempts to derive the arcsec-per-pixel scale and the solar-radius
    conversion factor from the supplied SunPy map.  Falls back to a nominal
    value of 950 arcsec/R_sun when map metadata is unavailable.

    Parameters
    ----------
    x0 : float
        Ellipse centre X coordinate in arcsec.
    y0 : float
        Ellipse centre Y coordinate in arcsec.
    a : float
        Semi-major axis in arcsec.
    b : float
        Semi-minor axis in arcsec.
    current_map : sunpy.map.GenericMap, optional
        The active map, used to derive scale metadata.
    current_filename : str, optional
        Basename of the current FITS file (for display and logging).

    Returns
    -------
    dict
        Keys: ``x0``, ``y0``, ``a_arcsec``, ``b_arcsec``, ``a_rsun``,
        ``b_rsun``, ``arcsec_per_pixel`` (float or None), ``notes`` (str),
        ``display_text`` (str).
    """
    arcsec_per_pixel = None
    a_rsun = None
    b_rsun = None
    notes = ''

    if current_map is not None:
        # --- arcsec per pixel ---
        try:
            arcsec_per_pixel = float(current_map.scale[0].to(u.arcsec).value)
        except Exception:
            try:
                cdelt1 = (current_map.meta.get('CDELT1') or
                          current_map.meta.get('cdelt1'))
                if cdelt1 is not None:
                    arcsec_per_pixel = float(cdelt1) * 3600.0
            except Exception:
                arcsec_per_pixel = None

        # --- solar radius in arcsec ---
        try:
            solar_r_arcsec = float(current_map.rsun_obs.to(u.arcsec).value)
        except Exception:
            solar_r_arcsec = None

        if solar_r_arcsec is not None:
            a_rsun = a / solar_r_arcsec
            b_rsun = b / solar_r_arcsec
        else:
            solar_r_arcsec = 950.0
            notes = '[Unverified] used 950 arcsec per solar radius fallback'
            a_rsun = a / solar_r_arcsec
            b_rsun = b / solar_r_arcsec
    else:
        notes = '[Unverified] no map loaded; conversions not verified'
        solar_r_arcsec = 950.0
        a_rsun = a / solar_r_arcsec
        b_rsun = b / solar_r_arcsec

    display_text = (
        f'Center: ({x0:.1f}, {y0:.1f}) arcsec\n'
        f'a = {a:.1f} arcsec ({a_rsun:.4f} R_sun), '
        f'b = {b:.1f} arcsec ({b_rsun:.4f} R_sun)'
    )
    if arcsec_per_pixel is not None:
        display_text += f'\nArcsec/pixel = {arcsec_per_pixel:.3f}'
    if notes:
        display_text += f'\n{notes}'

    return {
        'x0': x0,
        'y0': y0,
        'a_arcsec': a,
        'b_arcsec': b,
        'a_rsun': a_rsun,
        'b_rsun': b_rsun,
        'arcsec_per_pixel': arcsec_per_pixel,
        'notes': notes,
        'display_text': display_text,
    }


def append_fit_to_csv(filepath, params, current_filename='N/A'):
    """
    Append one row of ellipse fit results to the CSV log file.

    Parameters
    ----------
    filepath : str
        Path to the CSV file created by :func:`init_fit_results_csv`.
    params : dict
        Dictionary as returned by :func:`extract_ellipse_params`.
    current_filename : str, optional
        Basename of the FITS file for the log entry.
    """
    with open(filepath, 'a', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([
            datetime.utcnow().isoformat(),
            current_filename,
            params['x0'], params['y0'],
            params['a_arcsec'], params['b_arcsec'],
            params['a_rsun'], params['b_rsun'],
            params['arcsec_per_pixel'] if params['arcsec_per_pixel'] is not None else '',
            params['notes'],
        ])


def export_fit_and_theta(x0, y0, a, b, theta_angles, script_dir):
    """
    Write ellipse fit parameters to a CSV and theta angles to a text file.

    Parameters
    ----------
    x0 : float
        Ellipse centre X in arcsec.
    y0 : float
        Ellipse centre Y in arcsec.
    a : float
        Semi-major axis in arcsec.
    b : float
        Semi-minor axis in arcsec.
    theta_angles : numpy.ndarray
        Array of angles (degrees) between field-line directions and the
        ellipsoid surface normals, as returned by :func:`create_3d_ellipsoid`.
    script_dir : str
        Directory where the output files are written.

    Returns
    -------
    tuple of (str, str)
        Paths to the written CSV file and text file respectively.
    """
    out_csv = os.path.join(script_dir, 'ellipse_fit.csv')
    out_txt = os.path.join(script_dir, 'theta_angles.txt')

    with open(out_csv, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['x0_arcsec', 'y0_arcsec', 'a_arcsec', 'b_arcsec'])
        writer.writerow([x0, y0, a, b])

    np.savetxt(out_txt, theta_angles)
    return out_csv, out_txt

