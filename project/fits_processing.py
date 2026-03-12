"""
fits_processing.py
------------------
Standalone functions for loading, upgrading, and differencing AIA FITS files.
None of these functions depend on Qt or the viewer widget.
"""

import os
import sunpy.map
import numpy as np
from aiapy.calibrate import register, update_pointing
from scipy import ndimage
from matplotlib import colors
import astropy.units as u


def detect_data_level(filepath):
    """
    Determine the data level of an AIA FITS file.

    Checks the filename for level indicators ('lev1', 'lev1.5', 'lev15'),
    then falls back to reading the FITS header keyword 'lvl_num' if needed.

    Parameters
    ----------
    filepath : str
        Path to the FITS file.

    Returns
    -------
    str
        One of 'lev1', 'lev15', or 'unknown'.
    """
    name = os.path.basename(filepath).lower()
    if 'lev1.5' in name or 'lev15' in name:
        return 'lev15'
    elif 'lev1' in name:
        return 'lev1'
    else:
        try:
            hdr = sunpy.map.Map(filepath).meta
            level = hdr.get('lvl_num', None)
            if level == 1.5:
                return 'lev15'
            elif level == 1:
                return 'lev1'
        except Exception as e:
            print(f'Could not read header from {filepath}: {e}')
        return 'unknown'


def load_fits_files(filepaths, progress_callback=None):
    """
    Load a list of FITS files into SunPy Map objects, sorted chronologically.

    Parameters
    ----------
    filepaths : list of str
        Paths to the FITS files to load.
    progress_callback : callable, optional
        Called as ``progress_callback(current, total)`` after each file is processed.

    Returns
    -------
    maps : list of sunpy.map.GenericMap
        Loaded maps sorted by observation time.
    files : list of str
        Corresponding file paths in the same sorted order.
    """
    maps_and_files = []
    total = len(filepaths)
    for i, f in enumerate(filepaths):
        try:
            m = sunpy.map.Map(f)
            maps_and_files.append((m, f))
        except Exception as e:
            print(f'Error loading {f}: {e}')
        if progress_callback:
            progress_callback(i + 1, total)

    maps_and_files.sort(key=lambda mf: mf[0].date)
    maps = [mf[0] for mf in maps_and_files]
    files = [mf[1] for mf in maps_and_files]
    return maps, files


def upgrade_maps_to_lv15(files, channel, temp_dir, progress_callback=None):
    """
    Upgrade AIA FITS files from level 1 to level 1.5.

    For each file:
    - Level-1 files are pointing-corrected, registered, and normalised by
      exposure time, then saved to a subdirectory of ``temp_dir``.
    - Already level-1.5 files are loaded directly.
    - Files of unknown level are skipped with a warning.

    Parameters
    ----------
    files : list of str
        Paths to the FITS files to process.
    channel : str
        AIA wavelength channel (e.g. '193'), used to organise output paths.
    temp_dir : str
        Root temporary directory where processed files are cached.
    progress_callback : callable, optional
        Called as ``progress_callback(current, total)`` after each file.

    Returns
    -------
    list of sunpy.map.GenericMap
        Processed (level 1.5) maps.
    """
    processed_maps = []
    total = len(files)
    for i, file in enumerate(files):
        data_level = detect_data_level(file)
        try:
            if data_level == 'lev1':
                output_filename = os.path.basename(file).replace('lev1', 'lev15')
                file_path = os.path.join(temp_dir, 'AIA', f'{channel}A', 'processed', 'lv15')
                os.makedirs(file_path, exist_ok=True)
                full_path = os.path.join(file_path, output_filename)
                if not os.path.exists(full_path):
                    m = sunpy.map.Map(file)
                    m = update_pointing(m, pointing_table=None)
                    m = register(m)
                    m = m / m.exposure_time
                    m.save(full_path, filetype='auto')
                    processed_maps.append(m)
                else:
                    processed_maps.append(sunpy.map.Map(full_path))
            elif data_level == 'lev15':
                processed_maps.append(sunpy.map.Map(file))
            else:
                print(f'Skipped unknown-level file: {file}')
        except Exception as e:
            print(f'Error upgrading {file}: {e}')
        if progress_callback:
            progress_callback(i + 1, total)
    return processed_maps


def create_running_diff_maps(processed_maps, progress_callback=None):
    """
    Compute running-difference maps from a sequence of processed AIA maps.

    Each difference is taken between frame ``i`` and frame ``i - 5``, then
    smoothed with a Gaussian filter (sigma=3).  The resulting maps use a
    symmetric ±50 DN/s colour normalisation.

    Parameters
    ----------
    processed_maps : list of sunpy.map.GenericMap
        At least 6 level-1.5 maps in chronological order.
    progress_callback : callable, optional
        Called as ``progress_callback(current, total)`` after each diff map.

    Returns
    -------
    list of sunpy.map.GenericMap
        Running-difference maps ready for display.

    Raises
    ------
    ValueError
        If fewer than 6 maps are supplied.
    """
    if len(processed_maps) < 6:
        raise ValueError('Need at least 6 images for running difference.')

    running_diff_maps = []
    n = len(processed_maps) - 5
    for idx, i in enumerate(range(5, len(processed_maps))):
        try:
            m0 = processed_maps[i - 5]
            m1 = processed_maps[i]
            diff = m1.quantity - m0.quantity
            smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
            diff_map = sunpy.map.Map(smoothed, m1.meta)
            diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
            running_diff_maps.append(diff_map)
        except Exception as e:
            print(f'Error creating diff map {i}: {e}')
        if progress_callback:
            progress_callback(idx + 1, n)
    return running_diff_maps

