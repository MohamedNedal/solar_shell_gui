# pfss_utils.py
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from pfsspy import Input, pfss, tracing
from astropy.constants import R_sun as const

def calculate_pfss_model(viewer):
    """
    Compute PFSS model and trace magnetic field lines.
    Stores results in the viewer object.
    """
    if viewer.gong_filepath is None:
        from PyQt5.QtWidgets import QMessageBox
        QMessageBox.warning(viewer, 'Warning', 'Please select a GONG FITS file first.')
        return

    # Set UI to busy
    viewer.label.setText('Starting PFSS calculation...')
    viewer.progress.setVisible(True)
    viewer.progress.setRange(0, 0)  # indeterminate
    viewer.repaint()
    
    from PyQt5.QtWidgets import QApplication
    QApplication.processEvents()

    try:
        # Load GONG map
        gong_map = sunpy.map.Map(viewer.gong_filepath)

        if 'cunit1' not in gong_map.meta:
            gong_map.meta['cunit1'] = u.deg

        viewer.gong_map = gong_map

        # PFSS parameters
        nrho = 50
        rss = 3

        pfss_in = Input(viewer.gong_map, nrho, rss)
        viewer.label.setText('Calculating PFSS model...')
        QApplication.processEvents()
        viewer.pfss_out = pfss(pfss_in)

        # Trace magnetic field lines
        num_footpoints_lat = 40
        num_footpoints_lon = 60
        r_trace = 1.05 * const

        lat = np.linspace(np.radians(-90), np.radians(90), num_footpoints_lat, endpoint=False)
        lon = np.linspace(np.radians(-180), np.radians(180), num_footpoints_lon, endpoint=False)

        lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')
        lat_flat, lon_flat = lat_grid.ravel() * u.rad, lon_grid.ravel() * u.rad

        seeds = SkyCoord(lon_flat, lat_flat, r_trace, frame=viewer.pfss_out.coordinate_frame)

        tracer = tracing.FortranTracer()
        viewer.label.setText('Tracing magnetic field lines...')
        QApplication.processEvents()
        viewer.field_lines = tracer.trace(seeds, viewer.pfss_out)

        viewer.label.setText('PFSS model calculated and magnetic field data loaded.')
    except Exception as e:
        viewer.label.setText(f'Error calculating PFSS model: {e}')
        viewer.pfss_out = None
        viewer.field_lines = None
        viewer.gong_map = None
    finally:
        viewer.progress.setVisible(False)
        viewer.progress.setRange(0, 100)
        QApplication.processEvents()


def detect_data_level(file):
    """
    Detect whether the AIA FITS file is Level 1 or Level 1.5.
    """
    import os
    name = os.path.basename(file).lower()
    if 'lev1.5' in name or 'lev15' in name:
        return 'lev15'
    elif 'lev1' in name:
        return 'lev1'
    else:
        try:
            hdr = sunpy.map.Map(file).meta
            level = hdr.get('lvl_num', None)
            if level == 1.5:
                return 'lev15'
            elif level == 1:
                return 'lev1'
        except Exception as e:
            print(f'Could not read header from {file}: {e}')
        return 'unknown'
    
def select_gong_file(viewer):
        filepath, _ = QFileDialog.getOpenFileName(
            viewer, 'Select GONG FITS File', '', 'FITS files (*.fits *.fts)')
        if filepath:
            viewer.gong_filepath = filepath
            viewer.gong_file_input.setText(os.path.basename(filepath))
            viewer.pfss_button.setEnabled(True)
            viewer.label.setText(f'GONG file selected: {os.path.basename(filepath)}')
        else:
            viewer.gong_filepath = None
            viewer.gong_file_input.setText('No GONG file selected')
            viewer.pfss_button.setEnabled(False)
            viewer.label.setText('No GONG file selected.')