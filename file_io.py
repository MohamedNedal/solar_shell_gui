import os
import csv
from datetime import datetime
import numpy as np
from PyQt5.QtWidgets import QFileDialog, QMessageBox
import sunpy.map
import matplotlib.colors as colors

def detect_data_level(file):
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

def load_files(viewer):
    files, _ = QFileDialog.getOpenFileNames(viewer, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
    if not files:
        return
    viewer.raw_files = files
    n = len(files)
    viewer.progress.setVisible(True)
    viewer.progress.setRange(0,n)
    viewer.progress.setValue(0)
    viewer.label.setText('Loading FITS files...')

    viewer.maps_and_files = []
    maps_and_files = []
    for i,f in enumerate(files):
        try:
            m = sunpy.map.Map(f)
            viewer.maps_and_files.append((m,f))
            maps_and_files.append((m,f))
        except Exception as e:
            print(f'Error loading {f}: {e}')
        viewer.progress.setValue(i+1)
        viewer.progress.repaint()

    maps_and_files.sort(key=lambda mf: mf[0].date)
    vmin_text = viewer.vmin_input.text()
    vmax_text = viewer.vmax_input.text()
    v1 = float(vmin_text) if vmin_text else -50.0
    v2 = float(vmax_text) if vmax_text else 50.0

    viewer.maps = [mf[0] for mf in maps_and_files]
    viewer.files = [mf[1] for mf in maps_and_files]
    viewer.current_index = 0
    viewer.progress.setVisible(False)
    
    if viewer.maps:
        viewer.label.setText(f'{len(viewer.maps)} files loaded.')
        viewer.plot_map(viewer.maps[0])
    else:
        viewer.label.setText('No maps loaded.')

def export_fit_and_theta(viewer):
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        out_fit = os.path.join(script_dir, 'ellipse_fit.csv')
        out_theta = os.path.join(script_dir, 'theta_surface.csv')

        np.savetxt(out_theta, viewer.theta_surface)
        viewer.label.setText(f'Exported fit params to {out_fit}\nExported theta to {out_theta}')
    except Exception as e:
        viewer.label.setText(f'Export error: {e}')