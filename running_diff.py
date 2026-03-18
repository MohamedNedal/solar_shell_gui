import numpy as np
from scipy import ndimage
from PyQt5.QtWidgets import QMessageBox
from matplotlib import colors
import sunpy.map

def set_processed_maps_from_loaded(viewer):
    if not getattr(viewer, 'processed_maps', None) and viewer.maps:
        viewer.processed_maps = viewer.maps.copy()

def create_running_diff_maps(viewer):
    set_processed_maps_from_loaded(viewer)

    if len(viewer.processed_maps) < 6:
        QMessageBox.warning(viewer, 'Warning', 'Need at least 6 images for running difference.')
        return

    viewer.running_diff_maps = []
    n = len(viewer.processed_maps) - 5
    viewer.progress.setVisible(True)
    viewer.progress.setRange(0, n)
    viewer.progress.setValue(0)
    viewer.label.setText('Creating running-difference maps...')
    
    vmin_text = viewer.vmin_input.text()
    vmax_text = viewer.vmax_input.text()
    v1 = float(vmin_text) if vmin_text else -50.0
    v2 = float(vmax_text) if vmax_text else 50.0

    for idx, i in enumerate(range(5, len(viewer.processed_maps))):
        try:
            m0 = viewer.processed_maps[i-5]
            m1 = viewer.processed_maps[i]
            diff = m1.quantity - m0.quantity
            smoothed = ndimage.gaussian_filter(diff, sigma=[3,3])
            diff_map = sunpy.map.Map(smoothed, m1.meta)
            diff_map.plot_settings['norm'] = colors.Normalize(vmin=v1, vmax=v2)
            viewer.running_diff_maps.append(diff_map)
        except Exception as e:
            print(f'Error creating diff map {i}: {e}')
        viewer.progress.setValue(idx+1)
        viewer.progress.repaint()
    
    viewer.maps = viewer.running_diff_maps
    viewer.current_index = 0
    viewer.progress.setVisible(False)
    
    if viewer.maps:
        viewer.plot_map(viewer.maps[0])
    viewer.label.setText(f'Created {len(viewer.running_diff_maps)} running difference maps.')