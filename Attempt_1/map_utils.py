# map_utils.py
import os
from PyQt5.QtWidgets import QFileDialog, QApplication, QMessageBox
import numpy as np
import sunpy.map
from scipy import ndimage
from matplotlib import colors

def load_files(viewer):
    files, _ = QFileDialog.getOpenFileNames(viewer, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
    if not files:
        return
    viewer.raw_files = files
    n = len(files)
    viewer.progress.setVisible(True)
    viewer.progress.setRange(0, n)
    viewer.progress.setValue(0)
    viewer.label.setText('Loading FITS files...')
    QApplication.processEvents()

    viewer.maps_and_files = []
    for i, f in enumerate(files):
        try:
            m = sunpy.map.Map(f)
            viewer.maps_and_files.append((m, f))
        except Exception as e:
            print(f'Error loading {f}: {e}')
        viewer.progress.setValue(i + 1)
        QApplication.processEvents()

    viewer.maps_and_files.sort(key=lambda mf: mf[0].date)
    v1 = float(viewer.vmin_input.text()) if viewer.vmin_input.text() else -50.0
    v2 = float(viewer.vmax_input.text()) if viewer.vmax_input.text() else 50.0

    viewer.maps = [mf[0] for mf in viewer.maps_and_files]
    viewer.files = [mf[1] for mf in viewer.maps_and_files]
    viewer.current_index = 0
    viewer.progress.setVisible(False)
    if viewer.maps:
        from map_utils import plot_map
        plot_map(viewer, viewer.maps[0])
        viewer.label.setText(f'{len(viewer.maps)} files loaded.')
    else:
        viewer.label.setText('No maps loaded.')


def plot_map(viewer, amap):
    viewer.canvas.figure.clf()
    viewer.ellipse_artist = None
    ax = viewer.canvas.figure.add_subplot(111, projection=amap)
    if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
        amap.plot(axes=ax)
    else:
        amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

    ax.grid(False)
    if viewer.ellipse_center is not None:
        from ellipse_utils import draw_ellipse
        draw_ellipse(viewer, viewer.ellipse_center[0], viewer.ellipse_center[1], current_map=amap, redraw_canvas=False)

    viewer.canvas.draw()
    current_filename = os.path.basename(viewer.files[viewer.current_index]) if viewer.files else 'N/A'
    viewer.label.setText(f'Showing: {current_filename}')

def update_display(viewer):
    viewer.figure.clear()
    m = viewer.maps[viewer.current_index]
    ax = viewer.figure.add_subplot(111, projection=m)
    m.plot(axes=ax)
    viewer.canvas.draw_idle()


def show_prev_map(viewer):
    if viewer.maps:
        viewer.current_index = max(0, viewer.current_index - 1)
        plot_map(viewer, viewer.maps[viewer.current_index])


def show_next_map(viewer):
    if viewer.maps:
        viewer.current_index = min(len(viewer.maps) - 1, viewer.current_index + 1)
        plot_map(viewer, viewer.maps[viewer.current_index])


def display_first_map(viewer):
    if not viewer.maps:
        viewer.label.setText('No maps loaded.')
        return
    viewer.current_index = 0
    plot_map(viewer, viewer.maps[0])


def display_last_map(viewer):
    if not viewer.maps:
        viewer.label.setText('No maps loaded.')
        return
    viewer.current_index = len(viewer.maps) - 1
    plot_map(viewer, viewer.maps[viewer.current_index])


def create_running_diff_maps(viewer):
    from ellipse_utils import set_processed_maps_from_loaded
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
    QApplication.processEvents()

    v1 = float(viewer.vmin_input.text()) if viewer.vmin_input.text() else -50.0
    v2 = float(viewer.vmax_input.text()) if viewer.vmax_input.text() else 50.0

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
        viewer.progress.setValue(idx + 1)
        QApplication.processEvents()

    viewer.maps = viewer.running_diff_maps
    viewer.current_index = 0
    viewer.progress.setVisible(False)
    if viewer.maps:
        plot_map(viewer, viewer.maps[0])
    viewer.label.setText(f'Created {len(viewer.running_diff_maps)} running difference maps.')