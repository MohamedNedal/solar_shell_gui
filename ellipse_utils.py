# ellipse_utils.py
import numpy as np
from PyQt5.QtWidgets import QMessageBox
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv
from datetime import datetime

def set_processed_maps_from_loaded(viewer):
    if not viewer.processed_maps and viewer.maps:
        viewer.processed_maps = viewer.maps.copy()


def draw_ellipse(viewer, x0, y0, current_map=None, redraw_canvas=True):
    try:
        if not viewer.maps:
            viewer.label.setText('No map to draw ellipse on.')
            return

        if current_map is None:
            current_map = viewer.maps[viewer.current_index]

        if not viewer.canvas.figure.axes:
            from .map_utils import plot_map
            plot_map(viewer, current_map)
            ax = viewer.canvas.figure.axes[0]
        else:
            ax = viewer.canvas.figure.axes[0]

        if viewer.ellipse_artist is not None:
            if viewer.ellipse_artist in ax.lines:
                viewer.ellipse_artist.remove()
            viewer.ellipse_artist = None

        viewer.ellipse_center = (x0, y0)
        a = viewer.a_slider.value()
        b = viewer.b_slider.value()

        theta = np.linspace(0, 2*np.pi, 100)
        x_ell = x0 + a * np.cos(theta)
        y_ell = y0 + b * np.sin(theta)
        coords = SkyCoord(x_ell*u.arcsec, y_ell*u.arcsec, frame=current_map.coordinate_frame)

        viewer.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)

        if redraw_canvas:
            viewer.canvas.draw()
            viewer.label.setText('Ellipse updated.')

    except Exception as e:
        viewer.label.setText(f'Ellipse error: {str(e)}')


def draw_ellipse_from_input(viewer):
    try:
        x0 = float(viewer.x_input.text())
        y0 = float(viewer.y_input.text())
        draw_ellipse(viewer, x0, y0, redraw_canvas=True)
    except ValueError:
        viewer.label.setText('Invalid center coordinates for ellipse.')
    except IndexError:
        viewer.label.setText('Ellipse error: No map loaded to draw on.')


def update_ellipse(viewer):
    if viewer.ellipse_center is not None and viewer.maps:
        current_map = viewer.maps[viewer.current_index]
        draw_ellipse(viewer, viewer.ellipse_center[0], viewer.ellipse_center[1], current_map=current_map)
    elif not viewer.maps:
        viewer.label.setText('No map loaded to update ellipse.')


def extract_ellipse_params(viewer):
    notes = ''
    try:
        x0 = float(viewer.x_input.text())
        y0 = float(viewer.y_input.text())
        a = float(viewer.a_slider.value())
        b = float(viewer.b_slider.value())
    except ValueError:
        viewer.label.setText('Invalid ellipse parameters.')
        return

    current_filename = viewer.files[viewer.current_index] if viewer.files else 'N/A'
    arcsec_per_pixel = None
    a_rsun = None
    b_rsun = None

    try:
        current_map = viewer.maps[viewer.current_index]
    except Exception:
        current_map = None

    if current_map is not None:
        try:
            arcsec_per_pixel = float(current_map.scale[0].to(u.arcsec).value)
        except Exception:
            try:
                cdelt1 = current_map.meta.get('CDELT1') or current_map.meta.get('cdelt1')
                if cdelt1 is not None:
                    arcsec_per_pixel = float(cdelt1) * 3600.0
            except Exception:
                arcsec_per_pixel = None

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

    display_text = (f'Center: ({x0:.1f}, {y0:.1f}) arcsec\n'
                    f'a = {a:.1f} arcsec ({a_rsun:.4f} R_sun), '
                    f'b = {b:.1f} arcsec ({b_rsun:.4f} R_sun)')
    if arcsec_per_pixel is not None:
        display_text += f'\nArcsec/pixel = {arcsec_per_pixel:.3f}'
    if notes:
        display_text += f'\n{notes}'

    viewer.fit_result_label.setText(display_text)
    viewer.label.setText('Fit extracted and saved.')

    try:
        with open(viewer.fit_results_csv, 'a', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow([
                datetime.utcnow().isoformat(),
                current_filename,
                x0, y0,
                a, b,
                a_rsun, b_rsun,
                arcsec_per_pixel if arcsec_per_pixel is not None else '',
                notes
            ])
    except Exception as e:
        viewer.label.setText(f'Failed to save fit: {e}')