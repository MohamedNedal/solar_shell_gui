import os
import sys
import tempfile
import webbrowser
import csv
from datetime import datetime

from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider, QCheckBox, QProgressBar
)
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
import sunpy.map
import numpy as np
from aiapy.calibrate import register, update_pointing
from scipy import ndimage
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.sun import constants as const

import pfsspy
import pfsspy.tracing as tracing

import plotly.graph_objects as go
import plotly.io as pio


class AIAViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('AIA FITS Viewer')
        self.setGeometry(100, 100, 1200, 1000)

        self.ellipse_center = None
        self.ellipse_artist = None

        self.files = []
        self.maps = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0

        self.temp_data_dir = os.path.join(tempfile.gettempdir(), 'solar_shell_temp_data')
        os.makedirs(self.temp_data_dir, exist_ok=True)

        self.channel = '193'
        self.theta = np.linspace(0, 2*np.pi, 300)

        self.gong_filepath = None
        self.gong_map = None
        self.pfss_out = None
        self.field_lines = None

        # CSV file for fits results
        self.fit_results_csv = os.path.join(self.temp_data_dir, 'ellipse_fits.csv')
        if not os.path.exists(self.fit_results_csv):
            with open(self.fit_results_csv, 'w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow([
                    'timestamp', 'filename', 'x0_arcsec', 'y0_arcsec',
                    'a_arcsec', 'b_arcsec', 'a_rsun', 'b_rsun', 'arcsec_per_pixel', 'notes'
                ])

        self.init_ui()

    def init_ui(self):
        # Buttons
        self.load_button = QPushButton('Load FITS')
        self.display_button = QPushButton('Display')
        self.upgrade_button = QPushButton('Upgrade')
        self.run_diff_button = QPushButton('Run-Diff')
        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')
        self.save_button = QPushButton('Save')

        self.load_button.clicked.connect(self.load_files)
        self.display_button.clicked.connect(self.display_first_map)
        self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
        self.run_diff_button.clicked.connect(self.create_running_diff_maps)
        self.prev_button.clicked.connect(self.show_prev_map)
        self.next_button.clicked.connect(self.show_next_map)
        self.save_button.clicked.connect(self.save_current_map)

        self.ellipse_button = QPushButton('Ellipse')
        self.fit_button = QPushButton('Fit')

        # Inputs for ellipse center
        self.x_input = QLineEdit()
        self.x_input.setPlaceholderText('X center [arcsec]')
        self.y_input = QLineEdit()
        self.y_input.setPlaceholderText('Y center [arcsec]')

        # Sliders for axes (treated as arcsec values)
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setMinimum(10)
        self.a_slider.setMaximum(1000)
        self.a_slider.setValue(200)
        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setMinimum(10)
        self.b_slider.setMaximum(1000)
        self.b_slider.setValue(100)

        self.a_slider.valueChanged.connect(self.update_ellipse)
        self.b_slider.valueChanged.connect(self.update_ellipse)

        self.ellipse_button.clicked.connect(self.draw_ellipse_from_input)
        self.fit_button.clicked.connect(self.extract_ellipse_params)

        # 3D Ellipsoid functionality
        self.show_3d_button = QPushButton('Show 3D Ellipsoid')
        self.show_3d_button.clicked.connect(self.show_3d_ellipsoid_plot)

        self.show_shell_checkbox = QCheckBox('Show Ellipsoid Shell')
        self.show_shell_checkbox.setChecked(True)
        self.show_normals_checkbox = QCheckBox('Show Normals')
        self.show_normals_checkbox.setChecked(False)
        self.show_radials_checkbox = QCheckBox('Show Radial Lines')
        self.show_radials_checkbox.setChecked(False)
        self.show_field_lines_checkbox = QCheckBox('Show Magnetic Field Lines')
        self.show_field_lines_checkbox.setChecked(True)

        self.n_lat_radial_slider = QSlider(Qt.Horizontal)
        self.n_lat_radial_slider.setMinimum(5)
        self.n_lat_radial_slider.setMaximum(50)
        self.n_lat_radial_slider.setValue(10)
        self.n_lat_radial_slider.setToolTip('Number of radial lines (vertical)')

        self.n_lon_radial_slider = QSlider(Qt.Horizontal)
        self.n_lon_radial_slider.setMinimum(10)
        self.n_lon_radial_slider.setMaximum(100)
        self.n_lon_radial_slider.setValue(30)
        self.n_lon_radial_slider.setToolTip('Number of radial lines (horizontal)')

        # GONG and PFSS controls
        self.select_gong_button = QPushButton('Select GONG File')
        self.select_gong_button.clicked.connect(self.select_gong_file)
        self.gong_file_input = QLineEdit()
        self.gong_file_input.setPlaceholderText('No GONG file selected')
        self.gong_file_input.setReadOnly(True)

        self.pfss_button = QPushButton('Calculate PFSS')
        self.pfss_button.clicked.connect(self.calculate_pfss_model)
        self.pfss_button.setEnabled(False)

        # Status label and progress bar
        self.label = QLabel('Ready.')
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setMinimum(0)
        self.progress.setMaximum(100)

        # Fit result label (displayed under Fit button)
        self.fit_result_label = QLabel('')
        self.fit_result_label.setWordWrap(True)

        # === MAIN LAYOUTS =======================================================
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        # Adjust spacing/padding on the right side menu (user asked for this)
        right_layout.setSpacing(8)  # space between widgets
        right_layout.setContentsMargins(10, 10, 10, 10)  # left, top, right, bottom padding

        # --- top buttons row ----------------------------------------------------
        btn_layout = QHBoxLayout()
        for b in [self.load_button, self.display_button, self.upgrade_button,
                self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
            btn_layout.addWidget(b)

        left_layout.addLayout(btn_layout)

        # --- big canvas in the center ------------------------------------------
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        left_layout.addWidget(self.canvas, stretch=1)

        # status area under canvas
        status_layout = QHBoxLayout()
        status_layout.addWidget(self.label, stretch=1)
        status_layout.addWidget(self.progress, stretch=0)
        left_layout.addLayout(status_layout)

        # Final assembly
        main_layout.addLayout(left_layout, stretch=3)
        main_layout.addLayout(right_layout, stretch=1)
        self.setLayout(main_layout)

        # RIGHT side: ellipse controls
        ellipse_layout = QVBoxLayout()
        ellipse_layout.addWidget(QLabel('Ellipse Center:'))
        ellipse_layout.addWidget(self.x_input)
        ellipse_layout.addWidget(QLabel('Y center [arcsec]:'))
        ellipse_layout.addWidget(self.y_input)
        ellipse_layout.addWidget(QLabel('Semi-Major Axis [arcsec]:'))
        ellipse_layout.addWidget(self.a_slider)
        ellipse_layout.addWidget(QLabel('Semi-Minor Axis [arcsec]:'))
        ellipse_layout.addWidget(self.b_slider)
        ellipse_layout.addWidget(self.ellipse_button)
        ellipse_layout.addWidget(self.fit_button)
        ellipse_layout.addWidget(self.fit_result_label)  # shows fit values
        right_layout.addLayout(ellipse_layout)

        # GONG/PFSS controls
        gong_pfss_layout = QVBoxLayout()
        gong_pfss_layout.addWidget(QLabel('GONG Data for PFSS:'))
        gong_pfss_layout.addWidget(self.select_gong_button)
        gong_pfss_layout.addWidget(self.gong_file_input)
        gong_pfss_layout.addWidget(self.pfss_button)
        right_layout.addLayout(gong_pfss_layout)

        # 3D controls
        _3d_layout = QVBoxLayout()
        _3d_layout.addWidget(self.show_3d_button)
        _3d_layout.addWidget(self.show_shell_checkbox)
        _3d_layout.addWidget(self.show_normals_checkbox)
        _3d_layout.addWidget(self.show_radials_checkbox)
        _3d_layout.addWidget(QLabel('Radial Lines Vertical:'))
        _3d_layout.addWidget(self.n_lat_radial_slider)
        _3d_layout.addWidget(QLabel('Radial Lines Horizontal:'))
        _3d_layout.addWidget(self.n_lon_radial_slider)
        _3d_layout.addWidget(self.show_field_lines_checkbox)
        right_layout.addLayout(_3d_layout)

    def select_gong_file(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, 'Select GONG FITS File', '', 'FITS files (*.fits *.fts)')
        if filepath:
            self.gong_filepath = filepath
            self.gong_file_input.setText(os.path.basename(filepath))
            self.pfss_button.setEnabled(True)
            self.label.setText(f'GONG file selected: {os.path.basename(filepath)}')
        else:
            self.gong_filepath = None
            self.gong_file_input.setText('No GONG file selected')
            self.pfss_button.setEnabled(False)
            self.label.setText('No GONG file selected.')

    def calculate_pfss_model(self):
        if self.gong_filepath is None:
            QMessageBox.warning(self, 'Warning', 'Please select a GONG FITS file first.')
            return

        # set UI to busy
        self.label.setText('Starting PFSS calculation...')
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)  # indeterminate / busy
        self.repaint()
        QApplication.processEvents()

        try:
            gong_map = sunpy.map.Map(self.gong_filepath)

            if 'cunit1' not in gong_map.meta:
                gong_map.meta['cunit1'] = u.deg

            self.gong_map = gong_map

            nrho = 50
            rss = 3
            pfss_in = pfsspy.Input(self.gong_map, nrho, rss)
            self.label.setText('Calculating PFSS model...')
            QApplication.processEvents()
            self.pfss_out = pfsspy.pfss(pfss_in)

            num_footpoints_lat = 40
            num_footpoints_lon = 60
            r_trace = 1.05 * const.radius

            lat = np.linspace(np.radians(-90), np.radians(90), num_footpoints_lat, endpoint=False)
            lon = np.linspace(np.radians(-180), np.radians(180), num_footpoints_lon, endpoint=False)

            lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')
            lat_flat, lon_flat = lat_grid.ravel()*u.rad, lon_grid.ravel()*u.rad

            seeds = SkyCoord(lon_flat, lat_flat, r_trace, frame=self.pfss_out.coordinate_frame)
            tracer = tracing.FortranTracer()
            self.label.setText('Tracing magnetic field lines...')
            QApplication.processEvents()
            self.field_lines = tracer.trace(seeds, self.pfss_out)
            self.label.setText('PFSS model calculated and magnetic field data loaded.')
        except Exception as e:
            self.label.setText(f'Error calculating PFSS model: {e}')
            self.pfss_out = None
            self.field_lines = None
            self.gong_map = None
        finally:
            self.progress.setVisible(False)
            self.progress.setRange(0, 100)
            QApplication.processEvents()

    def detect_data_level(self, file):
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

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
        if not files:
            return

        self.files = files
        n = len(files)
        self.maps = []
        self.progress.setVisible(True)
        self.progress.setRange(0, n)
        self.progress.setValue(0)
        self.label.setText('Loading FITS files...')
        QApplication.processEvents()

        for i, f in enumerate(files):
            try:
                m = sunpy.map.Map(f)
                self.maps.append(m)
            except Exception as e:
                print(f'Error loading {f}: {e}')
            self.progress.setValue(i + 1)
            QApplication.processEvents()

        self.current_index = 0
        self.progress.setVisible(False)
        if self.maps:
            self.label.setText(f'{len(self.maps)} files loaded.')
            self.plot_map(self.maps[0])
        else:
            self.label.setText('No maps loaded.')

    def display_first_map(self):
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = 0
        self.plot_map(self.maps[0])

    def show_prev_map(self):
        if self.maps:
            self.current_index = max(0, self.current_index - 1)
            self.plot_map(self.maps[self.current_index])

    def show_next_map(self):
        if self.maps:
            self.current_index = min(len(self.maps) - 1, self.current_index + 1)
            self.plot_map(self.maps[self.current_index])

    def plot_map(self, amap):
        self.canvas.figure.clf()
        self.ellipse_artist = None

        ax = self.canvas.figure.add_subplot(111, projection=amap)

        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

        ax.grid(False)

        if self.ellipse_center is not None:
            self.draw_ellipse(self.ellipse_center[0], self.ellipse_center[1], current_map=amap, redraw_canvas=False)

        self.canvas.draw()
        # show only basename to avoid huge strings
        current_filename = os.path.basename(self.files[self.current_index]) if self.files else 'N/A'
        self.label.setText(f'Showing: {current_filename}')

    def upgrade_to_lv15(self):
        if not self.files:
            QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
            return

        self.processed_maps = []
        n = len(self.files)
        self.progress.setVisible(True)
        self.progress.setRange(0, n)
        self.progress.setValue(0)
        self.label.setText('Upgrading files to level 1.5...')
        QApplication.processEvents()

        for i, file in enumerate(self.files):
            data_level = self.detect_data_level(file)
            try:
                if data_level == 'lev1':
                    output_filename = os.path.basename(file).replace('lev1', 'lev15')
                    file_path = os.path.join(self.temp_data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
                    os.makedirs(file_path, exist_ok=True)
                    full_path = os.path.join(file_path, output_filename)

                    if not os.path.exists(full_path):
                        m = sunpy.map.Map(file)
                        m = update_pointing(m, pointing_table=None)
                        m = register(m)
                        m = m / m.exposure_time
                        m.save(full_path, filetype='auto')
                        self.processed_maps.append(m)
                    else:
                        self.processed_maps.append(sunpy.map.Map(full_path))

                elif data_level == 'lev15':
                    self.processed_maps.append(sunpy.map.Map(file))
                else:
                    print(f'Skipped unknown-level file: {file}')
            except Exception as e:
                print(f'Error upgrading {file}: {e}')

            self.progress.setValue(i + 1)
            QApplication.processEvents()

        self.maps = self.processed_maps
        self.current_index = 0
        self.progress.setVisible(False)
        if self.maps:
            self.plot_map(self.maps[0])
            self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
        else:
            self.label.setText('No maps were processed or loaded.')

    def set_processed_maps_from_loaded(self):
        if not self.processed_maps and self.maps:
            self.processed_maps = self.maps.copy()

    def create_running_diff_maps(self):
        self.set_processed_maps_from_loaded()

        if len(self.processed_maps) < 6:
            QMessageBox.warning(self, 'Warning',
                                'Need at least 6 images for running difference.')
            return

        self.running_diff_maps = []
        n = len(self.processed_maps) - 5
        self.progress.setVisible(True)
        self.progress.setRange(0, n)
        self.progress.setValue(0)
        self.label.setText('Creating running-difference maps...')
        QApplication.processEvents()

        for idx, i in enumerate(range(5, len(self.processed_maps))):
            try:
                m0 = self.processed_maps[i-5]
                m1 = self.processed_maps[i]
                diff = m1.quantity - m0.quantity
                smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
                diff_map = sunpy.map.Map(smoothed, m1.meta)
                diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
                self.running_diff_maps.append(diff_map)
            except Exception as e:
                print(f'Error creating diff map {i}: {e}')
            self.progress.setValue(idx + 1)
            QApplication.processEvents()

        self.maps = self.running_diff_maps
        self.current_index = 0
        self.progress.setVisible(False)
        if self.maps:
            self.plot_map(self.maps[0])
        self.label.setText(f'Created {len(self.running_diff_maps)} running difference maps.')

    def draw_ellipse(self, x0, y0, current_map=None, redraw_canvas=True):
        try:
            if not self.maps:
                self.label.setText('No map to draw ellipse on.')
                return

            if current_map is None:
                current_map = self.maps[self.current_index]

            if not self.canvas.figure.axes:
                self.plot_map(current_map)
                ax = self.canvas.figure.axes[0]
            else:
                ax = self.canvas.figure.axes[0]

            if self.ellipse_artist is not None and self.ellipse_artist in ax.lines:
                self.ellipse_artist.remove()
                self.ellipse_artist = None
            elif self.ellipse_artist is not None:
                self.ellipse_artist = None

            self.ellipse_center = (x0, y0)

            a = self.a_slider.value()
            b = self.b_slider.value()

            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)

            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec,
                              frame=current_map.coordinate_frame)

            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)

            if redraw_canvas:
                self.canvas.draw()
                self.label.setText('Ellipse updated.')
        except Exception as e:
            self.label.setText(f'Ellipse error: {str(e)}')

    def draw_ellipse_from_input(self):
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            self.draw_ellipse(x0, y0, redraw_canvas=True)
        except ValueError:
            self.label.setText('Invalid center coordinates for ellipse.')
        except IndexError:
            self.label.setText('Ellipse error: No map loaded to draw on.')

    def update_ellipse(self):
        if self.ellipse_center is not None and self.maps:
            current_map = self.maps[self.current_index]
            self.draw_ellipse(self.ellipse_center[0], self.ellipse_center[1], current_map=current_map)
        elif not self.maps:
            self.label.setText('No map loaded to update ellipse on.')

    def extract_ellipse_params(self):
        """
        Show fit values below the 'Fit' button and store them in a CSV.
        Also convert semi-major/minor to solar radii using map metadata if possible.
        """
        notes = ''
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            a = float(self.a_slider.value())  # arcsec (we assume slider in arcsec)
            b = float(self.b_slider.value())
        except ValueError:
            self.label.setText('Invalid ellipse parameters.')
            return

        current_filename = os.path.basename(self.files[self.current_index]) if self.files else 'N/A'
        arcsec_per_pixel = None

        # Try getting arcsec per pixel from map.scale if available
        a_rsun = None
        b_rsun = None
        try:
            current_map = self.maps[self.current_index]
        except Exception:
            current_map = None

        if current_map is not None:
            # arcsec per pixel
            try:
                # current_map.scale is a pair of Quantity objects in astropy units
                arcsec_per_pixel = float(current_map.scale[0].to(u.arcsec).value)
            except Exception:
                # try header keywords
                try:
                    cdelt1 = current_map.meta.get('CDELT1') or current_map.meta.get('cdelt1')
                    if cdelt1 is not None:
                        # cdelt1 usually in degrees/pixel
                        arcsec_per_pixel = float(cdelt1) * 3600.0
                except Exception:
                    arcsec_per_pixel = None

            # solar radius in arcsec from map metadata if present
            try:
                # many sunpy maps have rsun_obs in units of arcsec
                solar_r_arcsec = float(current_map.rsun_obs.to(u.arcsec).value)
            except Exception:
                solar_r_arcsec = None

            if solar_r_arcsec is not None:
                a_rsun = a / solar_r_arcsec
                b_rsun = b / solar_r_arcsec
            else:
                # [Unverified] fallback to 950 arcsec per solar radius (user suggested value)
                solar_r_arcsec = 950.0
                notes = '[Unverified] used 950 arcsec per solar radius fallback'
                a_rsun = a / solar_r_arcsec
                b_rsun = b / solar_r_arcsec
        else:
            # no map loaded: still write values but mark as unverified
            notes = '[Unverified] no map loaded; conversions not verified'
            # fallback conversions using 950 arcsec/rad
            solar_r_arcsec = 950.0
            a_rsun = a / solar_r_arcsec
            b_rsun = b / solar_r_arcsec

        # Prepare display text
        display_text = (
            f'Center: ({x0:.1f}, {y0:.1f}) arcsec\n'
            f'a = {a:.1f} arcsec ({a_rsun:.4f} R_sun), '
            f'b = {b:.1f} arcsec ({b_rsun:.4f} R_sun)'
        )
        if arcsec_per_pixel is not None:
            display_text += f'\nArcsec/pixel = {arcsec_per_pixel:.3f}'
        if notes:
            display_text += f'\n{notes}'

        # Show in GUI under Fit button
        self.fit_result_label.setText(display_text)
        self.label.setText('Fit extracted and saved.')

        # Append to CSV
        try:
            with open(self.fit_results_csv, 'a', newline='') as fh:
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
            self.label.setText(f'Failed to save fit: {e}')

    def save_current_map(self):
        if not self.maps:
            QMessageBox.warning(self, 'Warning', 'No map to save.')
            return

        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Save current map as image',
            f'{self.current_index:03d}.png',
            'PNG (*.png);;PDF (*.pdf)', options=options)

        if file_path:
            self.figure.tight_layout()
            self.figure.savefig(file_path, bbox_inches='tight', pad_inches=0.05, dpi=300)
            self.label.setText(f'Saved map to {file_path}')

    # create_3d_ellipsoid and show_3d_ellipsoid_plot unchanged (kept for brevity)
    # I'll include them unchanged from your original code, but with the same interface.
    def create_3d_ellipsoid(self, ellipse_params, sunpy_map, show_shell=True, show_normals=False, show_radials=False, n_lat_radial=10, n_lon_radial=30, show_field_lines=False):
        # identical implementation to your original function (kept as-is)
        # ... (copy the function body from your original script here)
        theta, phi = np.mgrid[0:np.pi:50j, 0:2*np.pi:50j]
        x_sun = np.sin(theta) * np.cos(phi)
        y_sun = np.sin(theta) * np.sin(phi)
        z_sun = np.cos(theta)
        fig = go.Figure()
        fig.add_trace(go.Surface(x=x_sun, y=y_sun, z=z_sun, colorscale=[[0, 'orange'], [1, 'orange']], showscale=False, name='Sun', hoverinfo='skip'))
        x0, y0 = ellipse_params['x0'], ellipse_params['y0']
        a_radius, b_radius = ellipse_params['a'], ellipse_params['b']
        hpc_coord = SkyCoord(x0 * u.arcsec, y0 * u.arcsec,
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
        nx_all = 2 * (x_src - xshift) / (b_minor**2)
        ny_all = 2 * (y_src - yshift) / (b_minor**2)
        nz_all = 2 * (z_src - zshift) / (b_major**2)
        norm_magnitude_all = np.sqrt(nx_all**2 + ny_all**2 + nz_all**2)
        nx_all /= norm_magnitude_all
        ny_all /= norm_magnitude_all
        nz_all /= norm_magnitude_all
        nx_flat = nx_all[mask].flatten()
        ny_flat = ny_all[mask].flatten()
        nz_flat = nz_all[mask].flatten()
        x_outer = x_src[mask].flatten()
        y_outer = y_src[mask].flatten()
        z_outer = z_src[mask].flatten()
        all_points_for_theta = []
        all_dirs_for_theta = []
        if show_radials:
            radial_length = 2.0
            theta_radial, phi_radial = np.mgrid[0:np.pi:complex(n_lat_radial), 0:2*np.pi:complex(n_lon_radial)]
            rdx = np.sin(theta_radial) * np.cos(phi_radial)
            rdy = np.sin(theta_radial) * np.sin(phi_radial)
            rdz = np.cos(theta_radial)
            rdx_flat = rdx.flatten()
            rdy_flat = rdy.flatten()
            rdz_flat = rdz.flatten()
            for i in range(len(rdx_flat)):
                r_vec = np.array([rdx_flat[i], rdy_flat[i], rdz_flat[i]])
                r_pts = np.linspace(0, radial_length, 20)
                lx = r_pts * r_vec[0]
                ly = r_pts * r_vec[1]
                lz = r_pts * r_vec[2]
                fig.add_trace(go.Scatter3d(x=lx, y=ly, z=lz, mode='lines', line=dict(color='gray', width=1), showlegend=False, hoverinfo='skip', opacity=0.3))
                if not show_field_lines or self.field_lines is None:
                    all_points_for_theta.append(r_vec * radial_length)
                    all_dirs_for_theta.append(r_vec)
        if show_field_lines and self.field_lines:
            all_points_for_theta = []
            all_dirs_for_theta = []
            for n, field_line in enumerate(self.field_lines):
                color = {0:'black', -1:'blue', 1:'red'}.get(field_line.polarity)
                coords = field_line.coords
                coords.representation_type = 'cartesian'
                x_field = coords.x / const.radius
                y_field = coords.y / const.radius
                z_field = coords.z / const.radius
                fig.add_trace(go.Scatter3d(x=x_field, y=y_field, z=z_field, mode='lines', line=dict(color=color, width=2), showlegend=False))
                for i in range(len(x_field)-1):
                    field_point = np.array([x_field[i].value, y_field[i].value, z_field[i].value])
                    all_points_for_theta.append(field_point)
                    field_vector = np.array([x_field[i+1].value - x_field[i].value, y_field[i+1].value - y_field[i].value, z_field[i+1].value - z_field[i].value])
                    field_vector = field_vector / np.linalg.norm(field_vector)
                    all_dirs_for_theta.append(field_vector)
        all_points_for_theta = np.array(all_points_for_theta)
        all_dirs_for_theta = np.array(all_dirs_for_theta)
        theta_angles = np.full_like(x_outer, np.nan)
        if len(all_points_for_theta) > 0 and len(x_outer) > 0:
            for i in range(len(x_outer)):
                surf_pt = np.array([x_outer[i], y_outer[i], z_outer[i]])
                normal_vec = np.array([nx_flat[i], ny_flat[i], nz_flat[i]])
                dists = np.linalg.norm(all_points_for_theta - surf_pt, axis=1)
                idx = np.argmin(dists)
                line_dir = all_dirs_for_theta[idx]
                cos_theta = np.dot(normal_vec, line_dir)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta_angles[i] = np.degrees(np.arccos(np.abs(cos_theta)))
        theta_surface = np.full_like(x_src, np.nan)
        theta_surface[mask] = theta_angles
        if show_shell:
            fig.add_trace(go.Scatter3d(x=[xshift], y=[yshift], z=[zshift], mode='markers', marker=dict(size=8, color='black'), showlegend=False, name='Shell Center'))
            fig.add_trace(go.Surface(x=x_display, y=y_display, z=z_display, surfacecolor=theta_surface, colorscale='Viridis', cmin=0, cmax=90, opacity=1, showscale=True, colorbar=dict(title=dict(text='Theta Angle (degrees)', side='right'), len=0.6), name='Shell'))
        if show_normals:
            step = max(1, len(x_outer) // 100)
            x_norm_sample = x_outer[::step]
            y_norm_sample = y_outer[::step]
            z_norm_sample = z_outer[::step]
            nx_sample = nx_all[mask][::step] * 0.2
            ny_sample = ny_all[mask][::step] * 0.2
            nz_sample = nz_all[mask][::step] * 0.2
            x_lines = []
            y_lines = []
            z_lines = []
            for i in range(len(x_norm_sample)):
                x_lines.extend([x_norm_sample[i], x_norm_sample[i] + nx_sample[i], None])
                y_lines.extend([y_norm_sample[i], y_norm_sample[i] + ny_sample[i], None])
                z_lines.extend([z_norm_sample[i], z_norm_sample[i] + nz_sample[i], None])
            fig.add_trace(go.Scatter3d(x=x_lines, y=y_lines, z=z_lines, mode='lines', line=dict(color='black', width=3), showlegend=False, name='Normal Vectors'))
        fig.update_layout(scene=dict(xaxis=dict(range=[-2, 2], title='X (Solar Radii)'), yaxis=dict(range=[-2, 2], title='Y (Solar Radii)'), zaxis=dict(range=[-2, 2], title='Z (Solar Radii)'), aspectmode='cube', camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))), width=1024, height=768, title='Interactive 3D Solar Shell Visualization', showlegend=True)
        return fig

    def show_3d_ellipsoid_plot(self):
        if self.ellipse_center is None or not self.maps:
            QMessageBox.warning(self, 'Warning', 'Please draw an ellipse and load a map first.')
            return

        if self.show_field_lines_checkbox.isChecked() and self.field_lines is None:
            QMessageBox.warning(self, 'Warning', 'Magnetic field data is not loaded. Please select a GONG file and click "PFSS" first, or uncheck "Show Magnetic Field Lines".')
            return

        try:
            x0, y0 = self.ellipse_center
            a_radius = self.a_slider.value()
            b_radius = self.b_slider.value()

            ellipse_params = {
                'x0': x0,
                'y0': y0,
                'a': a_radius,
                'b': b_radius
            }

            current_map = self.maps[self.current_index]
            show_shell = self.show_shell_checkbox.isChecked()
            show_normals = self.show_normals_checkbox.isChecked()
            show_radials = self.show_radials_checkbox.isChecked()
            n_lat_radial = self.n_lat_radial_slider.value()
            n_lon_radial = self.n_lon_radial_slider.value()
            show_field_lines = self.show_field_lines_checkbox.isChecked()

            fig = self.create_3d_ellipsoid(
                ellipse_params,
                current_map,
                show_shell,
                show_normals,
                show_radials,
                n_lat_radial,
                n_lon_radial,
                show_field_lines
            )

            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as f:
                temp_filepath = f.name
                pio.write_html(fig, file=f, auto_open=False, include_plotlyjs='cdn')

            webbrowser.open(f'file://{temp_filepath}')
            QMessageBox.information(self, '3D Plot Displayed',
                                    f'The 3D ellipsoid plot has been opened in your default web browser at:\n{temp_filepath}\n\n'
                                    'This is a workaround for the "Could not find QtWebEngineProcess" error.')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to create 3D plot: {str(e)}')


def main():
    app = QApplication(sys.argv)
    viewer = AIAViewer()
    viewer.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
    