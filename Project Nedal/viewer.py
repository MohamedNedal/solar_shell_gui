"""
viewer.py
---------
Main PyQt5 widget (AIAViewer) that wires together the GUI and delegates all
heavy computation to the helper modules.
"""

import os
import tempfile
import webbrowser
import sunpy.map
from sunpy.sun import constants as const
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider,
    QCheckBox, QProgressBar,
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
from scipy import ndimage
import plotly.io as pio

from fits_processing import load_fits_files, upgrade_maps_to_lv15, create_running_diff_maps
from ellipse_tools import (
    init_fit_results_csv,
    extract_ellipse_params,
    append_fit_to_csv,
    export_fit_and_theta,
)
from pfss_model import calculate_pfss_model
from visualization_3d import create_3d_ellipsoid


class AIAViewer(QWidget):
    """
    Main application window for viewing and analysing AIA FITS image sequences.

    Provides controls for:
    - Loading and navigating solar EUV image sequences.
    - Upgrading raw (level-1) FITS files to level 1.5.
    - Computing running-difference movies.
    - Overlaying and fitting ellipse regions of interest.
    - Computing PFSS magnetic field models from GONG magnetograms.
    - Rendering interactive 3-D CME shell visualisations.
    """

    def __init__(self):
        """Initialise application state, temp directories, the CSV log, and the UI."""
        super().__init__()
        self.setWindowTitle('AIA FITS Viewer')
        self.setGeometry(100, 100, 1200, 1000)

        self.ellipse_center = None
        self.ellipse_artist = None
        self.theta_angles = np.array([])
        self.files = []
        self.maps = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0
        self.channel = '193'
        self.theta = np.linspace(0, 2 * np.pi, 300)

        self.gong_filepath = None
        self.gong_map = None
        self.pfss_out = None
        self.field_lines = None

        self.temp_data_dir = os.path.join(tempfile.gettempdir(), 'solar_shell_temp_data')
        os.makedirs(self.temp_data_dir, exist_ok=True)
        print(tempfile.gettempdir())

        self.fit_results_csv = os.path.join(self.temp_data_dir, 'ellipse_fits.csv')
        init_fit_results_csv(self.fit_results_csv)

        self.init_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def init_ui(self):
        """Build and wire up all PyQt5 widgets and layout managers."""
        # --- Buttons ---
        self.load_button = QPushButton('Load FITS')
        self.display_button = QPushButton('Display')
        self.upgrade_button = QPushButton('Upgrade')
        self.run_diff_button = QPushButton('Run-Diff')
        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')
        self.save_button = QPushButton('Save')
        self.first_button = QPushButton('First Image')
        self.last_button = QPushButton('Last Image')
        self.ellipse_button = QPushButton('Ellipse')
        self.fit_button = QPushButton('Fit')
        self.export_params_button = QPushButton('Export Params')
        self.show_3d_button = QPushButton('Show 3D Ellipsoid')
        self.export_3d_png_button = QPushButton('Export 3D PNG')
        self.select_gong_button = QPushButton('Select GONG File')
        self.pfss_button = QPushButton('Calculate PFSS')
        self.pfss_button.setEnabled(False)

        # --- Text inputs ---
        self.vmin_input = QLineEdit()
        self.vmin_input.setPlaceholderText('Vmin')
        self.vmax_input = QLineEdit()
        self.vmax_input.setPlaceholderText('Vmax')
        self.x_input = QLineEdit()
        self.x_input.setPlaceholderText('X center [arcsec]')
        self.y_input = QLineEdit()
        self.y_input.setPlaceholderText('Y center [arcsec]')
        self.gong_file_input = QLineEdit()
        self.gong_file_input.setPlaceholderText('No GONG file selected')
        self.gong_file_input.setReadOnly(True)
        self.frame_input = QLineEdit()
        self.frame_input.setPlaceholderText('Image No.')

        # --- Sliders ---
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setMinimum(10)
        self.a_slider.setMaximum(1000)
        self.a_slider.setValue(200)

        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setMinimum(10)
        self.b_slider.setMaximum(1000)
        self.b_slider.setValue(100)

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

        # --- Checkboxes ---
        self.show_shell_checkbox = QCheckBox('Show Ellipsoid Shell')
        self.show_shell_checkbox.setChecked(True)
        self.show_normals_checkbox = QCheckBox('Show Normals')
        self.show_normals_checkbox.setChecked(False)
        self.show_radials_checkbox = QCheckBox('Show Radial Lines')
        self.show_radials_checkbox.setChecked(False)
        self.show_field_lines_checkbox = QCheckBox('Show Magnetic Field Lines')
        self.show_field_lines_checkbox.setChecked(True)

        # --- Status / progress ---
        self.label = QLabel('Ready.')
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setMinimum(0)
        self.progress.setMaximum(100)
        self.fit_result_label = QLabel('')
        self.fit_result_label.setWordWrap(True)

        # --- Signal connections ---
        self.load_button.clicked.connect(self.load_files)
        self.display_button.clicked.connect(self.display_first_map)
        self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
        self.run_diff_button.clicked.connect(self.run_running_diff)
        self.prev_button.clicked.connect(self.show_prev_map)
        self.next_button.clicked.connect(self.show_next_map)
        self.first_button.clicked.connect(self.display_first_map)
        self.last_button.clicked.connect(self.display_last_map)
        
        self.save_button.clicked.connect(self.save_current_map)
        self.ellipse_button.clicked.connect(self.draw_ellipse_from_input)
        self.fit_button.clicked.connect(self.on_fit_button)
        self.export_params_button.clicked.connect(self.on_export_params)
        self.show_3d_button.clicked.connect(self.show_3d_ellipsoid_plot)
        self.export_3d_png_button.clicked.connect(self.export_3d_png_silent)
        self.select_gong_button.clicked.connect(self.select_gong_file)
        self.pfss_button.clicked.connect(self.on_calculate_pfss)
        self.a_slider.valueChanged.connect(self.update_ellipse)
        self.b_slider.valueChanged.connect(self.update_ellipse)

        # --- Layouts ---
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        right_layout.setSpacing(8)
        right_layout.setContentsMargins(10, 10, 10, 10)
        frame_layout = QHBoxLayout()
        frame_layout.addWidget(self.frame_input)

        # Apply button for frame input
        apply_layout = QHBoxLayout()
        self.apply_button = QPushButton('Apply')
        self.apply_button.clicked.connect(self.apply_frameinput)
        apply_layout.addWidget(self.apply_button)

        btn_layout = QHBoxLayout()
        for btn in [
            self.load_button, self.display_button, self.upgrade_button,
            self.run_diff_button, self.prev_button, self.next_button,
            self.first_button, self.last_button, self.frame_input, self.apply_button, self.save_button,
        ]:
            btn_layout.addWidget(btn)
        left_layout.addLayout(btn_layout)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        left_layout.addWidget(self.canvas, stretch=1)

        status_layout = QHBoxLayout()
        status_layout.addWidget(self.label, stretch=1)
        status_layout.addWidget(self.progress, stretch=0)
        left_layout.addLayout(status_layout)

        main_layout.addLayout(left_layout, stretch=3)
        main_layout.addLayout(right_layout, stretch=1)
        self.setLayout(main_layout)

        # image layout
        image_layout = QHBoxLayout()
        image_layout.addWidget(QLabel('Vmin & Vmax'))
        image_layout.addWidget(self.vmin_input)
        right_layout.addLayout(image_layout)
        image_layout.addWidget(self.vmax_input)

        #Apply button for vmin/max
        apply_layout = QVBoxLayout()
        self.apply_button = QPushButton('Apply')
        self.apply_button.clicked.connect(self.apply_vminmax)
        apply_layout.addWidget(self.apply_button)
        right_layout.addLayout(apply_layout)

        # Right panel: ellipse controls
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
        ellipse_layout.addWidget(self.fit_result_label)
        ellipse_layout.addWidget(self.export_params_button)
        right_layout.addLayout(ellipse_layout)

        # Right panel: GONG / PFSS
        gong_layout = QVBoxLayout()
        gong_layout.addWidget(QLabel('GONG Data for PFSS:'))
        gong_layout.addWidget(self.select_gong_button)
        gong_layout.addWidget(self.gong_file_input)
        gong_layout.addWidget(self.pfss_button)
        right_layout.addLayout(gong_layout)

        # Right panel: 3-D visualisation
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
        _3d_layout.addWidget(self.export_3d_png_button)
        right_layout.addLayout(_3d_layout)

    # ------------------------------------------------------------------
    # File / map navigation
    # ------------------------------------------------------------------

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
        if not files:
            return
        
        self.raw_files = files

        n = len(files)
        self.progress.setVisible(True)
        self.progress.setRange(0, n)
        self.progress.setValue(0)
        self.label.setText('Loading FITS files...')
        QApplication.processEvents()

        # load maps once and keep the original filename with each map
        self.maps_and_files = []
        maps_and_files = []
        for i, f in enumerate(files):
            try:
                m = sunpy.map.Map(f)
                self.maps_and_files.append((m, f))
                maps_and_files.append((m, f))
            except Exception as e:
                print(f'Error loading {f}: {e}')
            self.progress.setValue(i + 1)
            QApplication.processEvents()

        # sort by the map observation time (chronological)
        maps_and_files.sort(key=lambda mf: mf[0].date)
        vmin_text = self.vmin_input.text()
        vmax_text = self.vmax_input.text()
        v1 = float(vmin_text) if vmin_text else -50.0
        v2 = float(vmax_text) if vmax_text else 50.0
        
        # unzip into two lists that are in the same chronological order
        self.maps = [mf[0] for mf in maps_and_files]
        self.files = [mf[1] for mf in maps_and_files]
        self.current_index = 0
        self.progress.setVisible(False)
        if self.maps:
            self.label.setText(f'{len(self.maps)} files loaded.')
            # show the first (chronological) map explicitly ->self.plot_map(self.maps[0], v1, v2)----
            self.plot_map(self.maps[0])
        else:
            self.label.setText('No maps loaded.')

    def display_first_map(self):
        """Jump to and display the first (earliest) map in the sequence."""
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = 0
        vmin_text = self.vmin_input.text()
        vmax_text = self.vmax_input.text()
        v1 = float(vmin_text) if vmin_text else -50.0
        v2 = float(vmax_text) if vmax_text else 50.0
        self.plot_map(self.maps[0])

    def display_last_map(self):
        """Jump to and display the last (latest) map in the sequence."""
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = len(self.maps) - 1
        self.plot_map(self.maps[self.current_index])

    def show_prev_map(self):
        """Step back one frame and display it, clamped at the first frame."""
        if self.maps:
            self.current_index = max(0, self.current_index - 1)
            self.plot_map(self.maps[self.current_index])

    def show_next_map(self):
        """Step forward one frame and display it, clamped at the last frame."""
        if self.maps:
            self.current_index = min(len(self.maps) - 1, self.current_index + 1)
            self.plot_map(self.maps[self.current_index])

    def plot_map(self, amap):
        """
        Render a SunPy map on the central Matplotlib canvas.

        Redraws any existing ellipse overlay on top of the new image.

        Parameters
        ----------
        amap : sunpy.map.GenericMap
            The map to render.
        """
        self.canvas.figure.clf()
        self.ellipse_artist = None
        ax = self.canvas.figure.add_subplot(111, projection=amap)
        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9) * u.percent)
        ax.grid(False)
        if self.ellipse_center is not None:
            self.draw_ellipse(
                self.ellipse_center[0], self.ellipse_center[1],
                current_map=amap, redraw_canvas=False,
            )
        self.canvas.draw()
        name = os.path.basename(self.files[self.current_index]) if self.files else 'N/A'
        self.label.setText(f'Showing: {name}')

    def save_current_map(self):
        """
        Prompt the user for a path and save the current canvas as a
        high-resolution PNG or PDF (300 dpi).
        """
        if not self.maps:
            QMessageBox.warning(self, 'Warning', 'No map to save.')
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Save current map as image',
            f'{self.current_index:03d}.png',
            'PNG (*.png);;PDF (*.pdf)',
        )
        if file_path:
            self.figure.tight_layout()
            self.figure.savefig(file_path, bbox_inches='tight', pad_inches=0.05, dpi=300)
            self.label.setText(f'Saved map to {file_path}')

    # ------------------------------------------------------------------
    # FITS processing
    # ------------------------------------------------------------------

    def upgrade_to_lv15(self):
        """
        Upgrade all loaded FITS files from AIA level 1 to level 1.5.

        Applies pointing correction, image registration, and exposure
        normalisation.  Results are cached to disk and replace the current
        map sequence.
        """
        if not self.files:
            QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
            return

        n = len(self.files)
        self.progress.setVisible(True)
        self.progress.setRange(0, n)
        self.progress.setValue(0)
        self.label.setText('Upgrading files to level 1.5...')
        QApplication.processEvents()

        def _cb(current, total):
            self.progress.setValue(current)
            QApplication.processEvents()

        self.processed_maps = upgrade_maps_to_lv15(
            self.files, self.channel, self.temp_data_dir, progress_callback=_cb)
        self.maps = self.processed_maps
        self.current_index = 0
        self.progress.setVisible(False)

        if self.maps:
            self.plot_map(self.maps[0])
            self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
        else:
            self.label.setText('No maps were processed or loaded.')

    def set_processed_maps_from_loaded(self):
        """
        Populate ``processed_maps`` from ``maps`` if the Upgrade step was skipped.

        Used as a fallback before running-difference computation.
        """
        if not self.processed_maps and self.maps:
            self.processed_maps = self.maps.copy()

    def update_display(self):
        self.figure.clear()
        m = self.maps[self.current_index]
        ax = self.figure.add_subplot(111, projection=m)
        m.plot(axes=ax)
        ax.grid(False)
        self.canvas.draw_idle()

    def apply_vminmax(self):
        v1 = float(self.vmin_input.text())
        v2 = float(self.vmax_input.text())
        m = self.maps[self.current_index]
        m.plot_settings['norm'] = colors.Normalize(vmin=v1, vmax=v2)
        self.maps[self.current_index] = m
        self.update_display()

    def apply_frameinput(self):
        f = int(self.frame_input.text())
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = f - 1
        max_maps = len(self.maps)

        if not (1 <= f <= max_maps):
            self.label.setText(f"Valid range: 1–{max_maps}")
            return

        self.current_index = f - 1
        m, filename = self.maps_and_files[self.current_index]
        self.plot_map(self.maps[self.current_index])

    def run_running_diff(self):
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

        vmin_text = self.vmin_input.text()
        vmax_text = self.vmax_input.text()
        v1 = float(vmin_text) if vmin_text else -50.0
        v2 = float(vmax_text) if vmax_text else 50.0
        
        for idx, i in enumerate(range(5, len(self.processed_maps))):
            try:
                m0 = self.processed_maps[i-5]
                m1 = self.processed_maps[i]
                diff = m1.quantity - m0.quantity
                smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
                diff_map = sunpy.map.Map(smoothed, m1.meta)
                diff_map.plot_settings['norm'] = colors.Normalize(vmin=v1, vmax=v2)
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

    # ------------------------------------------------------------------
    # Ellipse tools
    # ------------------------------------------------------------------

    def draw_ellipse(self, x0, y0, current_map=None, redraw_canvas=True):
        """
        Draw (or redraw) the ellipse overlay on the active Matplotlib axis.

        Parameters
        ----------
        x0 : float
            Ellipse centre X in arcsec.
        y0 : float
            Ellipse centre Y in arcsec.
        current_map : sunpy.map.GenericMap, optional
            Map whose coordinate frame is used.  Defaults to the current map.
        redraw_canvas : bool
            Whether to immediately refresh the Qt canvas.
        """
        try:
            if not self.maps:
                self.label.setText('No map to draw ellipse on.')
                return
            if current_map is None:
                current_map = self.maps[self.current_index]
            if not self.canvas.figure.axes:
                self.plot_map(current_map)
            ax = self.canvas.figure.axes[0]

            if self.ellipse_artist is not None and self.ellipse_artist in ax.lines:
                self.ellipse_artist.remove()
            self.ellipse_artist = None
            self.ellipse_center = (x0, y0)

            a = self.a_slider.value()
            b = self.b_slider.value()
            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)
            coords = SkyCoord(
                x_ell * u.arcsec, y_ell * u.arcsec,
                frame=current_map.coordinate_frame,
            )
            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)
            if redraw_canvas:
                self.canvas.draw()
                self.label.setText('Ellipse updated.')
        except Exception as e:
            self.label.setText(f'Ellipse error: {str(e)}')

    def draw_ellipse_from_input(self):
        """
        Read the ellipse centre from the X/Y text fields and draw the overlay.
        Shows an error message in the status bar for invalid inputs.
        """
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            self.draw_ellipse(x0, y0, redraw_canvas=True)
        except ValueError:
            self.label.setText('Invalid center coordinates for ellipse.')
        except IndexError:
            self.label.setText('Ellipse error: No map loaded to draw on.')

    def update_ellipse(self):
        """
        Redraw the ellipse after a slider value changes, if a centre is set.
        """
        if self.ellipse_center is not None and self.maps:
            self.draw_ellipse(
                self.ellipse_center[0], self.ellipse_center[1],
                current_map=self.maps[self.current_index],
            )
        elif not self.maps:
            self.label.setText('No map loaded to update ellipse on.')

    def on_fit_button(self):
        """
        Extract ellipse parameters in arcsec and R_sun, update the GUI label,
        and append the result to the persistent fit-results CSV log.
        """
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            a = float(self.a_slider.value())
            b = float(self.b_slider.value())
        except ValueError:
            self.label.setText('Invalid ellipse parameters.')
            return

        name = os.path.basename(self.files[self.current_index]) if self.files else 'N/A'
        try:
            current_map = self.maps[self.current_index]
        except Exception:
            current_map = None

        params = extract_ellipse_params(x0, y0, a, b, current_map, name)
        self.fit_result_label.setText(params['display_text'])
        self.label.setText('Fit extracted and saved.')

        try:
            append_fit_to_csv(self.fit_results_csv, params, name)
        except Exception as e:
            self.label.setText(f'Failed to save fit: {e}')

    def on_export_params(self):
        """
        Export the current ellipse fit to ``ellipse_fit.csv`` and theta angles
        to ``theta_angles.txt``, both written to the script directory.
        """
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            a = float(self.a_slider.value())
            b = float(self.b_slider.value())
        except ValueError:
            self.label.setText('Export error: invalid ellipse parameters.')
            return
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            out_csv, out_txt = export_fit_and_theta(
                x0, y0, a, b, self.theta_angles, script_dir)
            self.label.setText(f'Exported fit to {out_csv}')
            self.label.setText(f'Exported theta to {out_txt}')
        except Exception as e:
            self.label.setText(f'Export error: {e}')

    # ------------------------------------------------------------------
    # GONG / PFSS
    # ------------------------------------------------------------------

    def select_gong_file(self):
        """
        Open a file dialog to choose a GONG magnetogram FITS file and
        enable the PFSS calculation button once a valid file is selected.
        """
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

    def on_calculate_pfss(self):
        """
        Run PFSS model computation from the selected GONG file.

        Shows a busy (indeterminate) progress bar during the calculation and
        stores the resulting PFSS solution and field lines for 3-D visualisation.
        """
        if self.gong_filepath is None:
            QMessageBox.warning(self, 'Warning', 'Please select a GONG FITS file first.')
            return

        self.label.setText('Starting PFSS calculation...')
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)  # indeterminate / busy
        self.repaint()
        QApplication.processEvents()

        def _status(msg):
            self.label.setText(msg)
            QApplication.processEvents()

        try:
            self.gong_map, self.pfss_out, self.field_lines = calculate_pfss_model(
                self.gong_filepath, progress_callback=_status)
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

    # ------------------------------------------------------------------
    # 3-D visualisation
    # ------------------------------------------------------------------

    def _build_ellipse_params(self):
        """
        Assemble the ellipse parameter dictionary from the current UI state.

        Returns
        -------
        dict
            Keys: ``'x0'``, ``'y0'``, ``'a'``, ``'b'``.
        """
        x0, y0 = self.ellipse_center
        return {'x0': x0, 'y0': y0, 'a': self.a_slider.value(), 'b': self.b_slider.value()}

    def show_3d_ellipsoid_plot(self):
        """
        Build the 3-D ellipsoid visualisation and open it in the default web
        browser as an interactive Plotly HTML page.
        """
        if self.ellipse_center is None or not self.maps:
            QMessageBox.warning(
                self, 'Warning', 'Please draw an ellipse and load a map first.')
            return
        if self.show_field_lines_checkbox.isChecked() and self.field_lines is None:
            QMessageBox.warning(
                self, 'Warning',
                'Magnetic field data is not loaded. Please select a GONG file and '
                'click "Calculate PFSS" first, or uncheck "Show Magnetic Field Lines".',
            )
            return
        try:
            fig, self.theta_angles = create_3d_ellipsoid(
                self._build_ellipse_params(),
                self.maps[self.current_index],
                field_lines=self.field_lines,
                show_shell=self.show_shell_checkbox.isChecked(),
                show_normals=self.show_normals_checkbox.isChecked(),
                show_radials=self.show_radials_checkbox.isChecked(),
                n_lat_radial=self.n_lat_radial_slider.value(),
                n_lon_radial=self.n_lon_radial_slider.value(),
                show_field_lines=self.show_field_lines_checkbox.isChecked(),
            )
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as f:
                temp_filepath = f.name
                pio.write_html(fig, file=f, auto_open=False, include_plotlyjs='cdn')
            webbrowser.open(f'file://{temp_filepath}')
            QMessageBox.information(
                self, '3D Plot Displayed',
                f'Opened in your default browser:\n{temp_filepath}',
            )
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to create 3D plot: {str(e)}')

    def export_3d_png_silent(self):
        """
        Silently render the 3-D ellipsoid and save it as a zero-padded PNG
        (e.g. ``frame_001.png``) in the script directory.

        Requires the ``kaleido`` package for static image export.
        """
        if self.ellipse_center is None or not self.maps:
            self.label.setText('Draw an ellipse and load a map first.')
            return
        try:
            fig, self.theta_angles = create_3d_ellipsoid(
                self._build_ellipse_params(),
                self.maps[self.current_index],
                field_lines=self.field_lines,
                show_shell=self.show_shell_checkbox.isChecked(),
                show_normals=self.show_normals_checkbox.isChecked(),
                show_radials=self.show_radials_checkbox.isChecked(),
                n_lat_radial=self.n_lat_radial_slider.value(),
                n_lon_radial=self.n_lon_radial_slider.value(),
                show_field_lines=self.show_field_lines_checkbox.isChecked(),
            )
            script_dir = os.path.dirname(os.path.abspath(__file__))
            out_png = os.path.join(script_dir, f'frame_{self.current_index + 1:03d}.png')
            pio.write_image(fig, out_png, width=1600, height=1200, scale=2)
            self.label.setText(f'3D PNG saved to {out_png}')
        except Exception as e:
            self.label.setText(f'PNG export failed: {e}')

