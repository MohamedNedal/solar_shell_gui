import sys
import os
import csv
from datetime import datetime
import tempfile
import webbrowser
import numpy as np
import pfsspy
import pfsspy.tracing as tracing

import plotly.graph_objects as go
import plotly.io as pio

from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QSlider, QPushButton, QLineEdit,
    QFileDialog, QProgressBar, QMessageBox, QCheckBox, QVBoxLayout, QHBoxLayout
)
from PyQt5.QtCore import Qt
from matplotlib import colors
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# SunPy + Astropy imports
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy import units as u

# Absolute imports from your modules
from upgrades import upgrade_to_lv15
from pfss_utils import (calculate_pfss_model, detect_data_level, select_gong_file)
from map_utils import (
    load_files, display_first_map, display_last_map,
    show_prev_map, show_next_map, plot_map, create_running_diff_maps,
    update_display
)
from ellipse_utils import (
    draw_ellipse, draw_ellipse_from_input, update_ellipse,
    extract_ellipse_params, set_processed_maps_from_loaded
)
from three_d_utils import (
    create_3d_ellipsoid, show_3d_ellipsoid_plot, export_3d_png_silent
)


class AIAViewer(QWidget):

    def draw_ellipse_from_input(self):
        draw_ellipse_from_input(self)

    def draw_ellipse(self):
        draw_ellipse(self)

    def update_ellipse(self):
        update_ellipse(self)

    def update_display(self):
        update_display(self)

    def create_running_diff_maps(viewer):
        create_running_diff_maps

    def __init__(self):
        super().__init__()
        
    

        # -----------------------------
        # Data
        # -----------------------------
        self.maps = []
        self.files = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0
        self.ellipse_center = None
        self.theta_angles = None
        self.field_lines = None
        self.gong_filepath = None
        self.gong_map = None
        self.pfss_out = None
        # -----------------------------
        # UI elements
        # -----------------------------
        self.label = QLabel('Ready', self)
        self.progress = QProgressBar(self)
        self.vmin_input = QLineEdit(self)
        self.vmax_input = QLineEdit(self)
        self.frame_input = QLineEdit(self)
        self.x_input = QLineEdit(self)
        self.y_input = QLineEdit(self)
        self.a_slider = QSlider(Qt.Horizontal)
        self.b_slider = QSlider(Qt.Horizontal)
        self.n_lat_radial_slider = QSlider(Qt.Horizontal)
        self.n_lon_radial_slider = QSlider(Qt.Horizontal)
        self.show_shell_checkbox = QCheckBox('Show Shell')
        self.show_normals_checkbox = QCheckBox('Show Normals')
        self.show_radials_checkbox = QCheckBox('Show Radials')
        self.show_field_lines_checkbox = QCheckBox('Show Magnetic Field Lines')
        self.fit_result_label = QLabel('Fit results will appear here.')

        # Placeholder for matplotlib canvas
        from matplotlib.figure import Figure
        # === MAIN LAYOUTS =======================================================
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        # CSV for storing fits
        self.fit_results_csv = os.path.join(os.getcwd(), 'fit_results.csv')
        if not os.path.exists(self.fit_results_csv):
            with open(self.fit_results_csv, 'w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow([
                    'timestamp', 'filename', 'x0', 'y0', 'a', 'b', 'a_rsun', 'b_rsun',
                    'arcsec_per_pixel', 'notes'
                ])

        # -----------------------------
        # Layout
        # -----------------------------
        self.setWindowTitle("AIA Viewer")
        self.resize(1200, 800)

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

        # Inputs for vmin & vmax
        self.vmin_input.setPlaceholderText('Vmin')
        self.vmax_input.setPlaceholderText('Vmax')

        self.load_button = QPushButton('Load FITS')
        self.display_button = QPushButton('Display')
        self.upgrade_button = QPushButton('Upgrade')
        self.run_diff_button = QPushButton('Run-Diff')
        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')
        self.save_button = QPushButton('Save')
        
        self.first_button = QPushButton('First Image')
        self.last_button = QPushButton('Last Image')
        self.first_button.clicked.connect(self.display_first_map)
        self.last_button.clicked.connect(self.display_last_map)
        # self.frame_input = QLineEdit()
        self.frame_input.setPlaceholderText('Image No.')
        frame_layout = QHBoxLayout()
        frame_layout.addWidget(self.frame_input)

        # Apply button for frame input
        apply_layout = QHBoxLayout()
        self.apply_button = QPushButton('Apply')
        self.apply_button.clicked.connect(self.apply_frameinput)
        apply_layout.addWidget(self.apply_button)

        self.load_button.clicked.connect(self.load_files)
        self.display_button.clicked.connect(self.display_first_map)
        self.upgrade_button.clicked.connect(lambda: self.upgrade_to_lv15)
        self.run_diff_button.clicked.connect(self.create_running_diff_maps)
        self.prev_button.clicked.connect(self.show_prev_map)
        self.next_button.clicked.connect(self.show_next_map)
        self.save_button.clicked.connect(self.save_current_map)

        self.ellipse_button = QPushButton('Ellipse')
        self.fit_button = QPushButton('Fit')
        
        self.export_params_button = QPushButton('Export Params')        
        self.export_params_button.clicked.connect(self.export_fit_and_theta)

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
        self.fit_button.clicked.connect(lambda: self.extract_ellipse_params)

        self.export_3d_png_button = QPushButton('Export 3D PNG')
        self.export_3d_png_button.clicked.connect(lambda: self.export_3d_png_silent)
        
        # 3D Ellipsoid functionality
        self.show_3d_button = QPushButton('Show 3D Ellipsoid')
        self.show_3d_button.clicked.connect(lambda: self.show_3d_ellipsoid_plot)

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
        self.select_gong_button.clicked.connect(lambda: self.select_gong_file)
        self.gong_file_input = QLineEdit()
        self.gong_file_input.setPlaceholderText('No GONG file selected')
        self.gong_file_input.setReadOnly(True)

        self.pfss_button = QPushButton('Calculate PFSS')
        self.pfss_button.clicked.connect(lambda: self.calculate_pfss_model)
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

        # --- top buttons row ----------------------------------------------------
        btn_layout = QHBoxLayout()
        for b in [self.load_button, self.display_button, self.upgrade_button,
                self.run_diff_button, self.prev_button, self.next_button,
                  self.first_button, self.last_button, self.frame_input, self.apply_button, self.save_button]:
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
        main_layout.addLayout(image_layout)
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
        ellipse_layout.addWidget(self.export_params_button)
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
        _3d_layout.addWidget(self.export_3d_png_button)
        right_layout.addLayout(_3d_layout)

        # -----------------------------
        # Example: connect sliders to ellipse update
        # -----------------------------
        self.a_slider.valueChanged.connect(self.update_display)
        self.b_slider.valueChanged.connect(self.update_display)
        
    
    # -----------------------------
    # Wrapper functions calling module utilities
    # -----------------------------
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
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = len(self.maps) - 1
        self.plot_map(self.maps[self.current_index])

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

    def export_fit_and_theta(self):
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            out_fit = os.path.join(script_dir, 'ellipse_fit.csv')
            out_theta = os.path.join(script_dir, 'theta_surface.csv')
            
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            a = float(self.a_slider.value())
            b = float(self.b_slider.value())
            print("Export exists:", os.path.exists(out_fit))

            self.label.setText(f'Exported fit + XYZ + theta to {out_fit}')
                
            np.savetxt(out_theta, self.theta_surface)
                
            self.label.setText(
                f'Exported fit params to {out_fit}\n'
                f'Exported theta to {out_theta}'
            )

        except Exception as e:
            self.label.setText(f'Export error: {e}')
    
    def apply_vminmax(self):
        print("vmin text:", repr(self.vmin_input.text()))
        print("vmax text:", repr(self.vmax_input.text()))
        try:
            v1 = float(self.vmin_input.text())
            v2 = float(self.vmax_input.text())
        except ValueError:
            print("Invalid vmin/vmax input")
            return
        
        print("vmin text:", repr(self.vmin_input.text()))
        print("vmax text:", repr(self.vmax_input.text()))
        m = self.maps[self.current_index]
        m.plot_settings['norm'] = colors.Normalize(vmin=v1, vmax=v2)
        self.maps[self.current_index] = m
        self.update_display()