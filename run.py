import os
import sys
import tempfile # For creating temporary files
import webbrowser # For opening HTML files in a browser

from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider, QCheckBox, QDialog, QSizePolicy
)
from PyQt5.QtCore import Qt, QUrl

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
        self.ellipse_artist = None # This will hold the actual Line2D object

        self.files = []
        self.maps = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0
        
        # Use a temporary directory for GONG data processing outputs if needed,
        # but GONG file itself will be selected by user.
        self.temp_data_dir = os.path.join(tempfile.gettempdir(), 'solar_shell_temp_data')
        os.makedirs(self.temp_data_dir, exist_ok=True)

        self.channel = '193'
        self.theta = np.linspace(0, 2*np.pi, 300)

        self.gong_filepath = None # To store the path of the selected GONG file
        self.gong_map = None # Store the GONG map
        self.pfss_out = None # Store PFSS output
        self.field_lines = None # Store traced field lines

        self.init_ui()
    
    def init_ui(self):
        # --- Create Widgets ---

        # Buttons
        self.load_button = QPushButton('Load FITS Files')
        self.display_button = QPushButton('Display First Map')
        self.upgrade_button = QPushButton('Upgrade to Level 1.5')
        self.run_diff_button = QPushButton('Create Running Diff')
        self.prev_button = QPushButton('Previous Map')
        self.next_button = QPushButton('Next Map')
        self.save_button = QPushButton('Save Current Map')

        # Define btn_layout here, right after all buttons it contains are defined
        # Changed to QVBoxLayout for vertical stacking in the right menu
        btn_layout = QVBoxLayout() 
        for b in [self.load_button, self.display_button, self.upgrade_button,
                  self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
            btn_layout.addWidget(b)

        self.ellipse_button = QPushButton('Draw Ellipse') 
        self.fit_button = QPushButton('Extract Ellipse Params') 
        
        # Inputs for ellipse center
        self.x_input = QLineEdit()
        self.x_input.setPlaceholderText('X center [arcsec]')
        self.y_input = QLineEdit()
        self.y_input.setPlaceholderText('Y center [arcsec]')

        # Sliders for axes
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setMinimum(10)
        self.a_slider.setMaximum(1000)
        self.a_slider.setValue(200)
        self.a_value_label = QLabel(f'{self.a_slider.value()} arcsec') # Label for A slider value

        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setMinimum(10)
        self.b_slider.setMaximum(1000)
        self.b_slider.setValue(100)
        self.b_value_label = QLabel(f'{self.b_slider.value()} arcsec') # Label for B slider value

        # 3D Ellipsoid functionality
        self.show_3d_button = QPushButton('Show 3D Ellipsoid')

        self.show_shell_checkbox = QCheckBox('Show Ellipsoid Shell')
        self.show_shell_checkbox.setChecked(True) # Default to true
        self.show_normals_checkbox = QCheckBox('Show Normals')
        self.show_normals_checkbox.setChecked(False) # Default to false
        self.show_radials_checkbox = QCheckBox('Show Radial Lines')
        self.show_radials_checkbox.setChecked(False)
        self.show_field_lines_checkbox = QCheckBox('Show Magnetic Field Lines')
        self.show_field_lines_checkbox.setChecked(True) # Default to true for magnetic field

        self.n_lat_radial_slider = QSlider(Qt.Horizontal)
        self.n_lat_radial_slider.setMinimum(5)
        self.n_lat_radial_slider.setMaximum(50)
        self.n_lat_radial_slider.setValue(10)
        self.n_lat_radial_slider.setToolTip('Number of radial lines (vertical)')
        self.n_lat_radial_value_label = QLabel(f'{self.n_lat_radial_slider.value()}') # Label for N_lat slider

        self.n_lon_radial_slider = QSlider(Qt.Horizontal)
        self.n_lon_radial_slider.setMinimum(10)
        self.n_lon_radial_slider.setMaximum(100)
        self.n_lon_radial_slider.setValue(30)
        self.n_lon_radial_slider.setToolTip('Number of radial lines (horizontal)')
        self.n_lon_radial_value_label = QLabel(f'{self.n_lon_radial_slider.value()}') # Label for N_lon slider

        # GONG and PFSS controls
        self.select_gong_button = QPushButton('Select GONG File')
        self.gong_file_input = QLineEdit()
        self.gong_file_input.setPlaceholderText('No GONG file selected')
        self.gong_file_input.setReadOnly(True)

        self.pfss_button = QPushButton('Calculate PFSS')
        self.pfss_button.setEnabled(False) # Disable until a GONG file is selected

        # Matplotlib Figure and Canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.label = QLabel('Ready.')

        # --- Connect Signals ---
        self.load_button.clicked.connect(self.load_files)
        self.display_button.clicked.connect(self.display_first_map)
        self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
        self.run_diff_button.clicked.connect(self.create_running_diff_maps)
        self.prev_button.clicked.connect(self.show_prev_map)
        self.next_button.clicked.connect(self.show_next_map)
        self.save_button.clicked.connect(self.save_current_map)

        self.ellipse_button.clicked.connect(self.draw_ellipse_from_input)
        self.fit_button.clicked.connect(self.extract_ellipse_params)

        self.a_slider.valueChanged.connect(self.update_ellipse)
        self.a_slider.valueChanged.connect(lambda value: self.a_value_label.setText(f'{value} arcsec'))
        self.b_slider.valueChanged.connect(self.update_ellipse)
        self.b_slider.valueChanged.connect(lambda value: self.b_value_label.setText(f'{value} arcsec'))

        self.n_lat_radial_slider.valueChanged.connect(lambda value: self.n_lat_radial_value_label.setText(f'{value}'))
        self.n_lon_radial_slider.valueChanged.connect(lambda value: self.n_lon_radial_value_label.setText(f'{value}'))

        self.select_gong_button.clicked.connect(self.select_gong_file)
        self.pfss_button.clicked.connect(self.calculate_pfss_model)
        self.show_3d_button.clicked.connect(self.show_3d_ellipsoid_plot)


        # --- Layout Setup ---
        main_layout = QHBoxLayout() # Main layout to divide left and right sections

        # Left Section: Main Plot and Sliders
        left_section_layout = QVBoxLayout()

        # Top part of left section: Main plot and status label
        plot_and_status_layout = QVBoxLayout()
        plot_and_status_layout.addWidget(self.canvas, 5) # Increased stretch factor for canvas
        plot_and_status_layout.addWidget(self.label, 1) # Decreased stretch factor for label
        left_section_layout.addLayout(plot_and_status_layout)

        # Bottom part of left section: Sliders
        sliders_layout = QVBoxLayout()
        sliders_layout.addWidget(QLabel('Ellipse Parameters:'))
        
        a_slider_row = QHBoxLayout()
        a_slider_row.addWidget(QLabel('Semi-Major Axis (a):'))
        a_slider_row.addWidget(self.a_slider)
        a_slider_row.addWidget(self.a_value_label)
        sliders_layout.addLayout(a_slider_row)

        b_slider_row = QHBoxLayout()
        b_slider_row.addWidget(QLabel('Semi-Minor Axis (b):'))
        b_slider_row.addWidget(self.b_slider)
        b_slider_row.addWidget(self.b_value_label)
        sliders_layout.addLayout(b_slider_row)

        sliders_layout.addWidget(QLabel('Radial Line Density:'))

        n_lat_radial_slider_row = QHBoxLayout()
        n_lat_radial_slider_row.addWidget(QLabel('Vertical Lines:'))
        n_lat_radial_slider_row.addWidget(self.n_lat_radial_slider)
        n_lat_radial_slider_row.addWidget(self.n_lat_radial_value_label)
        sliders_layout.addLayout(n_lat_radial_slider_row)

        n_lon_radial_slider_row = QHBoxLayout()
        n_lon_radial_slider_row.addWidget(QLabel('Horizontal Lines:'))
        n_lon_radial_slider_row.addWidget(self.n_lon_radial_slider)
        n_lon_radial_slider_row.addWidget(self.n_lon_radial_value_label)
        sliders_layout.addLayout(n_lon_radial_slider_row)

        left_section_layout.addLayout(sliders_layout)
        main_layout.addLayout(left_section_layout)


        # Right Section: Buttons and Checkboxes (Menu)
        right_menu_layout = QVBoxLayout()
        right_menu_layout.setAlignment(Qt.AlignTop) # Align content to the top

        # File Operations
        right_menu_layout.addWidget(QLabel('--- File Operations ---'))
        right_menu_layout.addLayout(btn_layout) # Re-use existing btn_layout for main buttons

        # GONG/PFSS Section
        right_menu_layout.addWidget(QLabel('--- Magnetic Field ---'))
        right_menu_layout.addWidget(self.select_gong_button)
        right_menu_layout.addWidget(self.gong_file_input)
        right_menu_layout.addWidget(self.pfss_button)
        right_menu_layout.addWidget(self.show_field_lines_checkbox)

        # Ellipse Controls
        right_menu_layout.addWidget(QLabel('--- Ellipse Controls ---'))
        right_menu_layout.addWidget(QLabel('X center [arcsec]:'))
        right_menu_layout.addWidget(self.x_input)
        right_menu_layout.addWidget(QLabel('Y center [arcsec]:'))
        right_menu_layout.addWidget(self.y_input)
        right_menu_layout.addWidget(self.ellipse_button)
        right_menu_layout.addWidget(self.fit_button)

        # 3D Visualization Controls
        right_menu_layout.addWidget(QLabel('--- 3D Visualization ---'))
        right_menu_layout.addWidget(self.show_3d_button)
        right_menu_layout.addWidget(self.show_shell_checkbox)
        right_menu_layout.addWidget(self.show_normals_checkbox)
        right_menu_layout.addWidget(self.show_radials_checkbox)
        
        # Set a fixed width for the right menu
        right_menu_widget = QWidget()
        right_menu_widget.setLayout(right_menu_layout)
        right_menu_widget.setFixedWidth(250)

        main_layout.addWidget(right_menu_widget)
        
        self.setLayout(main_layout)

    def select_gong_file(self):
        """Opens a file dialog to select the GONG FITS file."""
        filepath, _ = QFileDialog.getOpenFileName(
            self, 'Select GONG FITS File', '', 'FITS files (*.fits *.fts)')
        if filepath:
            self.gong_filepath = filepath
            self.gong_file_input.setText(os.path.basename(filepath))
            self.pfss_button.setEnabled(True) # Enable PFSS button once file is selected
            self.label.setText(f'GONG file selected: {os.path.basename(filepath)}')
        else:
            self.gong_filepath = None
            self.gong_file_input.setText('No GONG file selected')
            self.pfss_button.setEnabled(False)
            self.label.setText('No GONG file selected.')

    def calculate_pfss_model(self):
        """
        Extrapolate the coronal magnetic field using the PFSS model
        and trace field lines using the selected GONG file.
        """
        if self.gong_filepath is None:
            QMessageBox.warning(self, 'Warning', 'Please select a GONG FITS file first.')
            return

        self.label.setText('Calculating PFSS model and tracing field lines...')
        QApplication.processEvents() # Update GUI immediately

        try:
            # Make a sunpy map of the file
            gong_map = sunpy.map.Map(self.gong_filepath)

            # Fix a bug in the GONG map file if 'cunit1' is missing
            if 'cunit1' not in gong_map.meta:
                gong_map.meta['cunit1'] = u.deg
            
            self.gong_map = gong_map # Store the loaded GONG map

            # Model the PFSS field lines
            nrho     = 70  # number of rho grid points
            rss      = 3   # source surface radius
            pfss_in  = pfsspy.Input(self.gong_map, nrho, rss)
            self.label.setText('Calculating PFSS model...')
            QApplication.processEvents()
            self.pfss_out = pfsspy.pfss(pfss_in)

            # Trace field lines
            num_footpoints_lat = 40
            num_footpoints_lon = 60
            # Start tracing from a radius slightly above the solar surface
            r_trace = 1.05 * const.radius

            lat = np.linspace(np.radians(-90), np.radians(90), num_footpoints_lat, endpoint=False)
            lon = np.linspace(np.radians(-180), np.radians(180), num_footpoints_lon, endpoint=False)

            lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')
            lat_flat, lon_flat = lat_grid.ravel()*u.rad, lon_grid.ravel()*u.rad

            # Make a 2D grid from these 1D points 
            seeds  = SkyCoord(lon_flat, lat_flat, r_trace, frame=self.pfss_out.coordinate_frame)
            tracer = tracing.FortranTracer()
            self.label.setText('Tracing magnetic field lines...')
            QApplication.processEvents()
            self.field_lines = tracer.trace(seeds, self.pfss_out)
            self.label.setText('PFSS model calculated and magnetic field data loaded.')

        except Exception as e:
            self.label.setText(f'Error calculating PFSS model: {e}')
            self.pfss_out = None
            self.field_lines = None
            self.gong_map = None # Clear gong_map if PFSS fails


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
        if files:
            self.files = files
            self.maps = [sunpy.map.Map(f) for f in files]
            self.current_index = 0
            self.label.setText(f'{len(files)} files loaded.')


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
        self.canvas.figure.clf() # Clear the figure
        self.ellipse_artist = None # Reset the ellipse artist reference
        
        ax = self.canvas.figure.add_subplot(111, projection=amap)

        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

        ax.grid(False)
        
        # If an ellipse was previously defined, redraw it on the new axes
        if self.ellipse_center is not None:
            # Pass the current map object for its coordinate frame
            self.draw_ellipse(self.ellipse_center[0], self.ellipse_center[1], current_map=amap, redraw_canvas=False)

        self.canvas.draw() # Draw the canvas once after all plotting
        self.label.setText(f'Showing: {self.files[self.current_index]}')


    def upgrade_to_lv15(self):
        if not self.files:
            QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
            return

        self.processed_maps = []

        for file in self.files:
            data_level = self.detect_data_level(file)
            try:
                if data_level == 'lev1':
                    output_filename = os.path.basename(file).replace('lev1', 'lev15')
                    file_path = os.path.join(self.temp_data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
                    os.makedirs(file_path, exist_ok=True)
                    full_path = os.path.join(file_path, output_filename)

                    if not os.path.exists(full_path):
                        m = sunpy.map.Map(file)
                        m = update_pointing(m)
                        m = register(m)
                        m = m / m.exposure_time
                        m.save(full_path, filetype='auto')
                        self.processed_maps.append(m)
                    else:
                        self.processed_maps.append(sunpy.map.Map(full_path))

                elif data_level == 'lev15':
                    self.processed_maps.append(sunpy.map.Map(file)) # Use original file if already level 1.5
                else:
                    print(f'Skipped unknown-level file: {file}')

            except Exception as e:
                print(f'Error upgrading {file}: {e}')

        self.maps = self.processed_maps
        self.current_index = 0
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

        for i in range(5, len(self.processed_maps)):
            try:
                m0 = self.processed_maps[i-5]
                m1 = self.processed_maps[i]
                diff = m1.quantity - m0.quantity
                smoothed = ndimage.gaussian_filter(diff, sigma=[3,3])
                diff_map = sunpy.map.Map(smoothed, m1.meta)
                diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
                self.running_diff_maps.append(diff_map)
            except Exception as e:
                print(f'Error creating diff map {i}: {e}')

        self.maps = self.running_diff_maps
        self.current_index = 0
        if self.maps:
            self.plot_map(self.maps[0])
        self.label.setText(f'Created {len(self.running_diff_maps)} running difference maps.')


    def draw_ellipse(self, x0, y0, current_map=None, redraw_canvas=True):
        """
        Draw or update the ellipse on the current plot using provided center coordinates.
        'current_map' can be optionally passed to avoid re-fetching if already known (e.g., from plot_map).
        """
        try:
            if not self.maps:
                self.label.setText('No map to draw ellipse on.')
                return

            if current_map is None:
                # If current_map is not passed, get it from the stored maps
                current_map = self.maps[self.current_index]

            # Ensure we have an active axes to draw on
            if not self.canvas.figure.axes:
                # This should ideally not happen if plot_map is always called first
                # but as a fallback, plot the current map to create axes.
                self.plot_map(current_map)
                ax = self.canvas.figure.axes[0]
            else:
                ax = self.canvas.figure.axes[0]  # Get current axes

            # Remove old ellipse if exists and is still attached to an axes
            if self.ellipse_artist is not None and self.ellipse_artist in ax.lines:
                self.ellipse_artist.remove()
                self.ellipse_artist = None # Clear reference after removal
            elif self.ellipse_artist is not None:
                 # If ellipse_artist exists but is not in current axes, it means
                 # the axes was cleared (e.g., by plot_map), so just clear the reference.
                 self.ellipse_artist = None

            # Update stored center
            self.ellipse_center = (x0, y0)
                
            a = self.a_slider.value()
            b = self.b_slider.value()
            
            # Calculate new ellipse coordinates
            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)
            
            # Use the passed or fetched current_map for the coordinate frame
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                             frame=current_map.coordinate_frame)
            
            # Draw new ellipse - correctly assign the Line2D object
            # ax.plot_coord returns a list, so unpack it
            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)
            
            if redraw_canvas:
                self.canvas.draw()
                self.label.setText('Ellipse updated.')
                
        except Exception as e:
            self.label.setText(f'Ellipse error: {str(e)}')


    def draw_ellipse_from_input(self):
        """Called when the 'Ellipse' button is clicked. Gets coordinates from input fields."""
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            # When drawing from input, we don't have the current map explicitly,
            # so draw_ellipse will fetch it from self.maps[self.current_index]
            self.draw_ellipse(x0, y0, redraw_canvas=True)
        except ValueError:
            self.label.setText('Invalid center coordinates for ellipse.')
        except IndexError:
            self.label.setText('Ellipse error: No map loaded to draw on.')


    def update_ellipse(self):
        """Update ellipse when sliders change"""
        # Only update if an ellipse has been drawn before and maps are loaded
        if self.ellipse_center is not None and self.maps:
            # Pass the current map explicitly to draw_ellipse for efficiency
            current_map = self.maps[self.current_index]
            self.draw_ellipse(self.ellipse_center[0], self.ellipse_center[1], current_map=current_map)
        elif not self.maps:
            self.label.setText('No map loaded to update ellipse on.')


    def extract_ellipse_params(self):
        try:
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            a = self.a_slider.value()
            b = self.b_slider.value()
            print(f'Fit extracted: center=({x0}, {y0}) arcsec, a={a}, b={b} arcsec')
            self.label.setText(f'Fit extracted: a={a}, b={b}')
        except ValueError:
            self.label.setText('Invalid ellipse parameters.')
    
    
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

    def create_3d_ellipsoid(self, ellipse_params, sunpy_map, show_shell=True, show_normals=False,
                            show_radials=False, n_lat_radial=10, n_lon_radial=30, show_field_lines=False):
        """Create 3D ellipsoid visualization using Plotly"""
        
        fig = go.Figure()
        
        # Create the Sun (sphere)
        theta, phi = np.mgrid[0:np.pi:50j, 0:2*np.pi:50j]
        x_sun = np.sin(theta) * np.cos(phi)
        y_sun = np.sin(theta) * np.sin(phi)
        z_sun = np.cos(theta)
        
        fig.add_trace(go.Surface(
            x=x_sun, y=y_sun, z=z_sun,
            colorscale=[[0, 'orange'], [1, 'orange']],
            showscale=False,
            name='Sun',
            hoverinfo='skip'
        ))
        
        # Get ellipse parameters
        x0, y0 = ellipse_params['x0'], ellipse_params['y0']
        a_radius, b_radius = ellipse_params['a'], ellipse_params['b']
        
        # Convert to heliographic coordinates (simplified)
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
        
        # Create ellipsoid surface
        # Normalizing by rsun_obs to scale to solar radii
        b_major = a_radius / sunpy_map.rsun_obs.value
        b_minor = b_radius / sunpy_map.rsun_obs.value
        
        shell_mesh_res = 50j # Resolution for ellipsoid surface mesh
        theta_src, phi_src = np.mgrid[0:np.pi:shell_mesh_res, 0:2*np.pi:shell_mesh_res]
        x_src = xshift + b_minor * np.sin(theta_src) * np.cos(phi_src)
        y_src = yshift + b_minor * np.sin(theta_src) * np.sin(phi_src)
        z_src = zshift + b_major * np.cos(theta_src)
        
        # Mask points inside the Sun
        r_shell = np.sqrt(x_src**2 + y_src**2 + z_src**2)
        mask = r_shell > 1 # Points outside the unit sphere (Sun)
        x_display = np.copy(x_src)
        y_display = np.copy(y_src)
        z_display = np.copy(z_src)
        x_display[~mask] = np.nan # Set points inside the Sun to NaN to make them invisible
        y_display[~mask] = np.nan
        z_display[~mask] = np.nan

        # Calculate normal vectors for the ellipsoid surface
        # Ellipsoid equation: ((x-x0)/b_minor)^2 + ((y-y0)/b_minor)^2 + ((z-z0)/b_major)^2 = 1
        # Gradient gives normal vector: (2(x-x0)/b_minor^2, 2(y-y0)/b_minor^2, 2(z-z0)/b_major^2)
        # We need normals only for the visible (masked) points
        nx_all = 2 * (x_src - xshift) / (b_minor**2)
        ny_all = 2 * (y_src - yshift) / (b_minor**2)
        nz_all = 2 * (z_src - zshift) / (b_major**2)
        
        norm_magnitude_all = np.sqrt(nx_all**2 + ny_all**2 + nz_all**2)
        nx_all /= norm_magnitude_all
        ny_all /= norm_magnitude_all
        nz_all /= norm_magnitude_all

        # Flatten the normal vectors for calculation with field lines
        nx_flat = nx_all[mask].flatten()
        ny_flat = ny_all[mask].flatten()
        nz_flat = nz_all[mask].flatten()
        
        # Flatten the visible shell points
        x_outer = x_src[mask].flatten()
        y_outer = y_src[mask].flatten()
        z_outer = z_src[mask].flatten()

        # Initialize all_points and all_dirs for theta calculation
        # This will hold the points and direction vectors of the lines we're comparing against normals
        all_points_for_theta = []
        all_dirs_for_theta = []

        # Radial directions for visualization and potential theta calculation
        if show_radials:
            radial_length = 2 # Length of radial lines in solar radii
            theta_radial, phi_radial = np.mgrid[0:np.pi:complex(n_lat_radial), 0:2*np.pi:complex(n_lon_radial)]
            rdx = np.sin(theta_radial) * np.cos(phi_radial)
            rdy = np.sin(theta_radial) * np.sin(phi_radial)
            rdz = np.cos(theta_radial)
            rdx_flat = rdx.flatten()
            rdy_flat = rdy.flatten()
            rdz_flat = rdz.flatten()
            
            for i in range(len(rdx_flat)):
                r_vec = np.array([rdx_flat[i], rdy_flat[i], rdz_flat[i]])
                r_pts = np.linspace(0, radial_length, 20) # Fewer points for plotting efficiency
                lx = r_pts * r_vec[0]
                ly = r_pts * r_vec[1]
                lz = r_pts * r_vec[2]
                
                fig.add_trace(go.Scatter3d(
                    x=lx, y=ly, z=lz,
                    mode='lines',
                    line=dict(color='gray', width=1),
                    showlegend=False,
                    hoverinfo='skip',
                    opacity=0.3
                ))
                
                # Store points and directions for theta calculation if no field lines
                # OR if both are shown, but field lines are prioritized for theta calculation
                if not show_field_lines or self.field_lines is None:
                    # Store the end point and direction of the radial line
                    all_points_for_theta.append(r_vec * radial_length) 
                    all_dirs_for_theta.append(r_vec)


        # Plot the magnetic field lines and populate all_points_for_theta/all_dirs_for_theta
        if show_field_lines and self.field_lines:
            # If magnetic field lines are shown, they take precedence for theta calculation
            all_points_for_theta = []
            all_dirs_for_theta = []

            for n, field_line in enumerate(self.field_lines):
                color = {0:'black', -1:'blue', 1:'red'}.get(field_line.polarity)
                coords = field_line.coords
                coords.representation_type = 'cartesian'
                x_field = coords.x / const.radius
                y_field = coords.y / const.radius
                z_field = coords.z / const.radius
                
                fig.add_trace(go.Scatter3d(
                    x=x_field, y=y_field, z=z_field,
                    mode='lines',
                    line=dict(color=color, width=2),
                    showlegend=False
                ))
                
                # Store points and directions for theta calculation
                for i in range(len(x_field)-1):
                    field_point = np.array([x_field[i].value, y_field[i].value, z_field[i].value])
                    all_points_for_theta.append(field_point)
                    
                    field_vector = np.array([x_field[i+1].value - x_field[i].value, 
                                                y_field[i+1].value - y_field[i].value, 
                                                z_field[i+1].value - z_field[i].value])
                    # Normalize field vector
                    field_vector = field_vector / np.linalg.norm(field_vector)
                    all_dirs_for_theta.append(field_vector)
        
        all_points_for_theta = np.array(all_points_for_theta)
        all_dirs_for_theta = np.array(all_dirs_for_theta)

        # Calculate theta angles between normal vectors and closest field/radial line
        theta_angles = np.full_like(x_outer, np.nan) # Initialize with NaN
        
        if len(all_points_for_theta) > 0 and len(x_outer) > 0:
            for i in range(len(x_outer)):
                surf_pt = np.array([x_outer[i], y_outer[i], z_outer[i]])
                normal_vec = np.array([nx_flat[i], ny_flat[i], nz_flat[i]])
                
                # Find the closest point in all_points_for_theta
                dists = np.linalg.norm(all_points_for_theta - surf_pt, axis=1)
                idx = np.argmin(dists)
                
                line_dir = all_dirs_for_theta[idx]
                
                # Calculate cosine of angle between normal and line direction
                cos_theta = np.dot(normal_vec, line_dir)
                cos_theta = np.clip(cos_theta, -1, 1) # Ensure within valid range for arccos
                theta_angles[i] = np.degrees(np.arccos(np.abs(cos_theta))) # Absolute value for angle 0-90

        # Map theta to surface
        theta_surface = np.full_like(x_src, np.nan)
        theta_surface[mask] = theta_angles
        
        # Show shell with theta angle coloring
        if show_shell:
            # Add shell center marker
            fig.add_trace(go.Scatter3d(
                x=[xshift], y=[yshift], z=[zshift],
                mode='markers',
                marker=dict(size=8, color='black'),
                showlegend=False,
                name='Shell Center'
            ))
            
            # Add shell surface with theta coloring
            fig.add_trace(go.Surface(
                x=x_display, y=y_display, z=z_display,
                surfacecolor=theta_surface,
                colorscale='Viridis', # Colormap for theta angles
                cmin=0, cmax=90, # Theta angle ranges from 0 to 90 degrees
                opacity=1,
                showscale=True,
                colorbar=dict(
                    title=dict(text='Theta Angle (degrees)', side='right'),
                    len=0.6
                ),
                name='Shell'
            ))
        
        # Show normal vectors
        if show_normals:
            # Sample normals to avoid too many arrows
            step = max(1, len(x_outer) // 100)  # Show max 100 arrows
            x_norm_sample = x_outer[::step]
            y_norm_sample = y_outer[::step]  
            z_norm_sample = z_outer[::step]
            nx_sample = nx_all[mask][::step] * 0.2 # Scale for visibility
            ny_sample = ny_all[mask][::step] * 0.2
            nz_sample = nz_all[mask][::step] * 0.2
            
            # Create line segments for normal vectors
            x_lines = []
            y_lines = []
            z_lines = []
            
            for i in range(len(x_norm_sample)):
                x_lines.extend([x_norm_sample[i], x_norm_sample[i] + nx_sample[i], None])
                y_lines.extend([y_norm_sample[i], y_norm_sample[i] + ny_sample[i], None])
                z_lines.extend([z_norm_sample[i], z_norm_sample[i] + nz_sample[i], None])
            
            fig.add_trace(go.Scatter3d(
                x=x_lines,
                y=y_lines,
                z=z_lines,
                mode='lines',
                line=dict(color='black', width=3),
                showlegend=False,
                name='Normal Vectors'
            ))

        # Layout settings
        fig.update_layout(
            scene=dict(
                xaxis=dict(range=[-2, 2], title='X (Solar Radii)'),
                yaxis=dict(range=[-2, 2], title='Y (Solar Radii)'),
                zaxis=dict(range=[-2, 2], title='Z (Solar Radii)'),
                aspectmode='cube',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            width=1024,
            height=768,
            title='Interactive 3D Solar Shell Visualization',
            showlegend=True # Show legend for Sun and Ellipsoid Shell
        )
        
        return fig

    def show_3d_ellipsoid_plot(self):
        """
        Gathers ellipse parameters, creates a 3D Plotly figure, and displays it
        by saving to a temporary HTML file and opening in the default web browser.
        This is a workaround for "Could not find QtWebEngineProcess" error.
        """
        if self.ellipse_center is None or not self.maps:
            QMessageBox.warning(self, 'Warning', 'Please draw an ellipse and load a map first.')
            return
        
        # Check if magnetic field lines are requested but not loaded
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

            # Create the Plotly figure
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

            # --- WORKAROUND: Save to temporary HTML and open in browser ---
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as f:
                temp_filepath = f.name
                pio.write_html(fig, file=f, auto_open=False, include_plotlyjs='cdn')
            
            webbrowser.open(f'file://{temp_filepath}')
            QMessageBox.information(self, '3D Plot Displayed', 
                                    f'The 3D ellipsoid plot has been opened in your default web browser at:\n{temp_filepath}\n\n'
                                    'This is a workaround for the "Could not find QtWebEngineProcess" error.')
            # --- END WORKAROUND ---

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to create 3D plot: {str(e)}')


def main():
    app = QApplication(sys.argv)
    viewer = AIAViewer()
    viewer.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
