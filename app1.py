import os
import sys
import tempfile # For creating temporary files
import webbrowser # For opening HTML files in a browser
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider, QCheckBox, QDialog
)
from PyQt5.QtCore import Qt, QUrl
# Commented out QWebEngineView as a workaround for "Could not find QtWebEngineProcess" error
# from PyQt5.QtWebEngineWidgets import QWebEngineView 

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
import sunpy.map
import numpy as np
from aiapy.calibrate import register, update_pointing
from scipy import ndimage
import astropy.units as u
from astropy.coordinates import SkyCoord

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
        self.data_dir = os.getcwd()
        self.channel = '193'
        
        self.theta = np.linspace(0, 2*np.pi, 300)

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

        # Sliders for axes
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setMinimum(10)
        self.a_slider.setMaximum(1000)
        self.a_slider.setValue(200)
        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setMinimum(10)
        self.b_slider.setMaximum(1000)
        self.b_slider.setValue(100)

        # Connect sliders to update_ellipse
        self.a_slider.valueChanged.connect(self.update_ellipse)
        self.b_slider.valueChanged.connect(self.update_ellipse)

        # Connect the buttons
        self.ellipse_button.clicked.connect(self.draw_ellipse_from_input)
        self.fit_button.clicked.connect(self.extract_ellipse_params)

        # 3D Ellipsoid functionality
        self.show_3d_button = QPushButton('Show 3D Ellipsoid')
        self.show_3d_button.clicked.connect(self.show_3d_ellipsoid_plot)

        self.show_shell_checkbox = QCheckBox('Show Ellipsoid Shell')
        self.show_shell_checkbox.setChecked(True) # Default to true
        self.show_normals_checkbox = QCheckBox('Show Normals (Future)')
        self.show_normals_checkbox.setChecked(False) # Default to false

        # Layouts
        btn_layout = QHBoxLayout()
        for b in [self.load_button, self.display_button, self.upgrade_button,
                  self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
            btn_layout.addWidget(b)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.label = QLabel('Ready.')

        layout = QVBoxLayout()
        layout.addLayout(btn_layout)
        layout.addWidget(self.canvas)
        layout.addWidget(self.label)
        self.setLayout(layout)

        ellipse_layout = QVBoxLayout()
        ellipse_layout.addWidget(QLabel('Ellipse Center:'))
        ellipse_layout.addWidget(self.x_input)
        ellipse_layout.addWidget(QLabel('Y center [arcsec]:')) # Added label for clarity
        ellipse_layout.addWidget(self.y_input)
        ellipse_layout.addWidget(QLabel('Semi-Major Axis:'))
        ellipse_layout.addWidget(self.a_slider)
        ellipse_layout.addWidget(QLabel('Semi-Minor Axis:'))
        ellipse_layout.addWidget(self.b_slider)
        ellipse_layout.addWidget(self.ellipse_button)
        ellipse_layout.addWidget(self.fit_button)
        layout.addLayout(ellipse_layout)

        # Add 3D ellipsoid controls to the layout
        _3d_layout = QVBoxLayout()
        _3d_layout.addWidget(self.show_3d_button)
        _3d_layout.addWidget(self.show_shell_checkbox)
        _3d_layout.addWidget(self.show_normals_checkbox) # Checkbox for future use
        layout.addLayout(_3d_layout)
    

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
        self.ellipse_artist = None # Crucially, reset the ellipse artist reference
        
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
                    file_path = os.path.join(self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
                    os.makedirs(file_path, exist_ok=True)
                    full_path = os.path.join(file_path, output_filename)

                    if not os.path.exists(full_path):
                        m = sunpy.map.Map(file)
                        # aiapy now requires `pointing_table` kwarg, so just skip this for now
                        m = update_pointing(m, pointing_table=None)
                        m = register(m)
                        m = m / m.exposure_time
                        m.save(full_path, filetype='auto')
                        self.processed_maps.append(m)
                    else:
                        self.processed_maps.append(sunpy.map.Map(full_path))

                elif data_level == 'lev15':
                    self.processed_maps.append(sunpy.map.Map(full_path)) # Fixed: use full_path here
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
                smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
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

    def create_3d_ellipsoid(self, ellipse_params, sunpy_map, show_shell=True, show_normals=True):
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
        
        shell_lon = np.deg2rad(hgs_coord.lon.value)
        shell_lat = np.deg2rad(hgs_coord.lat.value)
        
        xshift = np.cos(shell_lat) * np.cos(shell_lon)
        yshift = np.cos(shell_lat) * np.sin(shell_lon)
        zshift = np.sin(shell_lat)
        
        # Create ellipsoid surface
        # Normalizing by rsun_obs to scale to solar radii
        b_major = a_radius / sunpy_map.rsun_obs.value
        b_minor = b_radius / sunpy_map.rsun_obs.value
        
        theta_src, phi_src = np.mgrid[0:np.pi:50j, 0:2*np.pi:50j]
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
        
        if show_shell:
            fig.add_trace(go.Surface(
                x=x_display, y=y_display, z=z_display,
                colorscale='Viridis',
                opacity=0.8, # Slightly transparent for better visualization
                name='Ellipsoid Shell',
                showscale=False # No color scale for the ellipsoid
            ))
        
        # Normals (Future implementation)
        if show_normals:
            # This part is currently a placeholder.
            # You would add code here to calculate and plot normal vectors
            # For example, using go.Cone or go.Scatter3d with mode='lines'
            pass

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
            title='3D Solar Wave Ellipsoid',
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

            # Create the Plotly figure
            fig = self.create_3d_ellipsoid(ellipse_params, current_map, show_shell, show_normals)

            # --- WORKAROUND: Save to temporary HTML and open in browser ---
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as f:
                temp_filepath = f.name
                pio.write_html(fig, file=f, auto_open=False, include_plotlyjs='cdn')
            
            webbrowser.open(f'file://{temp_filepath}')
            QMessageBox.information(self, '3D Plot Displayed', 
                                    f'The 3D ellipsoid plot has been opened in your default web browser at:\n{temp_filepath}\n\n'
                                    'This is a workaround for the "Could not find QtWebEngineProcess" error.')
            # --- END WORKAROUND ---

            # Original QWebEngineView code (commented out for workaround)
            # plot_dialog = QDialog(self)
            # plot_dialog.setWindowTitle('3D Ellipsoid Visualization')
            # plot_dialog.setGeometry(150, 150, 1050, 800)
            # web_view = QWebEngineView()
            # html_content = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')
            # web_view.setHtml(html_content)
            # dialog_layout = QVBoxLayout()
            # dialog_layout.addWidget(web_view)
            # plot_dialog.setLayout(dialog_layout)
            # plot_dialog.exec_()

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to create 3D plot: {str(e)}')


def main():
    app = QApplication(sys.argv)
    viewer = AIAViewer()
    viewer.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
