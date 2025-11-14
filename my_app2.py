import os
import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider
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

class AIAViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('AIA FITS Viewer')
        self.setGeometry(100, 100, 1200, 1000)

        self.files = []
        self.maps = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0
        self.data_dir = os.getcwd()
        self.channel = '193'
        
        # Ellipse parameters
        self.ellipse_center = None
        self.ellipse_artist = None
        self.theta = np.linspace(0, 2*np.pi, 300)
        
        # Initialize UI
        self.init_ui()
    

    def init_ui(self):
        # Create widgets
        self.load_button = QPushButton('Load FITS')
        self.display_button = QPushButton('Display')
        self.upgrade_button = QPushButton('Upgrade')
        self.run_diff_button = QPushButton('Run-Diff')
        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')
        self.save_button = QPushButton('Save')
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

        # Connect signals
        self.load_button.clicked.connect(self.load_files)
        self.display_button.clicked.connect(self.display_first_map)
        self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
        self.run_diff_button.clicked.connect(self.create_running_diff_maps)
        self.prev_button.clicked.connect(self.show_prev_map)
        self.next_button.clicked.connect(self.show_next_map)
        self.save_button.clicked.connect(self.save_current_map)
        self.ellipse_button.clicked.connect(self.draw_ellipse)
        self.fit_button.clicked.connect(self.extract_ellipse_params)
        self.a_slider.valueChanged.connect(self.update_ellipse)
        self.b_slider.valueChanged.connect(self.update_ellipse)

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
        
        ellipse_layout = QVBoxLayout()
        ellipse_layout.addWidget(QLabel('Ellipse Center:'))
        ellipse_layout.addWidget(self.x_input)
        ellipse_layout.addWidget(self.y_input)
        ellipse_layout.addWidget(QLabel('Semi-Major Axis:'))
        ellipse_layout.addWidget(self.a_slider)
        ellipse_layout.addWidget(QLabel('Semi-Minor Axis:'))
        ellipse_layout.addWidget(self.b_slider)
        ellipse_btn_layout = QHBoxLayout()
        ellipse_btn_layout.addWidget(self.ellipse_button)
        ellipse_btn_layout.addWidget(self.fit_button)
        ellipse_layout.addLayout(ellipse_btn_layout)
        
        layout.addLayout(ellipse_layout)
        self.setLayout(layout)
    

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
        self.figure.clf()
        ax = self.figure.add_subplot(111, projection=amap)

        # Plot the map
        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)
        
        ax.grid(False)
        
        # Clear ellipse artist reference when changing maps
        self.ellipse_artist = None
        
        # Redraw ellipse if parameters exist
        if self.ellipse_center is not None:
            self.draw_ellipse(redraw=False)
        
        self.canvas.draw()
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

        self.maps = self.processed_maps
        self.current_index = 0
        self.plot_map(self.maps[0])
        self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
    

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


    def draw_ellipse(self, redraw=True):
        """Draw or update the ellipse on the current plot"""
        try:
            # Get current map and axes
            amap = self.maps[self.current_index]
            ax = self.figure.axes[0]
            
            # Get parameters from inputs
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            self.ellipse_center = (x0, y0)
            
            a = self.a_slider.value()
            b = self.b_slider.value()
            
            # Calculate new ellipse coordinates
            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                             frame=amap.coordinate_frame)
            
            # Remove old ellipse if it exists
            if self.ellipse_artist is not None:
                try:
                    self.ellipse_artist.remove()
                except:
                    pass  # Artist might already be removed
                self.ellipse_artist = None
            
            # Draw new ellipse
            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)
            
            if redraw:
                self.canvas.draw()
                self.label.setText('Ellipse drawn.')
                
        except ValueError:
            self.label.setText('Invalid center coordinates.')
        except Exception as e:
            self.label.setText(f'Ellipse error: {str(e)}')


    def update_ellipse(self):
        """Update ellipse when sliders change"""
        if self.ellipse_center is None:
            # Try to get center from inputs if it exists
            try:
                x0 = float(self.x_input.text())
                y0 = float(self.y_input.text())
                self.ellipse_center = (x0, y0)
            except ValueError:
                return  # No valid center yet
            
        if not self.maps or not self.figure.axes:
            return  # No map or axes available
            
        try:
            # Get current map and axes
            amap = self.maps[self.current_index]
            ax = self.figure.axes[0]
            
            # Get parameters
            x0, y0 = self.ellipse_center
            a = self.a_slider.value()
            b = self.b_slider.value()
            
            # Calculate new ellipse coordinates
            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                             frame=amap.coordinate_frame)
            
            # If artist exists, update it
            if self.ellipse_artist is not None:
                # Update existing artist
                self.ellipse_artist.set_data(coords.Tx.value, coords.Ty.value)
            else:
                # Create new artist
                self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)
            
            self.canvas.draw()
            
        except Exception as e:
            self.label.setText(f'Ellipse update error: {str(e)}')


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


def main():
    app = QApplication(sys.argv)
    viewer = AIAViewer()
    viewer.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()