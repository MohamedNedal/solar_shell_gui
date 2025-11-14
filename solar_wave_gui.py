"""
Solar Wave Analysis GUI
A comprehensive tool for analyzing AIA FITS files and fitting ellipses to coronal waves.
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.colors as colors
from scipy import ndimage
import plotly.graph_objects as go
import plotly.offline as pyo

from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                            QWidget, QPushButton, QLabel, QSlider, QFileDialog, 
                            QMessageBox, QProgressBar, QComboBox, QSpinBox,
                            QGroupBox, QGridLayout, QTextEdit, QSplitter)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QPixmap

import sunpy.map
from sunpy.net import Fido, attrs as a
from astropy.coordinates import SkyCoord
from aiapy.calibrate import register, update_pointing
import astropy.units as u
import astropy.constants as const


class DataProcessor(QThread):
    """Background thread for processing FITS files"""
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(list)
    
    def __init__(self, files, channel, data_dir):
        super().__init__()
        self.files = files
        self.channel = channel
        self.data_dir = data_dir


    def run(self):
        os.makedirs(f'{self.data_dir}/AIA/{self.channel}A/processed/lv15', exist_ok=True)
        
        aia_maps = []
        total_files = len(self.files)
        
        for i, file in enumerate(self.files):
            try:
                self.status.emit(f'Processing file {i+1}/{total_files}...')
                
                # Check data level from filename or header
                data_level = self.detect_data_level(file)
                
                if data_level == 'lev1':
                    # Process lev1 to lev1.5
                    output_filename = os.path.basename(file).replace('lev1', 'lev15')
                    file_path = f'{self.data_dir}/AIA/{self.channel}A/processed/lv15/{output_filename}'
                    
                    if not os.path.exists(file_path):
                        # Load and process the file
                        m = sunpy.map.Map(file)
                        m_updated = update_pointing(m)
                        m_registered = register(m_updated)
                        m_normalized = m_registered / m_registered.exposure_time
                        m_normalized.save(file_path, filetype='auto')
                        aia_maps.append(m_normalized)
                    else:
                        # Load existing processed file
                        aia_maps.append(sunpy.map.Map(file_path))
                        
                elif data_level == 'lev15':
                    # Data is already lev1.5, load directly
                    self.status.emit(f'File {i+1} is already lev1.5, loading directly...')
                    aia_maps.append(sunpy.map.Map(file))
                    
                else:
                    # Unknown data level
                    self.status.emit(f'Warning: Unknown data level for file {i+1}, attempting to load...')
                    aia_maps.append(sunpy.map.Map(file))
                
                progress_pct = int((i + 1) / total_files * 100)
                self.progress.emit(progress_pct)
                
            except Exception as e:
                self.status.emit(f'Error processing file {i+1}: {str(e)}')
                
        self.finished.emit(aia_maps)


    def detect_data_level(self, file_path):
        """
        Detect the data level of an AIA FITS file.
        Returns 'lev1', 'lev15', or 'unknown'
        """
        try:
            # First check filename for level indicators
            filename = os.path.basename(file_path).lower()
            
            if 'lev1' in filename and 'lev15' not in filename:
                return 'lev1'
            elif 'lev15' in filename or 'lev1.5' in filename:
                return 'lev15'
            
            # If filename doesn't indicate level, check FITS header
            from astropy.io import fits
            with fits.open(file_path) as hdul:
                header = hdul[0].header
                
                # Check for LVNUM keyword (common in AIA data)
                if 'LVNUM' in header:
                    lvnum = header['LVNUM']
                    if lvnum == 1.0:
                        return 'lev1'
                    elif lvnum == 1.5:
                        return 'lev15'
                
                # Check for other level indicators in header
                if 'LV_NUM' in header:
                    lv_num = header['LV_NUM']
                    if lv_num == 1.0:
                        return 'lev1'
                    elif lv_num == 1.5:
                        return 'lev15'
                
                # Check for processing history keywords
                if 'HISTORY' in header:
                    history = str(header.get('HISTORY', ''))
                    if 'register' in history.lower() or 'lev1.5' in history.lower():
                        return 'lev15'
                
                # Additional checks for calibration keywords that indicate lev1.5
                calibration_keywords = ['RSUN_REF', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']
                if all(keyword in header for keyword in calibration_keywords):
                    # Check if values look like they've been processed
                    if abs(header.get('CDELT1', 0) - 0.6) < 0.1:  # Typical lev1.5 pixel scale
                        return 'lev15'
            
            return 'unknown'
            
        except Exception as e:
            print(f"Error detecting data level for {file_path}: {e}")
            return 'unknown'


class MatplotlibWidget(QWidget):
    """Custom widget for embedding matplotlib plots"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=[10,8])
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        
    def clear(self):
        self.figure.clear()
        self.canvas.draw()


class EllipseFittingWidget(QWidget):
    """Widget for interactive ellipse fitting with sliders"""
    
    fitting_done = pyqtSignal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_map = None
        self.ellipse_params = {}
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Plot area
        self.plot_widget = MatplotlibWidget()
        layout.addWidget(self.plot_widget)
        
        # Controls
        controls_group = QGroupBox("Ellipse Parameters")
        controls_layout = QGridLayout()
        
        # X Center slider
        controls_layout.addWidget(QLabel("X Center:"), 0, 0)
        self.x_slider = QSlider(Qt.Horizontal)
        self.x_slider.setRange(-2000, 2000)
        self.x_slider.setValue(0)
        self.x_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.x_slider, 0, 1)
        self.x_label = QLabel("0")
        controls_layout.addWidget(self.x_label, 0, 2)
        
        # Y Center slider
        controls_layout.addWidget(QLabel("Y Center:"), 1, 0)
        self.y_slider = QSlider(Qt.Horizontal)
        self.y_slider.setRange(-2000, 2000)
        self.y_slider.setValue(0)
        self.y_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.y_slider, 1, 1)
        self.y_label = QLabel("0")
        controls_layout.addWidget(self.y_label, 1, 2)
        
        # Semi-major axis slider
        controls_layout.addWidget(QLabel("Semi-major:"), 2, 0)
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setRange(10, 1000)
        self.a_slider.setValue(200)
        self.a_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.a_slider, 2, 1)
        self.a_label = QLabel("200")
        controls_layout.addWidget(self.a_label, 2, 2)
        
        # Semi-minor axis slider
        controls_layout.addWidget(QLabel("Semi-minor:"), 3, 0)
        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setRange(10, 1000)
        self.b_slider.setValue(100)
        self.b_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.b_slider, 3, 1)
        self.b_label = QLabel("100")
        controls_layout.addWidget(self.b_label, 3, 2)
        
        controls_group.setLayout(controls_layout)
        layout.addWidget(controls_group)
        
        # Action buttons
        button_layout = QHBoxLayout()
        self.done_button = QPushButton("Done")
        self.done_button.clicked.connect(self.finish_fitting)
        button_layout.addWidget(self.done_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)
        
    def set_map(self, sunpy_map):
        """Set the current SunPy map for fitting"""
        self.current_map = sunpy_map
        self.plot_map()
        
    def plot_map(self):
        """Plot the current map with ellipse overlay"""
        if self.current_map is None:
            return
            
        self.plot_widget.figure.clear()
        ax = self.plot_widget.figure.add_subplot(111, projection=self.current_map)
        
        self.current_map.plot(axes=ax)
        self.update_ellipse()
        
        ax.set_title('Fit Ellipse to Coronal Wave')
        self.plot_widget.canvas.draw()
        
    def update_ellipse(self):
        """Update ellipse visualization based on slider values"""
        if self.current_map is None:
            return
            
        x0 = self.x_slider.value()
        y0 = self.y_slider.value()
        a = self.a_slider.value()
        b = self.b_slider.value()
        
        # Update labels
        self.x_label.setText(str(x0))
        self.y_label.setText(str(y0))
        self.a_label.setText(str(a))
        self.b_label.setText(str(b))
        
        # Clear previous ellipse and redraw
        ax = self.plot_widget.figure.axes[0] if self.plot_widget.figure.axes else None
        if ax is None:
            return
            
        # Remove previous ellipse lines
        for line in ax.lines:
            if hasattr(line, '_is_ellipse'):
                line.remove()
        
        # Draw new ellipse
        theta = np.linspace(0, 2*np.pi, 300)
        x_ell = x0 + a * np.cos(theta)
        y_ell = y0 + b * np.sin(theta)
        
        coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                        frame=self.current_map.coordinate_frame)
        line = ax.plot_coord(coords, color='red', lw=2)[0]
        line._is_ellipse = True
        
        ax.grid(False)

        self.plot_widget.canvas.draw()
        
    def finish_fitting(self):
        """Store fitting parameters and emit signal"""
        self.ellipse_params = {
            'x0': self.x_slider.value(),
            'y0': self.y_slider.value(),
            'a': self.a_slider.value(),
            'b': self.b_slider.value()
        }
        self.fitting_done.emit(self.ellipse_params)


class Visualization3D:
    """Class for creating 3D visualizations"""
    
    @staticmethod
    def create_3d_ellipsoid(ellipse_params, sunpy_map, show_shell=True, show_normals=True):
        """Create 3D ellipsoid visualization"""
        
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
        b_major = a_radius / sunpy_map.rsun_obs.value
        b_minor = b_radius / sunpy_map.rsun_obs.value
        
        theta_src, phi_src = np.mgrid[0:np.pi:50j, 0:2*np.pi:50j]
        x_src = xshift + b_minor * np.sin(theta_src) * np.cos(phi_src)
        y_src = yshift + b_minor * np.sin(theta_src) * np.sin(phi_src)
        z_src = zshift + b_major * np.cos(theta_src)
        
        # Mask points inside the Sun
        r_shell = np.sqrt(x_src**2 + y_src**2 + z_src**2)
        mask = r_shell > 1
        x_display = np.copy(x_src)
        y_display = np.copy(y_src)
        z_display = np.copy(z_src)
        x_display[~mask] = np.nan
        y_display[~mask] = np.nan
        z_display[~mask] = np.nan
        
        if show_shell:
            fig.add_trace(go.Surface(
                x=x_display, y=y_display, z=z_display,
                colorscale='Viridis',
                opacity=1,
                name='Ellipsoid Shell'
            ))
        
        # Layout settings
        fig.update_layout(
            scene=dict(
                xaxis=dict(range=[-2, 2]),
                yaxis=dict(range=[-2, 2]),
                zaxis=dict(range=[-2, 2]),
                aspectmode='cube',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            width=1024,
            height=768,
            title='3D Solar Wave Ellipsoid',
            showlegend=False
        )
        
        return fig


class MainWindow(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        self.files = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_image_index = 0
        self.ellipse_fits = []
        self.data_dir = os.getcwd()
        
        self.setup_ui()
        self.setWindowTitle('Solar Wave Analysis Tool')
        self.setGeometry(100, 100, 1400, 900)
        
    def setup_ui(self):
        """Setup the user interface"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout()
        
        # Top controls
        controls_layout = QHBoxLayout()
        
        self.select_files_btn = QPushButton('Select Files')
        self.select_files_btn.clicked.connect(self.select_files)
        controls_layout.addWidget(self.select_files_btn)
        
        # Channel selection
        controls_layout.addWidget(QLabel('Channel:'))
        self.channel_combo = QComboBox()
        self.channel_combo.addItems(['171', '193', '211', '304', '335', '94', '131'])
        self.channel_combo.setCurrentText('193')
        controls_layout.addWidget(self.channel_combo)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        controls_layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel('Ready')
        controls_layout.addWidget(self.status_label)
        
        controls_layout.addStretch()
        main_layout.addLayout(controls_layout)
        
        # Main content area with splitter
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - Image fitting
        self.fitting_widget = EllipseFittingWidget()
        self.fitting_widget.fitting_done.connect(self.on_fitting_done)
        splitter.addWidget(self.fitting_widget)
        
        # Right panel - Controls and info
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        
        # Image navigation
        nav_group = QGroupBox("Image Navigation")
        nav_layout = QVBoxLayout()
        
        nav_buttons = QHBoxLayout()
        self.done_btn = QPushButton('Done')
        self.done_btn.clicked.connect(self.on_fitting_done)
        self.done_btn.setEnabled(False)
        nav_buttons.addWidget(self.done_btn)
        
        self.next_btn = QPushButton('Next')
        self.next_btn.clicked.connect(self.next_image)
        self.next_btn.setEnabled(False)
        nav_buttons.addWidget(self.next_btn)
        
        self.finished_btn = QPushButton('Finished')
        self.finished_btn.clicked.connect(self.finish_fitting_session)
        self.finished_btn.setEnabled(False)
        nav_buttons.addWidget(self.finished_btn)
        
        nav_layout.addLayout(nav_buttons)
        
        self.image_info_label = QLabel('No images loaded')
        nav_layout.addWidget(self.image_info_label)
        
        nav_group.setLayout(nav_layout)
        right_layout.addWidget(nav_group)
        
        # 3D Visualization controls
        viz_group = QGroupBox("3D Visualization")
        viz_layout = QVBoxLayout()
        
        viz_buttons = QHBoxLayout()
        self.make_3d_btn = QPushButton('Make 3D')
        self.make_3d_btn.clicked.connect(self.create_3d_visualization)
        self.make_3d_btn.setEnabled(False)
        viz_buttons.addWidget(self.make_3d_btn)
        
        self.save_btn = QPushButton('Save')
        self.save_btn.clicked.connect(self.save_visualization)
        self.save_btn.setEnabled(False)
        viz_buttons.addWidget(self.save_btn)
        
        self.exit_btn = QPushButton('Exit')
        self.exit_btn.clicked.connect(self.close)
        viz_buttons.addWidget(self.exit_btn)
        
        viz_layout.addLayout(viz_buttons)
        viz_group.setLayout(viz_layout)
        right_layout.addWidget(viz_group)
        
        # Results display
        results_group = QGroupBox("Fitting Results")
        results_layout = QVBoxLayout()
        
        self.results_text = QTextEdit()
        self.results_text.setMaximumHeight(200)
        results_layout.addWidget(self.results_text)
        
        results_group.setLayout(results_layout)
        right_layout.addWidget(results_group)
        
        right_layout.addStretch()
        right_panel.setLayout(right_layout)
        splitter.addWidget(right_panel)
        
        # Set splitter proportions
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(splitter)
        central_widget.setLayout(main_layout)
        
    def select_files(self):
        """Open file dialog to select FITS files"""
        files, _ = QFileDialog.getOpenFileNames(
            self, 'Select AIA FITS Files', '',
            'FITS files (*.fits);;All files (*)')
        
        if files:
            self.files = sorted(files)
            self.status_label.setText(f'Selected {len(files)} files')
            self.process_files()
            
    def process_files(self):
        """Process selected FITS files in background thread"""
        if not self.files:
            return
            
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        
        channel = int(self.channel_combo.currentText())
        
        # Start background processing
        self.processor = DataProcessor(self.files, channel, self.data_dir)
        self.processor.progress.connect(self.progress_bar.setValue)
        self.processor.status.connect(self.status_label.setText)
        self.processor.finished.connect(self.on_processing_finished)
        self.processor.start()
        
    def on_processing_finished(self, processed_maps):
        """Handle completion of file processing"""
        self.processed_maps = processed_maps
        self.progress_bar.setVisible(False)
        
        if processed_maps:
            # Create running difference maps
            self.create_running_diff_maps()
            self.current_image_index = 0
            self.show_current_image()
            self.update_ui_state()
            self.status_label.setText(f'Loaded {len(processed_maps)} maps')
        else:
            self.status_label.setText('No maps could be processed')
            
    def create_running_diff_maps(self):
        """Create running difference maps from processed data"""
        if len(self.processed_maps) < 6:
            QMessageBox.warning(self, 'Warning', 
                              'Need at least 6 images for running difference')
            return
        
        self.running_diff_maps = []
        
        for i in range(5, len(self.processed_maps)):
            try:
                map_0 = self.processed_maps[i-5]
                map_1 = self.processed_maps[i]
                
                diff = map_1.quantity - map_0.quantity
                smoothed = ndimage.gaussian_filter(diff, sigma=[3,3])
                
                diff_map = sunpy.map.Map(smoothed, map_1.meta)
                diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
                
                self.running_diff_maps.append(diff_map)
            except Exception as e:
                print(f"Error creating running difference map {i}: {e}")
                
        self.status_label.setText(f'Created {len(self.running_diff_maps)} running difference maps')
        
    def show_current_image(self):
        """Display the current running difference image"""
        if not self.running_diff_maps or self.current_image_index >= len(self.running_diff_maps):
            return
            
        current_map = self.running_diff_maps[self.current_image_index]
        self.fitting_widget.set_map(current_map)
        
        self.image_info_label.setText(
            f'Image {self.current_image_index + 1} of {len(self.running_diff_maps)}')
        
    def on_fitting_done(self, params=None):
        """Handle completion of ellipse fitting for current image"""
        if params is None:
            params = self.fitting_widget.ellipse_params
            
        self.ellipse_fits.append({
            'image_index': self.current_image_index,
            'parameters': params
        })
        
        # Update results display
        results_text = f"Image {self.current_image_index + 1}:\n"
        results_text += f"  X Center: {params['x0']}\n"
        results_text += f"  Y Center: {params['y0']}\n" 
        results_text += f"  Semi-major: {params['a']}\n"
        results_text += f"  Semi-minor: {params['b']}\n\n"
        
        self.results_text.append(results_text)
        
        # Enable navigation buttons
        self.next_btn.setEnabled(True)
        self.finished_btn.setEnabled(True)
        
    def next_image(self):
        """Move to the next image for fitting"""
        if self.current_image_index < len(self.running_diff_maps) - 1:
            self.current_image_index += 1
            self.show_current_image()
            self.next_btn.setEnabled(False)
            self.done_btn.setEnabled(True)
        else:
            QMessageBox.information(self, 'Info', 'No more images to process')
            
    def finish_fitting_session(self):
        """Complete the fitting session and enable 3D visualization"""
        if self.ellipse_fits:
            self.make_3d_btn.setEnabled(True)
            self.status_label.setText('Fitting complete. Ready for 3D visualization.')
        else:
            QMessageBox.warning(self, 'Warning', 'No ellipse fits completed')
            
    def create_3d_visualization(self):
        """Create and display 3D ellipsoid visualization"""
        if not self.ellipse_fits:
            QMessageBox.warning(self, 'Warning', 'No fitting data available')
            return
            
        try:
            # Use the first fit for 3D visualization
            first_fit = self.ellipse_fits[0]
            params = first_fit['parameters']
            image_index = first_fit['image_index']
            
            if image_index < len(self.running_diff_maps):
                current_map = self.running_diff_maps[image_index]
                
                fig = Visualization3D.create_3d_ellipsoid(params, current_map)
                
                if fig:
                    # Save as HTML and open in browser
                    html_file = 'solar_wave_3d.html'
                    pyo.plot(fig, filename=html_file, auto_open=True)
                    self.save_btn.setEnabled(True)
                    self.status_label.setText('3D visualization created')
                else:
                    QMessageBox.warning(self, 'Error', 'Could not create 3D visualization')
                    
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Error creating 3D visualization: {str(e)}')
            
    def save_visualization(self):
        """Save the current visualization"""
        filename, _ = QFileDialog.getSaveFileName(
            self, 'Save Visualization', 'solar_wave_viz.png',
            'PNG files (*.png);;PDF files (*.pdf);;All files (*)')
        
        if filename:
            try:
                # Save the current matplotlib figure
                if self.fitting_widget.plot_widget.figure.axes:
                    self.fitting_widget.plot_widget.figure.savefig(filename, dpi=300, bbox_inches='tight')
                    self.status_label.setText(f'Saved to {filename}')
                else:
                    QMessageBox.warning(self, 'Warning', 'No visualization to save')
            except Exception as e:
                QMessageBox.critical(self, 'Error', f'Error saving file: {str(e)}')
                
    def update_ui_state(self):
        """Update UI button states based on current state"""
        has_maps = len(self.running_diff_maps) > 0
        self.done_btn.setEnabled(has_maps)
        
    def closeEvent(self, event):
        """Handle application close event"""
        reply = QMessageBox.question(self, 'Exit', 'Are you sure you want to exit?',
                                   QMessageBox.Yes | QMessageBox.No,
                                   QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


def main():
    """Main application entry point"""
    app = QApplication(sys.argv)
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

