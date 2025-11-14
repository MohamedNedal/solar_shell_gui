"""
Ellipse fitting widget for interactive wave front analysis
"""

import numpy as np
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, 
                           QGridLayout, QLabel, QSlider, QPushButton)
from PyQt5.QtCore import Qt, pyqtSignal
from astropy.coordinates import SkyCoord
import astropy.units as u

from .matplotlib_widget import MatplotlibWidget


class EllipseFittingWidget(QWidget):
    """Widget for interactive ellipse fitting with sliders"""
    
    fitting_done = pyqtSignal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_map = None
        self.ellipse_params = {}
        self.ellipse_line = None
        self.setup_ui()
        
    def setup_ui(self):
        """Setup the user interface"""
        layout = QVBoxLayout()
        
        # Plot area
        self.plot_widget = MatplotlibWidget()
        layout.addWidget(self.plot_widget)
        
        # Controls
        controls_group = QGroupBox("Ellipse Parameters")
        controls_layout = QGridLayout()
        
        # X Center slider
        controls_layout.addWidget(QLabel("X Center (arcsec):"), 0, 0)
        self.x_slider = QSlider(Qt.Horizontal)
        self.x_slider.setRange(-2000, 2000)
        self.x_slider.setValue(0)
        self.x_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.x_slider, 0, 1)
        self.x_label = QLabel("0")
        self.x_label.setMinimumWidth(60)
        controls_layout.addWidget(self.x_label, 0, 2)
        
        # Y Center slider
        controls_layout.addWidget(QLabel("Y Center (arcsec):"), 1, 0)
        self.y_slider = QSlider(Qt.Horizontal)
        self.y_slider.setRange(-2000, 2000)
        self.y_slider.setValue(0)
        self.y_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.y_slider, 1, 1)
        self.y_label = QLabel("0")
        self.y_label.setMinimumWidth(60)
        controls_layout.addWidget(self.y_label, 1, 2)
        
        # Semi-major axis slider
        controls_layout.addWidget(QLabel("Semi-major (arcsec):"), 2, 0)
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setRange(50, 1500)
        self.a_slider.setValue(300)
        self.a_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.a_slider, 2, 1)
        self.a_label = QLabel("300")
        self.a_label.setMinimumWidth(60)
        controls_layout.addWidget(self.a_label, 2, 2)
        
        # Semi-minor axis slider
        controls_layout.addWidget(QLabel("Semi-minor (arcsec):"), 3, 0)
        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setRange(50, 1500)
        self.b_slider.setValue(200)
        self.b_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.b_slider, 3, 1)
        self.b_label = QLabel("200")
        self.b_label.setMinimumWidth(60)
        controls_layout.addWidget(self.b_label, 3, 2)
        
        # Angle slider
        controls_layout.addWidget(QLabel("Angle (degrees):"), 4, 0)
        self.angle_slider = QSlider(Qt.Horizontal)
        self.angle_slider.setRange(0, 180)
        self.angle_slider.setValue(0)
        self.angle_slider.valueChanged.connect(self.update_ellipse)
        controls_layout.addWidget(self.angle_slider, 4, 1)
        self.angle_label = QLabel("0")
        self.angle_label.setMinimumWidth(60)
        controls_layout.addWidget(self.angle_label, 4, 2)
        
        controls_group.setLayout(controls_layout)
        layout.addWidget(controls_group)
        
        # Action buttons
        button_layout = QHBoxLayout()
        
        self.reset_button = QPushButton("Reset")
        self.reset_button.clicked.connect(self.reset_parameters)
        button_layout.addWidget(self.reset_button)
        
        self.auto_fit_button = QPushButton("Auto Fit")
        self.auto_fit_button.clicked.connect(self.auto_fit_ellipse)
        button_layout.addWidget(self.auto_fit_button)
        
        self.done_button = QPushButton("Accept Fit")
        self.done_button.clicked.connect(self.finish_fitting)
        button_layout.addWidget(self.done_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)
        
    def set_map(self, sunpy_map):
        """Set the current SunPy map for fitting"""
        self.current_map = sunpy_map
        self.plot_map()
        self.reset_parameters()
        
    def plot_map(self):
        """Plot the current map with ellipse overlay"""
        if self.current_map is None:
            return
            
        self.plot_widget.figure.clear()
        ax = self.plot_widget.figure.add_subplot(111, projection=self.current_map)
        
        # Plot the map
        self.current_map.plot(axes=ax)
        ax.set_title(f'Running Difference Map - {self.current_map.date.strftime("%Y-%m-%d %H:%M:%S")}')
        
        # Set appropriate limits
        ax.set_xlim(-1200, 1200)
        ax.set_ylim(-1200, 1200)
        
        self.plot_widget.canvas.draw()
        self.update_ellipse()
        
    def update_ellipse(self):
        """Update ellipse visualization based on slider values"""
        if self.current_map is None:
            return
            
        x0 = self.x_slider.value()
        y0 = self.y_slider.value()
        a = self.a_slider.value()
        b = self.b_slider.value()
        angle = self.angle_slider.value()
        
        # Update labels
        self.x_label.setText(str(x0))
        self.y_label.setText(str(y0))
        self.a_label.setText(str(a))
        self.b_label.setText(str(b))
        self.angle_label.setText(str(angle))
        
        # Get current axes
        ax = self.plot_widget.figure.axes[0] if self.plot_widget.figure.axes else None
        if ax is None:
            return
            
        # Remove previous ellipse
        if self.ellipse_line is not None:
            self.ellipse_line.remove()
            self.ellipse_line = None
        
        # Create ellipse points
        theta = np.linspace(0, 2*np.pi, 300)
        angle_rad = np.deg2rad(angle)
        
        # Parametric ellipse equations with rotation
        x_ell = x0 + a * np.cos(theta) * np.cos(angle_rad) - b * np.sin(theta) * np.sin(angle_rad)
        y_ell = y0 + a * np.cos(theta) * np.sin(angle_rad) + b * np.sin(theta) * np.cos(angle_rad)
        
        try:
            # Create coordinates
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                            frame=self.current_map.coordinate_frame)
            
            # Plot ellipse
            self.ellipse_line = ax.plot_coord(coords, color='red', linewidth=2, alpha=0.8)[0]
            
        except Exception as e:
            print(f"Error plotting ellipse: {e}")
            # Fallback to pixel coordinates
            try:
                x_pix, y_pix = self.current_map.world_to_pixel(
                    SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                           frame=self.current_map.coordinate_frame))
                self.ellipse_line = ax.plot(x_pix.value, y_pix.value, 'r-', linewidth=2, alpha=0.8)[0]
            except:
                pass
        
        self.plot_widget.canvas.draw()
        
    def reset_parameters(self):
        """Reset ellipse parameters to default values"""
        self.x_slider.setValue(0)
        self.y_slider.setValue(0)
        self.a_slider.setValue(300)
        self.b_slider.setValue(200)
        self.angle_slider.setValue(0)
        
    def auto_fit_ellipse(self):
        """Attempt automatic ellipse fitting (simplified version)"""
        if self.current_map is None:
            return
            
        try:
            # Get map data
            data = self.current_map.data
            
            # Find the center of mass of positive values
            positive_mask = data > np.percentile(data, 85)
            if np.sum(positive_mask) > 0:
                y_indices, x_indices = np.where(positive_mask)
                
                # Convert to world coordinates
                x_world = []
                y_world = []
                
                for i in range(0, len(x_indices), max(1, len(x_indices)//100)):  # Sample points
                    try:
                        coord = self.current_map.pixel_to_world(x_indices[i], y_indices[i])
                        x_world.append(coord.Tx.value)
                        y_world.append(coord.Ty.value)
                    except:
                        continue
                
                if x_world and y_world:
                    # Calculate center of mass
                    x_center = np.mean(x_world)
                    y_center = np.mean(y_world)
                    
                    # Estimate semi-axes from standard deviation
                    a_est = max(100, min(800, 2 * np.std(x_world)))
                    b_est = max(100, min(800, 2 * np.std(y_world)))
                    
                    # Update sliders
                    self.x_slider.setValue(int(x_center))
                    self.y_slider.setValue(int(y_center))
                    self.a_slider.setValue(int(a_est))
                    self.b_slider.setValue(int(b_est))
                    
        except Exception as e:
            print(f"Auto-fit error: {e}")
            # Fall back to center of image
            self.x_slider.setValue(0)
            self.y_slider.setValue(0)
            
    def finish_fitting(self):
        """Store fitting parameters and emit signal"""
        self.ellipse_params = {
            'x0': self.x_slider.value(),
            'y0': self.y_slider.value(),
            'a': self.a_slider.value(),
            'b': self.b_slider.value(),
            'angle': self.angle_slider.value(),
            'timestamp': self.current_map.date if self.current_map else None
        }
        self.fitting_done.emit(self.ellipse_params)
        
    def get_current_params(self):
        """Get current ellipse parameters"""
        return {
            'x0': self.x_slider.value(),
            'y0': self.y_slider.value(),
            'a': self.a_slider.value(),
            'b': self.b_slider.value(),
            'angle': self.angle_slider.value(),
            'timestamp': self.current_map.date if self.current_map else None
        }