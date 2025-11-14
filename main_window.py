import os
from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                             QFileDialog, QLabel, QComboBox, QProgressBar, QMessageBox, QTabWidget)
from PyQt5.QtCore import Qt
from ellipse_fitting import EllipseFittingWidget
from processor import DataProcessor, RunningDifferenceProcessor
from visualization import Visualization3D
import plotly
import plotly.express as px
from PyQt5.QtWebEngineWidgets import QWebEngineView

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Solar Wave Analysis Tool")
        self.setGeometry(100, 100, 1200, 800)
        
        # Data storage
        self.fits_files = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.ellipse_params_list = []
        self.current_index = 0
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Top controls
        top_controls = QHBoxLayout()
        
        self.load_button = QPushButton("Load FITS Files")
        self.load_button.clicked.connect(self.load_files)
        top_controls.addWidget(self.load_button)
        
        self.channel_combo = QComboBox()
        self.channel_combo.addItems(["94", "131", "171", "193", "211", "335"])
        top_controls.addWidget(QLabel("AIA Channel:"))
        top_controls.addWidget(self.channel_combo)
        
        self.process_button = QPushButton("Process Files")
        self.process_button.clicked.connect(self.process_files)
        self.process_button.setEnabled(False)
        top_controls.addWidget(self.process_button)
        
        main_layout.addLayout(top_controls)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main_layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel("Ready")
        main_layout.addWidget(self.status_label)
        
        # Tab widget
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Fitting tab
        self.fitting_tab = QWidget()
        self.fitting_tab_layout = QVBoxLayout(self.fitting_tab)
        self.tab_widget.addTab(self.fitting_tab, "Fitting")
        
        # Visualization tab
        self.visualization_tab = QWidget()
        self.visualization_tab_layout = QVBoxLayout(self.visualization_tab)
        self.tab_widget.addTab(self.visualization_tab, "Visualization")
        
        # Ellipse fitting widget
        self.ellipse_fitting_widget = EllipseFittingWidget()
        self.fitting_tab_layout.addWidget(self.ellipse_fitting_widget)
        
        # Navigation buttons
        nav_layout = QHBoxLayout()
        self.prev_button = QPushButton("Previous")
        self.prev_button.clicked.connect(self.prev_image)
        self.prev_button.setEnabled(False)
        nav_layout.addWidget(self.prev_button)
        
        self.next_button = QPushButton("Next")
        self.next_button.clicked.connect(self.next_image)
        self.next_button.setEnabled(False)
        nav_layout.addWidget(self.next_button)
        
        self.visualize_button = QPushButton("Generate 3D Visualization")
        self.visualize_button.clicked.connect(self.generate_visualization)
        self.visualize_button.setEnabled(False)
        nav_layout.addWidget(self.visualize_button)
        
        main_layout.addLayout(nav_layout)
        
        # Connect signals
        self.ellipse_fitting_widget.fitting_done.connect(self.handle_fitting_done)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select AIA FITS Files", "", "FITS Files (*.fits)"
        )
        if files:
            self.fits_files = files
            self.process_button.setEnabled(True)
            self.status_label.setText(f"Loaded {len(files)} files. Click 'Process Files' to continue.")

    def process_files(self):
        # Create and run processor thread
        self.processor = DataProcessor(
            self.fits_files, 
            self.channel_combo.currentText(),
            os.path.dirname(self.fits_files[0])
        )
        self.processor.progress.connect(self.progress_bar.setValue)
        self.processor.status.connect(self.status_label.setText)
        self.processor.finished.connect(self.handle_processing_finished)
        self.processor.error.connect(self.handle_processing_error)
        
        self.progress_bar.setVisible(True)
        self.status_label.setText("Processing files...")
        self.process_button.setEnabled(False)
        self.processor.start()

    def handle_processing_finished(self, processed_maps):
        self.progress_bar.setVisible(False)
        self.processed_maps = processed_maps
        self.status_label.setText(f"Processed {len(processed_maps)} maps.")
        
        # Create running difference maps
        try:
            self.running_diff_maps = RunningDifferenceProcessor.create_running_diff_maps(processed_maps)
            self.status_label.setText(f"Created {len(self.running_diff_maps)} running difference maps.")
            self.show_current_running_diff()
            self.next_button.setEnabled(True)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create running difference maps: {str(e)}")

    def show_current_running_diff(self):
        if self.current_index < len(self.running_diff_maps):
            self.ellipse_fitting_widget.set_map(self.running_diff_maps[self.current_index])
            self.status_label.setText(f"Fitting ellipse for frame {self.current_index+1}/{len(self.running_diff_maps)}")

    def handle_fitting_done(self, ellipse_params):
        self.ellipse_params_list.append(ellipse_params)
        self.visualize_button.setEnabled(True)
        
        if self.current_index < len(self.running_diff_maps) - 1:
            self.next_button.setEnabled(True)
        else:
            self.status_label.setText("All frames processed! Generate 3D visualization.")

    def next_image(self):
        if self.current_index < len(self.running_diff_maps) - 1:
            self.current_index += 1
            self.show_current_running_diff()
            self.prev_button.setEnabled(True)
            if self.current_index == len(self.running_diff_maps) - 1:
                self.next_button.setEnabled(False)

    def prev_image(self):
        if self.current_index > 0:
            self.current_index -= 1
            self.show_current_running_diff()
            self.next_button.setEnabled(True)
            if self.current_index == 0:
                self.prev_button.setEnabled(False)

    def generate_visualization(self):
        if not self.ellipse_params_list:
            return
            
        # Clear previous visualizations
        for i in reversed(range(self.visualization_tab_layout.count())):
            self.visualization_tab_layout.itemAt(i).widget().deleteLater()
        
        # Create 3D visualization for first frame
        fig_3d = Visualization3D.create_3d_ellipsoid(
            self.ellipse_params_list[0], 
            self.running_diff_maps[0]
        )
        html_3d = plotly.offline.plot(fig_3d, output_type='div', include_plotlyjs='cdn')
        web_view_3d = QWebEngineView()
        web_view_3d.setHtml(html_3d)
        self.visualization_tab_layout.addWidget(web_view_3d)
        
        # Create theta distribution
        fig_theta = Visualization3D.create_theta_distribution(self.ellipse_params_list)
        html_theta = plotly.offline.plot(fig_theta, output_type='div', include_plotlyjs='cdn')
        web_view_theta = QWebEngineView()
        web_view_theta.setHtml(html_theta)
        self.visualization_tab_layout.addWidget(web_view_theta)
        
        # Switch to visualization tab
        self.tab_widget.setCurrentIndex(1)