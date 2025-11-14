# # import sys
# # from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QVBoxLayout

# # def main():
# #     app = QApplication(sys.argv) # Handles the app lifecycle

# #     # Main window
# #     window = QWidget() # A blank canvas window
# #     window.setWindowTitle('My First GUI')
# #     window.setGeometry(100, 100, 300, 200)

# #     # Layout and widgets
# #     layout = QVBoxLayout() # Vertical stacking of widgets

# #     label = QLabel('Hello, click the button!') # For displaying text
# #     button = QPushButton('Click Me') # Our action trigger

# #     def on_click():
# #         label.setText('Button clicked!')

# #     button.clicked.connect(on_click)

# #     layout.addWidget(label)
# #     layout.addWidget(button)
# #     window.setLayout(layout)

# #     # Show the window
# #     window.show()
# #     sys.exit(app.exec_())

# # if __name__ == '__main__':
# #     main()

# # =========================================================================================

# # import sys
# # from PyQt5.QtWidgets import (
# #     QApplication, QWidget, QPushButton, QLabel,
# #     QVBoxLayout, QHBoxLayout, QFileDialog
# # )
# # from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# # from matplotlib.figure import Figure
# # import sunpy.map
# # import matplotlib.pyplot as plt

# # class AIAViewer(QWidget):
# #     def __init__(self):
# #         super().__init__()

# #         self.setWindowTitle('AIA FITS Viewer')
# #         self.setGeometry(200, 200, 800, 600)

# #         self.files = []
# #         self.maps = []
# #         self.current_index = 0

# #         self.init_ui()

# #     def init_ui(self):
# #         # Buttons
# #         self.load_button = QPushButton('Load FITS')
# #         self.display_button = QPushButton('Display')
# #         self.prev_button = QPushButton('Previous')
# #         self.next_button = QPushButton('Next')

# #         self.load_button.clicked.connect(self.load_files)
# #         self.display_button.clicked.connect(self.display_first_map)
# #         self.prev_button.clicked.connect(self.show_prev_map)
# #         self.next_button.clicked.connect(self.show_next_map)

# #         # Layouts
# #         btn_layout = QHBoxLayout()
# #         btn_layout.addWidget(self.load_button)
# #         btn_layout.addWidget(self.display_button)
# #         btn_layout.addWidget(self.prev_button)
# #         btn_layout.addWidget(self.next_button)

# #         # Matplotlib figure and canvas
# #         self.figure = Figure()
# #         self.canvas = FigureCanvas(self.figure)

# #         # Info label
# #         self.label = QLabel('No file loaded.')

# #         # Main layout
# #         layout = QVBoxLayout()
# #         layout.addLayout(btn_layout)
# #         layout.addWidget(self.canvas)
# #         layout.addWidget(self.label)

# #         self.setLayout(layout)

# #     def load_files(self):
# #         files, _ = QFileDialog.getOpenFileNames(
# #             self, 'Select FITS files', '', 'FITS files (*.fits *.fts)'
# #         )
# #         if files:
# #             self.files = files
# #             self.maps = [sunpy.map.Map(f) for f in files]
# #             self.current_index = 0
# #             self.label.setText(f'{len(files)} files loaded.')

# #     def display_first_map(self):
# #         if not self.maps:
# #             self.label.setText('No maps loaded.')
# #             return
# #         self.current_index = 0
# #         self.plot_map(self.maps[0])

# #     def show_prev_map(self):
# #         if not self.maps:
# #             return
# #         self.current_index = max(0, self.current_index - 1)
# #         self.plot_map(self.maps[self.current_index])

# #     def show_next_map(self):
# #         if not self.maps:
# #             return
# #         self.current_index = min(len(self.maps) - 1, self.current_index + 1)
# #         self.plot_map(self.maps[self.current_index])

# #     def plot_map(self, amap):
# #         self.figure.clf()  # Clear previous figure
# #         ax = self.figure.add_subplot(111, projection=amap.wcs)
# #         amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)
# #         ax.grid(False)
# #         self.canvas.draw()
# #         self.label.setText(f'Showing: {self.files[self.current_index]}')

# # def main():
# #     app = QApplication(sys.argv)
# #     viewer = AIAViewer()
# #     viewer.show()
# #     sys.exit(app.exec_())

# # if __name__ == '__main__':
# #     import astropy.units as u
# #     main()

# # =========================================================================================

# # import sys
# # import os
# # from PyQt5.QtWidgets import (
# #     QApplication, QWidget, QPushButton, QLabel, QFileDialog,
# #     QVBoxLayout, QHBoxLayout, QMessageBox
# # )
# # from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# # from matplotlib.figure import Figure
# # from matplotlib import colors
# # import sunpy.map
# # import numpy as np
# # from aiapy.calibrate import register, update_pointing
# # from scipy import ndimage
# # import astropy.units as u

# # class AIAViewer(QWidget):
# #     def __init__(self):
# #         super().__init__()
# #         self.setWindowTitle('AIA FITS Viewer')
# #         self.setGeometry(200, 200, 1000, 700)

# #         self.files = []
# #         self.maps = []
# #         self.processed_maps = []
# #         self.running_diff_maps = []
# #         self.current_index = 0
# #         self.data_dir = os.getcwd()
# #         self.channel = '171'

# #         self.init_ui()

# #     def init_ui(self):
# #         # Buttons
# #         self.load_button = QPushButton('Load FITS')
# #         self.display_button = QPushButton('Display')
# #         self.upgrade_button = QPushButton('Upgrade')
# #         self.run_diff_button = QPushButton('Run-Diff')
# #         self.prev_button = QPushButton('Previous')
# #         self.next_button = QPushButton('Next')
# #         self.save_button = QPushButton('Save')

# #         self.load_button.clicked.connect(self.load_files)
# #         self.display_button.clicked.connect(self.display_first_map)
# #         self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
# #         self.run_diff_button.clicked.connect(self.create_running_diff_maps)
# #         self.prev_button.clicked.connect(self.show_prev_map)
# #         self.next_button.clicked.connect(self.show_next_map)
# #         self.save_button.clicked.connect(self.save_current_map)

# #         # Layouts
# #         btn_layout = QHBoxLayout()
# #         for b in [self.load_button, self.display_button, self.upgrade_button,
# #                   self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
# #             btn_layout.addWidget(b)

# #         self.figure = Figure()
# #         self.canvas = FigureCanvas(self.figure)
# #         self.label = QLabel('Ready.')

# #         layout = QVBoxLayout()
# #         layout.addLayout(btn_layout)
# #         layout.addWidget(self.canvas)
# #         layout.addWidget(self.label)
# #         self.setLayout(layout)
    
# #     def detect_data_level(self, file):
# #         name = os.path.basename(file).lower()
# #         if 'lev1.5' in name or 'lev15' in name:
# #             return 'lev15'
# #         elif 'lev1' in name:
# #             return 'lev1'
# #         else:
# #             try:
# #                 hdr = sunpy.map.Map(file).meta
# #                 level = hdr.get('lvl_num', None)
# #                 if level == 1.5:
# #                     return 'lev15'
# #                 elif level == 1:
# #                     return 'lev1'
# #             except Exception as e:
# #                 print(f'Could not read header from {file}: {e}')
# #             return 'unknown'
    
# #     def load_files(self):
# #         files, _ = QFileDialog.getOpenFileNames(
# #             self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
# #         if files:
# #             self.files = files
# #             self.maps = [sunpy.map.Map(f) for f in files]
# #             self.current_index = 0
# #             self.label.setText(f'{len(files)} files loaded.')

# #     def display_first_map(self):
# #         if not self.maps:
# #             self.label.setText('No maps loaded.')
# #             return
# #         self.current_index = 0
# #         self.plot_map(self.maps[0])

# #     def show_prev_map(self):
# #         if self.maps:
# #             self.current_index = max(0, self.current_index - 1)
# #             self.plot_map(self.maps[self.current_index])

# #     def show_next_map(self):
# #         if self.maps:
# #             self.current_index = min(len(self.maps) - 1, self.current_index + 1)
# #             self.plot_map(self.maps[self.current_index])

# #     def plot_map(self, amap):
# #         self.canvas.figure.clf()
# #         ax = self.canvas.figure.add_subplot(111, projection=amap)

# #         # Only add clip_interval if norm is NOT already defined
# #         if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
# #             amap.plot(axes=ax)  # Don't pass clip_interval
# #         else:
# #             amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

# #         # ax.set_title(f'{amap.wavelength.value} Å  {amap.date.strftime("%Y-%m-%d %H:%M:%S")}')
# #         ax.grid(False)
# #         self.canvas.draw()
# #         # self.label.setText(f'Showing: {os.path.basename(amap.meta.get("filename", "map"))}')
# #         self.label.setText(f'Showing: {self.files[self.current_index]}')

# #     def upgrade_to_lv15(self):
# #         if not self.files:
# #             QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
# #             return

# #         self.processed_maps = []

# #         for file in self.files:
# #             data_level = self.detect_data_level(file)
# #             try:
# #                 if data_level == 'lev1':
# #                     output_filename = os.path.basename(file).replace('lev1', 'lev15')
# #                     file_path = os.path.join(self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
# #                     os.makedirs(file_path, exist_ok=True)
# #                     full_path = os.path.join(file_path, output_filename)

# #                     if not os.path.exists(full_path):
# #                         m = sunpy.map.Map(file)
# #                         # aiapy now requires `pointing_table` kwarg, so just skip this for now
# #                         m = update_pointing(m, pointing_table=None)
# #                         m = register(m)
# #                         m = m / m.exposure_time
# #                         m.save(full_path, filetype='auto')
# #                         self.processed_maps.append(m)
# #                     else:
# #                         self.processed_maps.append(sunpy.map.Map(full_path))

# #                 elif data_level == 'lev15':
# #                     self.processed_maps.append(sunpy.map.Map(file))

# #                 else:
# #                     print(f'Skipped unknown-level file: {file}')

# #             except Exception as e:
# #                 print(f'Error upgrading {file}: {e}')

# #         self.maps = self.processed_maps
# #         self.current_index = 0
# #         self.plot_map(self.maps[0])
# #         self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
    
# #         if self.processed_maps:
# #             self.maps = self.processed_maps
# #             self.current_index = 0
# #             self.plot_map(self.maps[0])
# #             self.label.setText(f'Loaded {len(self.maps)} maps.')
# #         else:
# #             self.label.setText('No maps were processed or loaded.')

# #     def set_processed_maps_from_loaded(self):
# #         if not self.processed_maps and self.maps:
# #             self.processed_maps = self.maps.copy()

# #     def create_running_diff_maps(self):
# #         self.set_processed_maps_from_loaded()

# #         if len(self.processed_maps) < 6:
# #             QMessageBox.warning(self, 'Warning',
# #                                 'Need at least 6 images for running difference.')
# #             return

# #         self.running_diff_maps = []

# #         for i in range(5, len(self.processed_maps)):
# #             try:
# #                 m0 = self.processed_maps[i-5]
# #                 m1 = self.processed_maps[i]
# #                 diff = m1.quantity - m0.quantity
# #                 smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
# #                 diff_map = sunpy.map.Map(smoothed, m1.meta)
# #                 diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
# #                 self.running_diff_maps.append(diff_map)
# #             except Exception as e:
# #                 print(f'Error creating diff map {i}: {e}')

# #         self.maps = self.running_diff_maps
# #         self.current_index = 0
# #         if self.maps:
# #             self.plot_map(self.maps[0])
# #         self.label.setText(f'Created {len(self.running_diff_maps)} running difference maps.')

# #     def save_current_map(self):
# #         if not self.maps:
# #             QMessageBox.warning(self, 'Warning', 'No map to save.')
# #             return

# #         options = QFileDialog.Options()
# #         file_path, _ = QFileDialog.getSaveFileName(
# #             self, 'Save current map as image',
# #             f'{self.current_index:03d}.png',
# #             'PNG (*.png);;PDF (*.pdf)', options=options)

# #         if file_path:
# #             self.figure.savefig(file_path)
# #             self.label.setText(f'Saved map to {file_path}')

# # def main():
# #     app = QApplication(sys.argv)
# #     viewer = AIAViewer()
# #     viewer.show()
# #     sys.exit(app.exec_())

# # if __name__ == '__main__':
# #     main()

# # =========================================================================================

# # import os
# # import sys
# # from PyQt5.QtWidgets import (
# #     QApplication, QWidget, QPushButton, QLabel, QFileDialog,
# #     QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider
# # )
# # from PyQt5.QtCore import Qt
# # from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# # from matplotlib.figure import Figure
# # from matplotlib import colors
# # import sunpy.map
# # import numpy as np
# # from aiapy.calibrate import register, update_pointing
# # from scipy import ndimage
# # import astropy.units as u
# # from astropy.coordinates import SkyCoord

# # class AIAViewer(QWidget):
# #     def __init__(self):
# #         super().__init__()
# #         self.setWindowTitle('AIA FITS Viewer')
# #         self.figure = Figure(figsize=[10,10])
# #         # self.setGeometry(100, 100, 1200, 1000)

# #         self.files = []
# #         self.maps = []
# #         self.processed_maps = []
# #         self.running_diff_maps = []
# #         self.current_index = 0
# #         self.data_dir = os.getcwd()
# #         self.channel = '193'
# #         self.ellipse_artist = None
        
# #         self.a_slider = QSlider(Qt.Horizontal)
# #         self.a_slider.setMinimum(10)
# #         self.a_slider.setMaximum(1000)
# #         self.a_slider.setValue(200)

# #         self.b_slider = QSlider(Qt.Horizontal)
# #         self.b_slider.setMinimum(10)
# #         self.b_slider.setMaximum(1000)
# #         self.b_slider.setValue(100)

# #         self.a_slider.valueChanged.connect(self.update_ellipse)
# #         self.b_slider.valueChanged.connect(self.update_ellipse)

# #         self.init_ui()
    

# #     def init_ui(self):
# #         # Buttons
# #         self.load_button = QPushButton('Load FITS')
# #         self.display_button = QPushButton('Display')
# #         self.upgrade_button = QPushButton('Upgrade')
# #         self.run_diff_button = QPushButton('Run-Diff')
# #         self.prev_button = QPushButton('Previous')
# #         self.next_button = QPushButton('Next')
# #         self.save_button = QPushButton('Save')

# #         self.load_button.clicked.connect(self.load_files)
# #         self.display_button.clicked.connect(self.display_first_map)
# #         self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
# #         self.run_diff_button.clicked.connect(self.create_running_diff_maps)
# #         self.prev_button.clicked.connect(self.show_prev_map)
# #         self.next_button.clicked.connect(self.show_next_map)
# #         self.save_button.clicked.connect(self.save_current_map)

# #         self.ellipse_button = QPushButton('Ellipse')
# #         self.fit_button = QPushButton('Fit')
        
# #         # Inputs for ellipse center
# #         self.x_input = QLineEdit()
# #         self.x_input.setPlaceholderText('X center [arcsec]')
# #         self.y_input = QLineEdit()
# #         self.y_input.setPlaceholderText('Y center [arcsec]')

# #         # Sliders for axes
# #         self.a_slider = QSlider(Qt.Horizontal)
# #         self.a_slider.setMinimum(10)
# #         self.a_slider.setMaximum(1000)
# #         self.a_slider.setValue(200)
# #         self.b_slider = QSlider(Qt.Horizontal)
# #         self.b_slider.setMinimum(10)
# #         self.b_slider.setMaximum(1000)
# #         self.b_slider.setValue(100)

# #         # Connect the buttons
# #         self.ellipse_button.clicked.connect(self.draw_ellipse)
# #         self.fit_button.clicked.connect(self.extract_ellipse_params)

# #         # Layouts
# #         btn_layout = QHBoxLayout()
# #         for b in [self.load_button, self.display_button, self.upgrade_button,
# #                   self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
# #             btn_layout.addWidget(b)

# #         self.figure = Figure()
# #         self.canvas = FigureCanvas(self.figure)
# #         self.label = QLabel('Ready.')

# #         layout = QVBoxLayout()
# #         layout.addLayout(btn_layout)
# #         layout.addWidget(self.canvas)
# #         layout.addWidget(self.label)
# #         self.setLayout(layout)

# #         ellipse_layout = QVBoxLayout()
# #         ellipse_layout.addWidget(QLabel('Ellipse Center:'))
# #         ellipse_layout.addWidget(self.x_input)
# #         ellipse_layout.addWidget(self.y_input)
# #         ellipse_layout.addWidget(QLabel('Semi-Major Axis:'))
# #         ellipse_layout.addWidget(self.a_slider)
# #         ellipse_layout.addWidget(QLabel('Semi-Minor Axis:'))
# #         ellipse_layout.addWidget(self.b_slider)
# #         ellipse_layout.addWidget(self.ellipse_button)
# #         ellipse_layout.addWidget(self.fit_button)
# #         layout.addLayout(ellipse_layout)
    

# #     def detect_data_level(self, file):
# #         name = os.path.basename(file).lower()
# #         if 'lev1.5' in name or 'lev15' in name:
# #             return 'lev15'
# #         elif 'lev1' in name:
# #             return 'lev1'
# #         else:
# #             try:
# #                 hdr = sunpy.map.Map(file).meta
# #                 level = hdr.get('lvl_num', None)
# #                 if level == 1.5:
# #                     return 'lev15'
# #                 elif level == 1:
# #                     return 'lev1'
# #             except Exception as e:
# #                 print(f'Could not read header from {file}: {e}')
# #             return 'unknown'
    
# #     def load_files(self):
# #         files, _ = QFileDialog.getOpenFileNames(
# #             self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
# #         if files:
# #             self.files = files
# #             self.maps = [sunpy.map.Map(f) for f in files]
# #             self.current_index = 0
# #             self.label.setText(f'{len(files)} files loaded.')

# #     def display_first_map(self):
# #         if not self.maps:
# #             self.label.setText('No maps loaded.')
# #             return
# #         self.current_index = 0
# #         self.plot_map(self.maps[0])

# #     def show_prev_map(self):
# #         if self.maps:
# #             self.current_index = max(0, self.current_index - 1)
# #             self.plot_map(self.maps[self.current_index])

# #     def show_next_map(self):
# #         if self.maps:
# #             self.current_index = min(len(self.maps) - 1, self.current_index + 1)
# #             self.plot_map(self.maps[self.current_index])

# #     def plot_map(self, amap):
# #         self.canvas.figure.clf()
# #         ax = self.canvas.figure.add_subplot(111, projection=amap)

# #         # Only add clip_interval if norm is NOT already defined
# #         if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
# #             amap.plot(axes=ax)  # Don't pass clip_interval
# #         else:
# #             amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

# #         # ax.set_title(f'{amap.wavelength.value} Å  {amap.date.strftime("%Y-%m-%d %H:%M:%S")}')
# #         ax.grid(False)
# #         self.canvas.draw()
# #         # self.label.setText(f'Showing: {os.path.basename(amap.meta.get("filename", "map"))}')
# #         self.label.setText(f'Showing: {self.files[self.current_index]}')

# #     def upgrade_to_lv15(self):
# #         if not self.files:
# #             QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
# #             return

# #         self.processed_maps = []

# #         for file in self.files:
# #             data_level = self.detect_data_level(file)
# #             try:
# #                 if data_level == 'lev1':
# #                     output_filename = os.path.basename(file).replace('lev1', 'lev15')
# #                     file_path = os.path.join(self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
# #                     os.makedirs(file_path, exist_ok=True)
# #                     full_path = os.path.join(file_path, output_filename)

# #                     if not os.path.exists(full_path):
# #                         m = sunpy.map.Map(file)
# #                         # aiapy now requires `pointing_table` kwarg, so just skip this for now
# #                         m = update_pointing(m, pointing_table=None)
# #                         m = register(m)
# #                         m = m / m.exposure_time
# #                         m.save(full_path, filetype='auto')
# #                         self.processed_maps.append(m)
# #                     else:
# #                         self.processed_maps.append(sunpy.map.Map(full_path))

# #                 elif data_level == 'lev15':
# #                     self.processed_maps.append(sunpy.map.Map(file))

# #                 else:
# #                     print(f'Skipped unknown-level file: {file}')

# #             except Exception as e:
# #                 print(f'Error upgrading {file}: {e}')

# #         self.maps = self.processed_maps
# #         self.current_index = 0
# #         self.plot_map(self.maps[0])
# #         self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
    
# #         if self.processed_maps:
# #             self.maps = self.processed_maps
# #             self.current_index = 0
# #             self.plot_map(self.maps[0])
# #             self.label.setText(f'Loaded {len(self.maps)} maps.')
# #         else:
# #             self.label.setText('No maps were processed or loaded.')

# #     def set_processed_maps_from_loaded(self):
# #         if not self.processed_maps and self.maps:
# #             self.processed_maps = self.maps.copy()

# #     def create_running_diff_maps(self):
# #         self.set_processed_maps_from_loaded()

# #         if len(self.processed_maps) < 6:
# #             QMessageBox.warning(self, 'Warning',
# #                                 'Need at least 6 images for running difference.')
# #             return

# #         self.running_diff_maps = []

# #         for i in range(5, len(self.processed_maps)):
# #             try:
# #                 m0 = self.processed_maps[i-5]
# #                 m1 = self.processed_maps[i]
# #                 diff = m1.quantity - m0.quantity
# #                 smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
# #                 diff_map = sunpy.map.Map(smoothed, m1.meta)
# #                 diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
# #                 self.running_diff_maps.append(diff_map)
# #             except Exception as e:
# #                 print(f'Error creating diff map {i}: {e}')

# #         self.maps = self.running_diff_maps
# #         self.current_index = 0
# #         if self.maps:
# #             self.plot_map(self.maps[0])
# #         self.label.setText(f'Created {len(self.running_diff_maps)} running difference maps.')


# #     def draw_ellipse(self):
# #         if not self.maps:
# #             self.label.setText('No maps loaded.')
# #             return

# #         try:
# #             x0 = float(self.x_input.text())
# #             y0 = float(self.y_input.text())
# #             a = self.a_slider.value()
# #             b = self.b_slider.value()
# #         except ValueError:
# #             self.label.setText('Invalid center coordinates.')
# #             return

# #         amap = self.maps[self.current_index]
# #         self.canvas.figure.clf()
# #         ax = self.canvas.figure.add_subplot(111, projection=amap)

# #         if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
# #             amap.plot(axes=ax)
# #         else:
# #             amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

# #         # Ellipse
# #         theta = np.linspace(0, 2*np.pi, 300)
# #         x_ell = x0 + a * np.cos(theta)
# #         y_ell = y0 + b * np.sin(theta)

# #         coords = SkyCoord(x_ell * u.arcsec,
# #                           y_ell * u.arcsec,
# #                           frame=amap.coordinate_frame)

# #         ax.plot_coord(coords, color='red', lw=2)
# #         ax.grid(False)
# #         self.canvas.draw()
# #         self.label.setText('Ellipse drawn.')


# #     def update_ellipse(self):
# #         if not self.maps:
# #             return

# #         amap = self.maps[self.current_index]
# #         x0 = self.x_input.text()
# #         y0 = self.y_input.text()
# #         try:
# #             x0 = float(x0)
# #             y0 = float(y0)
# #         except ValueError:
# #             return  # ignore until both inputs are valid floats

# #         a = self.a_slider.value()
# #         b = self.b_slider.value()

# #         x_ell = x0 + a * np.cos(self.theta)
# #         y_ell = y0 + b * np.sin(self.theta)

# #         coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, frame=amap.coordinate_frame)

# #         ax = self.canvas.figure.axes[0]
# #         if self.ellipse_artist:
# #             self.ellipse_artist.remove()
# #         self.ellipse_artist = ax.plot_coord(coords, color='red', lw=2)[0]
# #         self.canvas.draw()


# #     def extract_ellipse_params(self):
# #         try:
# #             x0 = float(self.x_input.text())
# #             y0 = float(self.y_input.text())
# #             a = self.a_slider.value()
# #             b = self.b_slider.value()
# #             print(f'Fit extracted: center=({x0}, {y0}) arcsec, a={a}, b={b} arcsec')
# #             self.label.setText(f'Fit extracted: a={a}, b={b}')
# #             # Optionally: save to CSV or log file here
# #         except ValueError:
# #             self.label.setText('Invalid ellipse parameters.')
    

# #     def save_current_map(self):
# #         if not self.maps:
# #             QMessageBox.warning(self, 'Warning', 'No map to save.')
# #             return

# #         options = QFileDialog.Options()
# #         file_path, _ = QFileDialog.getSaveFileName(
# #             self, 'Save current map as image',
# #             f'{self.current_index:03d}.png',
# #             'PNG (*.png);;PDF (*.pdf)', options=options)

# #         if file_path:
# #             self.figure.tight_layout()
# #             self.figure.savefig(file_path, bbox_inches='tight', pad_inches=0.05, dpi=300)
# #             self.label.setText(f'Saved map to {file_path}')

# # def main():
# #     app = QApplication(sys.argv)
# #     viewer = AIAViewer()
# #     viewer.show()
# #     sys.exit(app.exec_())

# # if __name__ == '__main__':
# #     main()

# # =========================================================================================

# import os
# import sys
# from PyQt5.QtWidgets import (
#     QApplication, QWidget, QPushButton, QLabel, QFileDialog,
#     QVBoxLayout, QHBoxLayout, QMessageBox, QLineEdit, QSlider
# )
# from PyQt5.QtCore import Qt
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.figure import Figure
# from matplotlib import colors
# import sunpy.map
# import numpy as np
# from aiapy.calibrate import register, update_pointing
# from scipy import ndimage
# import astropy.units as u
# from astropy.coordinates import SkyCoord

# class AIAViewer(QWidget):
#     def __init__(self):
#         # init stuff
#         super().__init__()
#         self.setWindowTitle('AIA FITS Viewer')
#         # self.figure = Figure(figsize=[10,10])
#         self.setGeometry(100, 100, 1200, 1000)

#         self.files = []
#         self.maps = []
#         self.processed_maps = []
#         self.running_diff_maps = []
#         self.current_index = 0
#         self.data_dir = os.getcwd()
#         self.channel = '193'
#         self.ellipse_artist = None
        
#         self.a_slider = QSlider(Qt.Horizontal)
#         self.a_slider.setMinimum(10)
#         self.a_slider.setMaximum(1000)
#         self.a_slider.setValue(200)

#         self.b_slider = QSlider(Qt.Horizontal)
#         self.b_slider.setMinimum(10)
#         self.b_slider.setMaximum(1000)
#         self.b_slider.setValue(100)

#         self.theta = np.linspace(0, 2*np.pi, 300)
#         self.ellipse_artist = None

#         self.init_ui()
    

#     def init_ui(self):
#         # create widgets
#         # connect signals

#         # Buttons
#         self.load_button = QPushButton('Load FITS')
#         self.display_button = QPushButton('Display')
#         self.upgrade_button = QPushButton('Upgrade')
#         self.run_diff_button = QPushButton('Run-Diff')
#         self.prev_button = QPushButton('Previous')
#         self.next_button = QPushButton('Next')
#         self.save_button = QPushButton('Save')

#         self.load_button.clicked.connect(self.load_files)
#         self.display_button.clicked.connect(self.display_first_map)
#         self.upgrade_button.clicked.connect(self.upgrade_to_lv15)
#         self.run_diff_button.clicked.connect(self.create_running_diff_maps)
#         self.prev_button.clicked.connect(self.show_prev_map)
#         self.next_button.clicked.connect(self.show_next_map)
#         self.save_button.clicked.connect(self.save_current_map)

#         self.ellipse_button = QPushButton('Ellipse')
#         self.fit_button = QPushButton('Fit')
        
#         # Inputs for ellipse center
#         self.x_input = QLineEdit()
#         self.x_input.setPlaceholderText('X center [arcsec]')
#         self.y_input = QLineEdit()
#         self.y_input.setPlaceholderText('Y center [arcsec]')

#         # Sliders for axes
#         self.a_slider = QSlider(Qt.Horizontal)
#         self.a_slider.setMinimum(10)
#         self.a_slider.setMaximum(1000)
#         self.a_slider.setValue(200)
#         self.b_slider = QSlider(Qt.Horizontal)
#         self.b_slider.setMinimum(10)
#         self.b_slider.setMaximum(1000)
#         self.b_slider.setValue(100)

#         # Connect the buttons
#         self.ellipse_button.clicked.connect(self.draw_ellipse)
#         self.fit_button.clicked.connect(self.extract_ellipse_params)
#         self.a_slider.valueChanged.connect(self.update_ellipse)
#         self.b_slider.valueChanged.connect(self.update_ellipse)

#         # Layouts
#         btn_layout = QHBoxLayout()
#         for b in [self.load_button, self.display_button, self.upgrade_button,
#                   self.run_diff_button, self.prev_button, self.next_button, self.save_button]:
#             btn_layout.addWidget(b)

#         self.figure = Figure()
#         self.canvas = FigureCanvas(self.figure)
#         self.label = QLabel('Ready.')

#         layout = QVBoxLayout()
#         layout.addLayout(btn_layout)
#         layout.addWidget(self.canvas)
#         layout.addWidget(self.label)
#         self.setLayout(layout)

#         ellipse_layout = QVBoxLayout()
#         ellipse_layout.addWidget(QLabel('Ellipse Center:'))
#         ellipse_layout.addWidget(self.x_input)
#         ellipse_layout.addWidget(self.y_input)
#         ellipse_layout.addWidget(QLabel('Semi-Major Axis:'))
#         ellipse_layout.addWidget(self.a_slider)
#         ellipse_layout.addWidget(QLabel('Semi-Minor Axis:'))
#         ellipse_layout.addWidget(self.b_slider)
#         ellipse_layout.addWidget(self.ellipse_button)
#         ellipse_layout.addWidget(self.fit_button)
#         layout.addLayout(ellipse_layout)
    

#     def detect_data_level(self, file):
#         # method code
#         name = os.path.basename(file).lower()
#         if 'lev1.5' in name or 'lev15' in name:
#             return 'lev15'
#         elif 'lev1' in name:
#             return 'lev1'
#         else:
#             try:
#                 hdr = sunpy.map.Map(file).meta
#                 level = hdr.get('lvl_num', None)
#                 if level == 1.5:
#                     return 'lev15'
#                 elif level == 1:
#                     return 'lev1'
#             except Exception as e:
#                 print(f'Could not read header from {file}: {e}')
#             return 'unknown'
    

#     def load_files(self):
#         # method code
#         files, _ = QFileDialog.getOpenFileNames(
#             self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
#         if files:
#             self.files = files
#             self.maps = [sunpy.map.Map(f) for f in files]
#             self.current_index = 0
#             self.label.setText(f'{len(files)} files loaded.')


#     def display_first_map(self):
#         # method code
#         if not self.maps:
#             self.label.setText('No maps loaded.')
#             return
#         self.current_index = 0
#         self.plot_map(self.maps[0])


#     def show_prev_map(self):
#         # method code
#         if self.maps:
#             self.current_index = max(0, self.current_index - 1)
#             self.plot_map(self.maps[self.current_index])


#     def show_next_map(self):
#         # method code
#         if self.maps:
#             self.current_index = min(len(self.maps) - 1, self.current_index + 1)
#             self.plot_map(self.maps[self.current_index])


#     def plot_map(self, amap):
#         # method code
#         self.canvas.figure.clf()
#         amap = self.maps[self.current_index]
#         ax = self.canvas.figure.add_subplot(111, projection=amap)

#         # Only add clip_interval if norm is NOT already defined
#         if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
#             amap.plot(axes=ax)
#         else:
#             amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

#         ax.grid(False)
#         self.canvas.draw()
#         self.label.setText(f'Showing: {self.files[self.current_index]}')


#     def plot_map_with_ellipse(self):
#         """Plot current map and ellipse if exists."""
#         self.canvas.figure.clf()
#         amap = self.maps[self.current_index]
#         ax = self.canvas.figure.add_subplot(111, projection=amap)

#         if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
#             amap.plot(axes=ax)
#         else:
#             amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

#         ax.grid(False)

#         if hasattr(self, 'ellipse_center'):
#             x0, y0 = self.ellipse_center
#             a = self.a_slider.value()
#             b = self.b_slider.value()
#             theta = np.linspace(0, 2*np.pi, 300)
#             x_ell = x0 + a * np.cos(theta)
#             y_ell = y0 + b * np.sin(theta)
#             coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, frame=amap.coordinate_frame)
#             # plot_coord returns a list of Line2D, grab first
#             self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)

#         self.canvas.draw()


#     def upgrade_to_lv15(self):
#         # method code
#         if not self.files:
#             QMessageBox.warning(self, 'Warning', 'No files to upgrade.')
#             return

#         self.processed_maps = []

#         for file in self.files:
#             data_level = self.detect_data_level(file)
#             try:
#                 if data_level == 'lev1':
#                     output_filename = os.path.basename(file).replace('lev1', 'lev15')
#                     file_path = os.path.join(self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
#                     os.makedirs(file_path, exist_ok=True)
#                     full_path = os.path.join(file_path, output_filename)

#                     if not os.path.exists(full_path):
#                         m = sunpy.map.Map(file)
#                         # aiapy now requires `pointing_table` kwarg, so just skip this for now
#                         m = update_pointing(m, pointing_table=None)
#                         m = register(m)
#                         m = m / m.exposure_time
#                         m.save(full_path, filetype='auto')
#                         self.processed_maps.append(m)
#                     else:
#                         self.processed_maps.append(sunpy.map.Map(full_path))

#                 elif data_level == 'lev15':
#                     self.processed_maps.append(sunpy.map.Map(file))

#                 else:
#                     print(f'Skipped unknown-level file: {file}')

#             except Exception as e:
#                 print(f'Error upgrading {file}: {e}')

#         self.maps = self.processed_maps
#         self.current_index = 0
#         self.plot_map(self.maps[0])
#         self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
    
#         if self.processed_maps:
#             self.maps = self.processed_maps
#             self.current_index = 0
#             self.plot_map(self.maps[0])
#             self.label.setText(f'Loaded {len(self.maps)} maps.')
#         else:
#             self.label.setText('No maps were processed or loaded.')


#     def set_processed_maps_from_loaded(self):
#         # method code
#         if not self.processed_maps and self.maps:
#             self.processed_maps = self.maps.copy()


#     def create_running_diff_maps(self):
#         # method code
#         self.set_processed_maps_from_loaded()

#         if len(self.processed_maps) < 6:
#             QMessageBox.warning(self, 'Warning',
#                                 'Need at least 6 images for running difference.')
#             return

#         self.running_diff_maps = []

#         for i in range(5, len(self.processed_maps)):
#             try:
#                 m0 = self.processed_maps[i-5]
#                 m1 = self.processed_maps[i]
#                 diff = m1.quantity - m0.quantity
#                 smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
#                 diff_map = sunpy.map.Map(smoothed, m1.meta)
#                 diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
#                 self.running_diff_maps.append(diff_map)
#             except Exception as e:
#                 print(f'Error creating diff map {i}: {e}')

#         self.maps = self.running_diff_maps
#         self.current_index = 0
#         if self.maps:
#             self.plot_map(self.maps[0])
#         self.label.setText(f'Created {len(self.running_diff_maps)} running difference maps.')


#     def draw_ellipse(self):
#         # Called when user clicks "Ellipse"
#         try:
#             x0 = float(self.x_input.text())
#             y0 = float(self.y_input.text())
#         except ValueError:
#             self.label.setText('Invalid center coordinates.')
#             return

#         self.ellipse_center = (x0, y0)

#         # redraw map and ellipse
#         self.plot_map_with_ellipse()
#         self.label.setText('Ellipse drawn.')


#     def update_ellipse(self):
#         # Only update ellipse shape (axes) when sliders move
#         if not hasattr(self, 'ellipse_center'):
#             return  # no ellipse to update

#         if not hasattr(self, 'ellipse_artist'):
#             # no existing ellipse_artist: draw the full plot again
#             self.plot_map_with_ellipse()
#             return

#         x0, y0 = self.ellipse_center
#         a = self.a_slider.value()
#         b = self.b_slider.value()
#         theta = np.linspace(0, 2*np.pi, 300)
#         x_ell = x0 + a * np.cos(theta)
#         y_ell = y0 + b * np.sin(theta)
#         amap = self.maps[self.current_index]
#         coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, frame=amap.coordinate_frame)

#         # Update ellipse_artist data with new coords
#         self.ellipse_artist.set_data(coords.Tx.value, coords.Ty.value)
#         self.canvas.draw()


#     def extract_ellipse_params(self):
#         try:
#             x0 = float(self.x_input.text())
#             y0 = float(self.y_input.text())
#             a = self.a_slider.value()
#             b = self.b_slider.value()
#             print(f'Fit extracted: center=({x0}, {y0}) arcsec, a={a}, b={b} arcsec')
#             self.label.setText(f'Fit extracted: a={a}, b={b}')
#         except ValueError:
#             self.label.setText('Invalid ellipse parameters.')


#     def save_current_map(self):
#         # method code
#         if not self.maps:
#             QMessageBox.warning(self, 'Warning', 'No map to save.')
#             return

#         options = QFileDialog.Options()
#         file_path, _ = QFileDialog.getSaveFileName(
#             self, 'Save current map as image',
#             f'{self.current_index:03d}.png',
#             'PNG (*.png);;PDF (*.pdf)', options=options)

#         if file_path:
#             self.figure.tight_layout()
#             self.figure.savefig(file_path, bbox_inches='tight', pad_inches=0.05, dpi=300)
#             self.label.setText(f'Saved map to {file_path}')

# def main():
#     app = QApplication(sys.argv)
#     viewer = AIAViewer()
#     viewer.show()
#     sys.exit(app.exec_())

# if __name__ == '__main__':
#     main()

# =========================================================================================

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
        # init stuff
        super().__init__()
        self.setWindowTitle('AIA FITS Viewer')
        # self.figure = Figure(figsize=[10,10])
        self.setGeometry(100, 100, 1200, 1000)

        self.ellipse_center = None  # Add this to store ellipse center
        self.ellipse_artist = None

        self.files = []
        self.maps = []
        self.processed_maps = []
        self.running_diff_maps = []
        self.current_index = 0
        self.data_dir = os.getcwd()
        self.channel = '193'
        self.ellipse_artist = None
        
        self.a_slider = QSlider(Qt.Horizontal)
        self.a_slider.setMinimum(10)
        self.a_slider.setMaximum(1000)
        self.a_slider.setValue(200)

        self.b_slider = QSlider(Qt.Horizontal)
        self.b_slider.setMinimum(10)
        self.b_slider.setMaximum(1000)
        self.b_slider.setValue(100)

        self.theta = np.linspace(0, 2*np.pi, 300)
        self.ellipse_artist = None

        self.init_ui()
    

    def init_ui(self):
        # create widgets
        # connect signals

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
        self.ellipse_button.clicked.connect(self.draw_ellipse)
        self.fit_button.clicked.connect(self.extract_ellipse_params)

        # self.a_slider.valueChanged.connect(self.update_ellipse)
        # self.b_slider.valueChanged.connect(self.update_ellipse)

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
        ellipse_layout.addWidget(self.y_input)
        ellipse_layout.addWidget(QLabel('Semi-Major Axis:'))
        ellipse_layout.addWidget(self.a_slider)
        ellipse_layout.addWidget(QLabel('Semi-Minor Axis:'))
        ellipse_layout.addWidget(self.b_slider)
        ellipse_layout.addWidget(self.ellipse_button)
        ellipse_layout.addWidget(self.fit_button)
        layout.addLayout(ellipse_layout)
    

    def detect_data_level(self, file):
        # method code
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
        # method code
        files, _ = QFileDialog.getOpenFileNames(
            self, 'Select FITS files', '', 'FITS files (*.fits *.fts)')
        if files:
            self.files = files
            self.maps = [sunpy.map.Map(f) for f in files]
            self.current_index = 0
            self.label.setText(f'{len(files)} files loaded.')


    def display_first_map(self):
        # method code
        if not self.maps:
            self.label.setText('No maps loaded.')
            return
        self.current_index = 0
        self.plot_map(self.maps[0])


    def show_prev_map(self):
        # method code
        if self.maps:
            self.current_index = max(0, self.current_index - 1)
            self.plot_map(self.maps[self.current_index])


    def show_next_map(self):
        # method code
        if self.maps:
            self.current_index = min(len(self.maps) - 1, self.current_index + 1)
            self.plot_map(self.maps[self.current_index])


    def plot_map(self, amap):
        # method code
        self.canvas.figure.clf()
        amap = self.maps[self.current_index]
        ax = self.canvas.figure.add_subplot(111, projection=amap)

        # Only add clip_interval if norm is NOT already defined
        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

        ax.grid(False)
        self.canvas.draw()
        self.label.setText(f'Showing: {self.files[self.current_index]}')

        # Add ellipse drawing to the main plot method
        self.draw_ellipse(redraw=False)  # Don't redraw canvas yet


    def plot_map_with_ellipse(self):
        """Plot current map and ellipse if exists."""
        self.canvas.figure.clf()
        amap = self.maps[self.current_index]
        ax = self.canvas.figure.add_subplot(111, projection=amap)

        if 'norm' in amap.plot_settings and amap.plot_settings['norm'] is not None:
            amap.plot(axes=ax)
        else:
            amap.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

        ax.grid(False)

        if hasattr(self, 'ellipse_center'):
            x0, y0 = self.ellipse_center
            a = self.a_slider.value()
            b = self.b_slider.value()
            theta = np.linspace(0, 2*np.pi, 300)
            x_ell = x0 + a * np.cos(theta)
            y_ell = y0 + b * np.sin(theta)
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, frame=amap.coordinate_frame)
            # plot_coord returns a list of Line2D, grab first
            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)

        self.canvas.draw()


    def upgrade_to_lv15(self):
        # method code
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
                    self.processed_maps.append(sunpy.map.Map(file))

                else:
                    print(f'Skipped unknown-level file: {file}')

            except Exception as e:
                print(f'Error upgrading {file}: {e}')

        self.maps = self.processed_maps
        self.current_index = 0
        self.plot_map(self.maps[0])
        self.label.setText(f'Upgraded and loaded {len(self.maps)} maps.')
    
        if self.processed_maps:
            self.maps = self.processed_maps
            self.current_index = 0
            self.plot_map(self.maps[0])
            self.label.setText(f'Loaded {len(self.maps)} maps.')
        else:
            self.label.setText('No maps were processed or loaded.')


    def set_processed_maps_from_loaded(self):
        # method code
        if not self.processed_maps and self.maps:
            self.processed_maps = self.maps.copy()


    def create_running_diff_maps(self):
        # method code
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
            # Get current map
            amap = self.maps[self.current_index]
            ax = self.canvas.figure.axes[0]  # Get current axes
            
            # Remove old ellipse if exists
            if self.ellipse_artist is not None:
                self.ellipse_artist.remove()
                self.ellipse_artist = None
                
            # Get parameters
            if self.ellipse_center is None:
                # First time drawing - get from inputs
                x0 = float(self.x_input.text())
                y0 = float(self.y_input.text())
                self.ellipse_center = (x0, y0)
            else:
                # Use stored center
                x0, y0 = self.ellipse_center
                
            a = self.a_slider.value()
            b = self.b_slider.value()
            
            # Calculate new ellipse coordinates
            x_ell = x0 + a * np.cos(self.theta)
            y_ell = y0 + b * np.sin(self.theta)
            coords = SkyCoord(x_ell * u.arcsec, y_ell * u.arcsec, 
                             frame=amap.coordinate_frame)
            
            # Draw new ellipse
            self.ellipse_artist, = ax.plot_coord(coords, color='red', lw=2)
            
            if redraw:
                self.canvas.draw()
                self.label.setText('Ellipse updated.')
                
        except Exception as e:
            self.label.setText(f'Ellipse error: {str(e)}')


    def update_ellipse(self):
        """Update ellipse when sliders change"""
        if self.ellipse_artist is not None:
            self.draw_ellipse()  # Redraw with new parameters


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


    def draw_ellipse_button(self):  # Renamed from draw_ellipse
        """Handle ellipse button click - force new ellipse"""
        try:
            # Get new center coordinates
            x0 = float(self.x_input.text())
            y0 = float(self.y_input.text())
            self.ellipse_center = (x0, y0)
            
            # Draw new ellipse
            self.draw_ellipse()
            self.label.setText('Ellipse drawn.')
        except ValueError:
            self.label.setText('Invalid center coordinates.')
    
    
    def save_current_map(self):
        # method code
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
















