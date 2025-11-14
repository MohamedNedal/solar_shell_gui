"""
Data processing module for AIA FITS files
Handles loading, calibration, and processing of solar data
"""

import os
import numpy as np
from scipy import ndimage
import matplotlib.colors as colors
from PyQt5.QtCore import QThread, pyqtSignal

import sunpy.map
from aiapy.calibrate import register, update_pointing
import astropy.units as u


class DataProcessor(QThread):
    """Background thread for processing FITS files"""
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    
    def __init__(self, files, channel, data_dir):
        super().__init__()
        self.files = files
        self.channel = channel
        self.data_dir = data_dir
        self.processed_maps = []
        
    def run(self):
        """Main processing loop"""
        try:
            self._create_output_directory()
            self._process_files()
            self.finished.emit(self.processed_maps)
        except Exception as e:
            self.error.emit(str(e))
            return
            
        self.finished.emit(self.processed_maps)
    
    def _create_output_directory(self):
        """Create output directory structure"""
        output_dir = os.path.join(self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15')
        os.makedirs(output_dir, exist_ok=True)
        
    def _process_files(self):
        """Process all FITS files"""
        total_files = len(self.files)
        
        for i, file_path in enumerate(self.files):
            try:
                self.status.emit(f'Processing file {i+1}/{total_files}: {os.path.basename(file_path)}')
                
                # Detect and process based on data level
                data_level = self._detect_data_level(file_path)
                processed_map = self._process_single_file(file_path, data_level)
                
                if processed_map is not None:
                    self.processed_maps.append(processed_map)
                else:
                    self.status.emit(f'Warning: Could not process file {i+1}')
                
                # Update progress
                progress_pct = int((i + 1) / total_files * 100)
                self.progress.emit(progress_pct)
                
            except Exception as e:
                self.status.emit(f'Error processing file {i+1}: {str(e)}')
                continue
    
    def _process_single_file(self, file_path, data_level):
        """Process a single FITS file based on its data level"""
        try:
            if data_level == 'lev1':
                return self._process_lev1_file(file_path)
            elif data_level == 'lev15':
                return self._load_lev15_file(file_path)
            else:
                # Try to load as-is for unknown level
                self.status.emit(f'Unknown data level, attempting direct load...')
                return sunpy.map.Map(file_path)
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            return None
    
    def _process_lev1_file(self, file_path):
        """Process Level 1 file to Level 1.5"""
        output_filename = os.path.basename(file_path).replace('lev1', 'lev15')
        output_path = os.path.join(
            self.data_dir, 'AIA', f'{self.channel}A', 'processed', 'lv15', output_filename
        )
        
        # Check if already processed
        if os.path.exists(output_path):
            self.status.emit(f'Loading existing processed file...')
            return sunpy.map.Map(output_path)
        
        # Load and process the file
        aia_map = sunpy.map.Map(file_path)
        
        # Update pointing information
        aia_map_updated = update_pointing(aia_map)
        
        # Register (align) the map
        aia_map_registered = register(aia_map_updated)
        
        # Normalize by exposure time
        aia_map_normalized = aia_map_registered / aia_map_registered.exposure_time
        
        # Save processed file
        aia_map_normalized.save(output_path, filetype='auto', overwrite=True)
        
        return aia_map_normalized
    
    def _load_lev15_file(self, file_path):
        """Load existing Level 1.5 file"""
        return sunpy.map.Map(file_path)
    
    def _detect_data_level(self, file_path):
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
                history_keys = [key for key in header.keys() if 'HISTORY' in key]
                for key in history_keys:
                    history = str(header.get(key, ''))
                    if 'register' in history.lower() or 'lev1.5' in history.lower():
                        return 'lev15'
                
                # Additional checks for calibration keywords that indicate lev1.5
                calibration_keywords = ['RSUN_REF', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']
                if all(keyword in header for keyword in calibration_keywords):
                    # Check if values look like they've been processed
                    cdelt1 = header.get('CDELT1', 0)
                    if abs(cdelt1 - 0.6) < 0.1:  # Typical lev1.5 pixel scale
                        return 'lev15'
            
            return 'unknown'
            
        except Exception as e:
            print(f"Error detecting data level for {file_path}: {e}")
            return 'unknown'


class RunningDifferenceProcessor:
    """Class for creating running difference maps"""
    
    @staticmethod
    def create_running_diff_maps(processed_maps, skip_frames=5, gaussian_sigma=3):
        """
        Create running difference maps from processed data
        
        Parameters:
        -----------
        processed_maps : list
            List of SunPy maps
        skip_frames : int
            Number of frames to skip for difference (default: 5)
        gaussian_sigma : float
            Sigma for Gaussian smoothing (default: 3)
            
        Returns:
        --------
        list : Running difference maps
        """
        if len(processed_maps) < skip_frames + 1:
            raise ValueError(f'Need at least {skip_frames + 1} images for running difference')
        
        running_diff_maps = []
        
        for i in range(skip_frames, len(processed_maps)):
            try:
                map_base = processed_maps[i - skip_frames]
                map_current = processed_maps[i]
                
                # Calculate difference
                diff_data = map_current.quantity - map_base.quantity
                
                # Apply Gaussian smoothing
                if gaussian_sigma > 0:
                    diff_data = ndimage.gaussian_filter(diff_data.value, sigma=[gaussian_sigma, gaussian_sigma]) * diff_data.unit
                
                # Create new map with difference data
                diff_map = sunpy.map.Map(diff_data, map_current.meta)
                
                # Set appropriate plot settings for difference map
                diff_map.plot_settings['norm'] = colors.Normalize(vmin=-50, vmax=50)
                diff_map.plot_settings['cmap'] = 'RdBu_r'
                
                running_diff_maps.append(diff_map)
                
            except Exception as e:
                print(f"Error creating running difference map {i}: {e}")
                continue
                
        return running_diff_maps