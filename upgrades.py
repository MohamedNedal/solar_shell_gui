import os
import sunpy.map
from aiapy.calibrate import register  # only register is needed now

def upgrade_to_lv15(viewer):
    """
    Upgrade loaded AIA maps from level 1 to level 1.5 using aiapy.register.
    The processed maps are stored in viewer.processed_maps and viewer.maps is updated.
    """
    if not viewer.files:
        from PyQt5.QtWidgets import QMessageBox
        QMessageBox.warning(viewer, 'Warning', 'No files to upgrade.')
        return

    viewer.processed_maps = []
    n = len(viewer.files)
    viewer.progress.setVisible(True)
    viewer.progress.setRange(0, n)
    viewer.progress.setValue(0)
    viewer.label.setText('Upgrading files to level 1.5...')
    from PyQt5.QtWidgets import QApplication
    QApplication.processEvents()

    for i, file in enumerate(viewer.files):
        try:
            # detect file level
            level = detect_data_level(file)
            if level == 'lev1':
                # Generate output path in temporary processed folder
                output_filename = os.path.basename(file).replace('lev1', 'lev15')
                file_path = os.path.join(viewer.temp_data_dir, 'AIA', f'{viewer.channel}A', 'processed', 'lv15')
                os.makedirs(file_path, exist_ok=True)
                full_path = os.path.join(file_path, output_filename)

                if not os.path.exists(full_path):
                    # load map and register to level 1.5
                    m = sunpy.map.Map(file)
                    m = register(m)  # only registration is applied
                    m = m / m.exposure_time  # normalize
                    m.save(full_path, filetype='auto')
                    viewer.processed_maps.append(m)
                else:
                    viewer.processed_maps.append(sunpy.map.Map(full_path))

            elif level == 'lev15':
                viewer.processed_maps.append(sunpy.map.Map(file))
            else:
                print(f'Skipped unknown-level file: {file}')

        except Exception as e:
            print(f'Error upgrading {file}: {e}')

        viewer.progress.setValue(i + 1)
        QApplication.processEvents()

    viewer.maps = viewer.processed_maps
    viewer.current_index = 0
    viewer.progress.setVisible(False)
    if viewer.maps:
        viewer.plot_map(viewer.maps[0])
        viewer.label.setText(f'Upgraded and loaded {len(viewer.maps)} maps.')
    else:
        viewer.label.setText('No maps were processed or loaded.')

def detect_data_level(file):
    """
    Determine the data level of a FITS file: 'lev1', 'lev15', or 'unknown'.
    """
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