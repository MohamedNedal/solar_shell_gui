import os
import numpy as np
import sunpy.map
from aiapy.calibrate import register, update_pointing
from scipy import ndimage
from matplotlib import colors
import astropy.units as u


class AIAProcessor:

    def detect_data_level(self, file):
        name = os.path.basename(file).lower()
        if "lev1.5" in name or "lev15" in name:
            return "lev15"
        elif "lev1" in name:
            return "lev1"
        else:
            try:
                hdr = sunpy.map.Map(file).meta
                level = hdr.get("lvl_num", None)
                if level == 1.5:
                    return "lev15"
                elif level == 1:
                    return "lev1"
            except Exception:
                pass
        return "unknown"

    def load_maps(self, files):
        maps = []
        for f in files:
            try:
                maps.append(sunpy.map.Map(f))
            except Exception as e:
                print(f"Error loading {f}: {e}")
        maps.sort(key=lambda m: m.date)
        return maps

    def upgrade_to_lv15(self, files, channel, temp_dir):
        from aiapy.calibrate.util import get_pointing_table

        pt_table = get_pointing_table("lmsal")
        processed_maps = []

        for file in files:
            level = self.detect_data_level(file)
            try:
                if level == "lev1":
                    m = sunpy.map.Map(file)
                    m = update_pointing(m, pointing_table=pt_table)
                    m = register(m)
                    m = m / m.exposure_time
                    processed_maps.append(m)
                elif level == "lev15":
                    processed_maps.append(sunpy.map.Map(file))
            except Exception as e:
                print(f"Upgrade error: {e}")

        return processed_maps

    def create_running_diff(self, maps, vmin=-50, vmax=50):
        if len(maps) < 6:
            return []

        diff_maps = []
        for i in range(5, len(maps)):
            m0 = maps[i - 5]
            m1 = maps[i]
            diff = m1.quantity - m0.quantity
            smoothed = ndimage.gaussian_filter(diff, sigma=[3, 3])
            diff_map = sunpy.map.Map(smoothed, m1.meta)
            diff_map.plot_settings["norm"] = colors.Normalize(
                vmin=vmin, vmax=vmax
            )
            diff_maps.append(diff_map)

        return diff_maps