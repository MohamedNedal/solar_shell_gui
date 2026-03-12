# AIA FITS Viewer — Solar CME Shell Analyser

An interactive desktop application for loading, processing, and visualising AIA (Atmospheric Imaging Assembly) solar EUV image sequences, fitting CME ellipse regions of interest, and rendering 3-D magnetic field–shell interaction plots.

---

## Features

- **Load & navigate** multi-frame AIA FITS sequences, sorted chronologically
- **Upgrade** raw level-1 files to level-1.5 (pointing correction, registration, exposure normalisation)
- **Running-difference movies** with Gaussian smoothing for CME detection
- **Ellipse overlay** with interactive centre input and semi-axis sliders
- **Fit extraction** — converts ellipse axes from arcsec to solar radii and logs results to CSV
- **PFSS magnetic field model** computed from GONG synoptic magnetograms
- **Interactive 3-D visualisation** of the CME ellipsoid shell coloured by the angle between outward surface normals and local magnetic field-line directions
- **Silent PNG export** of 3-D frames for batch animation workflows

---

## Project Structure

```
├── main.py               # Entry point — run this
├── viewer.py             # PyQt5 main window; wires UI to all modules
├── fits_processing.py    # Load, level-upgrade, and running-diff (pure functions)
├── ellipse_tools.py      # Parameter extraction, unit conversion, CSV/text export
├── pfss_model.py         # GONG load + PFSS solve + field-line tracing
├── visualization_3d.py   # Plotly 3-D figure builder
└── README.md
```

The application follows a clean separation of concerns: `viewer.py` owns the GUI state and delegates all computation to the four helper modules, none of which import Qt.

---

## Requirements

### Python version

Python **3.9 – 3.11** is recommended. Python 3.12+ may require patched versions of some solar physics packages.

### Dependencies

Install all dependencies with:

```bash
pip install -r requirements.txt
```

`requirements.txt`:

```
PyQt5>=5.15
matplotlib>=3.7
sunpy>=5.0
aiapy>=0.7
astropy>=5.3
scipy>=1.11
numpy>=1.24
pfsspy>=1.1
plotly>=5.18
kaleido>=0.2          # required only for silent PNG export (export_3d_png_silent)
```

> **Note — PyQt5 on macOS (Apple Silicon):** install via conda if pip fails:
> ```bash
> conda install pyqt
> ```

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/your-username/aia-fits-viewer.git
cd aia-fits-viewer

# 2. (Recommended) create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt
```

---

## Running the Application

```bash
python main.py
```

This opens the main GUI window. All other modules are loaded automatically.

---

## Walkthrough

### 1 — Load FITS files

Click **Load FITS** and select one or more AIA `.fits` / `.fts` files. Files are sorted by observation time and the first frame is displayed immediately.

### 2 — (Optional) Upgrade to level 1.5

If your files are raw level-1 data, click **Upgrade**. The app applies `aiapy` pointing correction and image registration, normalises by exposure time, and caches the results in your system temp directory so subsequent runs skip the processing step.

### 3 — (Optional) Running difference

Click **Run-Diff** to compute a running-difference movie (requires ≥ 6 frames). Each frame is the difference between frame *i* and frame *i − 5*, Gaussian-smoothed and displayed on a symmetric ±50 DN/s colour scale. Navigate with **Previous / Next / First Image / Last Image**.

### 4 — Draw an ellipse

Enter an X and Y centre in the right-hand panel (in arcsec in the helioprojective frame), adjust the **Semi-Major** and **Semi-Minor** axis sliders, then click **Ellipse**. The overlay updates live as you move the sliders.

### 5 — Extract fit parameters

Click **Fit** to compute the ellipse axes in arcsec and solar radii. Results appear under the button and are appended to `ellipse_fits.csv` in your system temp directory. Click **Export Params** to also write `ellipse_fit.csv` and `theta_angles.txt` to the script directory.

### 6 — Load a GONG magnetogram (optional)

Click **Select GONG File** and pick a GONG synoptic magnetogram FITS file. Then click **Calculate PFSS** to solve the Potential Field Source Surface model and trace a 40 × 60 grid of magnetic field lines. This step is required for angle-coloured 3-D shells.

> GONG data can be downloaded from the [NSO GONG archive](https://gong2.nso.edu/archive/patch.pl?menutype=s).

### 7 — 3-D visualisation

Click **Show 3D Ellipsoid**. The figure opens as an interactive HTML page in your default browser. Toggle the checkboxes on the right panel to show or hide the shell surface, outward normal vectors, radial grid lines, and magnetic field lines.

To export a static PNG of the current frame, click **Export 3D PNG** (requires `kaleido`). Files are saved as `frame_001.png`, `frame_002.png`, etc. in the script directory.

---

## Output Files

| File | Location | Contents |
|---|---|---|
| `ellipse_fits.csv` | System temp dir | Running log of every **Fit** button press |
| `ellipse_fit.csv` | Script directory | Last exported fit (x0, y0, a, b) |
| `theta_angles.txt` | Script directory | Flat array of shell–field angles (degrees) |
| `frame_NNN.png` | Script directory | Silent 3-D PNG exports |
| `lv15/` cache | System temp dir | Processed level-1.5 FITS files |

---

## Data Sources

| Data | Source |
|---|---|
| AIA FITS files | [SDO/AIA via JSOC](http://jsoc.stanford.edu) or [SunPy's `Fido`](https://docs.sunpy.org/en/stable/generated/api/sunpy.net.Fido.html) |
| GONG magnetograms | [NSO GONG archive](https://gong2.nso.edu/archive/patch.pl?menutype=s) |

---

## Known Limitations

- **PFSS tracing** uses `pfsspy`'s Fortran tracer and can take several minutes for a full 40 × 60 seed grid.
- **Silent PNG export** requires the `kaleido` package. If it is missing, Plotly will raise a `ValueError`; the interactive HTML export works without it.
- Running-difference requires at least **6 frames**.
- The level-1 → 1.5 upgrade calls `aiapy.calibrate.update_pointing` with `pointing_table=None`, which fetches the pointing table from the internet. An offline fallback is not currently implemented.

---

## License

MIT License. See `LICENSE` for details.
