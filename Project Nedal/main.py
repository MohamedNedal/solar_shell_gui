"""
main.py
-------
Entry point for the AIA FITS Viewer application.

Run this script directly:

    python main.py

All other modules (viewer, fits_processing, ellipse_tools, pfss_model,
visualization_3d) are imported transitively from here.
"""

import sys
from PyQt5.QtWidgets import QApplication
from viewer import AIAViewer

def main():
    """Initialise the Qt application, create the main window, and start the event loop."""
    app = QApplication(sys.argv)
    window = AIAViewer()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

