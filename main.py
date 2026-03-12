import sys
from PyQt5.QtWidgets import QApplication
from viewer_ui import AIAViewer   # GUI + core logic

def main():
    app = QApplication(sys.argv)
    viewer = AIAViewer()
    viewer.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()