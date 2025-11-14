"""
Custom matplotlib widget for embedding plots in PyQt5
"""

from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class MatplotlibWidget(QWidget):
    """Custom widget for embedding matplotlib plots"""
    
    def __init__(self, parent=None, figsize=(10, 8)):
        super().__init__(parent)
        
        # Create matplotlib figure and canvas
        self.figure = Figure(figsize=figsize, dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self)
        
        # Set up layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        
        # Configure matplotlib for better appearance
        self.figure.patch.set_facecolor('white')
        plt.style.use('default')
        
    def clear(self):
        """Clear the figure and redraw"""
        self.figure.clear()
        self.canvas.draw()
        
    def save_figure(self, filename, **kwargs):
        """Save the current figure to file"""
        default_kwargs = {
            'dpi': 300,
            'bbox_inches': 'tight',
            'facecolor': 'white',
            'edgecolor': 'none'
        }
        default_kwargs.update(kwargs)
        self.figure.savefig(filename, **default_kwargs)