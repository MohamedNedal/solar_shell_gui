from PyQt5.QtWebEngineWidgets import QWebEngineView

def show_plotly_figure(fig):
    """Display Plotly figure in a QWebEngineView"""
    html = fig.to_html(include_plotlyjs='cdn')
    web_view = QWebEngineView()
    web_view.setHtml(html)
    return web_view