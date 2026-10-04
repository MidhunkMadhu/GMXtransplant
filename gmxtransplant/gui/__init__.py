"""Desktop application. Importing this package never imports Qt, so the CLI never loads it."""


def display_available(platform, environ):
    # macOS uses Qt's native Cocoa backend and does not require X11 variables.
    return not platform.startswith('linux') or any(
        environ.get(key) for key in ('DISPLAY', 'WAYLAND_DISPLAY', 'QT_QPA_PLATFORM'))


def main():
    import os
    import sys
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print('The desktop application needs Qt (PySide6), which is not installed in this '
              'Python environment (for example after a command-line-only install). Add it with: '
              'python -m pip install "PySide6>=6.6,<7". The gmxtransplant command-line tool '
              'works without it.', file=sys.stderr)
        return 2
    if not display_available(sys.platform, os.environ):
        print('No graphical display found. Use a Linux desktop or WSL 2 with WSLg.', file=sys.stderr)
        return 2
    wayland = os.environ.get('WAYLAND_DISPLAY')
    if (sys.platform.startswith('linux') and wayland and os.environ.get('DISPLAY')
            and not os.environ.get('QT_QPA_PLATFORM')
            and not os.path.exists(os.path.join(os.environ.get('XDG_RUNTIME_DIR', ''), wayland))):
        # WSLg often advertises Wayland without a socket; use X11 directly
        # instead of failing over with a warning on every launch.
        os.environ['QT_QPA_PLATFORM'] = 'xcb'
    from .window import MainWindow
    app = QApplication(sys.argv)
    app.setApplicationName('GMXtransplant')
    app.setOrganizationName('GMXtransplant')
    app.setStyle('Fusion')
    window = MainWindow()
    window.show()
    return app.exec()
