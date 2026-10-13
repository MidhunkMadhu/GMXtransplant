"""GMXtransplant desktop workbench."""
from copy import deepcopy
import codecs
import json
import os
from pathlib import Path
import signal
import subprocess
import shutil
import sys
import tempfile

from PySide6.QtCore import Qt, QEvent, QProcess, QProcessEnvironment, QSettings, QTimer
from PySide6.QtGui import QKeySequence, QShortcut, QTextCursor
from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QLineEdit, QStackedWidget, QPlainTextEdit, QButtonGroup,
    QFrame, QProgressBar, QListWidget, QScrollArea)

from viewers import find_viewer, launch as launch_viewer
from mode_registry import EDIT_MODES
from .editor import ConfigEditor, choose_path
from .model import (MODES, FOLDERS, EXAMPLE_KEYS, MORE_EXAMPLE_KEYS, build_job, load_yaml, example_mode, example_folder,
                    example_path, destination, example_document,
                    template_document, open_documentation, open_with_default_app)
from .progress import RunProgress

APPLICATION_TITLE = 'GMXTRANSPLANT MD System Builder'

def report_files(root):
    """The .txt reports of a run: in its folder and one folder down (e.g. addbinder poses)."""
    root = Path(root)
    found = [p for p in root.glob('*.txt')] + [p for p in root.glob('*/*.txt')]
    return [str(p.relative_to(root)) for p in sorted(found)
            if not any(part.startswith('.') for part in p.relative_to(root).parts)]


# Shown under the mode selector: what the selected mode does.
MODE_DESCRIPTIONS = {
    'charmprot': 'Replace the proteins and ligands of a prepared CHARMM-GUI system with those of a second '
                 'CHARMM-GUI system, keeping its membrane, water and ions.',
    'protein': 'Insert a new protein, with any bound ligands, into a prepared membrane system in place of the old one.',
    'lig': 'Replace one bound ligand with another; the protein, membrane, water and ions stay as they are.',
    'chl': 'Restore experimentally resolved cholesterol into a prepared membrane system and adjust the '
           'lipid composition. The restored cholesterols start with few lipid neighbours, so the run also '
           'writes repack/, a system set up for a short equilibration (about 1.4 ns) that holds the protein, '
           'ligands and restored cholesterols while the other lipids fill the gaps; a sample run script is '
           'also saved.',
    'addbinder': 'Place a ligand or protein in the water above or below a membrane protein, a set distance from '
                 'its outermost atom, and build one complete system per pose.',
    'popc': 'Replace the membrane of a prepared system with a well-equilibrated pure POPC bilayer (500 ns, '
            'CHARMM36), or with any equilibrated bilayer you give as a GRO file. The protein keeps its OPM '
            'placement; overlapping lipids, water and ions are removed, and the run writes equilibration/, a '
            'GROMACS system set up for the CHARMM-GUI step6 equilibration with a sample run script.',
    **{key: mode.description for key, mode in EDIT_MODES.items()},
    'minimize': 'Prepare a restrained OpenMM minimization folder for a system you have already assembled.',
}

SANS = '"Inter", "SF Pro Text", "Helvetica Neue", "Segoe UI", "DejaVu Sans", sans-serif'
MONO = '"Menlo", "SF Mono", "DejaVu Sans Mono", "Consolas", monospace'

# Three themes: a clean light one, a dark one in black, graphite and brushed
# silver, and a colourful one.
THEMES = {
    'Light': dict(
        bg='#f4f6fa', surface='#ffffff', list_bg='#f7f8fb', text='#1b2230', muted='#5f6b7d',
        border='#e3e7ee', input_border='#cfd6e1', accent='#3b6cf6', accent_hover='#2f5be0',
        on_accent='#ffffff', accent_soft='#e8eeff', accent_soft_text='#2347b8',
        sidebar='#ffffff', sidebar_border='#e3e7ee', sidebar_text='#3a4456', sidebar_muted='#97a1b2',
        nav_hover='#f0f3f9', nav_active_bg='#e8eeff', nav_active_text='#2347b8',
        brand='#111827', badge='#3b6cf6', hero='#111827',
        button='#ffffff', button_border='#d6dce6', button_hover='#f3f5f9', button_pressed='#e8ecf3',
        primary='#3b6cf6', primary_hover='#2f5be0', disabled_bg='#eef1f5', disabled_text='#9aa3b2',
        console_bg='#0f172a', console_text='#d6e2ff', tile_bg='#ffffff', tile_selected_bg='#f1f5ff',
        notice='#b45309', danger='#c03a3a', danger_bg='#fdeeee', danger_border='#f2c4c4',
        info_bg='#eef4ff', info_border='#c9d8fb', scroll='#cdd4df',
        tooltip_bg='#111827', tooltip_text='#f9fafb', font=SANS, heading_font=SANS,
        icon_fill='#dbe5ff', swatch=('#ffffff', '#3b6cf6'), dark=False),
    'Dark': dict(
        bg='#141414', surface='#1d1d1e', list_bg='#18181a', text='#e6e6e6', muted='#9c9c9c',
        border='#343434', input_border='#4b4b4d', accent='#c3c7cd', accent_hover='#e2e4e8',
        on_accent='#121212', accent_soft='#2c2d30', accent_soft_text='#eceef1',
        sidebar='qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #2c2c2d, stop:1 #0c0c0c)',
        sidebar_border='#3b3b3c', sidebar_text='#d2d2d2', sidebar_muted='#7c7c7c',
        nav_hover='#262627',
        nav_active_bg='qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #4c4d50, stop:1 #303133)',
        nav_active_text='#ffffff', brand='#ededed', badge='#b9bdc3', hero='#f2f2f2',
        button='qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #3b3b3c, stop:1 #262627)',
        button_solid='#303031', button_border='#505052', button_hover='#424244', button_pressed='#2a2a2b',
        primary='qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #f4f5f7, stop:0.48 #c9cdd2, stop:0.52 #b8bcc2, stop:1 #9da2a9)',
        primary_hover='qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #ffffff, stop:0.48 #dcdfe3, stop:0.52 #caced3, stop:1 #b0b5bb)',
        disabled_bg='#1f1f20', disabled_text='#5c5c5c',
        console_bg='#0a0a0a', console_text='#d4d4d4', tile_bg='#1a1a1b', tile_selected_bg='#2a2b2d',
        notice='#d6b98c', danger='#dd7a72', danger_bg='#2b1a19', danger_border='#5c302d',
        info_bg='#26272a', info_border='#4b4d51', scroll='#4b4b4d',
        tooltip_bg='#000000', tooltip_text='#e2e2e2', font=SANS,
        heading_font='"Baskerville", "Georgia", "DejaVu Serif", "Times New Roman", serif',
        icon_fill='#3b3d41', swatch=('#1d1d1e', '#c3c7cd'), dark=True),
    'Aurora': dict(
        bg='#120f24', surface='#1b1733', list_bg='#16122b', text='#ece9f9', muted='#a59fc7',
        border='#2d2752', input_border='#3d3570', accent='#f472b6', accent_hover='#f9a8d4',
        on_accent='#ffffff', accent_soft='#3a1f45', accent_soft_text='#fbcfe8',
        sidebar='qlineargradient(x1:0, y1:0, x2:0.35, y2:1, stop:0 #3b1d6e, stop:0.55 #1e1a4d, stop:1 #0d3b4f)',
        sidebar_border='#2d2752', sidebar_text='#ddd8f5', sidebar_muted='#8f86c0',
        nav_hover='rgba(255, 255, 255, 18)', nav_active_bg='rgba(244, 114, 182, 60)', nav_active_text='#ffffff',
        brand='#ffffff', badge='#67e8f9', hero='#ffffff',
        button='#231e42', button_border='#3d3570', button_hover='#2b2550', button_pressed='#342d60',
        primary='qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #ec4899, stop:1 #8b5cf6)',
        primary_hover='qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #f472b6, stop:1 #a78bfa)',
        disabled_bg='#1f1b38', disabled_text='#5f5888',
        console_bg='#0b0918', console_text='#c4b5fd', tile_bg='#1f1a3a', tile_selected_bg='#2a1d48',
        notice='#fcd34d', danger='#fb7185', danger_bg='#3a1626', danger_border='#6b2442',
        info_bg='#1c2a4a', info_border='#35508a', scroll='#3d3570',
        tooltip_bg='#0b0918', tooltip_text='#f5f3ff', font=SANS, heading_font=SANS,
        icon_fill='#4a2353', swatch=('#3b1d6e', '#f472b6'), dark=True),
}
DEFAULT_THEME = 'Light'
# Themes of earlier versions, kept in user settings.
LEGACY_THEMES = {'Teal': 'Light', 'Sand': 'Light', 'Midnight': 'Dark', 'VMD Classic': 'Dark', 'Vintage': 'Dark'}

STYLE = """
QWidget {{ font-family: {font}; font-size: 13px; color: {text}; }}
QMainWindow, QDialog, QWidget#workspace, QStackedWidget#pages {{ background: {bg}; }}
QWidget#page {{ background: {surface}; }}
QWidget#transparent {{ background: transparent; }}
QFrame#rail {{ background: {sidebar}; border: 1px solid {sidebar_border}; border-radius: 18px; }}
QFrame#rail QLabel {{ background: transparent; }}
QLabel#brand {{ color: {brand}; font-family: {heading_font}; font-size: 19px; font-weight: 800; }}
QLabel#brandSub {{ color: {badge}; font-size: 10px; font-weight: 700; letter-spacing: 2px; }}
QLabel#railHeading {{ color: {sidebar_muted}; font-size: 10px; font-weight: 700; letter-spacing: 2px; }}
QPushButton#nav {{ background: transparent; color: {sidebar_text}; border: none; border-radius: 10px;
    padding: 11px 14px; text-align: left; font-size: 14px; font-weight: 600; }}
QPushButton#nav:hover {{ background: {nav_hover}; }}
QPushButton#nav:checked {{ background: {nav_active_bg}; color: {nav_active_text}; }}
QPushButton#swatch {{ background: transparent; color: {sidebar_text}; border: 1px solid transparent;
    border-radius: 9px; padding: 6px 10px; text-align: left; font-weight: 600; }}
QPushButton#swatch:hover {{ background: {nav_hover}; }}
QPushButton#swatch:checked {{ border: 1px solid {badge}; color: {brand}; }}
QPushButton#railButton {{ background: {nav_hover}; color: {sidebar_text}; border: 1px solid {sidebar_border};
    border-radius: 10px; padding: 9px; font-weight: 600; }}
QPushButton#railButton:hover {{ color: {brand}; border-color: {badge}; }}
QLabel#hero {{ font-family: {heading_font}; font-size: 23px; font-weight: 750; color: {hero}; }}
QLabel#subtitle {{ color: {muted}; font-size: 13px; }}
QFrame#outputBar {{ background: {surface}; border: 1px solid {border}; border-radius: 12px; }}
QFrame#outputBar QLineEdit {{ border: none; background: transparent; padding: 6px 2px; }}
QLabel#overline {{ color: {muted}; font-size: 10px; font-weight: 700; letter-spacing: 1.5px; }}
QFrame#panel, QFrame#card {{ background: {surface}; border: 1px solid {border}; border-radius: 14px; }}
QFrame#strip {{ background: {accent_soft}; border: none; border-radius: 10px; }}
QFrame#strip QLabel {{ color: {accent_soft_text}; background: transparent; }}
QFrame#tile {{ background: {tile_bg}; border: 1px solid {border}; border-radius: 12px; }}
QFrame#tile:hover {{ border: 1px solid {accent}; }}
QFrame#tile[selected="true"] {{ background: {tile_selected_bg}; border: 2px solid {accent}; }}
QFrame#tile QLabel {{ background: transparent; }}
QLabel#tileGlyph {{ font-size: 19px; color: {accent}; font-weight: 700; }}
QLabel#tileTitle {{ font-weight: 700; font-size: 13px; }}
QLabel#tileText {{ color: {muted}; font-size: 11px; }}
QFrame#exampleCard {{ background: {surface}; border: 1px solid {border}; border-radius: 14px; }}
QFrame#exampleCard:hover {{ border: 1px solid {accent}; }}
QFrame#exampleCard QLabel {{ background: transparent; }}
QLabel#cardTitle {{ font-family: {heading_font}; font-size: 15px; font-weight: 700; }}
QLabel#chip {{ background: {accent_soft}; color: {accent_soft_text}; border-radius: 7px; padding: 2px 8px;
    font-size: 10px; font-weight: 700; letter-spacing: 1px; }}
QLabel#path {{ color: {muted}; font-family: {heading_font}; font-size: 11px; }}
QLabel#sectionTitle {{ font-family: {heading_font}; font-size: 18px; font-weight: 700; margin: 4px 0; }}
QLabel#runTitle {{ font-family: {heading_font}; font-size: 21px; font-weight: 750; color: {hero}; }}
QLabel#fieldLabel {{ font-weight: 600; }}
QCheckBox#fieldToggle {{ font-weight: 600; }}
QLabel#muted {{ color: {muted}; }}
QLabel#notice {{ color: {notice}; }}
QFrame#fieldRow {{ background: transparent; border: none; }}
QFrame#banner {{ border-radius: 10px; padding: 2px; }}
QFrame#banner[level="error"] {{ background: {danger_bg}; border: 1px solid {danger_border}; }}
QFrame#banner[level="error"] QLabel {{ color: {danger}; }}
QFrame#banner[level="info"] {{ background: {info_bg}; border: 1px solid {info_border}; }}
QLineEdit, QComboBox, QPlainTextEdit {{ background: {surface}; border: 1px solid {input_border}; border-radius: 8px;
    padding: 7px 9px; selection-background-color: {accent_soft}; selection-color: {text}; }}
QLineEdit:focus, QComboBox:focus, QPlainTextEdit:focus {{ border: 1px solid {accent}; }}
QComboBox::drop-down {{ width: 26px; border: none; }}
QComboBox::down-arrow {{ image: url("{arrow_icon}"); width: 10px; height: 10px; }}
QComboBox QAbstractItemView {{ background: {surface}; color: {text}; selection-background-color: {accent_soft};
    selection-color: {text}; }}
QPushButton {{ background: {button}; border: 1px solid {button_border}; border-radius: 9px; padding: 8px 14px;
    font-weight: 600; }}
QPushButton:hover {{ background: {button_hover}; }}
QPushButton:pressed {{ background: {button_pressed}; }}
QPushButton#primary {{ background: {primary}; color: {on_accent}; border: 1px solid {accent}; padding: 8px 18px; }}
QPushButton#primary:hover {{ background: {primary_hover}; }}
QPushButton#danger {{ color: {danger}; background: {danger_bg}; border-color: {danger_border}; }}
QPushButton#quiet, QPushButton#disclosure {{ background: transparent; color: {muted}; border: none; padding: 4px 8px; }}
QPushButton#disclosure {{ color: {accent}; font-weight: 700; padding-left: 0; }}
QPushButton#quiet:hover, QPushButton#disclosure:hover {{ color: {accent}; }}
QPushButton:disabled, QLineEdit:disabled, QComboBox:disabled {{ color: {disabled_text}; background: {disabled_bg};
    border-color: {border}; }}
QPushButton#primary:disabled, QPushButton#danger:disabled {{ color: {disabled_text}; background: {disabled_bg};
    border-color: {border}; }}
QCheckBox {{ spacing: 8px; padding: 3px; background: transparent; }}
QCheckBox::indicator {{ width: 16px; height: 16px; border: 1px solid {input_border}; border-radius: 5px;
    background: {surface}; }}
QCheckBox::indicator:hover {{ border-color: {accent}; }}
QCheckBox::indicator:checked {{ background: {accent}; border-color: {accent}; image: url("{check_icon}"); }}
QCheckBox::indicator:disabled {{ background: {disabled_bg}; border-color: {border}; }}
QListWidget {{ background: {list_bg}; border: none; border-radius: 10px; padding: 6px; outline: none; }}
QListWidget::item {{ padding: 9px 10px; border-radius: 7px; margin-bottom: 2px; }}
QListWidget::item:hover {{ background: {button_hover}; }}
QListWidget::item:selected {{ background: {accent_soft}; color: {accent_soft_text}; }}
QScrollArea {{ border: none; background: transparent; }}
QScrollArea > QWidget > QWidget#transparent {{ background: transparent; }}
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: {scroll}; border-radius: 4px; min-height: 30px; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0px; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 2px; }}
QScrollBar::handle:horizontal {{ background: {scroll}; border-radius: 4px; min-width: 30px; }}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{ width: 0px; }}
QProgressBar {{ background: {border}; border: none; border-radius: 3px; }}
QProgressBar::chunk {{ background: {primary}; border-radius: 3px; }}
QPlainTextEdit#console {{ background: {console_bg}; color: {console_text}; border: none; border-radius: 12px;
    font-family: {mono}; font-size: 12px; padding: 14px; }}
QSplitter::handle {{ background: transparent; }}
QToolTip {{ background: {tooltip_bg}; color: {tooltip_text}; padding: 7px; border: none; border-radius: 6px; }}
QStatusBar {{ background: {bg}; color: {muted}; }}
QFileDialog QListView, QFileDialog QTreeView {{ background: {surface}; color: {text}; }}
"""


def theme_icons(name, theme):
    """Check-mark and drop-down arrow drawn in the theme's colours."""
    folder = Path(tempfile.gettempdir()) / f'gmxtransplant-{os.getuid() if hasattr(os, "getuid") else 0}' / name.replace(' ', '_')
    folder.mkdir(parents=True, exist_ok=True)
    icons = {
        'check_icon': ('check.svg', f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16"><path d="M3.5 8.5l3 3 6-7" '
                                    f'fill="none" stroke="{theme["on_accent"]}" stroke-width="2.2" stroke-linecap="round" '
                                    f'stroke-linejoin="round"/></svg>'),
        'arrow_icon': ('arrow.svg', f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"><path d="M1.5 3.5l3.5 3.5 3.5-3.5" '
                                    f'fill="none" stroke="{theme["muted"]}" stroke-width="1.6" stroke-linecap="round" '
                                    f'stroke-linejoin="round"/></svg>'),
    }
    paths = {}
    for key, (filename, svg) in icons.items():
        path = folder / filename
        if not path.is_file() or path.read_text() != svg:
            path.write_text(svg)
        paths[key] = path.as_posix()
    return paths


def apply_theme(name):
    """Style every window of the application (main window and dialogs alike)."""
    from PySide6.QtGui import QPalette, QColor
    from PySide6.QtWidgets import QApplication
    name = name if name in THEMES else DEFAULT_THEME
    theme = THEMES[name]
    try:
        icons = theme_icons(name, theme)
    except OSError:
        icons = {'check_icon': '', 'arrow_icon': ''}  # Plain indicators still show state.
    app = QApplication.instance()
    palette = QPalette()
    for role, key in [(QPalette.ColorRole.Window, 'bg'), (QPalette.ColorRole.WindowText, 'text'),
                      (QPalette.ColorRole.Base, 'surface'), (QPalette.ColorRole.AlternateBase, 'list_bg'),
                      (QPalette.ColorRole.Text, 'text'), (QPalette.ColorRole.Button, 'button_solid'),
                      (QPalette.ColorRole.ButtonText, 'text'), (QPalette.ColorRole.Highlight, 'accent'),
                      (QPalette.ColorRole.HighlightedText, 'on_accent'),
                      (QPalette.ColorRole.ToolTipBase, 'tooltip_bg'),
                      (QPalette.ColorRole.ToolTipText, 'tooltip_text'), (QPalette.ColorRole.PlaceholderText, 'muted')]:
        palette.setColor(role, QColor(theme.get(key, theme.get('button'))))
    app.setPalette(palette)
    app.setStyleSheet(STYLE.format(**theme, **icons, mono=MONO))


def label(text, name=None):
    result = QLabel(text)
    result.setWordWrap(True)
    if name:
        result.setObjectName(name)
    return result


def button(text, callback, primary=False):
    result = QPushButton(text)
    result.clicked.connect(callback)
    if primary:
        result.setObjectName('primary')
    return result


def repolish(widget):
    widget.style().unpolish(widget)
    widget.style().polish(widget)



# Example cards: the original titles and descriptions.
EXAMPLE_TITLES = {'popc': 'Equilibrated POPC membrane: D1R',
                  'popc_cholesterol': 'Equilibrated POPC/cholesterol membrane: D1R',
                  'addbinder': 'Add binder: dopamine (extracellular)',
                  'addbinder_gprotein': 'Add binder: G protein (intracellular)',
                  **{key: mode.example_title for key, mode in EDIT_MODES.items()}}
EXAMPLE_DETAILS = {
    'charmprot': 'Choose reference and transplant folders. Replace proteins and ligands automatically while retaining the reference environment.',
    'protein': 'Insert a protein and its bound ligands into a prepared membrane system using matching topology inputs.',
    'lig': 'Replace one bound ligand (LI1) with another (LI2) using explicit heavy-atom correspondence, its ITP and its force field.',
    'chl': 'Restore five experimental cholesterol molecules into the prepared D1R system. Uses Open Babel, installed with GMXtransplant.',
    'popc': 'Replace the POPC/POPE/cholesterol membrane of the apo D1 receptor build with equilibrated pure POPC; the receptor keeps its OPM placement. Writes the system and a GROMACS equilibration folder.',
    'popc_cholesterol': 'Your own bilayer GRO: replace the membrane of the apo D1 receptor build with an equilibrated POPC/cholesterol 70:30 bilayer; the receptor keeps its OPM placement.',
    'addbinder': 'Extracellular ligand: place dopamine above the apo D1 receptor, once straight above at 20 Å and three times from random directions at 20–25 Å; one run-ready folder each.',
    **{key: mode.example_description for key, mode in EDIT_MODES.items()},
    'addbinder_gprotein': 'Intracellular protein: place the Gs heterotrimer below the apo D1 receptor in its receptor-bound arrangement: pulled straight out to 20 Å (shorter if the box needs it) and towards two random sides at 15 Å, one of them also randomly rotated; one run-ready folder each.',
}
PAGES = ['Configuration', 'Examples', 'Run & results']


def disclosure(text, panel):
    """A '▸ text' button that shows and hides `panel` (collapsed at first)."""
    toggle = QPushButton('▸  ' + text)
    toggle.setObjectName('disclosure')
    toggle.setCheckable(True)
    panel.hide()

    def expand(shown):
        panel.setVisible(shown)
        toggle.setText(('▾  ' if shown else '▸  ') + text)
    toggle.toggled.connect(expand)
    return toggle


class Tile(QFrame):
    """A clickable card: a mode on the Build page."""

    def __init__(self, title, callback, parent=None):
        super().__init__(parent)
        self.setObjectName('tile')
        self.setAttribute(Qt.WidgetAttribute.WA_Hover)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.callback = callback
        box = QHBoxLayout(self)
        box.setContentsMargins(12, 10, 12, 10)
        box.setSpacing(10)
        self.icon = QLabel()
        self.icon.setFixedSize(40, 40)
        box.addWidget(self.icon)
        box.addWidget(label(title, 'tileTitle'), 1)
        self.set_selected(False)

    def set_selected(self, selected):
        self.setProperty('selected', 'true' if selected else 'false')
        repolish(self)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.callback()
        super().mousePressEvent(event)

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key.Key_Space, Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self.callback()
        else:
            super().keyPressEvent(event)


class ResponsiveGrid(QWidget):
    """Cards in as many columns as fit, so the examples need no scrolling."""

    def __init__(self, widgets, minimum_width, maximum_columns, viewport=None, parent=None):
        super().__init__(parent)
        # Columns follow the visible width (a scroll area's viewport), never the cards' own minimum.
        self.viewport = viewport
        if viewport is not None:
            viewport.installEventFilter(self)
        self.setObjectName('transparent')
        self.widgets, self.minimum_width, self.maximum_columns = widgets, minimum_width, maximum_columns
        self.grid = QGridLayout(self)
        self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.setSpacing(14)
        self.columns = 0
        self.arrange(maximum_columns)

    def arrange(self, columns):
        if columns == self.columns:
            return
        self.columns = columns
        for widget in self.widgets:
            self.grid.removeWidget(widget)
        for index, widget in enumerate(self.widgets):
            self.grid.addWidget(widget, index // columns, index % columns)
        for column in range(self.maximum_columns):
            self.grid.setColumnStretch(column, 1 if column < columns else 0)

    def fit(self, width):
        self.arrange(max(1, min(self.maximum_columns, width // self.minimum_width)))

    def eventFilter(self, watched, event):
        if watched is self.viewport and event.type() == QEvent.Type.Resize:
            self.fit(event.size().width())
        return False

    def resizeEvent(self, event):
        if self.viewport is None:
            self.fit(event.size().width())
        super().resizeEvent(event)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APPLICATION_TITLE)
        self.resize(1320, 900)
        self.setMinimumSize(980, 700)
        self.settings = QSettings()
        saved = str(self.settings.value('theme', DEFAULT_THEME))
        self.theme = LEGACY_THEMES.get(saved, saved)
        if self.theme not in THEMES:
            self.theme = DEFAULT_THEME
        apply_theme(self.theme)
        self.documents = {}
        self.mode = 'charmprot'
        self.editor = None
        self.process = None
        self.job_temp = None
        self.last_job = None
        self.cancelled = False
        self.close_after_stop = False
        self.active_pid = 0
        self.decoder = codecs.getincrementaldecoder('utf-8')(errors='replace')
        central = QWidget()
        central.setObjectName('workspace')
        self.setCentralWidget(central)
        shell = QHBoxLayout(central)
        shell.setContentsMargins(16, 16, 16, 8)
        shell.setSpacing(18)
        shell.addWidget(self.build_rail())
        body = QVBoxLayout()
        body.setSpacing(12)
        shell.addLayout(body, 1)
        body.addLayout(self.build_top_bar())
        # Messages appear here, inside the window: a modal pop-up can open
        # behind the main window under WSLg and leave it unresponsive.
        self.banner = QFrame()
        self.banner.setObjectName('banner')
        banner_row = QHBoxLayout(self.banner)
        banner_row.setContentsMargins(14, 8, 8, 8)
        self.banner_text = label('')
        self.banner_text.setTextFormat(Qt.TextFormat.PlainText)
        self.banner_text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        banner_row.addWidget(self.banner_text, 1)
        dismiss = button('Dismiss', self.banner.hide)
        dismiss.setObjectName('quiet')
        banner_row.addWidget(dismiss, 0, Qt.AlignmentFlag.AlignTop)
        self.banner.hide()
        body.addWidget(self.banner)
        self.tabs = QStackedWidget()
        self.tabs.setObjectName('pages')
        body.addWidget(self.tabs, 1)
        self.build_configuration_tab()
        self.build_examples_tab()
        self.build_run_tab()
        for number in range(len(PAGES)):
            shortcut = QShortcut(QKeySequence(f'Ctrl+{number + 1}'), self)
            shortcut.activated.connect(lambda n=number: self.show_page(n))
        self.paint_pictograms()
        self.install_document(template_document(self.mode))
        self.select_mode(self.mode)
        self.show_page(0)
        self.update_destination()
        self.statusBar().showMessage('Ready')

    # -- frame --------------------------------------------------------------
    def build_rail(self):
        rail = QFrame()
        rail.setObjectName('rail')
        rail.setFixedWidth(208)
        side = QVBoxLayout(rail)
        side.setContentsMargins(16, 24, 16, 18)
        side.setSpacing(4)
        side.addWidget(label('GMX\nTRANSPLANT', 'brand'))
        side.addWidget(label('MD SYSTEM BUILDER', 'brandSub'))
        side.addSpacing(26)
        self.nav_group = QButtonGroup(self)
        self.nav_group.setExclusive(True)
        self.nav_buttons = []
        for number, name in enumerate(PAGES):
            control = QPushButton(f'{number + 1:02d}   {name}'.replace('&', '&&'))
            control.setObjectName('nav')
            control.setCheckable(True)
            control.setCursor(Qt.CursorShape.PointingHandCursor)
            control.setToolTip(f'{name}  (Ctrl+{number + 1})')
            control.clicked.connect(lambda checked=False, n=number: self.show_page(n))
            self.nav_group.addButton(control, number)
            self.nav_buttons.append(control)
            side.addWidget(control)
        side.addStretch()
        side.addWidget(label('THEME', 'railHeading'))
        side.addSpacing(4)
        self.theme_group = QButtonGroup(self)
        self.theme_group.setExclusive(True)
        self.theme_buttons = {}
        for name, theme in THEMES.items():
            control = QPushButton(name)
            control.setObjectName('swatch')
            control.setCheckable(True)
            control.setCursor(Qt.CursorShape.PointingHandCursor)
            control.setIcon(self.swatch_icon(*theme['swatch']))
            control.setChecked(name == self.theme)
            control.clicked.connect(lambda checked=False, n=name: self.set_theme(n))
            self.theme_group.addButton(control)
            self.theme_buttons[name] = control
            side.addWidget(control)
        side.addSpacing(14)
        self.documentation_button = button('Documentation', self.open_documentation)
        self.documentation_button.setObjectName('railButton')
        self.documentation_button.setToolTip('Open the PDF documentation shipped with GMXtransplant.')
        side.addWidget(self.documentation_button)
        return rail

    def paint_pictograms(self):
        from .icons import pictogram
        theme = THEMES[self.theme]
        cache = {}
        for target, mode in self.pictograms:
            if mode not in cache:
                cache[mode] = pictogram(mode, theme['accent'], theme['muted'], theme['icon_fill'])
            target.setPixmap(cache[mode])

    @staticmethod
    def swatch_icon(fill, accent):
        from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap
        pixmap = QPixmap(36, 36)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QColor('#8a94a0'))
        painter.setBrush(QColor(fill))
        painter.drawEllipse(2, 2, 32, 32)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(accent))
        painter.drawEllipse(11, 11, 14, 14)
        painter.end()
        return QIcon(pixmap)

    def build_top_bar(self):
        bar = QHBoxLayout()
        bar.setSpacing(18)
        titles = QVBoxLayout()
        titles.setSpacing(2)
        titles.addWidget(label('Molecular System Editor', 'hero'))
        titles.addStretch()
        bar.addLayout(titles, 1)
        output = QFrame()
        output.setObjectName('outputBar')
        output.setMinimumWidth(430)
        box = QVBoxLayout(output)
        box.setContentsMargins(14, 8, 8, 8)
        box.setSpacing(0)
        box.addWidget(label('Output folder  ·  required before running', 'muted'))
        row = QHBoxLayout()
        # The output folder starts as the folder gmxtransplant-gui was launched
        # from, like the command-line tool; it is not carried over between sessions.
        self.settings.remove('output_root')
        self.output_root = QLineEdit(str(Path.cwd()))
        self.output_root.setPlaceholderText('Where generated systems and example results are written')
        self.output_root.textChanged.connect(self.update_destination)
        row.addWidget(self.output_root, 1)
        self.output_browse = button('Choose folder…', self.choose_output)
        row.addWidget(self.output_browse)
        box.addLayout(row)
        self.output_hint = label('', 'notice')
        self.output_hint.hide()
        box.addWidget(self.output_hint)
        bar.addWidget(output, 1)
        return bar

    def show_page(self, number):
        self.tabs.setCurrentIndex(number)
        self.nav_buttons[number].setChecked(True)

    # -- pages --------------------------------------------------------------
    def build_configuration_tab(self):
        self.config_tab = QWidget()
        layout = QVBoxLayout(self.config_tab)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)
        self.pictograms = []
        self.mode_buttons = {}
        # The main modes are always shown; the system-edit modes sit under "More modes".
        main = [m for m in MODES if m not in EDIT_MODES]
        more = [m for m in MODES if m in EDIT_MODES]
        layout.addLayout(self.tile_grid(main))
        self.more_modes = QWidget()
        self.more_modes.setObjectName('transparent')
        more_box = QVBoxLayout(self.more_modes)
        more_box.setContentsMargins(0, 0, 0, 0)
        more_box.addLayout(self.tile_grid(more))
        self.more_modes_toggle = disclosure('More modes', self.more_modes)
        layout.addWidget(self.more_modes_toggle, 0, Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(self.more_modes)
        strip = QFrame()
        strip.setObjectName('strip')
        strip_box = QHBoxLayout(strip)
        strip_box.setContentsMargins(14, 9, 14, 9)
        self.mode_info = label('')
        self.mode_info.setTextFormat(Qt.TextFormat.RichText)
        strip_box.addWidget(self.mode_info)
        layout.addWidget(strip)
        self.form_page = QFrame()
        self.form_page.setObjectName('panel')
        self.form_layout = QVBoxLayout(self.form_page)
        self.form_layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.form_page, 1)
        actions = QHBoxLayout()
        actions.addStretch()
        self.schema_button = button('Check configuration', lambda: self.start_current('schema'))
        self.paths_button = button('Check input paths', lambda: self.start_current('paths'))
        self.run_button = button('Run pipeline', lambda: self.start_current('run'), True)
        actions.addWidget(self.schema_button)
        actions.addWidget(self.paths_button)
        actions.addSpacing(6)
        actions.addWidget(self.run_button)
        layout.addLayout(actions)
        self.tabs.addWidget(self.config_tab)

    def tile_grid(self, modes):
        tiles = QGridLayout()
        tiles.setSpacing(10)
        for index, mode in enumerate(modes):
            tile = Tile(MODES[mode], lambda m=mode: self.select_mode(m))
            self.pictograms.append((tile.icon, mode))
            tile.setToolTip(MODE_DESCRIPTIONS[mode])
            self.mode_buttons[mode] = tile
            tiles.addWidget(tile, index // 4, index % 4)
        for column in range(4):
            tiles.setColumnStretch(column, 1)
        return tiles

    def build_examples_tab(self):
        self.examples_tab = QWidget()
        outer = QVBoxLayout(self.examples_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content = QWidget()
        content.setObjectName('transparent')
        scroll.setWidget(content)
        outer.addWidget(scroll)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 4, 0)
        layout.setSpacing(12)
        self.example_buttons = []
        cards, more_cards = [], []
        for key in EXAMPLE_KEYS:
            mode = example_mode(key)
            card = QFrame()
            card.setObjectName('exampleCard')
            card.setAttribute(Qt.WidgetAttribute.WA_Hover)
            box = QVBoxLayout(card)
            box.setContentsMargins(16, 15, 16, 13)
            box.setSpacing(7)
            heading = QHBoxLayout()
            heading.setSpacing(8)
            icon = QLabel()
            icon.setFixedSize(40, 40)
            self.pictograms.append((icon, mode))
            heading.addWidget(icon, 0, Qt.AlignmentFlag.AlignTop)
            heading.addWidget(label(EXAMPLE_TITLES.get(key, MODES[mode]), 'cardTitle'), 1)
            box.addLayout(heading)
            box.addWidget(label(EXAMPLE_DETAILS[key], 'muted'))
            box.addStretch()
            box.addWidget(label(f'Outputs: examples/{example_folder(key)}/', 'path'))
            row = QHBoxLayout()
            row.setSpacing(8)
            open_button = button('Open in editor', lambda checked=False, k=key: self.open_example(k))
            open_button.setObjectName('disclosure')
            row.addWidget(open_button)
            row.addStretch()
            run = button('Run example', lambda checked=False, k=key: self.run_example(k), True)
            self.example_buttons.append(run)
            row.addWidget(run)
            box.addLayout(row)
            (more_cards if mode in EDIT_MODES or key in MORE_EXAMPLE_KEYS else cards).append(card)
        layout.addWidget(ResponsiveGrid(cards, 235, 4, viewport=scroll.viewport()))
        # Examples of the system-edit modes, under "More examples".
        self.more_examples = ResponsiveGrid(more_cards, 235, 4, viewport=scroll.viewport())
        self.more_examples_toggle = disclosure('More examples', self.more_examples)
        layout.addWidget(self.more_examples_toggle, 0, Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(self.more_examples)
        layout.addWidget(label('Rerunning replaces generated files in that example’s folder. Source inputs stay unchanged.', 'muted'))
        layout.addStretch()
        self.tabs.addWidget(self.examples_tab)

    def build_run_tab(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)
        head = QFrame()
        head.setObjectName('panel')
        head_box = QVBoxLayout(head)
        head_box.setContentsMargins(20, 16, 20, 16)
        top = QHBoxLayout()
        words = QVBoxLayout()
        words.setSpacing(2)
        self.run_title = label('Ready when you are', 'runTitle')
        words.addWidget(self.run_title)
        self.run_detail = label('Configuration checks do not process coordinates. Full runs validate the '
                                'assembled system.', 'muted')
        self.run_detail.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        words.addWidget(self.run_detail)
        top.addLayout(words, 1)
        self.cancel_button = button('Cancel run', self.cancel_run)
        self.cancel_button.setObjectName('danger')
        self.cancel_button.setEnabled(False)
        self.cancel_button.setToolTip('Stops the run in progress. Available only while a run is going.')
        top.addWidget(self.cancel_button, 0, Qt.AlignmentFlag.AlignTop)
        head_box.addLayout(top)
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setFixedHeight(6)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        head_box.addSpacing(8)
        head_box.addWidget(self.progress)
        layout.addWidget(head)
        self.progress_model = RunProgress()
        columns = QHBoxLayout()
        columns.setSpacing(12)
        summary_card = QFrame()
        summary_card.setObjectName('panel')
        summary_layout = QVBoxLayout(summary_card)
        summary_layout.setContentsMargins(20, 16, 20, 16)
        summary_layout.addWidget(label('AT A GLANCE', 'overline'))
        self.summary_stage = label('Ready to start', 'cardTitle')
        self.summary_stage.setTextFormat(Qt.TextFormat.PlainText)
        summary_layout.addWidget(self.summary_stage)
        self.summary_facts = label('', 'muted')
        self.summary_facts.setTextFormat(Qt.TextFormat.PlainText)
        self.summary_facts.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        summary_layout.addWidget(self.summary_facts)
        self.summary_notices = label('')
        self.summary_notices.setTextFormat(Qt.TextFormat.PlainText)
        self.summary_notices.setObjectName('notice')
        summary_layout.addWidget(self.summary_notices)
        self.summary_notices.hide()
        summary_layout.addStretch()
        columns.addWidget(summary_card, 3)
        results_card = QFrame()
        results_card.setObjectName('panel')
        results_layout = QVBoxLayout(results_card)
        results_layout.setContentsMargins(20, 16, 20, 16)
        results_layout.addWidget(label('REPORTS  ·  DOUBLE-CLICK TO OPEN', 'overline'))
        self.output_location = label('', 'path')
        self.output_location.setTextFormat(Qt.TextFormat.PlainText)
        self.output_location.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        results_layout.addWidget(self.output_location)
        self.results = QListWidget()
        self.results.setMinimumHeight(95)
        self.results.itemDoubleClicked.connect(self.open_result)
        results_layout.addWidget(self.results, 1)
        row = QHBoxLayout()
        row.addStretch()
        self.viewer_buttons = {}
        self.viewers = {name: find_viewer(name) for name in ('pymol', 'vmd')}
        for viewer, text in [('pymol', 'Open in PyMOL'), ('vmd', 'Open in VMD')]:
            if self.viewers[viewer] is None:
                continue  # Only offer viewers that are installed.
            control = button(text, lambda checked=False, name=viewer: self.open_viewer(name))
            control.setEnabled(False)
            control.setToolTip('Available after a completed run.')
            self.viewer_buttons[viewer] = control
            row.addWidget(control)
        results_layout.addLayout(row)
        columns.addWidget(results_card, 2)
        layout.addLayout(columns, 2)
        self.command_toggle = QPushButton('View command progress')
        self.command_toggle.setObjectName('disclosure')
        self.command_toggle.setCheckable(True)
        self.command_toggle.toggled.connect(self.toggle_command_progress)
        layout.addWidget(self.command_toggle, 0, Qt.AlignmentFlag.AlignLeft)
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setObjectName('console')
        self.console.setMaximumBlockCount(20000)
        layout.addWidget(self.console, 3)
        self.console.hide()
        self.tabs.addWidget(page)

    def show_message(self, text, level='error'):
        self.banner.setProperty('level', level)
        self.banner.style().unpolish(self.banner)
        self.banner.style().polish(self.banner)
        self.banner_text.setText(str(text))
        self.banner.show()

    def error(self, exc):
        self.show_message(exc, 'error')

    def set_theme(self, name):
        if name in THEMES:
            self.theme = name
            apply_theme(name)
            self.theme_buttons[name].setChecked(True)
            self.paint_pictograms()
            for tile in self.mode_buttons.values():
                repolish(tile)
            self.settings.setValue('theme', name)

    def install_document(self, raw):
        if self.editor is not None:
            self.form_layout.removeWidget(self.editor)
            self.editor.deleteLater()
        self.editor = ConfigEditor(raw, self.mode)
        self.form_layout.addWidget(self.editor)

    def current_raw(self):
        return self.editor.value()

    def mark_mode(self):
        for mode, tile in self.mode_buttons.items():
            tile.set_selected(mode == self.mode)

    def select_mode(self, new_mode):
        """Switch the Build page to another mode, keeping unsaved edits of the current one."""
        if new_mode != self.mode:
            try:
                self.documents[self.mode] = self.current_raw()
                raw = self.documents.get(new_mode) or template_document(new_mode)
                self.mode = new_mode
                self.install_document(raw)
            except Exception as exc:
                self.error(exc)
        if new_mode in EDIT_MODES and not self.more_modes_toggle.isChecked():
            self.more_modes_toggle.setChecked(True)
        self.mark_mode()
        self.update_destination()

    def choose_output(self):
        choose_path(self, 'directory', self.output_root.text() or str(Path.home()),
                    lambda path: path and self.output_root.setText(path))

    def update_destination(self):
        selected = bool(self.output_root.text().strip())
        busy = self.process is not None
        # Always clickable; a missing output folder is explained in the banner.
        for btn in [self.run_button, self.schema_button, self.paths_button, *self.example_buttons]:
            btn.setEnabled(not busy)
        if selected and not busy and 'output folder' in self.banner_text.text():
            self.banner.hide()
        self.mode_info.setText(f'<b>Current mode: {MODES[self.mode]}</b> — {MODE_DESCRIPTIONS[self.mode]}')
        # The line under the output folder appears only for a problem with it.
        problem = ''
        if selected:
            try:
                destination(self.output_root.text(), self.mode)
            except ValueError as exc:
                problem = str(exc)
        self.output_hint.setText(problem)
        self.output_hint.setVisible(bool(problem))

    def open_example(self, key):
        try:
            mode = example_mode(key)
            self.documents[self.mode] = self.current_raw()
            self.mode = mode
            # The example's own inputs, written out as the full paths they are.
            self.install_document(example_document(key))
            self.mark_mode()
            self.show_page(0)
            self.banner.hide()
            self.update_destination()
            self.statusBar().showMessage('Example opened for customization. Editor runs use the normal mode output folder.')
        except Exception as exc:
            self.error(exc)

    def start_current(self, action):
        try:
            self.start_job(build_job(self.current_raw(), self.output_root.text(), self.mode), action)
        except Exception as exc:
            self.error(exc)

    def run_example(self, key):
        try:
            path = example_path(key)
            self.start_job(build_job(load_yaml(path), self.output_root.text(), example_mode(key),
                                     example_base=path.parent, example_folder=example_folder(key)), 'run')
        except Exception as exc:
            self.error(exc)

    def open_documentation(self):
        try:
            path = open_documentation()
            self.statusBar().showMessage(f'Opening the documentation: {path}')
        except Exception as exc:
            self.show_message(f'The documentation could not be opened. {exc}')

    def start_job(self, job, action):
        if self.process is not None:
            return
        self.banner.hide()
        self.last_job, self.action = job, action
        self.cancelled = False
        self.active_pid = 0
        self.decoder.reset()
        self.console.clear()
        self.command_toggle.setChecked(False)
        self.progress_model = RunProgress(job['mode'])
        self.progress_model.stage = 'Starting the pipeline' if action == 'run' else 'Checking your inputs'
        self.update_summary()
        self.results.clear()
        self.output_location.setText('')
        for control in self.viewer_buttons.values():
            control.setEnabled(False)
        self.run_title.setText(('Running ' if action == 'run' else 'Checking ') + MODES[job['mode']].lower())
        self.run_detail.setText(job['directory'] if action == 'run' else 'Validation only · previous outputs remain unchanged')
        self.progress.setRange(0, 0)
        self.job_temp = tempfile.TemporaryDirectory(prefix='gmxtransplant-gui-')
        request = Path(self.job_temp.name) / 'job.json'
        request.write_text(json.dumps(job), encoding='utf-8')
        process = QProcess(self)
        self.process = process
        env = QProcessEnvironment.systemEnvironment()
        module_root = str(Path(__file__).resolve().parents[2])
        env.insert('PYTHONPATH', module_root + os.pathsep + env.value('PYTHONPATH'))
        env.insert('PYTHONUNBUFFERED', '1')
        process.setProcessEnvironment(env)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.readyReadStandardOutput.connect(self.read_output)
        process.finished.connect(self.finished)
        process.errorOccurred.connect(self.process_error)
        process.started.connect(lambda: setattr(self, 'active_pid', int(process.processId())))
        self.config_tab.setEnabled(False)
        self.examples_tab.setEnabled(False)
        self.output_root.setEnabled(False)
        self.output_browse.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.update_destination()
        self.show_page(2)
        process.start(sys.executable, ['-u', '-B', '-m', 'gmxtransplant.gui.worker', str(request), action])

    def toggle_command_progress(self, expanded):
        self.console.setVisible(expanded)
        self.command_toggle.setText('Hide command progress' if expanded else 'View command progress')

    def update_summary(self):
        model = self.progress_model
        self.summary_stage.setText(model.stage)
        self.summary_facts.setText('\n'.join(f'{key}: {value}' for key, value in model.facts.items()))
        notices = ([model.error] if model.error else []) + model.notices
        self.summary_notices.setText('\n'.join(notices))
        self.summary_notices.setVisible(bool(notices))

    def append_log(self, text):
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(text)
        self.console.ensureCursorVisible()
        self.progress_model.feed(text)
        self.update_summary()

    def read_output(self):
        if self.process:
            self.append_log(self.decoder.decode(bytes(self.process.readAllStandardOutput())))

    def process_error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            self.append_log('Could not start the pipeline process: ' + self.process.errorString() + '\n')
            self.finished(1, QProcess.ExitStatus.CrashExit)

    def signal_worker(self, sig, pid):
        try:
            # setsid is called immediately by the worker. Never signal our own group.
            if os.getpgid(pid) == pid:
                os.killpg(pid, sig)
            else:
                os.kill(pid, sig)
        except ProcessLookupError:
            pass

    def cancel_run(self):
        if self.process is None:
            return
        self.cancelled = True
        self.run_title.setText('Stopping…')
        self.progress_model.stage = 'Stopping the pipeline safely…'
        self.update_summary()
        self.cancel_button.setEnabled(False)
        pid = self.active_pid or int(self.process.processId())
        if pid:
            self.signal_worker(signal.SIGINT, pid)
            QTimer.singleShot(5000, lambda: self.force_stop(pid))

    def force_stop(self, pid):
        if self.process is not None and self.active_pid == pid:
            self.signal_worker(signal.SIGKILL, pid)

    def finished(self, code, status):
        if self.process is None:
            return
        self.read_output()
        self.append_log(self.decoder.decode(b'', final=True))
        if self.cancelled and self.active_pid > 0:
            # Ensure an external converter cannot continue writing after cancellation.
            try:
                os.killpg(self.active_pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
        succeeded = code == 0 and status == QProcess.ExitStatus.NormalExit and not self.cancelled
        text = 'Cancelled' if self.cancelled else 'Completed' if succeeded and self.action == 'run' else 'Checks passed' if succeeded else 'Failed · review the log'
        self.run_title.setText(text)
        self.progress.setRange(0, 1)
        self.progress.setValue(1 if succeeded else 0)
        self.append_log(f'\n{text} (exit code {code})\n')
        self.progress_model.finish(succeeded, self.cancelled, self.action)
        self.update_summary()
        if self.cancelled and self.action == 'run':
            root = Path(self.last_job['directory'])
            status_file = root / 'run-status.json'
            if status_file.exists():
                saved_status = json.loads(status_file.read_text())
                if saved_status.get('id') == self.last_job.get('id'):
                    saved_status.update(status='cancelled', exit_code=code)
                    status_file.write_text(json.dumps(saved_status))
        self.process.deleteLater()
        self.process = None
        self.job_temp.cleanup()
        self.job_temp = None
        self.config_tab.setEnabled(True)
        self.examples_tab.setEnabled(True)
        self.output_root.setEnabled(True)
        self.output_browse.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.update_destination()
        root = Path(self.last_job['directory'])
        if self.action == 'run' and root.is_dir():
            message = f'All outputs were saved in: {root}'
            self.output_location.setText(message)
            self.append_log(message + '\n')
        for viewer, control in self.viewer_buttons.items():
            control.setEnabled(self.scene_file(viewer) is not None)
        if self.action == 'run' and root.is_dir():
            for name in report_files(root):
                self.results.addItem(name)
        self.statusBar().showMessage(text)
        if self.close_after_stop:
            self.close()

    def scene_file(self, viewer):
        """The viewer's scene from the last run: the assembled system, else the minimization folder."""
        if not self.last_job:
            return None
        name = 'view.pml' if viewer == 'pymol' else 'view.vmd'
        root = Path(self.last_job['directory'])
        for folder in (root, root / 'openmm_minimization'):
            if (folder / name).is_file():
                return folder / name
        return None

    def open_viewer(self, viewer):
        scene, found = self.scene_file(viewer), self.viewers.get(viewer)
        if scene is None or found is None:
            return
        root = scene.parent
        args = [scene] if viewer == 'pymol' else ['-e', scene]
        try:
            # Detached, so closing GMXtransplant leaves the viewer open.
            launch_viewer(found, args, root)
        except OSError as exc:
            self.error(f'Could not start {viewer} ({found.executable}): {exc}')

    def open_result(self, item):
        path = Path(self.last_job['directory']) / item.text()
        try:
            if not open_with_default_app(path):
                self.show_message(f'No application is available to open {path}.')
        except Exception as exc:
            self.error(exc)

    def closeEvent(self, event):
        if self.process is not None:
            self.close_after_stop = True
            self.cancel_run()
            event.ignore()
            return
        event.accept()
