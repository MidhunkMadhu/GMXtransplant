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

from PySide6.QtCore import Qt, QProcess, QProcessEnvironment, QSettings, QTimer
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QLineEdit, QComboBox, QTabWidget, QPlainTextEdit,
    QFrame, QProgressBar, QListWidget, QScrollArea, QSizePolicy)

from viewers import find_viewer
from .editor import ConfigEditor, choose_path
from .model import (MODES, FOLDERS, EXAMPLE_KEYS, build_job, load_yaml, example_mode, example_folder,
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
    'minimize': 'Prepare a restrained OpenMM minimization folder for a system you have already assembled.',
}

THEMES = {
    'Teal': dict(bg='#f4f7f8', surface='#ffffff', text='#203443', muted='#536d7a', border='#dce6ea',
                 input_border='#bdcdd4', accent='#087f72', accent_hover='#066b60', on_accent='#ffffff',
                 accent_soft='#d8eee8', accent_soft_text='#075f56', sidebar='#142e3b', sidebar_text='#bacdd4',
                 brand='#f5fbfc', badge='#7de0c5', hero='#173b48', button='#eaf1f4', button_border='#c7d7de',
                 button_hover='#dceaef', button_pressed='#c8dfe3', disabled_bg='#edf1f3', disabled_text='#83959e',
                 console_bg='#142e3b', console_text='#dbedea', list_bg='#f5f8fa', tab='#e7eef1',
                 notice='#975522', danger='#a43c40', danger_bg='#fff1f0', danger_border='#e4bab9',
                 info_bg='#e6f4f1', info_border='#a9d8cd', scroll='#b8cbd3', dark=False),
    'Midnight': dict(bg='#11161c', surface='#1a212a', text='#dbe4ea', muted='#8fa1ad', border='#2a3440',
                     input_border='#3a4756', accent='#3fb8a4', accent_hover='#5cc9b7', on_accent='#0b1f1c',
                     accent_soft='#1f3a37', accent_soft_text='#8fe3d3', sidebar='#0b0f14', sidebar_text='#93a4b0',
                     brand='#eef6f8', badge='#5fd3b8', hero='#eef3f6', button='#232c37', button_border='#354251',
                     button_hover='#2c3845', button_pressed='#34424f', disabled_bg='#1c232b', disabled_text='#5f6f7b',
                     console_bg='#0b0f14', console_text='#cfe3df', list_bg='#161c23', tab='#1f2731',
                     notice='#e0a458', danger='#f08a8d', danger_bg='#3a1f22', danger_border='#6b3437',
                     info_bg='#16302c', info_border='#2f5f57', scroll='#3a4756', dark=True),
    'Sand': dict(bg='#f7f3ec', surface='#fffdf9', text='#3b2f24', muted='#7a6a58', border='#e6dccd',
                 input_border='#d3c4ae', accent='#b5562b', accent_hover='#9a4722', on_accent='#ffffff',
                 accent_soft='#f3dfd1', accent_soft_text='#83391a', sidebar='#3a2a1f', sidebar_text='#d9c8b4',
                 brand='#fff6ea', badge='#f2b37c', hero='#4a3222', button='#f1e9dd', button_border='#dccdb8',
                 button_hover='#e8dccb', button_pressed='#dfcfb9', disabled_bg='#f1ece4', disabled_text='#a8998a',
                 console_bg='#3a2a1f', console_text='#f3e6d6', list_bg='#faf6f0', tab='#efe7da',
                 notice='#9a5a12', danger='#a43c40', danger_bg='#fbeceb', danger_border='#e4bab9',
                 info_bg='#f6ece0', info_border='#e2c9ab', scroll='#d3c4ae', dark=False),
}
DEFAULT_THEME = 'Teal'

STYLE = """
QWidget {{ font-family: "Inter", "Segoe UI", "DejaVu Sans", sans-serif; font-size: 13px; color: {text}; }}
QMainWindow, QDialog, QWidget#workspace {{ background: {bg}; }}
QWidget#page {{ background: {surface}; }}
QFrame#sidebar {{ background: {sidebar}; border-radius: 16px; }}
QFrame#sidebar QLabel {{ background: transparent; }}
QLabel#brand {{ color: {brand}; font-size: 23px; font-weight: 700; }}
QLabel#sidebarText {{ color: {sidebar_text}; font-size: 13px; }}
QLabel#badge {{ color: {badge}; font-size: 11px; font-weight: 700; }}
QLabel#hero {{ font-size: 27px; font-weight: 700; color: {hero}; }}
QLabel#sectionTitle {{ font-size: 19px; font-weight: 650; margin: 6px 0; }}
QLabel#fieldLabel {{ font-weight: 600; }}
QCheckBox#fieldToggle {{ font-weight: 600; }}
QLabel#muted {{ color: {muted}; }}
QLabel#notice {{ color: {notice}; }}
QFrame#card {{ background: {surface}; border: 1px solid {border}; border-radius: 10px; }}
QFrame#fieldRow {{ background: transparent; border: none; }}
QFrame#banner {{ border-radius: 8px; padding: 2px; }}
QFrame#banner[level="error"] {{ background: {danger_bg}; border: 1px solid {danger_border}; }}
QFrame#banner[level="error"] QLabel {{ color: {danger}; }}
QFrame#banner[level="info"] {{ background: {info_bg}; border: 1px solid {info_border}; }}
QLineEdit, QComboBox, QPlainTextEdit {{ background: {surface}; border: 1px solid {input_border}; border-radius: 6px; padding: 7px; selection-background-color: {accent_soft}; selection-color: {text}; }}
QLineEdit:focus, QComboBox:focus, QPlainTextEdit:focus {{ border: 1px solid {accent}; }}
QComboBox::drop-down {{ width: 26px; border: none; }}
QComboBox::down-arrow {{ image: url("{arrow_icon}"); width: 10px; height: 10px; }}
QComboBox QAbstractItemView {{ background: {surface}; color: {text}; selection-background-color: {accent_soft}; selection-color: {text}; }}
QPushButton {{ background: {button}; border: 1px solid {button_border}; border-radius: 7px; padding: 8px 13px; font-weight: 600; }}
QPushButton:hover {{ background: {button_hover}; }}
QPushButton:pressed {{ background: {button_pressed}; }}
QPushButton#primary {{ background: {accent}; color: {on_accent}; border-color: {accent}; }}
QPushButton#primary:hover {{ background: {accent_hover}; }}
QPushButton#danger {{ color: {danger}; background: {danger_bg}; border-color: {danger_border}; }}
QPushButton#quiet, QPushButton#disclosure {{ background: transparent; color: {muted}; border: none; padding: 4px 8px; }}
QPushButton#disclosure {{ color: {accent}; font-weight: 700; padding-left: 0; }}
QPushButton#quiet:hover, QPushButton#disclosure:hover {{ color: {accent}; }}
QPushButton:disabled, QLineEdit:disabled, QComboBox:disabled {{ color: {disabled_text}; background: {disabled_bg}; border-color: {border}; }}
QPushButton#primary:disabled, QPushButton#danger:disabled {{ color: {disabled_text}; background: {disabled_bg}; border-color: {border}; }}
QCheckBox {{ spacing: 8px; padding: 3px; background: transparent; }}
QCheckBox::indicator {{ width: 16px; height: 16px; border: 1px solid {input_border}; border-radius: 4px; background: {surface}; }}
QCheckBox::indicator:hover {{ border-color: {accent}; }}
QCheckBox::indicator:checked {{ background: {accent}; border-color: {accent}; image: url("{check_icon}"); }}
QCheckBox::indicator:disabled {{ background: {disabled_bg}; border-color: {border}; }}
QTabWidget::pane {{ border: 1px solid {border}; background: {surface}; border-radius: 10px; }}
QTabBar::tab {{ background: {tab}; padding: 12px 22px; margin-right: 5px; border-top-left-radius: 8px; border-top-right-radius: 8px; font-weight: 600; }}
QTabBar::tab:selected {{ background: {surface}; color: {accent}; border-bottom: 3px solid {accent}; }}
QListWidget {{ background: {list_bg}; border: none; border-radius: 8px; padding: 5px; outline: none; }}
QListWidget::item {{ padding: 10px 8px; border-radius: 6px; margin-bottom: 3px; }}
QListWidget::item:selected {{ background: {accent_soft}; color: {accent_soft_text}; }}
QScrollArea {{ border: none; background: {surface}; }}
QScrollBar:vertical {{ background: {list_bg}; width: 10px; }}
QScrollBar::handle:vertical {{ background: {scroll}; border-radius: 5px; min-height: 30px; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0px; }}
QProgressBar {{ background: {border}; border: none; border-radius: 4px; height: 7px; }}
QProgressBar::chunk {{ background: {accent}; border-radius: 4px; }}
QPlainTextEdit#console {{ background: {console_bg}; color: {console_text}; border: none; font-family: "DejaVu Sans Mono", monospace; font-size: 12px; padding: 14px; }}
QToolTip {{ background: {sidebar}; color: {brand}; padding: 6px; border: none; }}
QFileDialog QListView, QFileDialog QTreeView {{ background: {surface}; color: {text}; }}
"""


def theme_icons(name, theme):
    """Check-mark and drop-down arrow drawn in the theme's colours."""
    folder = Path(tempfile.gettempdir()) / f'gmxtransplant-{os.getuid() if hasattr(os, "getuid") else 0}' / name
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
    theme = THEMES.get(name, THEMES[DEFAULT_THEME])
    try:
        icons = theme_icons(name if name in THEMES else DEFAULT_THEME, theme)
    except OSError:
        icons = {'check_icon': '', 'arrow_icon': ''}  # Plain indicators still show state.
    app = QApplication.instance()
    palette = QPalette()
    for role, key in [(QPalette.ColorRole.Window, 'bg'), (QPalette.ColorRole.WindowText, 'text'),
                      (QPalette.ColorRole.Base, 'surface'), (QPalette.ColorRole.AlternateBase, 'list_bg'),
                      (QPalette.ColorRole.Text, 'text'), (QPalette.ColorRole.Button, 'button'),
                      (QPalette.ColorRole.ButtonText, 'text'), (QPalette.ColorRole.Highlight, 'accent'),
                      (QPalette.ColorRole.HighlightedText, 'on_accent'), (QPalette.ColorRole.ToolTipBase, 'sidebar'),
                      (QPalette.ColorRole.ToolTipText, 'brand'), (QPalette.ColorRole.PlaceholderText, 'muted')]:
        palette.setColor(role, QColor(theme[key]))
    app.setPalette(palette)
    app.setStyleSheet(STYLE.format(**theme, **icons))


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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APPLICATION_TITLE)
        self.resize(1320, 900)
        self.setMinimumSize(960, 680)
        self.settings = QSettings()
        self.theme = str(self.settings.value('theme', DEFAULT_THEME))
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
        shell.setContentsMargins(18, 18, 18, 18)
        shell.setSpacing(22)
        sidebar = QFrame()
        sidebar.setObjectName('sidebar')
        sidebar.setFixedWidth(210)
        side = QVBoxLayout(sidebar)
        side.setContentsMargins(21, 27, 21, 22)
        side.addWidget(label('GMX\nTRANSPLANT', 'brand'))
        side.addWidget(label('MD SYSTEM BUILDER', 'badge'))
        side.addSpacing(38)
        for heading, detail in [('01  Choose your inputs', 'Browse coordinates and matching parameters.'),
                                ('02  Shape your system', 'Edit selections and scientific settings.'),
                                ('03  Review & run', 'Check configuration, run, inspect results.')]:
            side.addWidget(label(heading, 'badge'))
            side.addWidget(label(detail, 'sidebarText'))
            side.addSpacing(18)
        side.addStretch()
        self.documentation_button = button('Documentation', self.open_documentation)
        self.documentation_button.setToolTip('Open the PDF documentation shipped with GMXtransplant.')
        side.addWidget(self.documentation_button)
        shell.addWidget(sidebar)
        body = QVBoxLayout()
        body.setSpacing(12)
        shell.addLayout(body, 1)
        header = QHBoxLayout()
        titles = QVBoxLayout()
        titles.addWidget(label('Molecular System Editor', 'hero'))
        header.addLayout(titles, 1)
        theme_box = QVBoxLayout()
        theme_box.addWidget(label('Theme', 'muted'))
        self.theme_combo = QComboBox()
        self.theme_combo.addItems(list(THEMES))
        self.theme_combo.setCurrentText(self.theme)
        self.theme_combo.setToolTip('Colour theme for every window')
        self.theme_combo.currentTextChanged.connect(self.set_theme)
        theme_box.addWidget(self.theme_combo)
        theme_box.addStretch()
        header.addLayout(theme_box)
        body.addLayout(header)
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
        output_card = QFrame()
        output_card.setObjectName('card')
        out_box = QVBoxLayout(output_card)
        out_box.addWidget(label('Output folder  ·  required before running'))
        out_row = QHBoxLayout()
        # The output folder starts as the folder gmxtransplant-gui was launched
        # from, like the command-line tool; it is not carried over between sessions.
        self.settings.remove('output_root')
        self.output_root = QLineEdit(str(Path.cwd()))
        self.output_root.setPlaceholderText('Choose where generated systems and example results will be written')
        out_row.addWidget(self.output_root, 1)
        self.output_browse = button('Choose folder…', self.choose_output)
        out_row.addWidget(self.output_browse)
        out_box.addLayout(out_row)
        self.output_hint = label('', 'muted')
        self.output_hint.hide()
        out_box.addWidget(self.output_hint)
        body.addWidget(output_card)
        self.output_root.textChanged.connect(self.update_destination)
        self.tabs = QTabWidget()
        body.addWidget(self.tabs, 1)
        self.build_configuration_tab()
        self.build_examples_tab()
        self.build_run_tab()
        self.install_document(template_document(self.mode))
        self.update_destination()
        self.statusBar().showMessage('Ready')

    def build_configuration_tab(self):
        self.config_tab = QWidget()
        layout = QVBoxLayout(self.config_tab)
        row = QHBoxLayout()
        self.mode_combo = QComboBox()
        for key, name in MODES.items():
            self.mode_combo.addItem(name, key)
        self.mode_combo.currentIndexChanged.connect(self.change_mode)
        row.addWidget(self.mode_combo, 1)
        layout.addLayout(row)
        self.mode_info = label('')
        self.mode_info.setTextFormat(Qt.TextFormat.RichText)
        layout.addSpacing(10)
        layout.addWidget(self.mode_info)
        layout.addWidget(label('Input editor', 'sectionTitle'))
        layout.addWidget(label('Every input is a complete path. Browse for each one; nothing is '
                               'resolved against a shared folder.', 'muted'))
        self.form_page = QWidget()
        self.form_layout = QVBoxLayout(self.form_page)
        self.form_layout.setContentsMargins(0, 4, 0, 0)
        layout.addWidget(self.form_page, 1)
        actions = QHBoxLayout()
        self.schema_button = button('Check configuration', lambda: self.start_current('schema'))
        self.paths_button = button('Check input paths', lambda: self.start_current('paths'))
        self.run_button = button('Run pipeline', lambda: self.start_current('run'), True)
        actions.addWidget(self.schema_button)
        actions.addWidget(self.paths_button)
        actions.addStretch()
        actions.addWidget(self.run_button)
        layout.addLayout(actions)
        self.tabs.addTab(self.config_tab, 'Configuration')

    def build_examples_tab(self):
        self.examples_tab = QWidget()
        outer = QVBoxLayout(self.examples_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        scroll.setWidget(content)
        outer.addWidget(scroll)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(24, 20, 24, 20)
        layout.addWidget(label('Start with a complete example', 'sectionTitle'))
        layout.addWidget(label('Ready-to-run inputs are included. Select your output folder above, then run an example in one click.', 'muted'))
        descriptions = {
            'charmprot': 'Choose reference and transplant folders. Replace proteins and ligands automatically while retaining the reference environment.',
            'protein': 'Insert a protein and its bound ligands into a prepared membrane system using matching topology inputs.',
            'lig': 'Replace one bound ligand (LI1) with another (LI2) using explicit heavy-atom correspondence, its ITP and its force field.',
            'chl': 'Restore five experimental cholesterol molecules into the prepared D1R system. Uses Open Babel, installed with GMXtransplant.',
            'addbinder': 'Extracellular ligand: place dopamine above the apo D1 receptor, once straight above at 20 Å and three times from random directions at 20–25 Å; one run-ready folder each.',
            'addbinder_gprotein': 'Intracellular protein: place the Gs heterotrimer below the apo D1 receptor in its receptor-bound arrangement: pulled straight out to 20 Å (shorter if the box needs it) and towards two random sides at 15 Å, one of them also randomly rotated; one run-ready folder each.',
        }
        titles = {'addbinder': 'Add binder: dopamine (extracellular)',
                  'addbinder_gprotein': 'Add binder: G protein (intracellular)'}
        self.example_buttons = []
        for key in EXAMPLE_KEYS:
            mode = example_mode(key)
            card = QFrame()
            card.setObjectName('card')
            box = QVBoxLayout(card)
            box.setContentsMargins(20, 15, 20, 15)
            box.addWidget(label(titles.get(key, MODES[mode]), 'sectionTitle'))
            box.addWidget(label(descriptions[key], 'muted'))
            box.addWidget(label(f'Outputs: examples/{example_folder(key)}/', 'muted'))
            row = QHBoxLayout()
            row.addStretch()
            row.addWidget(button('Open in editor', lambda checked=False, k=key: self.open_example(k)))
            run = button('Run example', lambda checked=False, k=key: self.run_example(k), True)
            self.example_buttons.append(run)
            row.addWidget(run)
            box.addLayout(row)
            layout.addWidget(card)
        layout.addStretch()
        layout.addWidget(label('Rerunning replaces generated files in that example’s folder. Source inputs stay unchanged.', 'muted'))
        self.tabs.addTab(self.examples_tab, 'Examples')

    def build_run_tab(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        self.run_title = label('Ready when you are', 'sectionTitle')
        layout.addWidget(self.run_title)
        self.run_detail = label('Configuration checks do not process coordinates. Full runs validate the assembled system.', 'muted')
        layout.addWidget(self.run_detail)
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setFixedHeight(8)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        layout.addWidget(self.progress)
        self.progress_model = RunProgress()
        summary_card = QFrame()
        summary_card.setObjectName('card')
        summary_card.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        summary_layout = QVBoxLayout(summary_card)
        summary_layout.setContentsMargins(20, 16, 20, 16)
        summary_layout.addWidget(label('Progress at a glance', 'sectionTitle'))
        self.summary_stage = label('Ready to start')
        self.summary_stage.setTextFormat(Qt.TextFormat.PlainText)
        summary_layout.addWidget(self.summary_stage)
        self.summary_facts = label('', 'muted')
        self.summary_facts.setTextFormat(Qt.TextFormat.PlainText)
        summary_layout.addWidget(self.summary_facts)
        self.summary_notices = label('')
        self.summary_notices.setTextFormat(Qt.TextFormat.PlainText)
        self.summary_notices.setObjectName('notice')
        summary_layout.addWidget(self.summary_notices)
        self.summary_notices.hide()
        layout.addWidget(summary_card)
        self.command_toggle = QPushButton('View Command Progress')
        self.command_toggle.setCheckable(True)
        self.command_toggle.toggled.connect(self.toggle_command_progress)
        layout.addWidget(self.command_toggle, 0, Qt.AlignmentFlag.AlignLeft)
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setObjectName('console')
        self.console.setMaximumBlockCount(20000)
        layout.addWidget(self.console, 3)
        self.console.hide()
        self.output_location = label('', 'muted')
        self.output_location.setTextFormat(Qt.TextFormat.PlainText)
        self.output_location.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.output_location)
        layout.addWidget(label('Reports · double-click to open', 'muted'))
        self.results = QListWidget()
        self.results.setMinimumHeight(95)
        self.results.itemDoubleClicked.connect(self.open_result)
        layout.addWidget(self.results, 1)
        row = QHBoxLayout()
        self.cancel_button = button('Cancel run', self.cancel_run)
        self.cancel_button.setObjectName('danger')
        self.cancel_button.setEnabled(False)
        self.cancel_button.setToolTip('Stops the run in progress. Available only while a run is going.')
        row.addWidget(self.cancel_button)
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
        layout.addLayout(row)
        self.tabs.addTab(page, 'Run && results')

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
            self.settings.setValue('theme', name)

    def install_document(self, raw):
        if self.editor is not None:
            self.form_layout.removeWidget(self.editor)
            self.editor.deleteLater()
        self.editor = ConfigEditor(raw, self.mode)
        self.form_layout.addWidget(self.editor)

    def current_raw(self):
        return self.editor.value()

    def change_mode(self):
        new_mode = self.mode_combo.currentData()
        if new_mode == self.mode:
            return
        try:
            self.documents[self.mode] = self.current_raw()
            raw = self.documents.get(new_mode) or template_document(new_mode)
            self.mode = new_mode
            self.install_document(raw)
            self.update_destination()
        except Exception as exc:
            self.mode_combo.blockSignals(True)
            self.mode_combo.setCurrentIndex(self.mode_combo.findData(self.mode))
            self.mode_combo.blockSignals(False)
            self.error(exc)

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
            self.mode_combo.blockSignals(True)
            self.mode_combo.setCurrentIndex(self.mode_combo.findData(mode))
            self.mode_combo.blockSignals(False)
            # The example's own inputs, written out as the full paths they are.
            self.install_document(example_document(key))
            self.tabs.setCurrentIndex(0)
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
        self.tabs.setCurrentIndex(2)
        process.start(sys.executable, ['-u', '-B', '-m', 'gmxtransplant.gui.worker', str(request), action])

    def toggle_command_progress(self, expanded):
        self.console.setVisible(expanded)
        self.command_toggle.setText('Hide Command Progress' if expanded else 'View Command Progress')

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
            subprocess.Popen(found.command(*args), cwd=str(root), env={**os.environ, **found.env},
                             stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL, start_new_session=True)
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
