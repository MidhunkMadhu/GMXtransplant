"""Locate the optional PyMOL and VMD viewers (no Qt dependency).

A viewer is looked for, in order: an explicit path in GMXTRANSPLANT_PYMOL /
GMXTRANSPLANT_VMD, the PATH, the folder of the running Python (a conda or pip
install whose environment was not activated), and on macOS the application
bundles in /Applications and ~/Applications, where PyMOL and VMD are normally
installed without any command on the PATH.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
import shutil
import sys
from typing import Dict, List, Optional

# Folders searched for macOS application bundles; a list so tests can redirect it.
MAC_APPLICATION_FOLDERS = [Path('/Applications'), Path.home() / 'Applications']
OVERRIDE = {'pymol': 'GMXTRANSPLANT_PYMOL', 'vmd': 'GMXTRANSPLANT_VMD'}


@dataclass
class Viewer:
    name: str
    executable: str
    env: Dict[str, str] = field(default_factory=dict)

    def command(self, *args) -> List[str]:
        return [self.executable, *map(str, args)]


def _runnable(path) -> bool:
    return path is not None and Path(path).is_file() and os.access(path, os.X_OK)


def _mac_bundle(name: str) -> Optional[Viewer]:
    pattern = 'PyMOL*.app' if name == 'pymol' else 'VMD*.app'
    for folder in MAC_APPLICATION_FOLDERS:
        if not folder.is_dir():
            continue
        # Newest-looking bundle first when several versions are installed.
        for bundle in sorted(folder.glob(pattern), reverse=True):
            if name == 'pymol':
                executable = bundle / 'Contents' / 'MacOS' / 'PyMOL'
                if _runnable(executable):
                    return Viewer(name, str(executable))
            else:
                # VMD's binary needs VMDDIR, which its own launcher normally sets.
                vmddir = bundle / 'Contents' / 'vmd'
                for executable in sorted(vmddir.glob('vmd_MACOSX*')) if vmddir.is_dir() else []:
                    if _runnable(executable):
                        return Viewer(name, str(executable), {'VMDDIR': str(vmddir)})
                launcher = bundle / 'Contents' / 'MacOS' / 'startup.command'
                if _runnable(launcher):
                    return Viewer(name, str(launcher))
    return None


def find_viewer(name: str) -> Optional[Viewer]:
    """The PyMOL ('pymol') or VMD ('vmd') program to run, or None if not installed."""
    override = os.environ.get(OVERRIDE[name], '').strip()
    if override:
        resolved = shutil.which(override) or (override if _runnable(override) else None)
        return Viewer(name, resolved) if resolved else None
    resolved = shutil.which(name)
    if resolved:
        return Viewer(name, resolved)
    sibling = Path(sys.executable).parent / name
    if _runnable(sibling):
        return Viewer(name, str(sibling))
    if sys.platform == 'darwin':
        return _mac_bundle(name)
    return None
