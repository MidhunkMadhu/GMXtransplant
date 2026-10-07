"""Include the curated example inputs in wheels as well as source distributions.

Dependencies are chosen here so that one installation provides the desktop
application, the command-line tool and Open Babel. GMXTRANSPLANT_CLI_ONLY=1 leaves
out Qt for machines where PySide6 cannot be installed, and
GMXTRANSPLANT_NO_OPENBABEL=1 leaves out Open Babel (install it another way).
"""
import os
from pathlib import Path
from setuptools import setup

CORE = ["numpy>=1.23", "scipy>=1.9", "PyYAML>=6.0", "MDAnalysis>=2.4"]
APPLICATION = ["PySide6>=6.6,<7"]
# Open Babel (the `obabel` command) converts experimental cholesterol to
# CHARMM36. Ready-built packages exist for Linux and macOS 14+ (Darwin 23+);
# older macOS would try to compile it from source, so it is skipped there.
OPENBABEL = ['openbabel-wheel>=3.1.1.20; sys_platform != "darwin" or platform_release >= "23"']


def _flag(name):
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes"}


cli_only = _flag("GMXTRANSPLANT_CLI_ONLY")
requirements = CORE + ([] if cli_only else APPLICATION) + ([] if _flag("GMXTRANSPLANT_NO_OPENBABEL") else OPENBABEL)

root = Path(__file__).parent
entries = [('share/gmxtransplant/examples', [
    'configs/charmprot.yaml', 'configs/addbinder.yaml', 'configs/protein_replace.yaml',
    'configs/ligand_replace.yaml', 'configs/cholesterol_restore.yaml', 'configs/minimization.yaml'])]
# The documentation the GUI's Documentation button opens, resolved relative to
# the installed package rather than any location on the user's system.
if (root / 'docs' / 'GMXtransplant.pdf').is_file():
    entries.append(('share/gmxtransplant/docs', ['docs/GMXtransplant.pdf']))
for name, config in [('charmprot', 'charmprot.yaml'),
                     ('addbinder', 'addbinder.yaml'),
                     ('protein_insertion', 'protein_replace.yaml'),
                     ('ligand_replacement', 'ligand_replace.yaml'),
                     ('cholesterol_restoration', 'cholesterol_restore.yaml')]:
    base = Path('examples') / name
    configs = [str(path.relative_to(root)) for path in sorted((root / base).glob('*.yaml'))] or [str(base / config)]
    entries.append((f'share/gmxtransplant/{base}', configs + [str(base / 'README.md')]))
    for role in ['environment', 'replacement', 'experimental', 'reference', 'transplant',
                 'host', 'binder', 'dopamine', 'gprotein']:
        directory = root / base / role
        if not directory.is_dir():
            continue
        groups = {}
        for path in sorted(directory.rglob('*')):
            # The addbinder host ships its CHARMM-GUI run inputs (copied into each
            # pose folder) but not its 30 MB step5 PDB, which addbinder never reads.
            if name == 'addbinder' and role == 'host' and path.suffix.lower() == '.pdb':
                continue
            if path.is_file() and (path.suffix.lower() in {'.pdb', '.gro', '.top', '.itp', '.mol2', '.md', '.mdp'}
                                   or path.name == 'README'):
                relative = path.relative_to(root)
                groups.setdefault(f'share/gmxtransplant/{relative.parent}', []).append(str(relative))
        entries.extend(groups.items())
setup(package_dir={"": "src"}, data_files=entries, install_requires=requirements)
