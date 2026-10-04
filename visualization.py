"""Portable comparison scenes. No viewer is required for the scientific workflow."""
from collections import Counter
from pathlib import Path
import json
import os
import shutil
import subprocess
import tempfile
import warnings
import numpy as np

# Shared palette across every mode; source overlays are initially hidden.
PALETTE = {
    'inserted_protein': ('Incoming protein', [0.15, .48, .90], True),
    'retained_protein': ('Retained protein', [.48, .55, .65], True),
    'inserted_ligand': ('Incoming ligand', [.84, .22, .67], True),
    'retained_ligand': ('Retained ligand', [.62, .44, .72], True),
    'inserted_cholesterol': ('Incoming cholesterol', [.95, .68, .12], True),
    'retained_cholesterol': ('Retained cholesterol', [.10, .70, .59], True),
    'retained_lipid': ('Retained membrane', [.69, .76, .80], True),
    'retained_water': ('Water', [.55, .75, .94], False),
    'retained_ion': ('Ions', [.36, .65, .85], False),
    'retained_other': ('Other retained molecules', [.60, .65, .70], True),
    'removed_protein': ('Original protein replaced', [.97, .50, .16], False),
    'removed_ligand': ('Original ligand replaced', [.97, .50, .16], True),
    'removed_environment': ('Removed environment molecules', [.91, .24, .24], False),
    'reference_input': ('Reference input in final frame', [.97, .60, .25], False),
    'transplant_input': ('Transplant input aligned to final frame', [.22, .77, .87], False),
}


def atom_keys(atoms):
    # Assembly only reorders these atoms: exact float coordinates, no fuzzy match.
    return [(str(a.name), str(a.resname), *map(float, a.position)) for a in atoms]


CHOLESTEROL_RESNAMES = {'CHL', 'CHL1', 'CHOL', 'CLR'}


def categories(atoms, cholesterol=False):
    """Classify atoms for coloring.

    Cholesterol is a class of its own only in cholesterol-restoration mode,
    which is the one mode that reasons about it. Every other mode renders
    cholesterol exactly like any other membrane lipid.
    """
    from classify import classify_resnames
    classes = classify_resnames(atoms.resnames)
    return ['cholesterol' if cholesterol and str(a.resname).upper() in CHOLESTEROL_RESNAMES
            else classes.get(a.resname, 'other') for a in atoms]


def comparison_groups(final, reference=None, inserted=None, replaced=None, transplant=None,
                      cholesterol=False):
    """Color final atoms by their insertion provenance; retain aligned input overlays."""
    final_atoms = final.atoms
    incoming = Counter(atom_keys(inserted)) if inserted is not None else Counter()
    labels = []
    for key, category in zip(atom_keys(final_atoms), categories(final_atoms, cholesterol)):
        is_new = incoming[key] > 0
        if is_new:
            incoming[key] -= 1
        if is_new and category == 'other':
            category = 'ligand'
        role = ('inserted_' if is_new else 'retained_') + category
        if role not in PALETTE:
            role = 'retained_other'
        labels.append(role)
    groups = {role: final_atoms[np.array(labels) == role] for role in dict.fromkeys(labels)}
    if replaced is not None:
        protein = np.array(categories(replaced, cholesterol)) == 'protein'
        groups['removed_protein'] = replaced[protein]
        groups['removed_ligand'] = replaced[~protein]
    if reference is not None:
        ref = reference.atoms
        groups['reference_input'] = ref
        retained = Counter(atom_keys(final_atoms))
        replaced_keys = Counter(atom_keys(replaced)) if replaced is not None else Counter()
        removed = []
        for index, key in enumerate(atom_keys(ref)):
            if retained[key]:
                retained[key] -= 1
            elif replaced_keys[key]:
                replaced_keys[key] -= 1
            else:
                removed.append(index)
        groups['removed_environment'] = ref[removed]
    if transplant is not None:
        groups['transplant_input'] = transplant.atoms
    return groups


# Proteins and peptides of at least this many residues are drawn as cartoon;
# shorter peptides, ligands and other small molecules as licorice.
MIN_CARTOON_RESIDUES = 3
# Lines for the membrane, removed molecules and whole-system overlays, so they do not crowd the view.
LINE_ROLES = {'retained_lipid', 'retained_cholesterol', 'removed_environment'}


def _style(role, atoms):
    """How a scene object is drawn, from what it contains.

    cartoon: its protein/peptide part (>= MIN_CARTOON_RESIDUES residues) as cartoon.
    rest: 'licorice' for ligands and small molecules, 'lines' for the membrane
    and overlays, None when nothing else is left. spheres: 'tip' or 'ion'.
    """
    from classify import classify_resnames
    if role.endswith('_tip'):
        return {'cartoon': False, 'rest': None, 'spheres': 'tip'}
    classes = classify_resnames(set(atoms.resnames))
    kinds = [classes.get(r.resname, 'other') for r in atoms.residues]
    protein = sum(k == 'protein' for k in kinds)
    cartoon = protein >= MIN_CARTOON_RESIDUES
    if kinds and all(k == 'ion' for k in kinds):
        return {'cartoon': False, 'rest': None, 'spheres': 'ion'}
    if kinds and all(k == 'water' for k in kinds):
        return {'cartoon': False, 'rest': 'lines', 'spheres': None, 'water': True}
    rest = [k for k in kinds if k not in ('water', 'ion') and not (cartoon and k == 'protein')]
    overlay = role.endswith('_input')
    lines = overlay or role in LINE_ROLES or (kinds and all(k == 'lipid' for k in kinds))
    return {'cartoon': cartoon, 'rest': ('lines' if lines else 'licorice') if rest else None, 'spheres': None}


def _write_object(atoms, path, style):
    """PDB of one scene object; protein drawn as cartoon is ATOM, everything else HETATM."""
    import MDAnalysis as mda
    from classify import classify_resnames
    copy = mda.Merge(atoms)
    classes = classify_resnames(set(copy.atoms.resnames))
    copy.add_TopologyAttr('record_types', ['ATOM' if style['cartoon'] and classes.get(r) == 'protein'
                                           else 'HETATM' for r in copy.atoms.resnames])
    if atoms.universe.dimensions is not None:
        copy.dimensions = atoms.universe.dimensions
    copy.atoms.write(str(path))


def write_scene(directory, groups, mode, alignment='Inputs use the assembly alignment in the final coordinate frame.',
                palette=None):
    """Write viewer scripts and coordinates; make a real PSE when PyMOL is installed.

    `palette` adds mode-specific objects (name -> (label, rgb, visible)) to PALETTE.
    """
    palette = {**PALETTE, **(palette or {})}
    directory = Path(directory).absolute()
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.view-', dir=directory) as tmp:
        staged = Path(tmp)
        assets = staged / 'visualization'
        assets.mkdir()
        objects = []
        for role, atoms in groups.items():
            if role not in palette or len(atoms) == 0:
                continue
            label, rgb, visible = palette[role]
            style = _style(role, atoms)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _write_object(atoms, assets / (role + '.pdb'), style)
            objects.append(dict(name=role, label=label, color=rgb, visible=visible, atoms=len(atoms), style=style))
        manifest = dict(mode=mode, alignment=alignment, objects=objects,
                        note='Source overlays and removed molecules are comparisons, not part of the final system. Water and ions are initially hidden.')
        (assets / 'scene.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (assets / 'LEGEND.txt').write_text(alignment + '\n\n' + manifest['note'] + '\n\n' + '\n'.join(
            f"{o['name']}: {o['label']}; {o['atoms']} atoms; {'shown' if o['visible'] else 'hidden'}; RGB {o['color']}" for o in objects))
        # This Python loader can also be run directly by PyMOL after moving the folder.
        (assets / 'load.py').write_text('''from pathlib import Path
import json
from pymol import cmd
base = Path(__file__).resolve().parent
scene = json.loads((base / 'scene.json').read_text())
cmd.bg_color('white')
for obj in scene['objects']:
    name = obj['name']
    cmd.load(str(base / (name + '.pdb')), name)
    cmd.hide('everything', name)
    cmd.set_color('gmx_' + name, obj['color'])
    cmd.color('gmx_' + name, name)
    style = obj['style']
    # Proteins and peptides as cartoon, ligands and small molecules as licorice.
    if style['cartoon']:
        cmd.show('cartoon', name + ' and not hetatm')
    if style['rest']:
        rest = (name + ' and hetatm' if style['cartoon'] else name) + ('' if style.get('water') else ' and not solvent')
        if style['rest'] == 'licorice':
            cmd.show('sticks', rest)
            cmd.set('stick_radius', 0.25, name)
        else:
            cmd.show('lines', rest)
    if style['spheres']:
        cmd.show('spheres', name)
        cmd.set('sphere_scale', 0.6 if style['spheres'] == 'tip' else 0.35, name)
    if name.startswith('removed_') or name.endswith('_input'):
        cmd.set('stick_transparency', 0.45, name)
        cmd.set('cartoon_transparency', 0.55, name)
    if not obj['visible']:
        cmd.disable(name)
cmd.set('stick_radius', 0.14)
cmd.orient('enabled')
cmd.save(str(base.parent / 'view.pse'))
''')
        (staged / 'view.pml').write_text('python\nfrom pathlib import Path\nfrom pymol import cmd\n'
            + f"_gmx_root = Path({str(directory)!r})\n"
            + "if (Path.cwd() / 'visualization/scene.json').is_file():\n    _gmx_root = Path.cwd()\n"
            + "cmd.run(str(_gmx_root / 'visualization/load.py'), namespace='module')\npython end\n")
        vmd = ['# GMXtransplant comparison scene; source overlays are initially hidden.',
               'set gmx_root [file dirname [file normalize [info script]]]', 'display backgroundgradient off',
               'color Display Background white', 'display projection Orthographic']
        for index, obj in enumerate(objects):
            name = obj['name']; rgb = ' '.join(map(str, obj['color'])); color = 17 + index
            vmd += [f'color change rgb {color} {rgb}',
                    f'set gmx_mol [mol new [file join $gmx_root visualization {name}.pdb] type pdb waitfor all]',
                    f'mol rename $gmx_mol {name}', 'mol delrep 0 $gmx_mol',
                    f'mol color ColorID {color}']
            style = obj['style']
            # Proteins and peptides as cartoon, ligands and small molecules as licorice.
            if style['cartoon']:
                vmd += ['mol representation NewCartoon', 'mol selection protein', 'mol addrep $gmx_mol']
            if style['rest']:
                selection = ('all' if style.get('water') else
                             '{not protein and not water}' if style['cartoon'] else '{not water}')
                representation = 'Licorice 0.3 12 12' if style['rest'] == 'licorice' else 'Lines 1.0'
                vmd += [f'mol representation {representation}', f'mol selection {selection}',
                        'mol addrep $gmx_mol']
            if style['spheres']:
                vmd += [f"mol representation VDW {'0.8' if style['spheres'] == 'tip' else '0.5'} 12",
                        'mol selection all', 'mol addrep $gmx_mol']
            if not obj['visible']:
                vmd += ['mol off $gmx_mol']
        vmd += ['display resetview']
        (staged / 'view.vmd').write_text('\n'.join(vmd) + '\n')
        # Generate only owned artifacts. Never replace the parent output directory.
        for name in ['visualization', 'view.pml', 'view.vmd']:
            target = directory / name
            if target.is_symlink():
                raise ValueError(f'Visualization output is a symlink: {target}')
            if target.is_dir():
                shutil.rmtree(target)
            (staged / name).replace(target)
    session = directory / 'view.pse'
    if session.exists() or session.is_symlink():
        session.unlink()
    from viewers import find_viewer
    viewer = find_viewer('pymol')
    if viewer:
        try:
            result = subprocess.run(viewer.command('-cq', directory / 'view.pml'), env={**os.environ, **viewer.env},
                                    capture_output=True, text=True, timeout=300, cwd=directory)
            if result.returncode or not session.is_file():
                print('[VISUALIZATION] PyMOL ran but wrote no session; open view.pml in PyMOL to save view.pse.')
        except (OSError, subprocess.TimeoutExpired):
            print('[VISUALIZATION] PyMOL session could not be built; open view.pml in PyMOL to save view.pse.')
    else:
        print('[VISUALIZATION] PyMOL was not found (PATH, this Python environment, or on macOS '
              '/Applications); open view.pml in PyMOL to save view.pse, or set GMXTRANSPLANT_PYMOL.')
    # Both session files stay in the output folder so the assembled system can be
    # reopened later without rerunning the pipeline.
    print(f'[VISUALIZATION] Session files written to {directory}')
    for name in ('view.vmd', 'view.pml') + (('view.pse',) if session.is_file() else ()):
        print(f'[VISUALIZATION]   {directory / name}')
    return manifest


def safe_scene(directory, groups, mode, **kwargs):
    try:
        return write_scene(directory, groups, mode, **kwargs)
    except Exception as exc:
        print(f'[VISUALIZATION WARNING] Scientific outputs are available; comparison views failed: {exc}')
        return None


def safe_comparison_scene(directory, mode, *args, **kwargs):
    try:
        kwargs.setdefault('cholesterol', mode == 'chl')
        return safe_scene(directory, comparison_groups(*args, **kwargs), mode)
    except Exception as exc:
        print(f'[VISUALIZATION WARNING] Scientific outputs are available; comparison views failed: {exc}')
        return None
