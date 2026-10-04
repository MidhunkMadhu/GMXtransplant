"""Provenance and portable viewer regressions, independent of installed viewers."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import MDAnalysis as mda
import numpy as np
from visualization import comparison_groups, write_scene


def atoms(names, coords):
    u = mda.Universe.empty(len(names), n_residues=len(names), atom_resindex=np.arange(len(names)), trajectory=True)
    u.add_TopologyAttr('names', ['CA' if n == 'ALA' else 'C1' for n in names])
    u.add_TopologyAttr('resnames', names)
    u.add_TopologyAttr('resids', np.arange(1, len(names)+1))
    u.atoms.positions = coords
    return u


class VisualizationTests(unittest.TestCase):
    def test_same_name_ligands_have_different_provenance_and_colors(self):
        reference = atoms(['ALA', 'LI1', 'CHL1', 'TIP3'], [[0,0,0],[1,1,1],[3,3,3],[4,4,4]])
        incoming = atoms(['LI1'], [[2,2,2]])
        final = mda.Merge(reference.atoms[[0,2]], incoming.atoms)
        before = reference.atoms.positions.copy()
        groups = comparison_groups(final, reference, incoming.atoms, reference.atoms[[1]], incoming)
        self.assertEqual(len(groups['inserted_ligand']), 1)
        self.assertEqual(len(groups['removed_ligand']), 1)
        # Outside cholesterol mode, cholesterol is drawn as an ordinary lipid.
        self.assertNotIn('retained_cholesterol', groups)
        self.assertEqual(len(groups['retained_lipid']), 1)
        self.assertEqual(len(groups['removed_environment']), 1)
        np.testing.assert_array_equal(reference.atoms.positions, before)
        with tempfile.TemporaryDirectory(prefix='views with spaces ') as temp, patch('visualization.shutil.which', return_value=None):
            scene = write_scene(temp, groups, 'lig')
            by_name = {o['name']: o for o in scene['objects']}
            self.assertNotEqual(by_name['inserted_ligand']['color'], by_name['removed_ligand']['color'])
            self.assertFalse(by_name['reference_input']['visible'])
            self.assertTrue((Path(temp)/'view.pml').is_file())
            self.assertIn('[info script]', (Path(temp)/'view.vmd').read_text())
            self.assertEqual(json.loads((Path(temp)/'visualization/scene.json').read_text())['mode'], 'lig')

    def test_cholesterol_is_only_a_separate_group_in_cholesterol_mode(self):
        reference = atoms(['ALA', 'CHL1', 'POPC'], [[0, 0, 0], [3, 3, 3], [5, 5, 5]])
        final = mda.Merge(reference.atoms)
        plain = comparison_groups(final)
        self.assertNotIn('retained_cholesterol', plain)
        self.assertEqual(len(plain['retained_lipid']), 2)
        restoration = comparison_groups(final, cholesterol=True)
        self.assertEqual(len(restoration['retained_cholesterol']), 1)
        self.assertEqual(len(restoration['retained_lipid']), 1)

    def test_two_folder_schema_and_registry(self):
        from charmprot import load_charmprot_config
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'job.yaml'
            p.write_text(f'paths:\n  ref: {Path(temp) / "reference"}\ncharmprot:\n'
                         f'  reference: ${{ref}}\n  transplant: {Path(temp) / "donor"}\n')
            cfg=load_charmprot_config(p,check_paths=False,output_root=temp)
            self.assertEqual(cfg.reference,str(Path(temp)/'reference'))
            self.assertEqual(cfg.ligands,'auto')
            self.assertIsNone(cfg.reference_proteins)

    def test_first_gui_mode_and_no_path_add_remove(self):
        from gmxtransplant.gui.model import MODES
        self.assertEqual(next(iter(MODES)), 'charmprot')
        try:
            from PySide6.QtWidgets import QApplication, QPushButton
            from gmxtransplant.gui.editor import ValueEditor
        except ImportError:
            self.skipTest('Optional Qt is not installed')
        app = QApplication.instance() or QApplication([])
        editor=ValueEditor({'source':'input.gro'},path=('paths',))
        buttons=[b.text() for b in editor.findChildren(QPushButton)]
        self.assertNotIn('+ Add entry', buttons)
        self.assertNotIn('Remove', buttons)
        editor.deleteLater()


class RepresentationTests(unittest.TestCase):
    """Proteins and peptides as cartoon, ligands and small molecules as licorice, lipids as lines."""

    def test_style_follows_what_an_object_contains(self):
        import warnings
        import MDAnalysis as mda
        from visualization import _style
        root = Path(__file__).resolve().parent / 'examples' / 'addbinder'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            gs = mda.Universe(str(root / 'gprotein' / 'gs_trimer.pdb')).atoms
            dopamine = mda.Universe(str(root / 'dopamine' / 'dopamine.pdb')).atoms
            mixed = mda.Merge(gs[-897:], dopamine).atoms
        self.assertEqual(_style('binder_binderpose1', gs), {'cartoon': True, 'rest': None, 'spheres': None})
        self.assertEqual(_style('binder_binderpose1', dopamine), {'cartoon': False, 'rest': 'licorice', 'spheres': None})
        self.assertEqual(_style('binder_x', mixed), {'cartoon': True, 'rest': 'licorice', 'spheres': None})
        self.assertEqual(_style('binder_x', gs.residues[:2].atoms)['rest'], 'licorice')  # a dipeptide
        self.assertEqual(_style('retained_lipid', dopamine)['rest'], 'lines')
        self.assertEqual(_style('host_tip', dopamine[:1])['spheres'], 'tip')

    def test_scripts_draw_cartoon_licorice_and_lines(self):
        load = Path(__file__).resolve().parent / 'visualization.py'
        text = load.read_text()
        for snippet in ("cmd.show('cartoon', name + ' and not hetatm')", "cmd.show('sticks', rest)",
                        "cmd.show('lines', rest)", "'Licorice 0.3 12 12'", "'Lines 1.0'",
                        "'mol representation NewCartoon'"):
            self.assertIn(snippet, text)
