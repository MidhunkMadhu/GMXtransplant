import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import MDAnalysis as mda
import numpy as np

from config import (
    Config,
    ConfigError,
    LigandFitSpec,
    LigandReplaceSpec,
    NewLigandSpec,
    OriginalLigandSpec,
)
from ligand_replace import (
    LigandReplaceError,
    _apply_itp_atom_names,
    assemble_ligand_replacement,
)


def _universe(names, coordinates, resnames, atom_resindex):
    universe = mda.Universe.empty(
        len(names),
        n_residues=len(resnames),
        atom_resindex=np.asarray(atom_resindex),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", names)
    universe.add_TopologyAttr("resnames", resnames)
    universe.add_TopologyAttr("resids", list(range(1, len(resnames) + 1)))
    universe.atoms.positions = np.asarray(coordinates, dtype=np.float32)
    universe.dimensions = np.asarray([50, 50, 50, 90, 90, 90], dtype=np.float32)
    return universe


class LigandFitTests(unittest.TestCase):
    def setUp(self):
        self.old_by_name = {
            "C1": np.array([1.0, 2.0, 3.0]),
            "C2": np.array([3.0, 2.0, 3.0]),
            "N1": np.array([1.0, 5.0, 3.0]),
            "O1": np.array([1.0, 2.0, 7.0]),
        }
        self.environment_position = np.array([30.0, 31.0, 32.0])
        self.structure = _universe(
            ["CA", "C1", "C2", "N1", "O1", "H1"],
            [
                self.environment_position,
                self.old_by_name["C1"],
                self.old_by_name["C2"],
                self.old_by_name["N1"],
                self.old_by_name["O1"],
                [2.0, 3.0, 3.0],
            ],
            ["ALA", "OLD"],
            [0, 1, 1, 1, 1, 1],
        )

    def _new_ligand(self):
        # A translated/rotated ligand with a different atom order and different
        # hydrogen name. Autofit must use the four heavy atoms and pair by name.
        rotation = np.array(
            [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        )
        translation = np.array([11.0, -4.0, 8.0])
        names = ["H9", "O1", "C2", "N1", "C1"]
        target = {
            **self.old_by_name,
            "H9": np.array([2.0, 3.0, 3.0]),
        }
        coordinates = [
            (target[name] - translation) @ rotation for name in names
        ]
        return _universe(names, coordinates, ["NEW"], [0] * len(names))

    def _spec(self, method, old_atoms=None, new_atoms=None):
        return LigandReplaceSpec(
            enabled=True,
            structure_path="structure.gro",
            original_ligand=OriginalLigandSpec(resname="OLD"),
            new_ligand=NewLigandSpec(coord_path="new.mol2", resname="NEW"),
            fit=LigandFitSpec(
                method=method,
                old_ligand_fit_atoms=old_atoms or [],
                new_ligand_fit_atoms=new_atoms or [],
            ),
        )

    def test_autofit_uses_all_heavy_atoms_by_name_and_preserves_environment(self):
        incoming = self._new_ligand()
        with patch("ligand_replace.load_new_ligand", return_value=incoming):
            result = assemble_ligand_replacement(
                self.structure, self._spec("autofit")
            )

        self.assertEqual(result.n_fit_atoms, 4)
        self.assertLess(result.fit_rmsd, 1e-5)
        self.assertEqual(
            result.fit_atom_pairs,
            [("C1", "C1"), ("C2", "C2"), ("N1", "N1"), ("O1", "O1")],
        )
        np.testing.assert_allclose(
            result.merged_universe.atoms.positions[0],
            self.environment_position,
            atol=1e-6,
        )
        positioned = {
            str(atom.name): atom.position.copy()
            for atom in result.new_ligand_ag_positioned
        }
        for name, expected in self.old_by_name.items():
            np.testing.assert_allclose(positioned[name], expected, atol=1e-5)

    def test_pairfit_remains_available(self):
        incoming = self._new_ligand()
        fit_names = ["C1", "C2", "N1"]
        with patch("ligand_replace.load_new_ligand", return_value=incoming):
            result = assemble_ligand_replacement(
                self.structure,
                self._spec("pairfit", fit_names, fit_names),
            )
        self.assertEqual(result.fit_method, "pairfit")
        self.assertEqual(result.n_fit_atoms, 3)
        self.assertLess(result.fit_rmsd, 1e-5)

    def test_nofit_keeps_incoming_coordinates_exactly(self):
        incoming = self._new_ligand()
        original = incoming.atoms.positions.copy()
        with patch("ligand_replace.load_new_ligand", return_value=incoming):
            result = assemble_ligand_replacement(
                self.structure, self._spec("nofit")
            )
        np.testing.assert_array_equal(
            result.new_ligand_ag_positioned.positions, original
        )
        self.assertEqual(result.n_fit_atoms, 0)
        self.assertIsNone(result.fit_rmsd)

    def test_autofit_rejects_different_heavy_atom_sets(self):
        incoming = self._new_ligand()
        incoming.atoms.names = ["H9", "S1", "C2", "N1", "C1"]
        with patch("ligand_replace.load_new_ligand", return_value=incoming):
            with self.assertRaisesRegex(
                LigandReplaceError, "identical.*heavy-atom sets"
            ):
                assemble_ligand_replacement(
                    self.structure, self._spec("autofit")
                )

    def test_itp_names_restore_generic_coordinate_names(self):
        incoming = _universe(
            ["C", "N", "H"],
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            ["LIG"],
            [0, 0, 0],
        )
        contents = """
[ moleculetype ]
LIG 3

[ atoms ]
1 CT 1 LIG C1 1 0.0
2 NT 1 LIG N1 2 0.0
3 HT 1 LIG H1 3 0.0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ligand.itp"
            path.write_text(contents)
            _apply_itp_atom_names(incoming, str(path), "LIG")
        self.assertEqual(list(incoming.atoms.names), ["C1", "N1", "H1"])

    def test_config_accepts_exactly_three_fit_modes(self):
        for method in ("pairfit", "autofit", "nofit"):
            old_atoms = ["C1", "C2", "N1"] if method == "pairfit" else []
            cfg = Config(mode="lig", ligand_replace=self._spec(method, old_atoms, old_atoms))
            cfg.validate(check_paths=False)

        cfg = Config(mode="lig", ligand_replace=self._spec("guessfit"))
        with self.assertRaisesRegex(
            ConfigError, "pairfit, autofit, nofit"
        ):
            cfg.validate(check_paths=False)


if __name__ == "__main__":
    unittest.main()
