import tempfile
import unittest
from pathlib import Path

import numpy as np

from charmm36_cholesterol import (
    DEFAULT_REFERENCE_PDB,
    PipelineError,
    audit_and_reconstruct_heavy_atoms,
    build_bonds,
    read_pdb_records,
    reference_atoms,
)


class CholesterolHeavyAtomAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temporary_directory = tempfile.TemporaryDirectory()
        reference_path = Path(cls._temporary_directory.name) / "reference.pdb"
        reference_path.write_text(DEFAULT_REFERENCE_PDB)
        cls.reference = reference_atoms(
            read_pdb_records(reference_path), reference_path
        )
        cls.heavy = [atom for atom in cls.reference if not atom["is_h"]]
        cls.key = ("X", "CHL1", 1, "")

    @classmethod
    def tearDownClass(cls):
        cls._temporary_directory.cleanup()

    def audit(self, atoms, *, reconstruct=True, max_missing=2):
        return audit_and_reconstruct_heavy_atoms(
            atoms,
            self.reference,
            self.key,
            reconstruct,
            max_missing,
            1.5,
        )

    def test_complete_reference_has_28_heavy_atoms(self):
        audited, report = self.audit(self.heavy)
        self.assertEqual(len(audited), 28)
        self.assertEqual(report["observed_heavy_atom_count"], 28)
        self.assertEqual(report["expected_heavy_atom_count"], 28)
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["modeled_heavy_atom_names"], [])

    def test_one_named_missing_atom_is_reconstructed(self):
        incomplete = [atom for atom in self.heavy if atom["name"] != "C27"]
        reconstructed, report = self.audit(incomplete)
        self.assertEqual(len(reconstructed), 28)
        self.assertEqual(report["observed_heavy_atom_count"], 27)
        self.assertEqual(report["missing_heavy_atom_names"], ["C27"])
        self.assertEqual(report["modeled_heavy_atom_names"], ["C27"])
        self.assertEqual(report["status"], "reconstructed")
        self.assertAlmostEqual(report["reference_fit_rmsd_angstrom"], 0.0)

    def test_flexible_tail_atom_uses_a_local_reference_frame(self):
        atoms = [dict(atom) for atom in self.heavy]
        names = {atom["name"]: index for index, atom in enumerate(atoms)}
        bonds = build_bonds(atoms)
        fixed_index = names["C17"]
        tail_root = names["C20"]

        # Rotate the complete side chain around C17-C20. This keeps its bond
        # graph intact while making a whole-cholesterol fit unsuitable for
        # placing terminal C27, which must follow the local tail conformation.
        tail_indices = set()
        frontier = [tail_root]
        while frontier:
            current = frontier.pop()
            if current in tail_indices:
                continue
            tail_indices.add(current)
            frontier.extend(
                neighbor for neighbor in bonds[current]
                if neighbor != fixed_index and neighbor not in tail_indices
            )
        pivot = np.array([
            atoms[fixed_index]["x"], atoms[fixed_index]["y"], atoms[fixed_index]["z"]
        ])
        root = np.array([
            atoms[tail_root]["x"], atoms[tail_root]["y"], atoms[tail_root]["z"]
        ])
        axis = (root - pivot) / np.linalg.norm(root - pivot)
        angle = np.deg2rad(45.0)
        cross = np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        rotation = (
            np.eye(3) * np.cos(angle)
            + (1.0 - np.cos(angle)) * np.outer(axis, axis)
            + np.sin(angle) * cross
        )
        for index in tail_indices:
            xyz = np.array([atoms[index]["x"], atoms[index]["y"], atoms[index]["z"]])
            xyz = pivot + rotation @ (xyz - pivot)
            atoms[index].update(x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))

        incomplete = [atom for atom in atoms if atom["name"] != "C27"]
        reconstructed, report = self.audit(incomplete)
        by_name = {atom["name"]: atom for atom in reconstructed}
        c25 = by_name["C25"]
        c27 = by_name["C27"]
        bond_length = np.linalg.norm(np.array([
            c27["x"] - c25["x"],
            c27["y"] - c25["y"],
            c27["z"] - c25["z"],
        ]))
        self.assertLess(bond_length, 1.85)
        self.assertEqual(report["status"], "reconstructed")
        self.assertIn("C27", report["modeled_atom_local_fit_rmsd_angstrom"])

    def test_more_than_configured_missing_atoms_is_rejected(self):
        incomplete = [
            atom for atom in self.heavy if atom["name"] not in {"C25", "C26", "C27"}
        ]
        with self.assertRaisesRegex(PipelineError, "exceeding max_missing_heavy_atoms"):
            self.audit(incomplete)

    def test_incomplete_unknown_naming_is_rejected_as_ambiguous(self):
        incomplete = [dict(atom) for atom in self.heavy if atom["name"] != "C27"]
        incomplete[0]["name"] = "CX"
        with self.assertRaisesRegex(PipelineError, "unambiguous reference-name subset"):
            self.audit(incomplete)

    def test_reconstruction_can_be_disabled(self):
        incomplete = [atom for atom in self.heavy if atom["name"] != "C27"]
        with self.assertRaisesRegex(PipelineError, "reconstruction is disabled"):
            self.audit(incomplete, reconstruct=False)


if __name__ == "__main__":
    unittest.main()
