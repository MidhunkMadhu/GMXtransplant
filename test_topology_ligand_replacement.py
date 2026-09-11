import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import MDAnalysis as mda
import numpy as np

from config import TopologySpec
from charges import compute_pre_neutralization_report
from config import ChargeSpec
from topology import (
    TopologyError,
    assemble_topology,
    audit_final_topology,
    build_topology_charge_model,
    collect_topology_definitions,
    parse_top_includes,
)


def _write_itp(path, moleculetype, atoms):
    rows = [
        "[ moleculetype ]",
        f"{moleculetype} 3",
        "",
        "[ atoms ]",
    ]
    for index, (atom_type, resname, atom_name, charge) in enumerate(atoms, 1):
        rows.append(
            f"{index} {atom_type} 1 {resname} {atom_name} {index} {charge}"
        )
    path.write_text("\n".join(rows) + "\n")


def _final_universe():
    universe = mda.Universe.empty(
        3,
        n_residues=2,
        atom_resindex=np.asarray([0, 1, 1]),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["CA", "N1", "C1"])
    universe.add_TopologyAttr("resnames", ["ALA", "UJU"])
    universe.add_TopologyAttr("resids", [1, 2])
    universe.atoms.positions = np.zeros((3, 3), dtype=np.float32)
    universe.dimensions = np.asarray([30, 30, 30, 90, 90, 90], dtype=np.float32)
    return universe


class LigandReplacementTopologyTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.template = self.root / "template"
        self.toppar = self.template / "toppar"
        self.output = self.root / "output"
        self.toppar.mkdir(parents=True)

        (self.toppar / "forcefield.itp").write_text("; parameter-only test file\n")
        _write_itp(
            self.toppar / "PROT.itp",
            "PROT",
            [("CT_OLD", "ALA", "CA", "0.0")],
        )
        self.old_ligand_itp = self.toppar / "UJU.itp"
        _write_itp(
            self.old_ligand_itp,
            "UJU",
            [
                ("C_OLD", "UJU", "C1", "-0.25"),
                ("N_OLD", "UJU", "N1", "0.25"),
            ],
        )
        (self.template / "topol.top").write_text(
            '#include "toppar/forcefield.itp"\n\n'
            "[ molecules ]\n"
            "PROT 1\n"
            "UJU 1\n"
        )

        self.new_ligand_itp = self.root / "new" / "UJU.itp"
        self.new_ligand_itp.parent.mkdir()
        self.new_forcefield_itp = self.new_ligand_itp.parent / "forcefield.itp"
        self.new_forcefield_itp.write_text("; replacement force field\n")
        _write_itp(
            self.new_ligand_itp,
            "UJU",
            [
                ("N_NEW", "UJU", "N1", "0.6"),
                ("C_NEW", "UJU", "C1", "0.4"),
            ],
        )
        self.topology = TopologySpec(
            enabled=True,
            receptor_toppar_dir=str(self.toppar),
            receptor_template_top=str(self.template / "topol.top"),
            environment_toppar_dir=str(self.toppar),
            ligand_itp_paths=[str(self.new_ligand_itp)],
            output_dir=str(self.output),
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_conflict_remains_an_error_without_replacement_precedence(self):
        with self.assertRaisesRegex(TopologyError, "conflicting definitions"):
            collect_topology_definitions(self.topology)

    def test_charge_model_uses_new_ligand_definition(self):
        model = build_topology_charge_model(
            self.topology,
            {"UJU"},
            replacement_ligand_itp_paths=[str(self.new_ligand_itp)],
        )

        self.assertEqual(model.files["UJU"], str(self.new_ligand_itp))
        self.assertEqual(model.definitions["UJU"].atom_names, ("N1", "C1"))
        self.assertAlmostEqual(model.definitions["UJU"].charge, 1.0)

    def test_generated_topology_copies_and_includes_only_new_ligand(self):
        self.output.mkdir()
        (self.output / "toppar").mkdir()
        (self.output / "toppar" / "stale.itp").write_text("; stale output\n")
        result = assemble_topology(
            _final_universe(),
            self.topology,
            replacement_ligand_itp_paths=[str(self.new_ligand_itp)],
        )

        copied_ligand = Path(result.toppar_dir) / "UJU.itp"
        self.assertEqual(copied_ligand.read_bytes(), self.new_ligand_itp.read_bytes())
        self.assertNotEqual(copied_ligand.read_bytes(), self.old_ligand_itp.read_bytes())
        copied_forcefield = Path(result.toppar_dir) / "forcefield.itp"
        self.assertEqual(
            copied_forcefield.read_bytes(), self.new_forcefield_itp.read_bytes()
        )
        self.assertFalse((Path(result.toppar_dir) / "stale.itp").exists())
        self.assertEqual(parse_top_includes(result.top_path)[0], "toppar/forcefield.itp")
        self.assertIn("toppar/UJU.itp", parse_top_includes(result.top_path))
        self.assertEqual(result.molecules_written, [("PROT", 1), ("UJU", 1)])
        audit = audit_final_topology(_final_universe(), result)
        self.assertEqual(audit["coordinate_atoms"], 3)
        self.assertAlmostEqual(audit["net_charge"], 1.0)

    def test_failed_staged_audit_preserves_existing_topology_bundle(self):
        self.output.mkdir()
        old_top = self.output / "topol.top"
        old_top.write_text("old topology\n")
        old_toppar = self.output / "toppar"
        old_toppar.mkdir()
        old_itp = old_toppar / "old.itp"
        old_itp.write_text("old parameters\n")

        with patch(
            "topology.audit_final_topology",
            side_effect=TopologyError("simulated staged audit failure"),
        ):
            with self.assertRaisesRegex(TopologyError, "simulated staged audit failure"):
                assemble_topology(
                    _final_universe(),
                    self.topology,
                    replacement_ligand_itp_paths=[str(self.new_ligand_itp)],
                )

        self.assertEqual(old_top.read_text(), "old topology\n")
        self.assertEqual(old_itp.read_text(), "old parameters\n")
        self.assertFalse(any(self.output.glob(".gmxtransplant-topology-*")))

    def test_removed_old_ligand_does_not_have_to_match_new_itp_order(self):
        model = build_topology_charge_model(
            self.topology,
            {"UJU"},
            replacement_ligand_itp_paths=[str(self.new_ligand_itp)],
        )
        final = _final_universe()
        old_ligand = mda.Merge(final.atoms[[2, 1]])
        new_ligand = final.atoms[1:]
        empty = final.atoms[0:0]

        report = compute_pre_neutralization_report(
            old_ligand,
            new_ligand,
            empty,
            empty,
            ChargeSpec(),
            topology_charge_model=model,
            resolve_original_charge=False,
        )

        self.assertTrue(np.isnan(report.charge_original_receptor))
        self.assertAlmostEqual(report.charge_replacement_receptor, 1.0)
        self.assertAlmostEqual(report.net_charge_after_clash_removal, 1.0)


if __name__ == "__main__":
    unittest.main()
