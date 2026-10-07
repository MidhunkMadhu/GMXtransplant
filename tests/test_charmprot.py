"""Topology/coordinate integration tests for automatic CHARMM-GUI transplants."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import MDAnalysis as mda
import numpy as np
import yaml

from charmprot import (inspect_system, load_charmprot_config, run_charmprot,
                       _sequence_pairs, align_systems)
from charmm_environment import CATALOGUE, environment_name
from config import ConfigError


def build_fixture(folder, donor=False, ligand=True, multi_lipid=True):
    """Tiny real GRO/ITP systems: donor is rotated and translated, ligand differs."""
    (folder / "toppar").mkdir(parents=True)
    definitions = {
        "PROA": [(i+1, rn, "CA", 0) for i, rn in enumerate(("ALA", "GLY", "SER", "LEU"))],
        "TIP3": [(1, "TIP3", "OH2", 0)],
        "CLA": [(1, "CLA", "CLA", -1)],
        "SOD": [(1, "SOD", "SOD", 1)],
    }
    coordinates = {"PROA": [(10, 10, 10), (14, 10, 10), (14, 14, 10), (14, 14, 14)],
                   "TIP3": [(40, 40, 40)], "CLA": [(60, 60, 60)], "SOD": [(70, 70, 70)]}
    if ligand:
        name = "NEW" if donor else "OLD"
        definitions[name] = [(1, name, "C1", -1 if donor else 0), (1, name, "C2", 0)]
        coordinates[name] = [(18, 14, 14), (19, 14, 14)]
    if multi_lipid:
        definitions["POPC"] = [(1, "HEAD", "C1", 0), (2, "TAIL", "C2", 0)]
        # A clash in just the head must remove BOTH residues.
        coordinates["POPC"] = [(10.1, 10.1, 10.1), (20, 20, 20)]
    order = ["PROA"] + (["NEW" if donor else "OLD"] if ligand else [])
    order += (["POPC"] if multi_lipid else []) + ["SOD", "CLA", "TIP3"]
    ff = "[ defaults ]\n1 2 yes 1.0 1.0\n[ atomtypes ]\nC 6 12.01 0 A 0.3 0.1\n"
    (folder / "toppar" / "forcefield.itp").write_text(ff)
    gro_lines, residue_number, last_residue = [], 0, None
    atom_number = 0
    for name in order:
        rows = definitions[name]
        itp = f"[ moleculetype ]\n{name} 3\n[ atoms ]\n"
        for i, ((resnr, resname, atom, charge), position) in enumerate(zip(rows, coordinates[name]), 1):
            itp += f"{i} C {resnr} {resname} {atom} {i} {charge} 12.01\n"
            identity = (name, resnr)
            if identity != last_residue:
                residue_number += 1
                last_residue = identity
            atom_number += 1
            pos = np.array(position, dtype=float)
            if donor:
                pos = pos @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [35, 40, 5]
            gro_lines.append(f"{residue_number:5d}{resname:<5}{atom:>5}{atom_number:5d}" +
                             "".join(f"{p/10:8.3f}" for p in pos))
        (folder / "toppar" / f"{name}.itp").write_text(itp)
    (folder / "step5_input.gro").write_text("fixture\n" + str(atom_number) + "\n" +
        "\n".join(gro_lines) + "\n 10.0 10.0 10.0\n")
    (folder / "topol.top").write_text('#include "toppar/forcefield.itp"\n' +
        "\n".join(f'#include "toppar/{name}.itp"' for name in order) +
        "\n[ system ]\nfixture\n[ molecules ]\n" + "\n".join(f"{name} 1" for name in order) + "\n")


class CharmProtTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref, self.donor = self.root / "reference", self.root / "donor"
        build_fixture(self.ref)
        build_fixture(self.donor, donor=True)
        self.config = self.root / "input.yaml"
        # Input folders are always full paths; the output root replaces the
        # former "relative to the configuration file" behaviour.
        self.folders = f"reference: {self.ref}\ntransplant: {self.donor}\n"
        self.config.write_text(self.folders)

    def run_config(self, **kwargs):
        return run_charmprot(self.config, output_root=self.root, **kwargs)

    def test_catalogue_covers_primary_sources_and_aliases_without_broad_prefix(self):
        self.assertGreater(len(CATALOGUE["names"]), 900)
        for name in "SOD POT CLA NA+ K+ CL- WAT TIP3 TIP3P OPC SWM4 POP POPE POPS DOPE DOPC CHL1 CLR".split():
            self.assertIsNotNone(environment_name(name), name)
        for name in ("LI1", "LI2", "ATP", "POPCUSTOM", "POISON"):
            self.assertIsNone(environment_name(name), name)
        self.assertEqual(environment_name("FE2P")[0], "ion")

    def test_detection_and_override(self):
        system = inspect_system(self.ref)
        self.assertEqual(system.report["selected_proteins"], ["PROA"])
        self.assertEqual(system.report["selected_ligands"], ["OLD"])
        self.assertEqual([m.category for m in system.molecules if m.name == "POPC"], ["lipid"])
        system = inspect_system(self.ref, ligands=["POPC"])
        self.assertEqual(system.report["selected_ligands"], ["POPC"])
        self.assertFalse(next(m for m in system.molecules if m.name == "OLD").selected)

    def test_sequence_gaps_and_rigid_ligand_transform(self):
        self.assertEqual(_sequence_pairs(["ALA", "GLY", "SER"], ["ALA", "SER"]), [(0, 0), (2, 1)])
        ref, donor = inspect_system(self.ref), inspect_system(self.donor)
        report = align_systems(ref, donor)
        self.assertLess(report["rmsd_angstrom"], 1e-5)
        mol = next(m for m in donor.molecules if m.name == "NEW")
        np.testing.assert_allclose(donor.universe.atoms[mol.start].position, [18, 14, 14], atol=1e-4)

    def test_end_to_end_removes_whole_multiresidue_molecule_and_neutralizes(self):
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            self.assertEqual(self.run_config(), 0)
        folder = self.root / "charmprot_output"
        report = json.loads((folder / "charmprot_report.json").read_text())
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["clash_removals"][0]["atoms"], 2)
        self.assertEqual(report["charge"]["removed_ions"][0]["moleculetype"], "CLA")
        self.assertAlmostEqual(report["charge"]["after"], 0)
        result = mda.Universe(str(folder / "step5_input.gro"))
        self.assertEqual(len(result.atoms), 8)
        self.assertIn("NEW", result.residues.resnames)
        self.assertNotIn("OLD", result.residues.resnames)
        self.assertNotIn("HEAD", result.residues.resnames)
        self.assertNotIn("TAIL", result.residues.resnames)
        self.assertEqual(report["topology_audit"]["topology_atoms"], 8)

    def test_ignore_removes_reference_ligand_and_omits_transplant_ligand(self):
        self.config.write_text(self.config.read_text() + "ligands: ignore\n")
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            self.run_config()
        result = mda.Universe(str(self.root / "charmprot_output/step5_input.gro"))
        # Reference protein and ligand are removed; only the transplant protein is inserted.
        self.assertNotIn("OLD", result.residues.resnames)
        self.assertNotIn("NEW", result.residues.resnames)
        report = json.loads((self.root / "charmprot_output/charmprot_report.json").read_text())
        self.assertEqual(report["reference"]["selected_ligands"], ["OLD"])
        self.assertEqual(report["transplant"]["selected_ligands"], [])

    def test_itp_order_mismatch_fails_before_output(self):
        gro = self.ref / "step5_input.gro"
        gro.write_text(gro.read_text().replace("   CA", "   XX", 1))
        with self.assertRaisesRegex(ConfigError, "atom-order mismatch"):
            inspect_system(self.ref)

    def test_forcefield_conflict_publishes_only_failure_report(self):
        ff = self.donor / "toppar/forcefield.itp"
        ff.write_text(ff.read_text().replace("0.3 0.1", "0.3 0.9"))
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            with self.assertRaisesRegex(ConfigError, "Incompatible force-field"):
                self.run_config()
        files = {p.name for p in (self.root / "charmprot_output").iterdir()}
        self.assertEqual(files, {"charmprot_report.json", "charmprot_report.txt"})

    def test_invalid_configs_and_output_protection(self):
        for suffix in ("neutralize: yesplease\n", "clash_distance: .nan\n", "reference_proteins: [MISSING]\n",
                       "output_dir: reference/subfolder\n", "ligands: ignore\ntransplant_ligands: [NEW]\n"):
            self.config.write_text(self.folders + suffix)
            with self.subTest(suffix=suffix):
                with self.assertRaises(ConfigError):
                    if "MISSING" in suffix:
                        inspect_system(self.ref, proteins=["MISSING"])
                    else:
                        load_charmprot_config(self.config, output_root=self.root)

    def test_dry_run_does_not_open_inputs(self):
        self.config.write_text(f"reference: {self.root / 'missing_ref'}\n"
                               f"transplant: {self.root / 'missing_donor'}\n")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.run_config(dry_run=True), 0)
        self.assertFalse((self.root / "charmprot_output").exists())

    def test_zero_ligands_is_explicit_and_unreferenced_itps_are_ignored(self):
        (self.ref / "toppar/UNUSED.itp").write_text("this is not a molecule\n")
        system = inspect_system(self.ref, ligands=[])
        self.assertEqual(system.report["selected_ligands"], [])
        self.assertFalse(next(m for m in system.molecules if m.name == "OLD").selected)

    def test_conditional_include_is_rejected(self):
        top = self.ref / "topol.top"
        top.write_text(top.read_text().replace('#include "toppar/PROA.itp"',
                       '#ifdef ALTERNATE\n#include "toppar/PROA.itp"\n#endif'))
        with self.assertRaisesRegex(ConfigError, "Conditional include"):
            inspect_system(self.ref)

    def test_existing_output_directory_is_replaced(self):
        folder = self.root / "charmprot_output"
        folder.mkdir()
        (folder / "keep.txt").write_text("stale run")
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            self.assertEqual(self.run_config(), 0)
        self.assertTrue((folder / "keep.txt").exists())
        self.assertTrue((folder / "step5_input.gro").exists())
        self.assertFalse((self.root / "charmprot_output.backup").exists())

    def test_failed_rerun_replaces_prior_successful_output(self):
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            self.assertEqual(self.run_config(), 0)
        ff = self.donor / "toppar/forcefield.itp"
        ff.write_text(ff.read_text().replace("0.3 0.1", "0.3 0.9"))
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            with self.assertRaisesRegex(ConfigError, "Incompatible force-field"):
                self.run_config()
        folder = self.root / "charmprot_output"
        files = {p.name for p in folder.iterdir()}
        self.assertEqual(files, {"charmprot_report.json", "charmprot_report.txt"})
        report = json.loads((folder / "charmprot_report.json").read_text())
        self.assertEqual(report["status"], "failed")
        self.assertFalse((self.root / "charmprot_output.backup").exists())

    def test_multi_chain_alignment_uses_all_chains(self):
        for folder in (self.ref, self.donor):
            top = folder / "topol.top"
            text = top.read_text().replace('#include "toppar/PROA.itp"',
                '#include "toppar/PROA.itp"\n#include "toppar/PROB.itp"')
            top.write_text(text.replace("PROA 1", "PROA 1\nPROB 1"))
            source = folder / "toppar/PROA.itp"
            (folder / "toppar/PROB.itp").write_text(source.read_text().replace("PROA", "PROB"))
            gro = folder / "step5_input.gro"
            lines = gro.read_text().splitlines()
            lines[1] = str(int(lines[1]) + 4)
            # A distinct second chain, with topology residue IDs independent of GRO IDs.
            extra = []
            for line in lines[2:6]:
                extra.append(f"{int(line[:5])+4:5d}" + line[5:20] +
                             f"{float(line[20:28])+2:8.3f}" + line[28:])
            gro.write_text("\n".join(lines[:6] + extra + lines[6:]) + "\n")
        ref, donor = inspect_system(self.ref), inspect_system(self.donor)
        self.assertEqual(ref.report["selected_proteins"], ["PROA", "PROB"])
        report = align_systems(ref, donor)
        self.assertEqual(report["fit_atoms"], 8)
        self.assertEqual(len(report["chains"]), 2)

    def test_periodic_clash_removes_complete_molecule(self):
        from charmprot import remove_clashes, CharmProtSpec
        ref, donor = inspect_system(self.ref), inspect_system(self.donor)
        align_systems(ref, donor)
        lipid = next(m for m in ref.molecules if m.name == "POPC")
        ref.universe.atoms[lipid.start].position = [110.1, 10.1, 10.1]
        selected = [m for m in donor.molecules if m.selected]
        indices = np.concatenate([np.arange(m.start, m.stop) for m in selected])
        kept, removed, retained = remove_clashes(ref, donor.universe.atoms[indices], CharmProtSpec(str(self.ref), str(self.donor)))
        self.assertFalse(any(m.name == "POPC" for m in kept))
        self.assertEqual(next(m for m in removed if m["moleculetype"] == "POPC")["atoms"], 2)
        self.assertEqual(retained, [])

    def test_keep_lipids_retains_clashing_lipids_only(self):
        from charmprot import remove_clashes, CharmProtSpec
        ref, donor = inspect_system(self.ref), inspect_system(self.donor)
        align_systems(ref, donor)
        lipid = next(m for m in ref.molecules if m.name == "POPC")
        ref.universe.atoms[lipid.start].position = [110.1, 10.1, 10.1]
        selected = [m for m in donor.molecules if m.selected]
        indices = np.concatenate([np.arange(m.start, m.stop) for m in selected])
        spec = CharmProtSpec(str(self.ref), str(self.donor), keep_lipids=True)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            kept, removed, retained = remove_clashes(ref, donor.universe.atoms[indices], spec)
        self.assertTrue(any(m.name == "POPC" for m in kept))
        self.assertFalse(any(m["category"] in ("lipid", "sterol") for m in removed))
        self.assertEqual([m["moleculetype"] for m in retained], ["POPC"])
        self.assertIn("WARNING: keep_lipids retained 1", out.getvalue())


EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "charmprot"


@unittest.skipUnless((EXAMPLE / "reference" / "toppar" / "forcefield.itp").is_file(), "example data not present")
class CharmprotExampleTopologyTests(unittest.TestCase):
    def test_cmap_tables_stay_within_gromacs_line_limit(self):
        from charmprot import _parameter_lines
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "charmprot.yaml"
            config.write_text(yaml.safe_dump({"charmprot": {
                "reference": str(EXAMPLE / "reference"), "transplant": str(EXAMPLE / "transplant")}}))
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                self.assertEqual(run_charmprot(config, output_root=tmp), 0)
            written = Path(tmp) / "charmprot_output" / "toppar" / "forcefield.itp"
            # GROMACS grompp rejects any line of 4095 characters or more.
            self.assertLess(max(map(len, written.read_text().splitlines())), 4095)
            parsed = _parameter_lines(type("System", (), {"files": [written]})())["cmaptypes"]
            source = _parameter_lines(type("System", (), {"files": [EXAMPLE / "reference/toppar/forcefield.itp"]})())
            self.assertEqual(parsed, source["cmaptypes"])


if __name__ == "__main__":
    unittest.main()
