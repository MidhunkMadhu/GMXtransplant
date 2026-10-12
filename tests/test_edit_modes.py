"""Tests of the system-edit modes (membrane, salt, protonate, mutate, hmr, resize, flood, oligomer).

The end-to-end classes build each bundled D1R example (about 10-40 s each); they are
marked slow and skipped with GMXT_FAST=1 (tools/run_tests.py --fast).
"""
import contextlib
import io
import json
import os
from pathlib import Path
import stat
import tempfile
import unittest
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import yaml

from tests import slow
from config import ConfigError
from mode_registry import EDIT_MODES
from protein_edit import ProteinItp, replace_residue, residue_state
from topology import parse_top_molecules
import charmm36_library as c36

ROOT = Path(__file__).resolve().parent.parent
HOST = ROOT / "examples" / "addbinder" / "host"
PROA = HOST / "toppar" / "PROA.itp"


def _module(mode):
    import importlib
    return importlib.import_module(EDIT_MODES[mode].module)


def _universe(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return mda.Universe(str(path))


def _build(mode, tmp, edits=None):
    """Build the bundled example of `mode` (optionally with changed settings) into tmp."""
    example = ROOT / "examples" / EDIT_MODES[mode].example_dir / EDIT_MODES[mode].example_yaml
    path = example
    if edits:
        document = yaml.safe_load(example.read_text())
        document[mode].update(edits)
        path = example.parent / f".test_{mode}_{os.getpid()}.yaml"
        path.write_text(yaml.safe_dump(document))
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            code = _module(mode).run(str(path), output_root=tmp)
    finally:
        if path != example:
            path.unlink()
    folder = Path(tmp) / f"{mode}_output"
    return code, folder, json.loads((folder / f"{mode}_report.json").read_text())


class RegistryTests(unittest.TestCase):
    def test_every_mode_has_example_template_and_loader(self):
        from run_pipeline import _read_example
        for key, mode in EDIT_MODES.items():
            with self.subTest(mode=key):
                self.assertTrue((ROOT / "examples" / mode.example_dir / mode.example_yaml).is_file())
                text = _read_example(key)
                self.assertIn(f"{key}:", text)
                spec = getattr(_module(key), mode.spec)
                self.assertEqual(spec.__name__, mode.spec)

    def test_examples_validate_in_dry_run(self):
        for key, mode in EDIT_MODES.items():
            with self.subTest(mode=key), tempfile.TemporaryDirectory() as tmp:
                path = ROOT / "examples" / mode.example_dir / mode.example_yaml
                cfg = _module(key).load_config(str(path), check_paths=True, output_root=tmp)
                self.assertEqual(Path(cfg.output_dir).parent, Path(tmp))
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(_module(key).run(str(path), dry_run=True, output_root=tmp), 0)

    def test_unknown_and_bad_options_are_rejected(self):
        cases = {"salt": {"concentration": -1.0}, "hmr": {"hydrogen_mass": 9.0}, "resize": {"add_water": 0.0},
                 "flood": {"side": "middle"}, "oligomer": {"copies": []}, "protonate": {"residues": [
                     {"molecule": "PROA", "resid": 70, "state": "XYZ"}]},
                 "mutate": {"mutations": [{"molecule": "PROA", "resid": 103, "to": "UNK"}]},
                 "membrane": {"max_tilt": 60.0}}
        for mode, edit in cases.items():
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "c.yaml"
                path.write_text(yaml.safe_dump({mode: {"host": str(HOST), **edit}}))
                with self.assertRaises(ConfigError):
                    _module(mode).load_config(str(path), check_paths=False, output_root=tmp)
                path.write_text(yaml.safe_dump({mode: {"host": str(HOST), "no_such_option": 1}}))
                with self.assertRaises(ConfigError):
                    _module(mode).load_config(str(path), check_paths=False, output_root=tmp)

    def test_relative_host_outside_the_package_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.yaml"
            path.write_text(yaml.safe_dump({"salt": {"host": "relative/host"}}))
            with self.assertRaises(ConfigError):
                _module("salt").load_config(str(path), check_paths=False, output_root=tmp)


class LibraryTests(unittest.TestCase):
    def test_host_protein_needs_no_parameters(self):
        text, report = c36.supplement([HOST / "toppar" / "forcefield.itp"], [PROA.read_text()])
        self.assertEqual(report, {})

    def test_new_ion_brings_its_type_and_nbfix(self):
        text, report = c36.supplement([HOST / "toppar" / "forcefield.itp"], [c36.ion_itp("POT")], {"CLA", "OC"})
        self.assertEqual(report["atomtypes"], 1)
        self.assertIn("CLA POT", text.replace("POT CLA", "CLA POT"))

    def test_missing_parameters_raise(self):
        bogus = c36.ion_itp("POT").replace(" POT      1 ", " XXQ      1 ")
        with self.assertRaises(ConfigError):
            c36.supplement([HOST / "toppar" / "forcefield.itp"], [bogus])


class ProteinEditTests(unittest.TestCase):
    """Rebuilding a residue with its own template must reproduce CHARMM-GUI's topology exactly."""

    @classmethod
    def setUpClass(cls):
        cls.text = PROA.read_text()
        cls.positions = _universe(HOST / "step5_input.gro").atoms[:5206].positions.astype(float)

    @staticmethod
    def terms(itp):
        from collections import Counter
        from protein_edit import INDEX_COLUMNS
        out = Counter()
        for item in itp.items:
            width = INDEX_COLUMNS.get(item.section, 0) if item.kind == "row" else 0
            if not width:
                continue
            ids = [(itp.atoms[int(t) - 1].resnr, itp.atoms[int(t) - 1].name) for t in item.tokens[:width]]
            improper = item.section == "dihedrals" and item.tokens[4] == "2"
            if item.section != "cmap" and not improper:
                ids = min(ids, ids[::-1])
            out[(item.section, item.tokens[4] if item.section == "dihedrals" else "", tuple(ids),
                 tuple(item.tokens[width:]) if item.section == "position_restraints" else ())] += 1
        return out

    def test_identity_rebuild_matches_charmm_gui(self):
        base = ProteinItp(self.text)
        reference = self.terms(base)
        first = {}
        for resnr, atoms in base.residues().items():
            first.setdefault(base.atoms[atoms[0]].resname, resnr)
        tested = 0
        for resname, resnr in sorted(first.items()):
            names = [base.atoms[i].name for i in base.residues()[resnr]]
            state = residue_state(resname, names)
            if state is None:
                continue  # the capped N-terminal residue
            with self.subTest(residue=f"{resname}{resnr}"):
                itp = ProteinItp(self.text)
                positions, _ = replace_residue(itp, resnr, state, self.positions)
                self.assertEqual(self.terms(itp), reference)
                np.testing.assert_allclose(positions, self.positions)
                tested += 1
        self.assertGreaterEqual(tested, 19)

    def test_rebuilt_hydrogens_match_the_input(self):
        from scipy.optimize import linear_sum_assignment
        base = ProteinItp(self.text)
        for resnr in (29, 53, 79, 80):    # PHE, HSD, PRO, TRP
            atoms = base.residues()[resnr]
            names = [base.atoms[i].name for i in atoms]
            heavy = [n for i, n in zip(atoms, names) if not base.atoms[i].type.startswith("H")]
            itp = ProteinItp(self.text)
            new, _ = replace_residue(itp, resnr, residue_state(base.atoms[atoms[0]].resname, names),
                                     self.positions, keep=heavy)
            hydrogens = [i for i in itp.residues()[resnr] if itp.atoms[i].name not in heavy]
            original = [atoms[names.index(itp.atoms[i].name)] for i in hydrogens]
            distance = np.linalg.norm(new[hydrogens][:, None] - self.positions[original][None], axis=2)
            rows, cols = linear_sum_assignment(distance)
            self.assertLess(float(np.median(distance[rows, cols])), 0.15)

    def test_state_is_read_from_atoms(self):
        base = ProteinItp(self.text)
        names = [base.atoms[i].name for i in base.residues()[70]]
        self.assertEqual(base.atoms[base.residues()[70][0]].resname, "ASP")
        self.assertEqual(residue_state("ASP", names), "ASPP")

    def test_terminal_and_disulfide_residues_are_refused(self):
        for resnr in (21, 96):   # capped N-terminus; C96-C186 disulfide
            with self.subTest(resnr=resnr), self.assertRaises(ConfigError):
                replace_residue(ProteinItp(self.text), resnr, "ALA", self.positions)


class GromacsCheckTests(unittest.TestCase):
    def _fake(self, folder, code, output):
        script = Path(folder) / "fakegmx"
        script.write_text(f"#!/bin/sh\ncat <<'EOF'\n{output}\nEOF\nexit {code}\n")
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        return str(script)

    def test_pass_fail_and_skip(self):
        from gromacs_check import check_folder
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "topol.top").write_text("")
            (Path(tmp) / "step5_input.gro").write_text("")
            passed = check_folder(tmp, self._fake(tmp, 0, "NOTE 1 [file a]:\n  x\n"))
            self.assertEqual((passed["status"], passed["notes"]), ("passed", 1))
            failed = check_folder(tmp, self._fake(tmp, 1, "Fatal error:\nNo default Bond types\n\n"))
            self.assertEqual(failed["status"], "failed")
            self.assertIn("No default Bond types", failed["errors"][0])
            self.assertEqual(check_folder(Path(tmp) / "missing")["status"], "skipped")

    def test_slurm_script(self):
        from slurm_script import script
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "sample_run.sh").write_text("#!/bin/bash\n")
            text = script(tmp, "dardel", account="naiss0000-0-0").read_text()
            self.assertIn("#SBATCH -A naiss0000-0-0", text)
            self.assertIn("module load PDC/26.03", text)
            self.assertIn("export MDRUN='srun gmx_mpi mdrun -pin on'", text)


def _check_common(test, folder, report):
    test.assertEqual(report["status"], "complete")
    test.assertAlmostEqual(report["final_net_charge"], 0.0, places=3)
    universe = _universe(folder / "step5_input.gro")
    counts = dict(parse_top_molecules(str(folder / "topol.top")))
    test.assertEqual(report["final_molecules"], [list(x) for x in parse_top_molecules(str(folder / "topol.top"))])
    test.assertTrue((folder / "index.ndx").is_file())
    return universe, counts


@slow
class HmrRunTests(unittest.TestCase):
    def test_masses(self):
        from itp import parse_itp
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("hmr", tmp)
            self.assertEqual(code, 0)
            _check_common(self, folder, report)
            protein = next(p for p in (folder / "toppar").glob("molecule_*.itp")
                           if "PROA" in parse_itp(str(p)))
            masses = [float(line.split()[7]) for line in protein.read_text().splitlines()
                      if line.strip() and line.split()[0].isdigit() and len(line.split()) >= 8
                      and line.split()[3].isalpha()]
            original = [float(line.split()[7]) for line in PROA.read_text().splitlines()
                        if line.strip() and line.split()[0].isdigit() and len(line.split()) >= 8
                        and line.split()[3].isalpha()]
            self.assertAlmostEqual(sum(masses), sum(original), places=2)
            self.assertIn(3.024, masses)
            self.assertIn("dt                      = 0.004", (folder / "step7_production.mdp").read_text())


@slow
class SaltRunTests(unittest.TestCase):
    def test_kcl(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("salt", tmp)
            _, counts = _check_common(self, folder, report)
            self.assertNotIn("SOD", counts)
            self.assertGreater(counts["POT"], 250)
            self.assertAlmostEqual(report["salt"]["final"]["concentration_molar"], 0.15, delta=0.005)
            self.assertIn("POT", (folder / "toppar" / "forcefield.itp").read_text())


@slow
class ResizeRunTests(unittest.TestCase):
    def test_taller_box(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("resize", tmp)
            universe, counts = _check_common(self, folder, report)
            host = _universe(HOST / "step5_input.gro")
            self.assertAlmostEqual(float(universe.dimensions[2] - host.dimensions[2]), 20.0, places=2)
            self.assertGreater(counts["TIP3"], 108361 + 8000)
            oxygens = universe.select_atoms("resname TIP3 and name OH2")
            pairs = capped_distance(oxygens.positions, oxygens.positions, 2.2, box=universe.dimensions,
                                    return_distances=False)
            # No two waters on top of each other at the seams.
            self.assertEqual([(a, b) for a, b in pairs if a != b], [])


@slow
class FloodRunTests(unittest.TestCase):
    def test_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("flood", tmp)
            universe, counts = _check_common(self, folder, report)
            self.assertEqual(counts["LDP"], 10)
            copies = universe.select_atoms("resname LDP and not name H*")
            protein = universe.atoms[:5206]
            protein = protein[[not n.startswith("H") for n in protein.names]]
            _, distances = capped_distance(copies.positions, protein.positions, 10.0, box=universe.dimensions)
            self.assertEqual(len(distances), 0)


@slow
class OligomerRunTests(unittest.TestCase):
    def test_dimer(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("oligomer", tmp)
            universe, counts = _check_common(self, folder, report)
            self.assertEqual(counts["PROA"], 2)
            self.assertGreaterEqual(report["oligomer"]["closest_approach_angstrom"], 10.0)
            first, second = universe.atoms[:5206], universe.atoms[5206:10412]
            np.testing.assert_allclose(first.positions, _universe(HOST / "step5_input.gro").atoms[:5206].positions,
                                       atol=1e-3)
            self.assertTrue(all(second.names == first.names))


@slow
class MembraneRunTests(unittest.TestCase):
    def test_popc_cholesterol(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("membrane", tmp)
            universe, counts = _check_common(self, folder, report)
            self.assertNotIn("POPE", counts)
            self.assertAlmostEqual(counts["CHL1"] / (counts["CHL1"] + counts["POPC"]), 0.3, delta=0.05)
            host = _universe(HOST / "step5_input.gro")
            np.testing.assert_allclose(universe.atoms[:5206].positions, host.atoms[:5206].positions, atol=1e-3)
            self.assertTrue((folder / "equilibration" / "analyze_equilibration.py").is_file())


@slow
class ProtonateRunTests(unittest.TestCase):
    def test_states(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("protonate", tmp)
            universe, counts = _check_common(self, folder, report)
            from itp import parse_itp
            protein = next(d for p in (folder / "toppar").glob("molecule_*.itp")
                           for n, d in parse_itp(str(p)).items() if n == "PROA")
            names = {r: [a.atom_name for a in protein.atoms if a.resnr == r] for r in (70, 164)}
            self.assertNotIn("HD2", names[70])
            self.assertIn("HE2", names[164])
            self.assertEqual({a.resname for a in protein.atoms if a.resnr == 164}, {"HSP"})
            self.assertAlmostEqual(protein.charge, parse_itp(str(PROA))["PROA"].charge, places=3)
            he2 = universe.atoms[[a.nr - 1 for a in protein.atoms if a.resnr == 164 and a.atom_name in ("NE2", "HE2")]]
            self.assertAlmostEqual(float(np.linalg.norm(he2.positions[0] - he2.positions[1])), 1.0, delta=0.05)


@slow
class MutateRunTests(unittest.TestCase):
    def test_alanines_and_a_rotamer(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("mutate", tmp)
            universe, counts = _check_common(self, folder, report)
            from itp import parse_itp
            protein = next(d for p in (folder / "toppar").glob("molecule_*.itp")
                           for n, d in parse_itp(str(p)).items() if n == "PROA")
            self.assertEqual({a.resname for a in protein.atoms if a.resnr in (103, 198)}, {"ALA"})
            self.assertAlmostEqual(protein.charge, parse_itp(str(PROA))["PROA"].charge + 1, places=3)
        with tempfile.TemporaryDirectory() as tmp:
            code, folder, report = _build("mutate", tmp, {"mutations": [
                {"molecule": "PROA", "resid": 198, "to": "PHE"}]})
            change = report["mutate"]["changes"][0]
            self.assertEqual(change["rotamer"]["tried"], 12)
            self.assertLess(change["rotamer"]["score"], 5.0)


if __name__ == "__main__":
    unittest.main()
