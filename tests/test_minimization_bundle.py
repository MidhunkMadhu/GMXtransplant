import contextlib
from dataclasses import replace
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from config import ConfigError, MinimizationSpec, load_minimization_config
from minimization_bundle import BundleError, prepare_minimization, validate_preparation
from tests.test_minimization import TOP, GRO


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "input.gro").write_text(GRO.replace("1.00000", "4.00000"))
        (self.root / "topol.top").write_text(TOP)
        self.spec = MinimizationSpec(coordinates_path=str(self.root / "input.gro"),
            topology_path=str(self.root / "topol.top"), output_dir=str(self.root / "bundle"),
            constraints="none", rigid_water=False, platform="CPU", tolerance_kj_mol_nm=1)

    def prepare(self, **kwargs):
        return prepare_minimization(replace(self.spec, **kwargs))

    def test_independent_inputs_manifest_hashes_and_no_execution(self):
        with patch("subprocess.run", side_effect=AssertionError("Engine execution forbidden")):
            result = self.prepare()
        self.assertEqual(result["status"], "prepared_not_run")
        bundle = Path(result["output_dir"])
        manifest = json.loads((bundle / "manifest.json").read_text())
        for name, checksum in manifest["files_sha256"].items():
            self.assertEqual(hashlib.sha256((bundle / name).read_bytes()).hexdigest(), checksum)
        self.assertFalse((bundle / "gromacs").exists())
        self.assertFalse((bundle / "openmm").exists())
        folder = bundle
        self.assertEqual((folder / "input.gro").read_bytes(), (self.root / "input.gro").read_bytes())
        self.assertTrue((folder / "topol.top").is_file())
        self.assertTrue((folder / "toppar").is_dir())
        self.assertFalse((folder / "results").exists())
        for script in list(folder.glob("*.sh")) + list(folder.glob("*.slurm")):
            subprocess.run(["bash", "-n", str(script)], check=True)
        slurm = (folder / "submit_cpu.slurm").read_text()
        self.assertIn("--cpus-per-task=8", slurm)
        self.assertNotIn("--account", slurm)
        self.assertNotIn("--partition", slurm)

    def test_no_overwrite(self):
        self.prepare()
        with self.assertRaisesRegex(BundleError, "overwrite"):
            self.prepare()

    def test_mismatch_leaves_no_bundle_or_lock(self):
        (self.root / "input.gro").write_text(GRO.replace("C1", "XX"))
        with self.assertRaisesRegex(BundleError, "identity mismatch"):
            self.prepare()
        self.assertFalse((self.root / "bundle").exists())
        self.assertFalse((self.root / ".bundle.lock").exists())

    def test_count_mismatch(self):
        (self.root / "topol.top").write_text(TOP.replace("MOL  1", "MOL  2"))
        with self.assertRaisesRegex(BundleError, "more atoms"):
            self.prepare()

    def test_nested_relative_includes_duplicate_basenames_and_relocation(self):
        ff = self.root / "external"
        (ff / "a").mkdir(parents=True)
        (ff / "b").mkdir()
        (ff / "a/common.itp").write_text("; a\n")
        (ff / "b/common.itp").write_text("; b\n")
        (ff / "forcefield.itp").write_text('#include "a/common.itp"\n#include "b/common.itp"\n')
        (self.root / "topol.top").write_text('#include "forcefield.itp"\n' + TOP)
        self.prepare(include_dir=str(ff))
        copied = self.root / "copied"
        shutil.copytree(self.root / "bundle", copied)
        from minimization_bundle import _audit
        self.assertEqual(_audit(copied, {}, "relocated")["atom_count"], 2)
        self.assertEqual(len(list((copied / "toppar").glob("*common.itp"))), 2)
        self.assertNotIn(str(self.root), (copied / "topol.top").read_text())

    def test_missing_even_inactive_include_rejected(self):
        (self.root / "topol.top").write_text('#ifdef UNUSED\n#include "missing.itp"\n#endif\n' + TOP)
        with self.assertRaisesRegex(BundleError, "Unresolved include"):
            self.prepare()

    def test_cycle_rejected(self):
        (self.root / "topol.top").write_text('#include "topol.top"\n' + TOP)
        with self.assertRaisesRegex(BundleError, "Cyclic"):
            self.prepare()

    def test_active_defines_control_atom_identity(self):
        (self.root / "topol.top").write_text(TOP.replace("[ atoms ]", "#ifdef WRONG\n#define ATOM XX\n#else\n#define ATOM C1\n#endif\n[ atoms ]").replace("MOL  C1", "MOL  ATOM"))
        self.prepare()
        with self.assertRaisesRegex(BundleError, "identity mismatch"):
            self.prepare(output_dir=str(self.root / "bad"), defines={"WRONG": 1})

    def test_topology_position_restraints_are_rejected(self):
        top = TOP.replace("[ system ]", "#ifdef POSRES\n[ position_restraints ]\n1 1 100 100 100\n#endif\n[ system ]")
        (self.root / "topol.top").write_text(top)
        self.prepare()
        with self.assertRaisesRegex(BundleError, "position_restraints"):
            self.prepare(output_dir=str(self.root / "bad"), defines={"POSRES": 1})

    def test_openmm_implicit_flexible_water_definition_is_preserved(self):
        top = TOP.replace("[ system ]", "#ifdef FLEXIBLE\n; flexible water branch\n[ angles ]\n#else\n[ settles ]\n#endif\n[ system ]")
        (self.root / "topol.top").write_text(top)
        self.prepare()
        expanded = (self.root / "bundle/audit_openmm.top").read_text()
        self.assertIn("[ angles ]", expanded)
        self.assertNotIn("[ settles ]", expanded)

    def test_openmm_identity_uses_raw_source_names(self):
        try:
            from openmm import app
        except ImportError:
            self.skipTest("OpenMM not installed")
        from minimization import _validate_coordinate_topology_identity
        top = TOP.replace("MOL  C1", "SER  HN").replace("MOL  C2", "SER  CA")
        gro = GRO.replace("MOL", "SER").replace("C1", "HN").replace("C2", "CA")
        (self.root / "input.gro").write_text(gro)
        (self.root / "topol.top").write_text(top)
        parsed = app.GromacsTopFile(str(self.root / "topol.top"), includeDir=str(self.root))
        self.assertEqual(list(parsed.topology.atoms())[0].name, "H")
        _validate_coordinate_topology_identity(str(self.root / "input.gro"), parsed)

    def test_unsupported_preprocessor_rejected(self):
        (self.root / "topol.top").write_text("#if 1\n" + TOP + "#endif\n")
        with self.assertRaisesRegex(BundleError, "Unsupported #if"):
            self.prepare()

    def test_resources_and_injection_are_validated(self):
        for changes in ({"resources": {"cpus_per_task": 0}},
                        {"resources": {"unknown": 1}},
                        {"resources": {"memory": "8G\nmalicious"}},
                        {"defines": {"X": "$(command)"}},
                        {"restraint_force_constant_kj_mol_nm2": float("nan")},
                        {"restraint_force_constant_kj_mol_nm2": 0},
                        {"restraint_residue_classes": {"CUSTOM": "typo"}},
                        {"max_iterations": 0}):
            with self.subTest(changes=changes), self.assertRaises((BundleError, ConfigError)):
                validate_preparation(replace(self.spec, **changes))

    def test_old_execution_cli_is_removed(self):
        from run_pipeline import cli
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
            cli(["--minimize"])
        self.assertEqual(caught.exception.code, 2)

    def test_existing_cli_and_dry_run(self):
        from run_pipeline import cli
        with contextlib.redirect_stdout(io.StringIO()):
            result = cli(["prepare-minimization", "--coordinates", self.spec.coordinates_path,
                          "--topology", self.spec.topology_path, "--output", self.spec.output_dir])
        self.assertEqual(result, 0)
        with contextlib.redirect_stdout(io.StringIO()):
            result = cli(["prepare-minimization", "--coordinates", "/missing/input.gro",
                          "--topology", "/missing/topol.top", "--dry-run"])
        self.assertEqual(result, 0)

    def test_integrated_prepare_uses_fresh_outputs_only(self):
        from run_pipeline import _post_pipeline_minimization_spec
        from config import Config, TopologySpec, OutputSpec
        cfg = Config(topology=TopologySpec(enabled=True, output_dir="new"),
                     output=OutputSpec(gro_path="new/final.gro"),
                     minimization=replace(self.spec, topology_path="old.top"))
        result = _post_pipeline_minimization_spec(cfg)
        self.assertEqual(result.coordinates_path, "new/final.gro")
        self.assertEqual(result.topology_path, "new/topol.top")
        self.assertEqual(result.include_dir, "new")

    def test_generated_openmm_runner_after_relocation(self):
        try:
            import openmm
        except ImportError:
            self.skipTest("OpenMM not installed")
        self.prepare(nonbonded_method="CutoffPeriodic", switch_distance_nm=None)
        copied = self.root / "copied openmm"
        shutil.copytree(self.root / "bundle", copied)
        env = dict(os.environ, PYTHONPATH="", OPENMM_CPU_THREADS="2")
        result = subprocess.run([sys.executable, str(copied / "minimize_openmm.py"), "--platform", "CPU"],
                                cwd=self.root, env=env, text=True, capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        report = json.loads((copied / "results/minimization_report.json").read_text())
        self.assertTrue(report["converged"])
        self.assertEqual(report["openmm"]["platform"], "CPU")
        self.assertEqual(report["openmm"]["properties"]["Threads"], "2")
        self.assertTrue((copied / "results/minimized.gro").is_file())
        rerun = subprocess.run([sys.executable, str(copied / "minimize_openmm.py")],
                               capture_output=True, timeout=60)
        self.assertNotEqual(rerun.returncode, 0)



if __name__ == "__main__":
    unittest.main()
