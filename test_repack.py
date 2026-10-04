"""Tests for the GROMACS repacking bundle written by cholesterol restoration."""
import contextlib
import io
from pathlib import Path
import shutil
import tempfile
import unittest

import yaml

from config import ConfigError, load_config
from itp import parse_itp
from repack import make_restrained_cholesterol_itp, restraint_coverage
from topology import parse_top_molecules

ROOT = Path(__file__).resolve().parent
EXAMPLE = ROOT / "examples" / "cholesterol_restoration"
TOPPAR = EXAMPLE / "environment" / "toppar"


class RestraintTests(unittest.TestCase):
    def test_restrained_cholesterol_keeps_parameters_and_holds_every_heavy_atom(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "CHLR.itp"
            path.write_text(make_restrained_cholesterol_itp(str(TOPPAR / "CHL1.itp"), "CHL1", "CHLR"))
            original, restrained = parse_itp(str(TOPPAR / "CHL1.itp"))["CHL1"], parse_itp(str(path))["CHLR"]
            self.assertEqual(original.atoms, restrained.atoms)
            text = path.read_text()
            self.assertNotIn("POSRES_FC_LIPID", text)
            coverage = restraint_coverage(restrained)
            self.assertEqual((coverage["heavy_atoms"], coverage["restrained_heavy_atoms"]), (28, 28))
            self.assertEqual(coverage["restraint_values"], ["POSRES_FC_CHLR"])
            # Bonded terms are untouched: only the header, the name and the POSRES block change.
            source = (TOPPAR / "CHL1.itp").read_text()
            for section in ("[ bonds ]", "[ angles ]", "[ dihedrals ]"):
                self.assertEqual(text.split(section)[1].split("#ifdef")[0].split("[")[0].strip(),
                                 source.split(section)[1].split("#ifdef")[0].split("[")[0].strip())

    def test_coverage_ignores_lone_pairs(self):
        # G4C carries two CGenFF halogen lone pairs (LPH); they are not atoms to restrain.
        coverage = restraint_coverage(parse_itp(str(TOPPAR / "G4C.itp"))["G4C"])
        self.assertEqual((coverage["heavy_atoms"], coverage["restrained_heavy_atoms"]), (30, 30))
        protein = restraint_coverage(parse_itp(str(TOPPAR / "PROA.itp"))["PROA"])
        self.assertEqual(protein["missing_heavy_atoms"], [])
        self.assertEqual(protein["restraint_values"], ["POSRES_FC_BB", "POSRES_FC_SC"])


class RepackConfigTests(unittest.TestCase):
    def _load(self, repack):
        raw = yaml.safe_load((EXAMPLE / "cholesterol_restore.yaml").read_text())
        raw["cholesterol"]["repack"] = repack
        for key, value in raw["paths"].items():
            if isinstance(value, str) and (EXAMPLE / value).exists():
                raw["paths"][key] = str(EXAMPLE / value)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.yaml"
            path.write_text(yaml.safe_dump(raw))
            return load_config(str(path), mode="chl", check_paths=False, output_root=tmp)

    def test_defaults_hold_at_1000_and_end_with_one_ns(self):
        repack = self._load({}).cholesterol.repack
        self.assertTrue(repack.enabled)
        self.assertEqual(repack.lipid_restraint, 0.0)
        self.assertEqual([s.restraint for s in repack.stages], [4000.0, 4000.0, 2000.0, 1000.0, 1000.0])
        last = repack.stages[-1]
        self.assertEqual((last.integrator, last.nsteps * last.dt), ("npt", 1000.0))
        self.assertAlmostEqual(sum(s.nsteps * s.dt for s in repack.stages if s.integrator != "steep"), 1375.0)

    def test_rejects_bad_values(self):
        for repack in ({"unknown": 1}, {"moleculetype": "CHL R"}, {"lipid_restraint": -1},
                       {"stages": []}, {"stages": [{"name": "a", "integrator": "md", "nsteps": 10}]},
                       {"stages": [{"name": "a", "integrator": "npt", "nsteps": 10, "dt": 0.01}]},
                       {"stages": [{"name": "a", "integrator": "npt", "nsteps": 10, "extra": 1}]}):
            with self.subTest(repack=repack), self.assertRaises(ConfigError):
                self._load(repack)


class RepackBundleTests(unittest.TestCase):
    """The cholesterol example, end to end: the bundle next to the restored system."""

    @classmethod
    def setUpClass(cls):
        if shutil.which("obabel") is None:
            raise unittest.SkipTest("Open Babel is needed for the cholesterol example")
        from run_pipeline import main
        cls.tmp = tempfile.TemporaryDirectory()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            cls.code = main(str(EXAMPLE / "cholesterol_restore.yaml"), mode="chl", output_root=cls.tmp.name)
        cls.out = Path(cls.tmp.name)
        cls.bundle = cls.out / "repack"

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_topology_lists_the_restored_cholesterols_as_their_own_type(self):
        self.assertEqual(self.code, 0)
        main = parse_top_molecules(str(self.out / "topol.top"))
        repack = parse_top_molecules(str(self.bundle / "topol.top"))
        self.assertEqual(dict(repack)["CHLR"], 5)
        # Same molecules in the same order; only the five restored CHL1 are renamed.
        self.assertEqual([(("CHL1" if n == "CHLR" else n), c) for n, c in repack], main)
        text = (self.bundle / "topol.top").read_text()
        self.assertIn('#include "toppar/CHLR.itp"', text)
        # Self-contained: every include, the GRO and the index are inside the folder.
        for line in text.splitlines():
            if line.startswith("#include"):
                self.assertTrue((self.bundle / line.split('"')[1]).is_file(), line)
        self.assertEqual((self.bundle / "step5_input.gro").read_bytes(), (self.out / "step5_input.gro").read_bytes())
        self.assertEqual((self.bundle / "index.ndx").read_bytes(), (self.out / "index.ndx").read_bytes())
        self.assertEqual({p.name for p in self.bundle.iterdir()} - {f"{s}.mdp" for s in (
            "step6.0_minimization", "step6.1_equilibration", "step6.2_equilibration", "step6.3_equilibration",
            "step6.4_equilibration")}, {"toppar", "topol.top", "step5_input.gro", "index.ndx", "sample_run.sh"})

    def test_stages_hold_solute_and_free_the_lipids(self):
        last = (self.bundle / "step6.4_equilibration.mdp").read_text()
        self.assertIn("-DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=1000.0 -DPOSRES_FC_LIPID=0.0", last)
        self.assertIn("-DPOSRES_FC_CHLR=1000.0", last)
        self.assertIn("nsteps                  = 500000", last)
        self.assertIn("tc_grps                 = SOLU MEMB SOLV", last)
        self.assertIn("pcoupltype              = semiisotropic", last)
        first_md = (self.bundle / "step6.1_equilibration.mdp").read_text()
        self.assertIn("gen-vel                 = yes", first_md)
        self.assertIn("ref_t                   = 310 310 310", first_md)
        script = (self.bundle / "sample_run.sh").read_text()
        self.assertIn("-r step5_input.gro -p topol.top -n index.ndx", script)
        self.assertEqual(script.count("$GMX grompp"), 5)

    def test_report_mentions_the_sample_run_script(self):
        report = (self.out / "cholesterol_restoration_report.txt").read_text()
        self.assertIn("Repacking bundle:", report)
        self.assertNotIn("WARNING: G4C", report)
        self.assertIn("sample run script is also saved", report)
        self.assertIn("1.375 ns of MD", (self.bundle / "sample_run.sh").read_text())


    def test_report_locates_restored_cholesterols_and_removed_molecules(self):
        import json
        import MDAnalysis as mda
        data = json.loads((self.out / "cholesterol_restoration_report.json").read_text())
        summary = data["restoration_summary"]
        self.assertEqual(len(summary["restored"]), 5)
        gro = mda.Universe(str(self.out / "step5_input.gro"))
        for item in summary["restored"]:
            residue = gro.residues[item["final_residue_number"] - 1]
            self.assertEqual(str(residue.resname), "CHL1")
            self.assertEqual([residue.atoms[0].index + 1, residue.atoms[-1].index + 1], item["atoms"])
        self.assertEqual(sum(len(r["removed_nearby"]) for r in summary["restored"]),
                         summary["removed_nearby"]["count"])
        far = data["composition"]["lipid_residues_removed"]
        self.assertEqual(summary["removed_far"]["count"], len(far))
        self.assertTrue(all(x["distance_to_protein_angstrom"] > 20.0 for x in far))
        text = (self.out / "cholesterol_restoration_report.txt").read_text()
        self.assertIn("Restored cholesterols: 5, final residues", text)
        self.assertIn("removed next to it:", text)


if __name__ == "__main__":
    unittest.main()
