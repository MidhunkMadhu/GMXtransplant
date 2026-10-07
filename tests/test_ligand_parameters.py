"""A replacement ligand whose bonded parameters are missing must stop the run."""
from pathlib import Path
import tempfile
import unittest

from topology import TopologyError, check_ligand_parameters

FORCEFIELD = """[ bondtypes ]
CT1 CT2 1 0.153 186000.0
[ angletypes ]
CT1 CT2 CT3 5 114.0 488.0 0.0 0.0
[ dihedraltypes ]
X CT2 CT2 X 9 0.0 0.8 3
"""
LIGAND = """[ moleculetype ]
LIG 3
[ atoms ]
1 CT1 1 LIG C1 1 0.0 12.011
2 CT2 1 LIG C2 2 0.0 12.011
3 CT3 1 LIG C3 3 0.0 12.011
4 CT2 1 LIG C4 4 0.0 12.011
[ bonds ]
1 2 1
2 3 1
[ angles ]
1 2 3 5
[ dihedrals ]
1 2 4 3 9
"""


class LigandParameterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "LIG.itp").write_text(LIGAND)

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_bond_is_named(self):
        (self.root / "forcefield.itp").write_text(FORCEFIELD)
        with self.assertRaisesRegex(TopologyError, r"bonds: CT2-CT3.*forcefield.itp"):
            check_ligand_parameters([str(self.root / "LIG.itp")], [str(self.root / "forcefield.itp")])

    def test_reversed_order_and_wildcards_are_accepted(self):
        (self.root / "forcefield.itp").write_text(FORCEFIELD.replace(
            "[ angletypes ]", "CT3 CT2 1 0.153 186000.0\n[ angletypes ]"))
        check_ligand_parameters([str(self.root / "LIG.itp")], [str(self.root / "forcefield.itp")])


EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "ligand_replacement"


@unittest.skipUnless((EXAMPLE / "replacement" / "forcefield.itp").is_file(), "example data not present")
class ExplicitForcefieldTests(unittest.TestCase):
    """forcefield_path selects the ligand force field wherever it is stored."""

    def run_example(self, forcefield):
        import contextlib, io, shutil, warnings, yaml
        from run_pipeline import cli
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(EXAMPLE, root / "src")
            moved = root / "parameters" / "LI2_charmm_gui.itp"
            moved.parent.mkdir()
            shutil.move(root / "src/replacement/forcefield.itp", moved)
            raw = yaml.safe_load((root / "src/ligand_replace.yaml").read_text())
            for key, value in raw["paths"].items():
                if isinstance(value, str) and "/" in value and not value.startswith("/"):
                    raw["paths"][key] = str(root / "src" / value)
            raw["paths"]["new_ligand_forcefield"] = str(moved) if forcefield else ""
            if not forcefield:
                del raw["ligand_replace"]["new_ligand"]["forcefield_path"]
            config = root / "run.yaml"
            config.write_text(yaml.safe_dump(raw))
            out = io.StringIO()
            with warnings.catch_warnings(), contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
                warnings.simplefilter("ignore")
                code = cli(["--mode", "lig", "-i", str(config), "-o", str(root / "out")])
            return code, out.getvalue()

    def test_forcefield_in_another_folder_is_used(self):
        code, log = self.run_example(forcefield=True)
        self.assertEqual(code, 0, log[-2000:])
        self.assertIn("LI2_charmm_gui.itp", log)

    def test_without_ligand_forcefield_the_run_stops(self):
        code, log = self.run_example(forcefield=False)
        self.assertNotEqual(code, 0)
        self.assertIn("uses bonded terms with no parameters", log)


if __name__ == "__main__":
    unittest.main()
