import json
import tempfile
import unittest
from pathlib import Path

from config import MinimizationSpec
from minimization import run_minimization


TOP = """
[ defaults ]
1 2 yes 0.5 0.833333

[ atomtypes ]
CT  12.011  0.0  A  0.340  0.276144

[ moleculetype ]
MOL  3

[ atoms ]
1  CT  1  MOL  C1  1  0.0  12.011
2  CT  1  MOL  C2  1  0.0  12.011

[ bonds ]
1  2  1  0.100  100000.0

[ system ]
Two atom minimization smoke test

[ molecules ]
MOL  1
"""

GRO = """Two atom minimization smoke test
    2
    1MOL     C1    1   0.400   0.500   0.500
    1MOL     C2    2   0.600   0.500   0.500
   1.00000   1.00000   1.00000
"""


class OpenMMMinimizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import openmm  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("OpenMM optional dependency is not installed")

    def test_minimizes_gromacs_bundle_and_writes_audited_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gro = root / "input.gro"
            top = root / "topol.top"
            output = root / "minimized.gro"
            report = root / "minimization_report"
            gro.write_text(GRO, encoding="utf-8")
            top.write_text(TOP, encoding="utf-8")
            spec = MinimizationSpec(
                enabled=True,
                coordinates_path=str(gro),
                topology_path=str(top),
                output_gro_path=str(output),
                report_path=str(report),
                nonbonded_method="NoCutoff",
                switch_distance_nm=None,
                constraints="none",
                rigid_water=False,
                platform="Reference",
                precision="double",
                tolerance_kj_mol_nm=1.0,
                max_iterations=5000,
            )

            result = run_minimization(spec)

            self.assertTrue(result["converged"])
            self.assertTrue(output.is_file())
            self.assertTrue((root / "minimization_report.txt").is_file())
            json_path = root / "minimization_report.json"
            self.assertTrue(json_path.is_file())
            written = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertEqual(written["system"]["particles"], 2)
            self.assertLessEqual(
                written["minimizer_progress"]["last_objective_rms_kj_mol_nm"],
                written["settings"]["tolerance_kj_mol_nm"],
            )
            self.assertLess(
                written["after"]["potential_energy_kj_mol"],
                written["before"]["potential_energy_kj_mol"],
            )


if __name__ == "__main__":
    unittest.main()
