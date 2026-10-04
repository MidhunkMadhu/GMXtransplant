"""Automatic maximum-common-substructure placement (fit.method: mcsfit)."""
from pathlib import Path
import unittest

import numpy as np

from substructure import LigandGraph, SubstructureError, maximum_common_substructures, best_substructure_fit


def graph(elements, bonds, names=None):
    names = names or [f"{e}{i + 1}" for i, e in enumerate(elements)]
    return LigandGraph(names, elements, set(bonds), "test")


RING = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]


class Atoms:
    def __init__(self, names, positions):
        self.names, self.positions = np.array(names), np.asarray(positions, dtype=float)


class SubstructureTests(unittest.TestCase):
    def test_largest_shared_core_and_symmetric_alternatives(self):
        # Toluene-like (ring + methyl) against a ring with an amine chain: the
        # shared part is the ring plus the attached carbon, with two ring directions.
        tol = graph(["C"] * 7, RING + [(0, 6)])
        pea = graph(["C"] * 8 + ["N"], RING + [(0, 6), (6, 7), (7, 8)])
        mappings = maximum_common_substructures(tol, pea)
        self.assertTrue(all(len(m) == 7 for m in mappings))
        self.assertEqual(len(mappings), 2)

    def test_bonds_must_match_not_just_elements(self):
        chain = graph(["C"] * 6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
        ring = graph(["C"] * 6, RING)
        # A ring shares only an open 5-atom path with a chain, never all 6 atoms.
        self.assertEqual(len(maximum_common_substructures(chain, ring)[0]), 5)

    def test_lowest_rmsd_mapping_is_chosen(self):
        angles = np.arange(6) * np.pi / 3
        ring = np.c_[np.cos(angles), np.sin(angles), np.zeros(6)] * 1.4
        old_xyz = np.vstack([ring, [[2.9, 0, 0]]])
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        new_xyz = old_xyz @ rotation.T + 5.0
        old = graph(["C"] * 7, RING + [(0, 6)], [f"A{i}" for i in range(7)])
        new = graph(["C"] * 7, RING + [(0, 6)], [f"B{i}" for i in range(7)])
        old_names, new_names, rmsd, tried = best_substructure_fit(
            Atoms(old.names, old_xyz), Atoms(new.names, new_xyz), old, new)
        self.assertLess(rmsd, 1e-6)
        self.assertEqual(new_names, [n.replace("A", "B") for n in old_names])
        self.assertGreaterEqual(tried, 2)

    def test_no_shared_core_is_reported(self):
        old = graph(["S", "S"], [(0, 1)])
        new = graph(["C", "C", "C"], [(0, 1), (1, 2)])
        with self.assertRaisesRegex(SubstructureError, "fewer than three"):
            best_substructure_fit(Atoms(old.names, np.eye(2, 3)), Atoms(new.names, np.eye(3)), old, new)


EXAMPLE = Path(__file__).resolve().parent / "examples" / "ligand_replacement"


@unittest.skipUnless((EXAMPLE / "ligand_replace.yaml").is_file(), "example data not present")
class ExampleMcsfitTests(unittest.TestCase):
    def test_example_finds_the_same_core_as_the_explicit_pairs(self):
        import contextlib, io, tempfile, warnings, yaml
        from run_pipeline import cli
        with tempfile.TemporaryDirectory() as tmp:
            raw = yaml.safe_load((EXAMPLE / "ligand_replace.yaml").read_text())
            explicit = set(zip(raw["ligand_replace"]["fit"]["old_ligand_fit_atoms"],
                               raw["ligand_replace"]["fit"]["new_ligand_fit_atoms"]))
            for key, value in raw["paths"].items():
                if isinstance(value, str) and "/" in value and not value.startswith("/"):
                    raw["paths"][key] = str(EXAMPLE / value)
            raw["ligand_replace"]["fit"] = {"method": "mcsfit"}
            config = Path(tmp) / "run.yaml"
            config.write_text(yaml.safe_dump(raw))
            out = io.StringIO()
            with warnings.catch_warnings(), contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
                warnings.simplefilter("ignore")
                code = cli(["--mode", "lig", "-i", str(config), "-o", str(Path(tmp) / "out")])
            log = out.getvalue()
            self.assertEqual(code, 0, log[-2000:])
            self.assertIn("mcsfit: the largest common substructure has 9 heavy atoms", log)
            self.assertIn("RMSD over 9 heavy atom(s): 0.399 A", log)
            import json
            report = json.loads(next((Path(tmp) / "out").glob("*_report.json")).read_text())
            pairs = {tuple(p) for p in report["ligand_replacement"]["fit_atom_pairs"]}
            self.assertEqual(pairs, explicit)


if __name__ == "__main__":
    unittest.main()
