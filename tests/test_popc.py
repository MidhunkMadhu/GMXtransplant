"""Geometry, configuration and end-to-end tests for the popc mode."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import yaml

from charmprot import inspect_system
from config import ConfigError
from popc_membrane import (PopcSpec, _midplane_shift, _min_image_xy, default_reference_path,
                           load_popc_config, load_reference, measure_host, run_popc, tile)
from topology import parse_top_molecules
from tests import slow

HOST = Path(__file__).resolve().parent.parent / "examples" / "addbinder" / "host"


def _write(folder, document):
    path = Path(folder) / "popc.yaml"
    path.write_text(yaml.safe_dump(document))
    return path


class GeometryTests(unittest.TestCase):
    def test_midplane_from_tail_ends_when_planes_are_half_a_box_apart(self):
        # Leaflet P planes 38 A apart in a 76 A box: their circular mean is undefined,
        # the tail ends at the centre are not.
        height, centre = 76.0, 70.0
        tails = np.mod(centre + np.random.default_rng(1).normal(0, 2.0, 200), height)
        shift = _midplane_shift(tails, height)
        self.assertAlmostEqual(float(np.mod(centre + shift, height)), height / 2, delta=0.5)

    def test_minimum_image_in_a_square_cell(self):
        cell = np.array([[100.0, 0.0], [0.0, 100.0]])
        # Two points 98 A apart in x are 2 A apart through the boundary.
        delta = np.array([1.0, 50.0]) - np.array([99.0, 50.0])
        self.assertAlmostEqual(float(np.linalg.norm(_min_image_xy(delta, cell))), 2.0)

    def test_hexagonal_minimum_image(self):
        cell = np.array([[100.0, 0.0], [50.0, 100.0 * np.sqrt(3) / 2]])
        delta = cell[0] + cell[1] - np.array([1.0, 0.0])
        self.assertAlmostEqual(float(np.linalg.norm(_min_image_xy(delta, cell))), 1.0, places=6)

    def test_tiling_covers_the_host_cell_once(self):
        rng = np.random.default_rng(0)
        ref_cell = np.array([[30.0, 0.0], [0.0, 30.0]])
        points = rng.uniform(0, 30, size=(50, 2))
        for host_cell in (np.array([[30.0, 0.0], [0.0, 30.0]]), np.array([[71.0, 0.0], [0.0, 64.0]]),
                          np.array([[60.0, 0.0], [30.0, 60.0 * np.sqrt(3) / 2]])):
            placed = tile(points, ref_cell, host_cell, np.array([7.0, 3.0]))
            xy = np.array([points[k] + offset for k, offset in placed])
            fractional = xy @ np.linalg.inv(host_cell)
            self.assertTrue(np.all((fractional >= 0) & (fractional < 1)))
            density = len(placed) / abs(np.linalg.det(host_cell))
            self.assertAlmostEqual(density, 50 / 900, delta=0.25 * 50 / 900)
            self.assertEqual(len({tuple(np.round(p, 6)) for p in xy}), len(xy))
        same = tile(points, ref_cell, ref_cell, np.zeros(2))
        self.assertEqual(sorted(k for k, _ in same), list(range(50)))


class ConfigTests(unittest.TestCase):
    def test_dry_run_and_bundled_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, {"popc": {"host": str(HOST)}})
            cfg = load_popc_config(path, check_paths=False, output_root=tmp)
            self.assertEqual(cfg.output_dir, str(Path(tmp).absolute() / "popc_output"))
            self.assertIsNone(cfg.reference_gro)
            self.assertEqual(len(cfg.equilibration.stages), 7)
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(run_popc(path, dry_run=True, output_root=tmp), 0)
            self.assertIn("validated", out.getvalue())
        self.assertTrue(default_reference_path().is_file())

    def test_bad_values_are_rejected(self):
        cases = [({"unknown": 1}, "Unknown"), ({"resname_map": {"CLOL": 1}}, "resname_map"),
                 ({"protein_clash_distance": -1}, "protein_clash_distance"),
                 ({"equilibration": {"stages": [{"name": "a", "integrator": "md", "nsteps": 1}]}}, "integrator"),
                 ({"lateral_search": 0}, "lateral_search"), ({"restore_salt": "yes"}, "restore_salt")]
        with tempfile.TemporaryDirectory() as tmp:
            for extra, message in cases:
                with self.subTest(extra=extra):
                    path = _write(tmp, {"popc": {"host": str(HOST), **extra}})
                    with self.assertRaisesRegex(ConfigError, message):
                        load_popc_config(path, check_paths=False, output_root=tmp)

    def test_defaults_and_retired_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, {"popc": {"host": str(HOST), "anchor_count": 3, "anchor_separation": 20.0}})
            cfg = load_popc_config(path, check_paths=False, output_root=tmp)
            self.assertEqual(cfg.retired_settings, ["anchor_count", "anchor_separation"])
            self.assertEqual((cfg.protein_clash_distance, cfg.water_clash_distance), (2.0, 2.0))

    def test_relative_host_is_rejected_outside_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, {"popc": {"host": "host"}})
            with self.assertRaisesRegex(ConfigError, "in full"):
                load_popc_config(path, check_paths=False, output_root=tmp)


class ReferenceTests(unittest.TestCase):
    def test_bundled_bilayer(self):
        host = inspect_system(HOST)
        reference = load_reference(default_reference_path(), host.definitions, PopcSpec(host=str(HOST)), ["POPC"])
        report = reference.report
        self.assertEqual((report["lipids"], report["waters"]), (512, 20480))
        self.assertEqual(report["leaflets"]["upper"]["lipids"], 256)
        self.assertAlmostEqual(report["phosphate_thickness_angstrom"], 39.0, delta=1.0)
        self.assertAlmostEqual(report["area_per_lipid_angstrom2"]["upper"], 64.5, delta=0.5)
        heads = reference.heads
        self.assertTrue(np.all(np.abs(heads[:, 2]) < 30))
        for (_, lipid), side in zip(reference.lipids, reference.leaflets):
            self.assertLess(np.ptp(lipid, axis=0).max(), 35.0)          # whole
            self.assertLess(np.abs(lipid[:, 2]).max(), 32.0)           # inside its bilayer
            self.assertEqual(side, "upper" if lipid[19, 2] > 0 else "lower")


    def test_cholesterol_is_matched_by_its_atom_names(self):
        host = inspect_system(HOST)
        path = Path(__file__).resolve().parent.parent / "examples" / "popc" / "reference" / \
            "popc_chol_70-30_charmm36_303K_100ns.gro.gz"
        reference = load_reference(path, host.definitions, PopcSpec(host=str(HOST)), ["CHL1", "POPC", "POPE"])
        self.assertEqual(reference.report["leaflets"]["upper"]["composition"], {"POPC": 42, "CHL1": 18})
        self.assertIn("CLOL", reference.report["residues_matched"])
        self.assertIn("SOL", reference.report["ignored_residues"])


class HostTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.host = inspect_system(HOST)

    def test_midplane_and_planes(self):
        frame = measure_host(self.host, PopcSpec(host=str(HOST)))
        upper, lower = frame.planes["upper"]["mean_z"], frame.planes["lower"]["mean_z"]
        self.assertAlmostEqual(upper - lower, 37.3, delta=0.5)
        self.assertAlmostEqual(frame.midplane, (upper + lower) / 2, delta=1.0)
        for side in ("upper", "lower"):
            self.assertLess(frame.planes[side]["tilt_degrees"], 1.0)
        removed = sum(sum(c.values()) for c in frame.counts.values())
        self.assertEqual(removed, 105 + 261 + 156)
        given = measure_host(self.host, PopcSpec(host=str(HOST), membrane_center_z=120.0))
        self.assertEqual(given.midplane, 120.0)

    def test_tilted_membrane_fails(self):
        with self.assertRaisesRegex(ConfigError, "tilted"):
            measure_host(self.host, PopcSpec(host=str(HOST), max_tilt=0.01))


@slow
class EndToEndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        # A full host path: an installed package resolves relative paths only in its own examples.
        config = _write(cls.tmp.name, {"popc": {"host": str(HOST)}})
        with contextlib.redirect_stdout(io.StringIO()):
            cls.status = run_popc(config, output_root=cls.tmp.name)
        cls.out = Path(cls.tmp.name) / "popc_output"
        cls.report = json.loads((cls.out / "popc_report.json").read_text())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.final = mda.Universe(str(cls.out / "step5_input.gro"))
            cls.host = mda.Universe(str(HOST / "step5_input.gro"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_completes_with_pure_popc(self):
        self.assertEqual(self.status, 0)
        self.assertEqual(self.report["status"], "complete")
        names = [name for name, _ in parse_top_molecules(str(self.out / "topol.top"))]
        self.assertNotIn("POPE", names)
        self.assertNotIn("CHL1", names)
        self.assertEqual(names[:2], ["PROA", "POPC"])
        self.assertAlmostEqual(self.report["final_net_charge"], 0.0, places=3)

    def test_protein_keeps_its_coordinates(self):
        before = self.host.select_atoms("segid PROA or resname ALA ARG ASN ASP CYS GLN GLU GLY HSD HSE HSP "
                                        "ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL")
        n = len(before)
        np.testing.assert_allclose(self.final.atoms[:n].positions, before.positions, atol=1e-3)

    def test_new_popc_clears_the_protein_and_sits_on_the_host_midplane(self):
        protein = self.final.atoms[:len(self.host.select_atoms("resname ALA ARG ASN ASP CYS GLN GLU GLY HSD "
                                                               "HSE HSP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"))]
        heavy = protein[[not n.startswith("H") for n in protein.names]]
        popc = self.final.select_atoms("resname POPC")
        popc_heavy = popc[[not n.startswith("H") for n in popc.names]]
        pairs = capped_distance(popc_heavy.positions, heavy.positions, 1.95, box=self.final.dimensions,
                                return_distances=False)
        self.assertEqual(len(pairs), 0)
        p = self.final.select_atoms("resname POPC and name P").positions[:, 2]
        midplane = self.report["host_membrane"]["midplane_z"]
        upper, lower = p[p > midplane], p[p <= midplane]
        self.assertAlmostEqual((upper.mean() + lower.mean()) / 2, midplane, delta=1.0)
        placed = self.report["placement"]["placed"]
        self.assertEqual((len(upper), len(lower)), (placed["upper POPC"], placed["lower POPC"]))

    def test_salt_and_water(self):
        salt = self.report["salt"]
        self.assertAlmostEqual(salt["final"]["concentration_molar"], salt["host"]["concentration_molar"],
                               delta=0.005)
        self.assertGreater(self.report["water"]["reference_waters"], 5000)
        water = self.final.select_atoms("resname TIP3 and name OH2")
        popc = self.final.select_atoms("resname POPC")
        heavy = popc[[not n.startswith("H") for n in popc.names]]
        close = capped_distance(water.positions, heavy.positions, 1.5, box=self.final.dimensions,
                                return_distances=False)
        self.assertEqual(len(close), 0)

    def test_equilibration_folder(self):
        folder = self.out / "equilibration"
        for name in ("topol.top", "step5_input.gro", "index.ndx", "sample_run.sh",
                     "step6.0_minimization.mdp", "step6.6_equilibration.mdp"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertTrue((folder / "toppar").is_dir())
        first = (folder / "step6.1_equilibration.mdp").read_text()
        self.assertIn("-DPOSRES_FC_BB=4000.0", first)
        self.assertIn("gen-vel                 = yes", first)
        self.assertIn("pcoupltype              = semiisotropic", (folder / "step6.6_equilibration.mdp").read_text())
        self.assertEqual(self.report["equilibration"]["warnings"], [])
        index = (self.out / "index.ndx").read_text()
        for group in ("[ SOLU_MEMB ]", "[ SOLV ]", "[ POPC_upper ]", "[ POPC_lower ]", "[ Upper_leaflet ]"):
            self.assertIn(group, index)


if __name__ == "__main__":
    unittest.main()
