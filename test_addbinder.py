"""Geometry, configuration and end-to-end tests for the addbinder mode."""
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

from addbinder import (AddBinderSpec, MembraneFrame, PoseSpec, _orientation_matrix,
                       _unwrap_sequential, approach_direction, check_pose, load_addbinder_config,
                       place_binder, run_addbinder)
from config import ConfigError
from scipy.spatial.transform import Rotation
from topology import parse_top_molecules

EXAMPLE = Path(__file__).resolve().parent / "examples" / "addbinder"
GPROTEIN = Path(__file__).resolve().parent / "examples" / "cholesterol_restoration" / "environment"


def _frame(lengths=(80.0, 80.0, 100.0)):
    return MembraneFrame(np.array(lengths), 0.0, 70.0, 30.0, "test", (0.0, 0.0))


class GeometryTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        # A flat slab-like binder: long in one direction, thin in another.
        self.binder = rng.normal(size=(30, 3)) * [6.0, 2.0, 0.5] @ np.array(
            [[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
        self.host = np.array([[40.0, 40.0, 50.0], [40.0, 40.0, 80.0], [45.0, 40.0, 75.0]])
        self.tip = {"index": 1, "position": self.host[1]}
        self.cfg = AddBinderSpec(host="h", binder_coordinates="b", binder_itp="i")

    def test_flat_puts_thinnest_axis_on_normal(self):
        pose = PoseSpec("p", 10.0, orientation="flat")
        rotation = _orientation_matrix(pose, self.binder - self.binder.mean(0))
        extent = np.ptp((self.binder - self.binder.mean(0)) @ rotation.T, axis=0)
        self.assertEqual(int(np.argmax(extent)), 0)
        self.assertEqual(int(np.argmin(extent)), 2)
        self.assertAlmostEqual(np.linalg.det(rotation), 1.0)

    def test_end_on_points_the_far_end_at_the_host(self):
        # A rod with a heavy blob at one end: the bare tail is the far-reaching end.
        rod = np.array([[x, 0.0, 0.0] for x in (-1.0, -0.8, -0.6, -0.4, 0.0, 1.5, 3.0, 4.5)])
        rod = rod @ Rotation.from_euler("zyz", [30, 60, 10], degrees=True).as_matrix().T
        centred = rod - rod.mean(0)
        for side, sign in (("upper", -1), ("lower", 1)):
            rotation = _orientation_matrix(PoseSpec("p", 10.0, orientation="end_on"), centred, side)
            z = (centred @ rotation.T)[:, 2]
            self.assertAlmostEqual(np.ptp(z), np.ptp(centred @ np.linalg.svd(centred)[2][0]), places=6)
            self.assertEqual(int(np.argmax(sign * z)), len(rod) - 1)
            flipped = _orientation_matrix(PoseSpec("p", 10.0, orientation="end_on", flip=True), centred, side)
            self.assertEqual(int(np.argmin(sign * (centred @ flipped.T)[:, 2])), len(rod) - 1)
            self.assertAlmostEqual(np.linalg.det(flipped), 1.0)

    def test_approach_stops_at_distance_from_nearest_host_atom(self):
        from scipy.spatial import cKDTree
        for side in ("upper", "lower"):
            self.cfg.side = side
            pose = PoseSpec("p", 12.0, orientation="flat", approach=[40.0, 135.0])
            placed, _ = place_binder(pose, self.cfg, _frame(), self.host, self.binder, np.ones(30, bool), self.tip)
            self.assertAlmostEqual(float(cKDTree(self.host).query(placed)[0].min()), 12.0, places=4)
            # The binder centre lies on the ray from the tip along the approach direction.
            offset = placed.mean(0) - self.tip["position"]
            direction = approach_direction(pose, side)
            self.assertAlmostEqual(float(np.linalg.norm(np.cross(offset, direction))), 0.0, places=6)
            self.assertGreater(float(offset @ direction), 0.0)
            self.assertAlmostEqual(float(np.degrees(np.arccos(abs(direction[2])))), 40.0)

    def test_centroid_distance_is_spherical_coordinates_from_the_tip(self):
        for side, sign in (("upper", 1), ("lower", -1)):
            self.cfg.side = side
            pose = PoseSpec("p", 25.0, orientation="flat", approach=[30.0, 90.0], distance_to="centroid")
            placed, _ = place_binder(pose, self.cfg, _frame(), self.host, self.binder, np.ones(30, bool), self.tip)
            expected = self.tip["position"] + 25.0 * np.array([0.0, np.sin(np.radians(30)), sign * np.cos(np.radians(30))])
            self.assertTrue(np.allclose(placed.mean(0), expected))

    def test_centroid_puts_the_heavy_atom_centre_there(self):
        frame = MembraneFrame(np.array([80.0, 80.0, 100.0]), 7.0, 70.0, 30.0, "test", (0.0, 0.0))
        pose = PoseSpec("p", 10.0, orientation="as_is", centroid=[12.0, 30.0, 65.0])
        placed, _ = place_binder(pose, self.cfg, frame, self.host, self.binder, np.ones(30, bool), self.tip)
        # place_binder works in the membrane-centred frame (z + shift).
        self.assertTrue(np.allclose(placed.mean(0), [12.0, 30.0, 72.0]))

    def test_distance_is_measured_from_tip_plane(self):
        for side, sign in (("upper", 1), ("lower", -1)):
            self.cfg.side = side
            tip = {"index": 0, "position": np.array([40.0, 40.0, 80.0 if side == "upper" else 20.0])}
            placed, _ = place_binder(PoseSpec("p", 10.0), self.cfg, _frame(), self.host,
                                     self.binder, np.ones(30, bool), tip)
            edge = placed[:, 2].min() if side == "upper" else placed[:, 2].max()
            self.assertAlmostEqual(edge, tip["position"][2] + sign * 10.0, places=6)
            self.assertTrue(np.allclose(placed[:, :2].mean(0), [40.0, 40.0]))

    def test_periodic_image_of_host_rejects_pose(self):
        # Box only 100 A tall: the host bottom (z = 50) repeats at z = 150.
        placed = np.array([[40.0, 40.0, 140.0]])
        checks, failures = check_pose(PoseSpec("p", 10.0), self.cfg, _frame(), self.host,
                                      np.zeros((0, 3)), placed)
        self.assertAlmostEqual(checks["gap_to_host_image_angstrom"], 10.0, places=3)
        self.assertEqual(failures, [])
        self.cfg.min_image_gap = 12.0
        _, failures = check_pose(PoseSpec("p", 10.0), self.cfg, _frame(), self.host,
                                 np.zeros((0, 3)), placed)
        self.assertEqual(len(failures), 1)
        self.assertIn("more box height", failures[0])

    def test_unwrap_joins_molecule_split_by_box(self):
        lengths = np.array([20.0, 20.0, 20.0])
        chain = np.array([[18.5, 5.0, 5.0], [19.8, 5.0, 5.0], [1.1, 5.0, 5.0], [2.4, 5.0, 5.0]])
        whole = _unwrap_sequential(chain, lengths)
        self.assertTrue(np.allclose(np.diff(whole[:, 0]), 1.3))


class ConfigTests(unittest.TestCase):
    def _load(self, extra, poses=None):
        data = {"addbinder": {"host": "/abs/host", "binder_coordinates": "/abs/b.pdb",
                              "binder_itp": "/abs/b.itp", **extra}}
        if poses is not None:
            data["addbinder"]["poses"] = poses
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.yaml"
            path.write_text(yaml.safe_dump(data))
            return load_addbinder_config(str(path), check_paths=False, output_root=tmp)

    def test_defaults_give_one_flat_pose(self):
        cfg = self._load({})
        self.assertEqual([(p.name, p.distance, p.orientation) for p in cfg.poses],
                         [("binderpose1", 10.0, "flat")])

    def test_auto_poses_get_different_orientations(self):
        cfg = self._load({}, [{"distance": 10}, {"distance": 15}, {}, {}])
        self.assertEqual([(p.orientation, p.flip) for p in cfg.poses],
                         [("flat", False), ("end_on", False), ("end_on", True), ("edge", False)])
        self.assertEqual([p.distance for p in cfg.poses], [10.0, 15.0, 10.0, 10.0])

    def test_same_orientation_gives_every_auto_pose_the_first(self):
        cfg = self._load({"same_orientation": True}, [{}, {"distance": 20}, {}])
        self.assertEqual({(p.orientation, p.flip, p.spin) for p in cfg.poses}, {("flat", False, 0.0)})

    def test_auto_skips_orientations_already_chosen_and_cycles_with_spin(self):
        cfg = self._load({}, [{"orientation": "flat"}, {}] + [{}] * 5)
        self.assertEqual([(p.orientation, p.flip, p.spin) for p in cfg.poses[:3]],
                         [("flat", False, 0.0), ("end_on", False, 0.0), ("end_on", True, 0.0)])
        self.assertEqual((cfg.poses[-1].orientation, cfg.poses[-1].spin), ("end_on", 45.0))

    def test_random_pose_settings(self):
        cfg = self._load({"distance": 15})
        self.assertEqual((cfg.random_poses, cfg.random_distance, cfg.random_seed), (3, [15.0, 20.0], 1))
        for extra in ({"random_poses": -1}, {"random_poses": 1.5}, {"random_distance": [25, 20]},
                      {"random_distance": [20]}, {"random_max_tilt": 90}, {"random_min_angle": -1}):
            with self.subTest(extra=extra), self.assertRaises(ConfigError):
                self._load(extra)
        with self.assertRaisesRegex(ConfigError, "random poses"):
            self._load({}, [{"name": "binderpose3"}])

    def test_centroid_and_approach_poses(self):
        cfg = self._load({}, [{"centroid": [60, 57, 180]}, {"approach": [30, 90], "distance": 22}])
        self.assertEqual(cfg.poses[0].centroid, [60.0, 57.0, 180.0])
        self.assertEqual((cfg.poses[1].approach, cfg.poses[1].distance), ([30.0, 90.0], 22.0))
        for pose in ({"centroid": [1, 2]}, {"centroid": [1, 2, 3], "distance": 20},
                     {"centroid": [1, 2, 3], "approach": [10, 0]}, {"approach": [95, 0]},
                     {"centroid": [1, 2, 3], "lateral_offset": [1, 0]}):
            with self.subTest(pose=pose), self.assertRaises(ConfigError):
                self._load({}, [pose])

    def test_distance_options(self):
        cfg = self._load({"random_rotated": 1, "random_from_input_position": True},
                         [{"distance": 20, "reduce_distance_by": 2, "min_distance": 10},
                          {"distance": 40, "distance_to": "centroid", "approach": [20, 0],
                           "orientation": "euler", "angles": [10, 20, 30]}])
        self.assertEqual((cfg.poses[0].reduce_distance_by, cfg.poses[0].min_distance), (2.0, 10.0))
        self.assertEqual(cfg.poses[1].distance_to, "centroid")
        self.assertEqual(cfg.random_rotated, 1)
        for extra in ({"random_rotated": 4}, {"random_rotated": -1}, {"random_rotated": 1.5}):
            with self.subTest(extra=extra), self.assertRaises(ConfigError):
                self._load(extra)
        for pose in ({"distance_to": "middle"}, {"reduce_distance_by": -1},
                     {"from_input_position": True, "distance_to": "centroid"},
                     {"centroid": [1, 2, 3], "reduce_distance_by": 2}):
            with self.subTest(pose=pose), self.assertRaises(ConfigError):
                self._load({}, [pose])

    def test_from_input_position_keeps_the_input_pose(self):
        cfg = self._load({}, [{"from_input_position": True, "distance": 12}])
        self.assertEqual((cfg.poses[0].from_input_position, cfg.poses[0].orientation), (True, "as_is"))
        for pose in ({"from_input_position": "yes"}, {"from_input_position": True, "centroid": [1, 2, 3]},
                     {"from_input_position": True, "lateral_offset": [1, 0]}):
            with self.subTest(pose=pose), self.assertRaises(ConfigError):
                self._load({}, [pose])

    def test_rejects_bad_values(self):
        for extra, poses in (({"side": "top"}, None), ({"unknown": 1}, None),
                             ({"distance": -1}, None), ({}, [{"name": "a"}, {"name": "a"}]),
                             ({}, [{"angles": [0, 90, 0]}]), ({}, [{"orientation": "random"}]),
                             ({}, [{"flip": True}]), ({}, [{"orientation": "flat", "flip": 1}]),
                             ({"same_orientation": "yes"}, None), ({"salt_concentration": "lots"}, None),
                             ({"salt_concentration": -0.1}, None)):
            with self.subTest(extra=extra, poses=poses), self.assertRaises(ConfigError):
                self._load(extra, poses)

    def test_salt_concentration(self):
        self.assertEqual(self._load({}).salt_concentration, "host")
        self.assertEqual(self._load({"salt_concentration": "0.15"}).salt_concentration, 0.15)
        self.assertEqual(self._load({}).ion_side, "any")

    def test_binder_itp_accepts_a_list(self):
        self.assertEqual(self._load({}).binder_itp, ["/abs/b.itp"])
        cfg = self._load({"binder_itp": ["/abs/a.itp", "/abs/b.itp"]})
        self.assertEqual(cfg.binder_itp, ["/abs/a.itp", "/abs/b.itp"])
        for value in ([], ["/abs/a.itp", "/abs/a.itp"], ["/abs/a.itp", ""], 3):
            with self.subTest(value=value), self.assertRaises(ConfigError):
                self._load({"binder_itp": value})

    def test_relative_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.yaml"
            path.write_text(yaml.safe_dump({"addbinder": {"host": "host", "binder_coordinates": "b.pdb",
                                                          "binder_itp": "b.itp"}}))
            with self.assertRaises(ConfigError):
                load_addbinder_config(str(path), check_paths=False, output_root=tmp)


class ExampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        with contextlib.redirect_stdout(io.StringIO()):
            cls.code = run_addbinder(str(EXAMPLE / "addbinder.yaml"), output_root=cls.tmp.name)
        cls.out = Path(cls.tmp.name) / "addbinder_output"

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_summary(self):
        self.assertEqual(self.code, 0)
        status = {row["pose"]: row["status"] for row in json.loads((self.out / "summary.json").read_text())}
        self.assertEqual(status, {name: "accepted" for name in ("binderpose1", "binderpose2", "binderpose3", "binderpose4")})

    def test_random_poses_differ_in_direction_position_and_distance(self):
        poses, centres = [], []
        for name in ("binderpose1", "binderpose2", "binderpose3", "binderpose4"):
            report = json.loads((self.out / name / "addbinder_report.json").read_text())
            pose = PoseSpec(**report["pose"])
            poses.append(pose)
            ligand = self._universe(self.out / name).select_atoms("resname LDP and not name H*")
            centres.append(ligand.positions.mean(0))
            if name != "binderpose1":
                self.assertEqual(report["random"]["seed"], 1)
                self.assertTrue(20.0 <= pose.distance <= 25.0)
                self.assertLessEqual(pose.approach[0], 60.0)
                self.assertGreaterEqual(report["checks"]["gap_to_host_angstrom"], pose.distance - 0.02)
        lengths = self._universe(self.out / "binderpose1").dimensions[:3]
        for i in range(4):
            for j in range(i):
                a, b = approach_direction(poses[i], "upper"), approach_direction(poses[j], "upper")
                self.assertGreaterEqual(np.degrees(np.arccos(np.clip(a @ b, -1, 1))), 25.0 - 1e-6)
                delta = centres[i] - centres[j]
                delta -= lengths * np.round(delta / lengths)
                self.assertGreaterEqual(np.linalg.norm(delta), 10.0 - 0.05)

    def test_viewer_scene_has_every_pose(self):
        for name in ("view.pml", "view.vmd"):
            self.assertTrue((self.out / name).is_file())
        scene = json.loads((self.out / "visualization" / "scene.json").read_text())
        names = {o["name"] for o in scene["objects"]}
        self.assertLessEqual({"retained_protein", "host_tip", "binder_binderpose1", "binder_binderpose2",
                              "binder_binderpose3", "binder_binderpose4"}, names)

    def test_rejected_pose_writes_only_report_and_centroid_pose_is_exact(self):
        raw = yaml.safe_load((EXAMPLE / "addbinder.yaml").read_text())["addbinder"]
        for key in ("host", "binder_coordinates", "binder_itp", "binder_forcefield"):
            raw[key] = str(EXAMPLE / raw[key])
        # The free water column above D1R is about 158 A; 170 A above the tip
        # leaves < 10 A to the periodic image of the receptor.
        raw["poses"] = [{"name": "too_far", "distance": 170.0},
                        {"name": "fixed", "orientation": "as_is", "centroid": [70.0, 50.0, 185.0]},
                        {"name": "shrinks", "distance": 170.0, "reduce_distance_by": 20.0, "min_distance": 10.0}]
        raw["random_poses"] = 0
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.yaml"
            path.write_text(yaml.safe_dump({"addbinder": raw}))
            with contextlib.redirect_stdout(io.StringIO()):
                code = run_addbinder(str(path), output_root=tmp)
            out = Path(tmp) / "addbinder_output"
            self.assertEqual(code, 0)
            shrunk = json.loads((out / "shrinks" / "addbinder_report.json").read_text())
            self.assertEqual(shrunk["status"], "accepted")
            self.assertEqual([a["distance"] for a in shrunk["distance_attempts"]], [170.0, 150.0])
            self.assertTrue(shrunk["distance_attempts"][0]["failures"])
            self.assertEqual(shrunk["pose"]["distance"], 150.0)
            ligand = self._universe(out / "fixed").select_atoms("resname LDP and not name H*")
            self.assertTrue(np.allclose(ligand.positions.mean(0), [70.0, 50.0, 185.0], atol=0.01))
            files = {p.name for p in (out / "too_far").iterdir()}
            self.assertEqual(files, {"addbinder_report.json", "addbinder_report.txt"})
            report = json.loads((out / "too_far" / "addbinder_report.json").read_text())
            self.assertIn("periodic image", report["failures"][0])

    def test_accepted_pose_is_consistent(self):
        folder = self.out / "binderpose1"
        report = json.loads((folder / "addbinder_report.json").read_text())
        self.assertEqual(report["final_net_charge"], 0.0)
        self.assertGreaterEqual(report["checks"]["gap_to_host_angstrom"], 20.0 - 0.02)
        self.assertGreaterEqual(report["checks"]["gap_to_host_image_angstrom"], 10.0 - 0.02)
        molecules = dict(parse_top_molecules(str(folder / "topol.top")))
        self.assertEqual(molecules["LDP"], 1)
        # Salt stays at the host's concentration; Cl- outnumbers Na+ by D1R (+14) + dopamine (+1).
        self.assertEqual(molecules["CLA"] - molecules["SOD"], 15)
        salt = report["charge"]
        self.assertAlmostEqual(salt["final_concentration_molar"], salt["host"]["concentration_molar"], delta=0.001)
        groups = self._index(folder)
        self.assertEqual(len(groups["Binder"]), 23)
        self.assertEqual(len(groups["Host_tip"]), 1)
        u = self._universe(folder)
        self.assertEqual(set(u.atoms[np.array(groups["Binder"]) - 1].resnames), {"LDP"})
        ligand = u.select_atoms("resname LDP and not name H*").positions
        tip = u.atoms[groups["Host_tip"][0] - 1].position
        self.assertGreater(ligand[:, 2].min() - tip[2], 20.0 - 0.02)

    def test_pose_folder_is_ready_for_the_charmm_gui_run_script(self):
        folder = self.out / "binderpose2"
        for name in ["README", "step6.0_minimization.mdp", "step7_production.mdp"]:
            self.assertEqual((folder / name).read_bytes(), (EXAMPLE / "host" / name).read_bytes())
        groups = self._index(folder)
        self.assertEqual(len(groups["SYSTEM"]), sum(len(groups[g]) for g in ("SOLU", "MEMB", "SOLV")))
        self.assertEqual(sorted(groups["SOLU_MEMB"]), sorted(groups["SOLU"] + groups["MEMB"]))
        self.assertLessEqual(set(groups["Binder"]), set(groups["SOLU"]))

    @staticmethod
    def _universe(folder):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mda.Universe(str(folder / "step5_input.gro"))

    @staticmethod
    def _index(folder):
        groups = {}
        for block in (folder / "index.ndx").read_text().split("[ ")[1:]:
            name, *rows = block.split("\n")
            groups[name.strip(" ]")] = [int(x) for row in rows for x in row.split()]
        return groups


class MultiMoleculeBinderTests(unittest.TestCase):
    """Binders made of several molecules, taken from the D1R-Gs example system."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            system = mda.Universe(str(GPROTEIN / "step5_input.gro"))
            # PROB (Galpha), PROC (Gbeta), PROD (Ggamma) follow the 5206-atom receptor.
            trimer = system.atoms[5206:5206 + 3944 + 5113 + 897]
            positions = _unwrap_sequential(trimer.positions.astype(float), system.dimensions[:3])
            trimer.positions = positions
            trimer.write(str(tmp / "trimer.gro"))
            gamma = trimer[3944 + 5113:]
            mixed = mda.Merge(gamma, mda.Universe(str(EXAMPLE / "dopamine" / "dopamine.pdb")).atoms)
            mixed.dimensions = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
            mixed.atoms.write(str(tmp / "gamma_ldp.gro"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _run(self, name, coordinates, itps, **extra):
        out = Path(self.tmp.name) / name
        raw = {"host": str(EXAMPLE / "host"), "binder_coordinates": str(Path(self.tmp.name) / coordinates),
               "binder_itp": [str(p) for p in itps], "random_poses": 0, **extra}
        path = Path(self.tmp.name) / f"{name}.yaml"
        path.write_text(yaml.safe_dump({"addbinder": raw}))
        with contextlib.redirect_stdout(io.StringIO()):
            code = run_addbinder(str(path), output_root=str(out))
        return code, out / "addbinder_output"

    def test_protein_and_ligand_binder_is_built_in_topology_order(self):
        code, out = self._run("mixed", "gamma_ldp.gro",
                              [GPROTEIN / "toppar" / "PROD.itp", EXAMPLE / "dopamine" / "LDP.itp"],
                              binder_forcefield=str(EXAMPLE / "dopamine" / "forcefield.itp"),
                              distance=3.0, min_image_gap=3.0, min_membrane_gap=3.0)
        self.assertEqual(code, 0)
        folder = out / "binderpose1"
        report = json.loads((folder / "addbinder_report.json").read_text())
        self.assertEqual(report["binder"]["moleculetype"], "PROD+LDP")
        self.assertEqual(report["binder"]["net_charge"], 1.0)
        self.assertEqual(report["final_net_charge"], 0.0)
        names = [name for name, _ in parse_top_molecules(str(folder / "topol.top"))]
        self.assertEqual(names[:3], ["PROA", "PROD", "LDP"])
        index = (folder / "index.ndx").read_text()
        binder = index.split("[ Binder ]")[1].split("[")[0].split()
        self.assertEqual(len(binder), 897 + 23)

    def test_trimer_too_tall_for_the_box_is_rejected_with_height_needed(self):
        code, out = self._run("trimer", "trimer.gro",
                              [GPROTEIN / "toppar" / f"{n}.itp" for n in ("PROB", "PROC", "PROD")],
                              side="lower", distance=120.0, min_image_gap=4.0, min_membrane_gap=4.0,
                              binder_forcefield=str(GPROTEIN / "toppar" / "forcefield.itp"))
        self.assertEqual(code, 3)
        report = json.loads((out / "binderpose1" / "addbinder_report.json").read_text())
        self.assertEqual(report["binder"]["moleculetype"], "PROB+PROC+PROD")
        self.assertEqual(report["binder"]["atoms"], 9954)
        self.assertEqual(report["binder"]["net_charge"], -7.0)
        hint = next(f for f in report["failures"] if "more box height" in f)
        needed = float(hint.split("About ")[1].split(" A")[0])
        # 120 A below D1R leaves ~34 A of the ~158 A water column; the flat trimer is ~57 A tall.
        self.assertGreater(needed, 15.0)

    def test_coordinate_mismatch_names_the_molecule(self):
        with self.assertRaisesRegex(ConfigError, "PROD atom 1"):
            self._run("swapped", "trimer.gro",
                      [GPROTEIN / "toppar" / f"{n}.itp" for n in ("PROB", "PROD", "PROC")])



class ParameterCoverageTests(unittest.TestCase):
    """Every binder molecule is checked: atom types, bonded terms and CMAP."""

    def setUp(self):
        from charmprot import _included_files
        toppar = (EXAMPLE / "host" / "toppar").resolve()
        self.host_ff = [str(p) for p in _included_files(toppar / "forcefield.itp", toppar)]
        self.gs = [str(EXAMPLE / "gprotein" / f"PRO{c}.itp") for c in "BCD"]

    def test_gs_needs_the_complex_force_field(self):
        from topology import TopologyError, check_parameters
        with self.assertRaisesRegex(TopologyError, r"bonds with no parameters: CC-CT1 \(residue LEU 246\)"):
            check_parameters(self.gs, self.host_ff, "advice")
        check_parameters(self.gs, self.host_ff + [str(EXAMPLE / "gprotein" / "forcefield.itp")], "advice")

    def test_non_standard_residue_is_named(self):
        import re
        from topology import TopologyError, check_parameters
        lines = Path(self.gs[2]).read_text().splitlines()
        row = next(i for i, line in enumerate(lines) if re.match(r"\s*20\s+\S+\s+\d+\s+\S+\s+\S+", line))
        fields = lines[row].split()
        lines[row] = lines[row].replace(fields[1], "ZPX1", 1)  # an atom type no force field here defines
        with tempfile.TemporaryDirectory() as tmp:
            modified = Path(tmp) / "PROD.itp"
            modified.write_text("\n".join(lines) + "\n")
            with self.assertRaisesRegex(TopologyError, rf"atom types not in \[ atomtypes \]: ZPX1 "
                                                        rf"\(residue {fields[3]} {fields[2]}\).*binder_forcefield"):
                from addbinder import PARAMETER_ADVICE
                check_parameters([str(modified)], self.host_ff + [str(EXAMPLE / "gprotein" / "forcefield.itp")],
                                 PARAMETER_ADVICE)


class GProteinExampleTests(unittest.TestCase):
    """The intracellular example: the Gs trimer pulled straight out from its bound arrangement."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        with contextlib.redirect_stdout(io.StringIO()):
            cls.code = run_addbinder(str(EXAMPLE / "addbinder_gprotein.yaml"), output_root=cls.tmp.name)
        cls.out = Path(cls.tmp.name) / "addbinder_output"

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_three_poses_built_in_topology_order(self):
        self.assertEqual(self.code, 0)
        for name, distance in (("binderpose1", 20.0), ("binderpose2", 15.0), ("binderpose3", 15.0)):
            report = json.loads((self.out / name / "addbinder_report.json").read_text())
            self.assertEqual(report["status"], "accepted")
            self.assertAlmostEqual(report["checks"]["gap_to_host_angstrom"], distance, delta=0.05)
            self.assertEqual(report["final_net_charge"], 0.0)
            names = [n for n, _ in parse_top_molecules(str(self.out / name / "topol.top"))]
            self.assertEqual(names[:4], ["PROA", "PROB", "PROC", "PROD"])

    def test_salt_keeps_the_host_concentration(self):
        for name in ("binderpose1", "binderpose2", "binderpose3"):
            salt = json.loads((self.out / name / "addbinder_report.json").read_text())["charge"]
            self.assertAlmostEqual(salt["host"]["concentration_molar"], 0.1537, places=3)
            self.assertAlmostEqual(salt["final_concentration_molar"], salt["host"]["concentration_molar"],
                                   delta=0.001)
            counts = salt["final_counts"]
            # D1R +14, Gs -7: seven more Cl- than Na+, and the pairs fit the remaining water.
            self.assertEqual(counts["CLA"] - counts["SOD"], 7)
            self.assertEqual(counts["SOD"], round(salt["target_concentration_molar"] * salt["final_waters"] / 55.51))
            self.assertEqual(salt["after"], 0.0)
        # binderpose1 needs Na+ added: they replace bulk waters away from protein and other ions.
        added = json.loads((self.out / "binderpose1" / "addbinder_report.json").read_text())["charge"]["added_ions"]
        self.assertEqual([a["moleculetype"] for a in added], ["SOD"] * 3)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            built = mda.Universe(str(self.out / "binderpose1" / "step5_input.gro"))
            host = mda.Universe(str(EXAMPLE / "host" / "step5_input.gro"))
        from MDAnalysis.lib.distances import distance_array
        ions = built.select_atoms("resname SOD CLA")
        new = distance_array(ions.positions, host.select_atoms("resname SOD CLA").positions,
                             box=built.dimensions).min(axis=1) > 0.01
        self.assertEqual(int(new.sum()), 3)
        protein = built.atoms[:5206 + 9954]
        protein = protein[[not n.startswith("H") for n in protein.names]]
        self.assertGreaterEqual(distance_array(ions[new].positions, protein.positions,
                                               box=built.dimensions).min(), 10.0 - 0.01)
        others = distance_array(ions[new].positions, ions.positions, box=built.dimensions)
        others[others < 0.01] = 99
        self.assertGreaterEqual(others.min(), 5.0 - 0.01)

    def test_random_poses_move_out_sideways_and_the_last_is_rotated(self):
        rotated = json.loads((self.out / "binderpose3" / "addbinder_report.json").read_text())
        self.assertEqual(rotated["pose"]["orientation"], "euler")
        self.assertFalse(np.allclose(rotated["rotation_matrix"], np.eye(3), atol=0.05))
        self.assertTrue(rotated["pose"]["from_input_position"])
        directions = []
        for name in ("binderpose2",):
            report = json.loads((self.out / name / "addbinder_report.json").read_text())
            pose = report["pose"]
            self.assertEqual((pose["orientation"], pose["from_input_position"]), ("as_is", True))
            self.assertTrue(np.allclose(report["rotation_matrix"], np.eye(3)))
            self.assertTrue(0 < pose["approach"][0] <= 45.0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                source = mda.Universe(str(EXAMPLE / "gprotein" / "gs_trimer.pdb")).atoms.positions
                built = mda.Universe(str(self.out / name / "step5_input.gro"))
            lengths = built.dimensions[:3]
            shift = built.atoms[5206:5206 + 9954].positions - source
            shift -= lengths * np.round((shift - shift[0]) / lengths)
            self.assertLess(float(shift.std(axis=0).max()), 0.01)  # a pure translation
            move = shift.mean(axis=0)
            directions.append(move / np.linalg.norm(move))
            tilt = np.degrees(np.arccos(-directions[-1][2]))
            self.assertAlmostEqual(tilt, pose["approach"][0], delta=0.1)
        # Directions differ by at least random_min_angle.
        from addbinder import PoseSpec, approach_direction
        a, b = (approach_direction(PoseSpec(approach=json.loads(
            (self.out / n / "addbinder_report.json").read_text())["pose"]["approach"]), "lower")
            for n in ("binderpose2", "binderpose3"))
        self.assertGreater(np.degrees(np.arccos(a @ b)), 25.0)

    def test_trimer_only_moves_along_the_membrane_normal(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            source = mda.Universe(str(EXAMPLE / "gprotein" / "gs_trimer.pdb")).atoms.positions
            built = mda.Universe(str(self.out / "binderpose1" / "step5_input.gro"))
        lengths = built.dimensions[:3]
        shift = built.atoms[5206:5206 + 9954].positions - source
        shift -= lengths * np.round((shift - shift[0]) / lengths)
        self.assertLess(float(shift.std(axis=0).max()), 0.01)
        self.assertLess(float(np.abs(shift.mean(axis=0)[:2]).max()), 0.01)
        self.assertLess(float(shift.mean(axis=0)[2]), -10.0)  # away from the membrane, to -z


if __name__ == "__main__":
    unittest.main()
