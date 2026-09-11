import tempfile
import unittest
from pathlib import Path

from config import ConfigError, load_config, load_minimization_config


MINIMAL_LIGAND_CONFIG = """
paths:
  environment_coordinates: &environment_coordinates environment.gro
  incoming_coordinates: &incoming_coordinates ligand.mol2

ligand_replace:
  structure_path: *environment_coordinates
  format: auto
  original_ligand:
    resname: OLD
  new_ligand:
    coord_path: *incoming_coordinates
    format: auto
    resname: NEW
  fit:
    method: autofit
"""


class ConfigurationSafetyTests(unittest.TestCase):
    def _write(self, text):
        directory = tempfile.TemporaryDirectory()
        path = Path(directory.name) / "config.yaml"
        path.write_text(text, encoding="utf-8")
        self.addCleanup(directory.cleanup)
        return path

    def test_path_registry_aliases_feed_mode_sections(self):
        config = load_config(
            str(self._write(MINIMAL_LIGAND_CONFIG)),
            mode="lig",
            check_paths=False,
        )
        self.assertEqual(config.ligand_replace.structure_path, "environment.gro")
        self.assertEqual(config.ligand_replace.new_ligand.coord_path, "ligand.mol2")

    def test_duplicate_yaml_key_is_rejected(self):
        text = MINIMAL_LIGAND_CONFIG.replace(
            "  structure_path: *environment_coordinates\n",
            "  structure_path: first.gro\n  structure_path: second.gro\n",
        )
        with self.assertRaisesRegex(ConfigError, "Duplicate YAML key.*structure_path"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_missing_config_is_reported_as_config_error(self):
        with self.assertRaisesRegex(ConfigError, "Cannot read config file"):
            load_config("definitely-missing.yaml", mode="lig", check_paths=False)

    def test_invalid_numeric_type_is_reported_without_conversion_traceback(self):
        text = MINIMAL_LIGAND_CONFIG + """
clash_detection:
  threshold: not-a-number
"""
        with self.assertRaisesRegex(ConfigError, "clash_detection.threshold"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_inspection_output_cannot_overwrite_an_input(self):
        text = MINIMAL_LIGAND_CONFIG + """
output:
  pdb_path: final.pdb
  gro_path: final.gro
  report_path: report
  inspection_pdb_path: environment.gro
"""
        with self.assertRaisesRegex(ConfigError, "must not overwrite an input"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_standalone_minimization_config(self):
        text = """
paths:
  coordinates: &coordinates assembled.gro
  topology: &topology topol.top
minimization:
  enabled: true
  coordinates_path: *coordinates
  topology_path: *topology
  output_gro_path: minimized.gro
  report_path: minimization_report
  algorithm: lbfgs
  nonbonded_method: PME
"""
        spec = load_minimization_config(
            str(self._write(text)), check_paths=False
        )
        self.assertEqual(spec.coordinates_path, "assembled.gro")
        self.assertEqual(spec.topology_path, "topol.top")
        self.assertEqual(spec.algorithm, "lbfgs")

    def test_minimization_refuses_to_overwrite_coordinates(self):
        text = """
minimization:
  enabled: true
  coordinates_path: assembled.gro
  topology_path: topol.top
  output_gro_path: assembled.gro
  report_path: minimization_report
"""
        with self.assertRaisesRegex(ConfigError, "paths collide"):
            load_minimization_config(str(self._write(text)), check_paths=False)

    def test_integrated_minimization_requires_generated_topology(self):
        text = MINIMAL_LIGAND_CONFIG + """
minimization:
  enabled: true
"""
        with self.assertRaisesRegex(ConfigError, "requires topology.enabled"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)


if __name__ == "__main__":
    unittest.main()
