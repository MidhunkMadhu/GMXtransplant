import tempfile
import unittest
import io
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

import yaml

from config import ConfigError, load_config, load_minimization_config


# Input paths are always given in full: nothing resolves against a shared root.
MINIMAL_LIGAND_CONFIG = """
paths:
  environment_coordinates: &environment_coordinates /inputs/environment.gro
  incoming_coordinates: &incoming_coordinates /inputs/ligand.mol2

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
        self.assertEqual(config.ligand_replace.structure_path, "/inputs/environment.gro")
        self.assertEqual(config.ligand_replace.new_ligand.coord_path, "/inputs/ligand.mol2")

    def test_relative_input_paths_are_rejected(self):
        text = MINIMAL_LIGAND_CONFIG.replace("/inputs/environment.gro", "environment.gro")
        with self.assertRaisesRegex(ConfigError, "must be given in full"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_output_filenames_land_in_the_output_folder(self):
        text = MINIMAL_LIGAND_CONFIG + """
output:
  pdb_path: final.pdb
  gro_path: final.gro
  report_path: report
  inspection_pdb_path: inspection.pdb
"""
        with tempfile.TemporaryDirectory() as out:
            config = load_config(str(self._write(text)), mode="lig",
                                 check_paths=False, output_root=out)
            self.assertEqual(config.output.gro_path, str(Path(out) / "final.gro"))
            self.assertEqual(config.topology.output_dir, str(Path(out).absolute()))

    def test_duplicate_yaml_key_is_rejected(self):
        text = MINIMAL_LIGAND_CONFIG.replace(
            "  structure_path: *environment_coordinates\n",
            "  structure_path: first.gro\n  structure_path: second.gro\n",
        )
        with self.assertRaisesRegex(ConfigError, "Duplicate YAML key.*structure_path"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_plain_path_references_and_null(self):
        text = MINIMAL_LIGAND_CONFIG.replace(
            "&environment_coordinates ", ""
        ).replace("&incoming_coordinates ", "").replace(
            "*environment_coordinates", '"${environment_coordinates}"'
        ).replace("*incoming_coordinates", '"${incoming_coordinates}"')
        text = text.replace("paths:\n", "paths:\n  reference: null\n")
        text += '\nrefgro: "${reference}"\n'
        config = load_config(str(self._write(text)), mode="lig", check_paths=False)
        self.assertEqual(config.ligand_replace.structure_path, "/inputs/environment.gro")
        self.assertEqual(config.ligand_replace.new_ligand.coord_path, "/inputs/ligand.mol2")
        self.assertIsNone(config.refgro)

    def test_invalid_path_references_fail_clearly(self):
        for reference, error in [
            ("${missing}", "Unknown path reference"),
            ("${environment_coordinates}/extra", "Invalid path reference"),
        ]:
            with self.subTest(reference=reference):
                text = MINIMAL_LIGAND_CONFIG.replace(
                    "*environment_coordinates", f'"{reference}"'
                )
                with self.assertRaisesRegex(ConfigError, error):
                    load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_examples_with_and_without_comments_are_equivalent(self):
        from run_pipeline import cli

        for mode in ("protein", "lig", "chl", "minimize"):
            with self.subTest(mode=mode):
                versions = []
                for option in ("--comments", "--no-comments"):
                    output = io.StringIO()
                    with redirect_stdout(output):
                        self.assertEqual(cli(["--show-example", mode, option]), 0)
                    text = output.getvalue()
                    self.assertEqual("#" in text, option == "--comments")
                    versions.append(yaml.safe_load(text))
                    path = str(self._write(text))
                    if mode == "minimize":
                        config = load_minimization_config(path, check_paths=False)
                        self.assertEqual(config.coordinates_path,
                                         "/path/to/assembled/step5_input.gro")
                    else:
                        config = load_config(path, mode=mode, check_paths=False)
                        if mode == "protein":
                            self.assertEqual(len(config.replacement_ligands), 2)
                            self.assertEqual(config.replacement_structure.protein_mask, ":1-963")
                            self.assertEqual(config.replacement_ligands[1].itp_path,
                                             "/path/to/replacement/toppar/LIG2.itp")
                self.assertEqual(*versions)

    def test_comment_options_require_example_generation(self):
        from run_pipeline import cli

        for option in ("--comments", "--no-comments"):
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                cli(["--mode", "lig", "-i", "config.yaml", option])
            self.assertEqual(error.exception.code, 2)

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
  inspection_pdb_path: /inputs/environment.gro
"""
        with self.assertRaisesRegex(ConfigError, "must not overwrite an input"):
            load_config(str(self._write(text)), mode="lig", check_paths=False)

    def test_standalone_minimization_config(self):
        text = """
paths:
  coordinates: &coordinates /inputs/assembled.gro
  topology: &topology /inputs/topol.top
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
        self.assertEqual(spec.coordinates_path, "/inputs/assembled.gro")
        self.assertEqual(spec.topology_path, "/inputs/topol.top")
        self.assertEqual(spec.algorithm, "lbfgs")

    def test_minimization_refuses_to_overwrite_coordinates(self):
        text = """
minimization:
  enabled: true
  coordinates_path: /inputs/assembled.gro
  topology_path: /inputs/topol.top
  output_gro_path: /inputs/assembled.gro
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
