"""GUI model and offscreen UI regression tests. Qt tests skip if PySide6 cannot be imported."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from gmxtransplant.gui.model import (MODES, EXAMPLES, defaults, build_job, validate_job,
    destination, example_path, template_path, load_yaml, read_yaml, dump_yaml,
    example_document, template_document)
from gmxtransplant.gui.worker import cleanup, execute


class GuiModelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.source = self.root / 'source'
        self.source.mkdir()
        # A user document names every input in full; nothing is anchored for it.
        self.raw = template_document('protein')
        self.raw['paths'].update(
            target_coordinates=str(self.source / 'inputs/environment/equilibrated.gro'),
            replacement_coordinates=str(self.source / 'inputs/replacement/step5_input.gro'),
            protein_toppar=str(self.source / 'inputs/replacement/toppar'),
            environment_toppar=str(self.source / 'inputs/environment/toppar'))
        for index, name in enumerate(('LIG1.itp', 'LIG2.itp')):
            self.raw['replacement_ligands'][index]['itp_path'] = str(
                self.source / 'inputs/replacement/toppar' / name)

    def tearDown(self):
        self.tmp.cleanup()

    def job(self, raw=None, mode='protein'):
        return build_job(raw or self.raw, self.root / 'output', mode)

    def test_fixed_mode_and_example_layout(self):
        for mode, folder in [('protein', 'protein'), ('lig', 'ligand'), ('chl', 'cholesterol')]:
            self.assertEqual(destination(self.root, mode), self.root / folder)
            self.assertEqual(destination(self.root, mode, True), self.root / 'examples' / folder)
        with self.assertRaises(ValueError):
            destination('', 'protein')

    def test_addbinder_example_writes_pose_folders_in_its_mode_folder(self):
        p = example_path('addbinder')
        job = build_job(load_yaml(p), self.root, 'addbinder', example_base=p.parent)
        self.assertEqual(job['directory'], str(self.root / 'examples' / 'addbinder'))
        self.assertEqual(job['config']['addbinder']['output_dir'], job['directory'])
        self.assertEqual(job['artifacts'][-4:], ['binderpose1', 'binderpose2', 'binderpose3', 'binderpose4'])
        self.assertLessEqual({'summary.txt', 'view.pml', 'view.vmd'}, set(job['artifacts']))
        raw = template_document('addbinder')
        raw['addbinder']['poses'] = [{'name': '../escape'}]
        with self.assertRaises(ValueError):
            build_job(raw, self.root, 'addbinder')

    def test_addbinder_pose_cannot_use_a_run_artifact_name(self):
        raw = template_document('addbinder')
        for name in ('run.yaml', 'summary.txt', 'visualization'):
            with self.subTest(name=name):
                raw['addbinder']['poses'] = [{'name': name}]
                with self.assertRaisesRegex(ValueError, 'reserved for run output'):
                    build_job(raw, self.root, 'addbinder')

    def test_all_example_inputs_validate_without_mutating_inputs(self):
        for mode in EXAMPLES:
            p = example_path(mode)
            raw = load_yaml(p)
            before = dump_yaml(raw)
            job = build_job(raw, self.root, mode, example_base=p.parent)
            validate_job(job, check_paths=True)
            self.assertEqual(before, dump_yaml(raw))
            self.assertFalse(Path(job['directory']).exists())

    def test_outputs_are_contained_and_inputs_stay_explicit(self):
        self.raw['output']['gro_path'] = '../../outside.gro'
        job = self.job()
        self.assertEqual(job['config']['output']['gro_path'], str(self.root / 'output/protein/outside.gro'))
        self.assertEqual(job['config']['target_box']['path'], str(self.source / 'inputs/environment/equilibrated.gro'))
        self.assertNotIn('paths', job['config'])

    def test_relative_inputs_are_rejected_without_a_base_folder(self):
        self.raw['paths']['target_coordinates'] = 'environment/equilibrated.gro'
        with self.assertRaisesRegex(ValueError, 'complete path'):
            self.job()

    def test_custom_references_and_itp_lists(self):
        self.raw['paths']['custom'] = str(self.source / 'my ligand.itp')
        self.raw['topology']['ligand_itp_paths'] = ['${custom}', str(self.source / 'second.itp')]
        job = self.job()
        self.assertEqual(job['config']['topology']['ligand_itp_paths'],
                         [str(self.source / 'my ligand.itp'), str(self.source / 'second.itp')])

    def test_all_diagnostic_outputs_redirected(self):
        raw = example_document('chl')
        raw['cholesterol']['merged_pdb_path'] = '/tmp/elsewhere/merged.pdb'
        job = self.job(raw, 'chl')
        self.assertEqual(job['config']['cholesterol']['merged_pdb_path'], str(self.root / 'output/cholesterol/merged.pdb'))

    def test_collision_and_source_overlap_rejected(self):
        self.raw['output']['pdb_path'] = 'run.yaml'
        with self.assertRaises(ValueError):
            self.job()
        self.raw['output']['pdb_path'] = 'step5_input.pdb'
        self.raw['paths']['target_coordinates'] = str(self.root / 'output/protein/input.gro')
        with self.assertRaises(ValueError):
            self.job()

    def test_symlink_output_rejected(self):
        (self.root / 'protein').symlink_to(self.source, target_is_directory=True)
        with self.assertRaises(ValueError):
            destination(self.root, 'protein')

    def test_duplicate_yaml_rejected(self):
        with self.assertRaises(ValueError):
            read_yaml('output: {}\noutput: {}')

    def test_cleanup_replaces_toppar_and_old_outputs_only(self):
        job = self.job()
        out = Path(job['directory'])
        out.mkdir(parents=True)
        (out / 'toppar').mkdir()
        (out / 'toppar/stale.itp').write_text('old')
        (out / 'notes.txt').write_text('keep')
        (out / 'old_report.txt').write_text('old')
        (out / '.gui-manifest.json').write_text(json.dumps({'artifacts': ['toppar', 'old_report.txt']}))
        cleanup(job)
        self.assertFalse((out / 'toppar').exists())
        self.assertFalse((out / 'old_report.txt').exists())
        self.assertEqual((out / 'notes.txt').read_text(), 'keep')

    def test_cleanup_preserves_existing_file_without_a_manifest(self):
        job = self.job()
        out = Path(job['directory'])
        out.mkdir(parents=True)
        existing = out / 'step5_input.gro'
        existing.write_text('user data')
        with self.assertRaisesRegex(ValueError, 'was not created by a previous GUI run'):
            cleanup(job)
        self.assertEqual(existing.read_text(), 'user data')
        self.assertFalse((out / '.gui-manifest.json').exists())

    def test_cleanup_preserves_unowned_file_after_output_name_changes(self):
        job = self.job()
        out = Path(job['directory'])
        out.mkdir(parents=True)
        old = out / 'step5_input.gro'
        old.write_text('previous GUI output')
        unowned = out / 'new_result.gro'
        unowned.write_text('user data')
        (out / '.gui-manifest.json').write_text(json.dumps({'artifacts': ['step5_input.gro']}))
        job['artifacts'].append('new_result.gro')
        with self.assertRaisesRegex(ValueError, 'was not created by a previous GUI run'):
            cleanup(job)
        self.assertEqual(old.read_text(), 'previous GUI output')
        self.assertEqual(unowned.read_text(), 'user data')

    def test_manifest_owns_only_outputs_created_by_the_run(self):
        job = self.job()

        def run(_args):
            Path('step5_input.gro').write_text('generated')
            return 0

        with patch('gmxtransplant.gui.worker.validate_job'), patch('run_pipeline.cli', side_effect=run):
            self.assertEqual(execute(job), 0)
        out = Path(job['directory'])
        owned = json.loads((out / '.gui-manifest.json').read_text())['artifacts']
        self.assertIn('step5_input.gro', owned)
        self.assertNotIn('view.pse', owned)
        (out / 'view.pse').write_text('user file')
        with self.assertRaisesRegex(ValueError, 'was not created by a previous GUI run'):
            cleanup(job)
        self.assertEqual((out / 'step5_input.gro').read_text(), 'generated')
        self.assertEqual((out / 'view.pse').read_text(), 'user file')

    def test_cleanup_rejects_symlink_before_deleting_anything(self):
        job = self.job()
        out = Path(job['directory'])
        out.mkdir(parents=True)
        (out / 'step5_input.gro').write_text('previous')
        (out / 'toppar').symlink_to(self.source, target_is_directory=True)
        with self.assertRaises(ValueError):
            cleanup(job)
        self.assertEqual((out / 'step5_input.gro').read_text(), 'previous')

    def test_failed_validation_does_not_clear_previous_outputs(self):
        job = self.job()
        out = Path(job['directory'])
        out.mkdir(parents=True)
        (out / 'step5_input.gro').write_text('previous')
        with self.assertRaises(ValueError):
            execute(job)
        self.assertEqual((out / 'step5_input.gro').read_text(), 'previous')

    def test_null_output_sections_use_pipeline_defaults(self):
        for section in ['output', 'ndx', 'topology', 'minimization']:
            self.raw[section] = None
        job = self.job()
        validate_job(job)
        self.assertEqual(job['config']['output']['gro_path'], str(self.root / 'output/protein/step5_input.gro'))

    def test_running_destination_is_locked_before_cleanup(self):
        import fcntl
        job = self.job()
        root = Path(job['directory'])
        root.mkdir(parents=True)
        with (root / '.gui.lock').open('w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with patch('gmxtransplant.gui.worker.validate_job'), patch('gmxtransplant.gui.worker.cleanup') as remove:
                with self.assertRaisesRegex(ValueError, 'Another run'):
                    execute(job)
                remove.assert_not_called()

    def test_standalone_minimization(self):
        raw = load_yaml(template_path('minimize'))
        raw['minimization']['coordinates_path'] = str(self.source / 'assembled.gro')
        raw['minimization']['topology_path'] = str(self.source / 'topol.top')
        job = self.job(raw, 'minimize')
        validate_job(job)
        self.assertEqual(job['config']['minimization']['output_dir'], str(self.root / 'output/minimization/openmm_minimization'))

    def test_documentation_pdf_ships_with_the_package(self):
        from gmxtransplant.gui.model import documentation_path
        path = documentation_path()
        self.assertTrue(path.is_file())
        self.assertEqual(path.read_bytes()[:5], b'%PDF-')
        self.assertEqual(path.name, 'GMXtransplant.pdf')

    def test_documentation_falls_back_from_viewer_to_browser_then_reports(self):
        from gmxtransplant.gui.model import open_documentation
        opened = []

        def browser(path):
            opened.append(path)
            return True

        self.assertTrue(open_documentation(viewer=lambda path: False, browser=browser))
        self.assertEqual(len(opened), 1)
        with self.assertRaisesRegex(RuntimeError, 'Could not open the documentation'):
            open_documentation(viewer=lambda path: False, browser=lambda path: False)

        def broken(path):
            raise OSError('no display')

        with self.assertRaisesRegex(RuntimeError, 'no display'):
            open_documentation(viewer=broken, browser=lambda path: False)

    def test_schema_matches_accepted_sections(self):
        import config
        for mode in EXAMPLES:
            schema = defaults(mode)
            for section, value in schema.items():
                if section in config._SECTION_KEYS:
                    self.assertLessEqual(set(value), config._SECTION_KEYS[section])


class ProgressTests(unittest.TestCase):
    def test_chunks_are_buffered_and_technical_output_is_excluded(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress()
        progress.feed('[4-6] Aligning replacement')
        self.assertEqual(progress.stage, 'Ready to start')
        progress.feed(' protein ...\n[config] path=/private/input.pdb\nRotation: [[1, 0, 0]]\n')
        self.assertEqual(progress.stage, 'Aligning the incoming protein to the target')
        self.assertFalse(progress.facts)
        progress.feed('RMSD (positionally-matched): 0.289 A (3842 atom pairs)\n')
        self.assertIn('0.289 Å', progress.facts['Alignment'])

    def test_addbinder_progress_reports_tip_and_poses(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress('addbinder')
        progress.feed('[frame] headgroup planes 70.4 / 32.8 A; tip HSD159:NE2 11.5 A beyond the upper plane\n')
        self.assertIn('HSD159:NE2, 11.5 Å', progress.facts['Protein tip'])
        progress.feed('[pose01_flat] accepted: host 10.2 A\n[pose04] rejected: host 22.1 A\n')
        self.assertEqual(progress.facts['Poses'], 'pose01_flat accepted; pose04 rejected')
        progress.feed('[complete] 1/2 pose(s) built; see summary.txt\n')
        self.assertEqual(progress.facts['Result'], '1 of 2 poses built; see summary.txt.')

    def test_cholesterol_progress_points_to_the_repacking_step(self):
        from gmxtransplant.gui.progress import RunProgress
        from gmxtransplant.gui.window import MODE_DESCRIPTIONS
        progress = RunProgress('chl')
        progress.feed('      Wrote /x/out/repack/: repacking equilibration (1.375 ns, 5 restored cholesterol(s) '
                      'held) with topol.top, toppar/, the GRO, index.ndx and mdp files; a sample run script is '
                      'also saved.\n')
        self.assertIn('repack/ holds a 1.375 ns repacking equilibration', progress.facts['Repacking'])
        self.assertIn('repack/', MODE_DESCRIPTIONS['chl'])

    def test_addbinder_glance_reads_the_current_pose_and_salt_lines(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress('addbinder')
        progress.feed('[binderpose1] as_is, moved out from its input position at 10 A: accepted: host 10.0 A, '
                      'host image 54.0 A\n      salt 0.154 M (host 0.154 M): 18 ion(s) removed, 0 added; net charge +0 e\n'
                      '[binderpose2] euler [1.0, 2.0, 3.0], approach tilt 40 azimuth 9 at 22.43 A: rejected\n'
                      '      salt 0.154 M (host 0.154 M): 16 ion(s) removed, 3 added; net charge +0 e\n')
        self.assertEqual(progress.facts['Poses'], 'binderpose1 accepted; binderpose2 rejected')
        self.assertEqual(progress.facts['Salt'], '0.154 M (host 0.154 M), ions per pose: binderpose1 18 removed / '
                                                 '0 added; binderpose2 16 removed / 3 added.')

    def test_cholesterol_glance_shows_restored_and_removed_molecules(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress('chl')
        progress.feed("Removed by class: {'lipid': 7} (lipids 7/500)\n")
        progress.feed('      Restored cholesterols: 5, final residues 1356-1360 (cholesterol 76-80 of 80 in topol.top).\n'
                      '      Removed next to them: 7 molecule(s) (POPC 6, POPE 1), 0.09-0.62 A from a restored cholesterol.\n'
                      '      Removed far away: 24 molecule(s) (POPE 10, CHL1 8, POPC 6), 20.2-40.8 A from the protein '
                      '(only beyond 20 A).\n')
        self.assertEqual(progress.facts['Restored cholesterol'],
                         '5 restored: final residues 1356-1360 (cholesterol 76-80 of 80).')
        self.assertTrue(progress.facts['Removed nearby'].startswith('7 molecule(s) (POPC 6, POPE 1)'))
        self.assertIn('20.2-40.8 A from the protein', progress.facts['Removed far away'])
        self.assertNotIn('Clash removal', progress.facts)

    def test_gprotein_example_card_has_its_own_folder(self):
        from gmxtransplant.gui.model import EXAMPLE_KEYS, example_folder, example_mode
        self.assertIn('addbinder_gprotein', EXAMPLE_KEYS)
        path = example_path('addbinder_gprotein')
        self.assertEqual((path.name, example_mode('addbinder_gprotein')), ('addbinder_gprotein.yaml', 'addbinder'))
        with tempfile.TemporaryDirectory() as tmp:
            job = build_job(load_yaml(path), tmp, 'addbinder', example_base=path.parent,
                            example_folder=example_folder('addbinder_gprotein'))
            self.assertEqual(Path(job['directory']).name, 'addbinder_gprotein')
            self.assertEqual(job['artifacts'][-3:], ['binderpose1', 'binderpose2', 'binderpose3'])
            cfg = validate_job(job, check_paths=True)
        self.assertTrue(cfg.poses[0].from_input_position and cfg.random_from_input_position)
        self.assertEqual(example_document('addbinder_gprotein')['addbinder']['side'], 'lower')

    def test_results_list_only_reports(self):
        from gmxtransplant.gui.window import report_files
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ('summary.txt', 'step5_input.gro', 'report.json', 'binderpose1/addbinder_report.txt',
                         'toppar/x.itp', '.hidden/a.txt'):
                (root / name).parent.mkdir(parents=True, exist_ok=True)
                (root / name).write_text('x')
            self.assertEqual(report_files(root), ['binderpose1/addbinder_report.txt', 'summary.txt'])

    def test_scientific_warnings_and_failures_remain_visible(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress()
        progress.feed('WARNING: PDB: 3 contacts accepted by configured cutoffs.\n')
        progress.feed('[OUTPUT VALIDATION ERROR] Coordinate check failed.')
        progress.finish(False)
        self.assertEqual(progress.error, 'Coordinate check failed.')
        self.assertEqual(len(progress.notices), 1)
        self.assertIn('could not finish', progress.stage)

    def test_cholesterol_announcements_do_not_imply_completed_stages(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress('chl')
        for i in range(1, 5):
            progress.feed(f'[{i}] Announcing operation ...\n')
        self.assertIn('Restoring cholesterol', progress.stage)
        self.assertFalse(progress.facts)
        progress.feed('Summary: 5 PASS, 0 FAIL, 5 checked\nInserted 5 cholesterol residue(s) (370 atoms).\n')
        self.assertIn('5 passed', progress.facts['Cholesterol conversion'])
        self.assertEqual(progress.facts['Cholesterol'], '5 molecules inserted.')
        progress.finish(True)
        self.assertEqual(progress.stage, 'Your outputs are ready.')

    def test_cancelled_status_does_not_claim_success(self):
        from gmxtransplant.gui.progress import RunProgress
        progress = RunProgress()
        progress.finish(True, cancelled=True)
        self.assertIn('cancelled', progress.stage)

    def test_macos_uses_native_display_without_x11_variables(self):
        from gmxtransplant.gui import display_available
        self.assertTrue(display_available('darwin', {}))
        self.assertFalse(display_available('linux', {}))
        self.assertTrue(display_available('linux', {'WAYLAND_DISPLAY': 'wayland-0'}))


@unittest.skipUnless(importlib.util.find_spec('PySide6'), 'Optional GUI dependencies are not installed')
class GuiEditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QApplication
        cls.app = QApplication.instance() or QApplication([])
        # Remembered output folder and theme go to a throwaway location, never
        # the user's own settings.
        cls.settings_dir = tempfile.TemporaryDirectory()
        cls.previous_settings_format = QSettings.defaultFormat()
        QSettings.setDefaultFormat(QSettings.Format.IniFormat)
        for fmt in (QSettings.Format.NativeFormat, QSettings.Format.IniFormat):
            QSettings.setPath(fmt, QSettings.Scope.UserScope, cls.settings_dir.name)

    @classmethod
    def tearDownClass(cls):
        from PySide6.QtCore import QSettings
        QSettings.setDefaultFormat(cls.previous_settings_format)
        cls.settings_dir.cleanup()

    def test_settings_are_isolated_from_the_user(self):
        from PySide6.QtCore import QSettings
        self.assertTrue(Path(QSettings().fileName()).resolve().is_relative_to(
            Path(self.settings_dir.name).resolve()))

    def test_form_roundtrips_every_template_and_example(self):
        from gmxtransplant.gui.editor import ConfigEditor
        for mode in MODES:
            paths = [template_path(mode)] + ([example_path(mode)] if mode in EXAMPLES else [])
            for path in paths:
                raw = load_yaml(path)
                editor = ConfigEditor(raw, mode)
                self.assertEqual(editor.value(), raw, str(path))
                editor.deleteLater()
        self.app.processEvents()

    def test_null_zero_false_and_text_remain_distinct(self):
        from gmxtransplant.gui.editor import ValueEditor
        for value in [0, 0.0, False, '', '0', 'from_itp', [], {}]:
            editor = ValueEditor(value, path=('test',))
            self.assertEqual(editor.value(), value)
            self.assertIs(type(editor.value()), type(value))
            editor.deleteLater()

    def test_concrete_fields_do_not_offer_type_selectors(self):
        from gmxtransplant.gui.editor import ValueEditor
        from PySide6.QtWidgets import QComboBox, QCheckBox, QLineEdit
        for value, path, expected in [(True, ('clash_detection', 'use_pbc'), QCheckBox),
                (1.2, ('clash_detection', 'threshold'), QLineEdit),
                ('/tmp/input.gro', ('target_box', 'path'), QLineEdit)]:
            editor = ValueEditor(value, path=path)
            self.assertIsInstance(editor.scalar, expected)
            self.assertFalse(editor.findChildren(QComboBox))
            self.assertEqual(editor.value(), value)
            editor.deleteLater()

    def test_ligand_charge_has_only_meaningful_alternatives(self):
        from gmxtransplant.gui.editor import ValueEditor
        editor = ValueEditor('from_itp', path=('replacement_ligands', 0, 'charge'))
        self.assertEqual(editor.value(), 'from_itp')
        self.assertEqual(editor.charge_source.count(), 2)
        editor.charge_source.setCurrentIndex(1)
        editor.scalar.setText('-1.5')
        self.assertEqual(editor.value(), -1.5)
        editor.charge_source.setCurrentIndex(0)
        self.assertEqual(editor.value(), 'from_itp')
        editor.deleteLater()

    def test_mapping_residue_labels_preserve_uppercase(self):
        from gmxtransplant.gui.editor import ValueEditor
        editor = ValueEditor({'CHL': ['CHL1'], 'POP': ['POPC', 'POPE']},
                             path=('name_restoration', 'pdb_to_full_resname'))
        self.assertEqual([r.key for r in editor.rows], ['CHL', 'POP'])
        # Entries are added inline, never through a pop-up dialog.
        self.assertFalse(editor.add_button.isEnabled())
        editor.new_key.setText('tip')
        self.assertTrue(editor.add_button.isEnabled())
        editor.add_button.click()
        self.assertEqual(editor.rows[-1].key, 'TIP')
        editor.rows[-1].editor.scalar.setText('TIP3')
        self.assertEqual(editor.value()['TIP'], ['TIP3'])
        editor.new_key.setText('TIP')
        editor.add_button.click()
        self.assertEqual(len(editor.rows), 3)
        self.assertIn('already', editor.add_note.text())
        editor.deleteLater()

    def test_name_lists_are_one_text_field(self):
        from gmxtransplant.gui.editor import ValueEditor
        from PySide6.QtWidgets import QPushButton
        editor = ValueEditor(['CHL1', 'CHL', 'CHOL'], path=('cholesterol', 'cholesterol_resnames'))
        self.assertEqual(editor.scalar.text(), 'CHL1, CHL, CHOL')
        self.assertFalse(editor.findChildren(QPushButton))
        editor.scalar.setText('CHL1 CLR,  CHOL')
        self.assertEqual(editor.value(), ['CHL1', 'CLR', 'CHOL'])
        editor.deleteLater()

    def test_new_lipid_targets_and_optional_numbers_are_numeric(self):
        from gmxtransplant.gui.editor import ValueEditor
        editor = ValueEditor({}, path=('cholesterol', 'composition', 'lipid_targets'))
        editor.add_item('popc')
        lipid = editor.rows[0].editor
        self.assertEqual(lipid.value(), {'upper': 0, 'lower': 0})
        lipid.rows[0].editor.scalar.setText('100')
        self.assertEqual(editor.value()['POPC']['upper'], 100)
        editor.deleteLater()

    def test_no_yaml_editor_and_command_log_collapses_on_rerun(self):
        from PySide6.QtCore import QEventLoop, QTimer
        from gmxtransplant.gui.window import MainWindow
        window = MainWindow()
        self.assertFalse(hasattr(window, 'yaml_text'))
        self.assertTrue(window.console.isHidden())
        window.command_toggle.click()
        self.assertFalse(window.console.isHidden())
        with tempfile.TemporaryDirectory() as output:
            # A new document starts with empty input fields, so check a filled one.
            window.install_document(example_document('charmprot'))
            window.output_root.setText(output)
            window.start_current('schema')
            self.assertTrue(window.console.isHidden())
            self.assertFalse(window.command_toggle.isChecked())
            loop = QEventLoop()
            timer = QTimer()
            timer.timeout.connect(lambda: loop.quit() if window.process is None else None)
            timer.start(20)
            QTimer.singleShot(10000, loop.quit)
            loop.exec()
            self.assertIsNone(window.process)
            self.assertEqual(window.run_title.text(), 'Checks passed')
            self.assertIn('Configuration check passed', window.console.toPlainText())
            self.assertIn('Configuration is valid', window.summary_facts.text())
            self.assertTrue(window.console.isHidden())
        window.deleteLater()

    def test_browsing_folder_updates_field(self):
        from gmxtransplant.gui.editor import ValueEditor
        editor = ValueEditor('', path=('topology', 'protein_toppar_dir'))
        with patch('gmxtransplant.gui.editor.choose_path',
                   side_effect=lambda parent, kind, initial, callback, save=False: callback('/tmp/with spaces/toppar')):
            editor.browse(True)
        self.assertEqual(editor.value(), '/tmp/with spaces/toppar')
        editor.deleteLater()

    def test_file_chooser_never_blocks_the_window(self):
        from gmxtransplant.gui.editor import choose_path
        from PySide6.QtCore import Qt
        from PySide6.QtWidgets import QWidget
        parent = QWidget()
        dialog = choose_path(parent, 'directory', '', lambda path: None)
        self.assertEqual(dialog.windowModality(), Qt.WindowModality.NonModal)
        self.assertFalse(dialog.isModal())
        dialog.close()
        parent.deleteLater()

    def test_addbinder_pose_can_be_added_in_the_form(self):
        from gmxtransplant.gui.editor import ConfigEditor
        raw = example_document('addbinder')
        editor = ConfigEditor(raw, 'addbinder')
        poses = next(row for row in editor.fields['addbinder'].editor.rows if row.key == 'poses')
        poses.editor.add_item()
        value = editor.value()
        self.assertEqual(len(value['addbinder']['poses']), 2)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = validate_job(build_job(value, tmp, 'addbinder'), check_paths=True)
        self.assertEqual((cfg.poses[-1].name, cfg.poses[-1].distance, cfg.poses[-1].orientation),
                         ('binderpose2', 20.0, 'end_on'))  # the second auto orientation
        editor.deleteLater()

    def test_optional_field_can_be_enabled(self):
        from gmxtransplant.gui.editor import ConfigEditor
        editor = ConfigEditor({'target_box': {'path': '/tmp/test.gro'}}, 'protein')
        section = editor.fields['target_box'].editor
        field = next(row for row in section.rows if row.key == 'box_dimensions')
        self.assertFalse(field.editor.isVisibleTo(editor))
        self.assertNotIn('box_dimensions', editor.value()['target_box'])
        field.toggle.setChecked(True)
        self.assertEqual(editor.value()['target_box']['box_dimensions'], [0.0, 0.0, 0.0, 90.0, 90.0, 90.0])
        field.toggle.setChecked(False)
        self.assertNotIn('box_dimensions', editor.value()['target_box'])
        editor.deleteLater()

    def test_switched_sections_have_one_checkbox_and_no_managed_folders(self):
        from gmxtransplant.gui.editor import ConfigEditor
        from PySide6.QtWidgets import QCheckBox, QLabel
        editor = ConfigEditor(example_document('lig'), 'lig')
        self.assertEqual(editor.fields['ndx'].kind, 'fixed')
        self.assertIsNone(editor.fields['ndx'].toggle)
        self.assertNotIn('enabled', [r.key for r in editor.fields['ndx'].editor.rows])
        for key in ('topology', 'minimization'):
            field = editor.fields[key]
            self.assertEqual(field.kind, 'toggle')
            boxes = [b for b in field.findChildren(QCheckBox) if b.text().lower() in ('enabled', 'include this section')]
            self.assertFalse(boxes, key)
            self.assertNotIn('enabled', [r.key for r in field.editor.rows])
        self.assertNotIn('output_dir', [r.key for r in editor.fields['topology'].editor.rows])
        texts = ' '.join(label.text() for label in editor.findChildren(QLabel))
        self.assertNotIn('Managed by the output folder', texts)
        self.assertNotIn('Leave unset', ' '.join(b.text() for b in editor.findChildren(QCheckBox)))
        editor.fields['topology'].toggle.setChecked(False)
        self.assertFalse(editor.value()['topology']['enabled'])
        editor.deleteLater()

    def test_minimization_shows_only_essential_settings(self):
        from gmxtransplant.gui.editor import ConfigEditor
        editor = ConfigEditor({'charmprot': {'reference': '/a', 'transplant': '/b'}}, 'charmprot')
        minimization = editor.fields['minimization']
        self.assertEqual(minimization.kind, 'toggle')
        self.assertNotIn('minimization', editor.value())
        minimization.toggle.setChecked(True)
        self.assertEqual({r.key for r in minimization.editor.rows},
                         {'restraint_force_constant_kj_mol_nm2', 'max_iterations', 'platform',
                          'restraint_residue_classes'})
        self.assertEqual(editor.value()['minimization'], {'enabled': True})
        editor.deleteLater()
        standalone = ConfigEditor(load_yaml(template_path('minimize')), 'minimize')
        keys = {r.key for r in standalone.fields['minimization'].editor.rows}
        self.assertLessEqual({'coordinates_path', 'topology_path'}, keys)
        self.assertNotIn('nonbonded_method', keys)
        self.assertNotIn('paths', standalone.fields)
        standalone.deleteLater()

    def test_input_file_list_has_no_checkbox(self):
        from gmxtransplant.gui.editor import ConfigEditor
        for mode in ('protein', 'lig', 'chl'):
            editor = ConfigEditor(example_document(mode), mode)
            self.assertEqual(editor.fields['paths'].kind, 'required')
            self.assertIsNone(editor.fields['paths'].toggle)
            self.assertIn('paths', editor.value())
            editor.deleteLater()

    def test_charmprot_keep_lipids_option(self):
        from gmxtransplant.gui.editor import ConfigEditor
        editor = ConfigEditor({'charmprot': {'reference': '/a', 'transplant': '/b'}}, 'charmprot')
        keep = next(r for r in editor.fields['charmprot'].editor.rows if r.key == 'keep_lipids')
        self.assertEqual(keep.kind, 'boolean')
        self.assertNotIn('keep_lipids', editor.value()['charmprot'])
        keep.editor.scalar.setChecked(True)
        self.assertTrue(editor.value()['charmprot']['keep_lipids'])
        editor.deleteLater()

    def test_mode_switch_preserves_unsaved_edits_and_errors_stay_in_window(self):
        from gmxtransplant.gui.window import MainWindow
        window = MainWindow()
        window.output_root.clear()
        # Buttons stay clickable; the reason a run cannot start is shown inline.
        self.assertTrue(window.run_button.isEnabled())
        self.assertTrue(all(b.isEnabled() for b in window.example_buttons))
        with patch('PySide6.QtWidgets.QMessageBox.warning') as popup:
            window.run_example('charmprot')
            popup.assert_not_called()
        self.assertTrue(window.banner.isVisibleTo(window))
        self.assertIn('output folder', window.banner_text.text())
        self.assertFalse(hasattr(window, 'base'))
        before = window.current_raw()
        window.mode_combo.setCurrentIndex(window.mode_combo.findData('chl'))
        window.mode_combo.setCurrentIndex(window.mode_combo.findData('charmprot'))
        self.assertEqual(window.current_raw(), before)
        window.deleteLater()

    def test_output_folder_starts_at_the_launch_directory(self):
        from PySide6.QtCore import QSettings
        from gmxtransplant.gui.window import MainWindow
        QSettings().setValue('output_root', '/somewhere/remembered')
        with tempfile.TemporaryDirectory() as launch:
            previous = os.getcwd()
            os.chdir(launch)
            try:
                window = MainWindow()
                self.assertEqual(Path(window.output_root.text()).resolve(), Path(launch).resolve())
                window.close()
                self.assertFalse(QSettings().contains('output_root'))
                window.deleteLater()
            finally:
                os.chdir(previous)

    def test_current_mode_is_described_under_the_mode_selector(self):
        from gmxtransplant.gui.window import MainWindow, MODE_DESCRIPTIONS
        window = MainWindow()
        for mode, name in MODES.items():
            window.mode_combo.setCurrentIndex(window.mode_combo.findData(mode))
            self.assertIn(name, window.mode_info.text())
            self.assertIn(MODE_DESCRIPTIONS[mode], window.mode_info.text())
            self.assertFalse(window.output_hint.isVisibleTo(window))
        window.deleteLater()

    def test_viewer_button_launches_the_found_program_with_the_scene(self):
        from gmxtransplant.gui import window as window_module
        from viewers import Viewer
        vmd = Viewer('vmd', '/Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64', {'VMDDIR': '/Applications/VMD.app/Contents/vmd'})
        with patch.object(window_module, 'find_viewer', side_effect=lambda name: vmd if name == 'vmd' else None):
            window = window_module.MainWindow()
        self.assertEqual(list(window.viewer_buttons), ['vmd'])
        with tempfile.TemporaryDirectory() as out:
            (Path(out) / 'view.vmd').write_text('')
            window.last_job = {'directory': out}
            with patch.object(window_module.subprocess, 'Popen') as popen:
                window.open_viewer('vmd')
            # The first process is VMD; on POSIX a second one holds VMD's input open.
            started = popen.call_args_list[0]
            self.assertEqual(started.args[0], [vmd.executable, '-e', str(Path(out) / 'view.vmd')])
            self.assertEqual(started.kwargs['env']['VMDDIR'], vmd.env['VMDDIR'])
            self.assertIsNot(started.kwargs['stdin'], window_module.subprocess.DEVNULL)
        window.deleteLater()

    def test_disabled_run_and_cancel_buttons_look_disabled(self):
        from gmxtransplant.gui.window import MainWindow, THEMES
        window = MainWindow()
        window.resize(1200, 800)
        window.show()
        for name in THEMES:
            window.theme_combo.setCurrentText(name)
            for control in (window.cancel_button, window.run_button):
                colours = []
                for enabled in (True, False):
                    control.setEnabled(enabled)
                    self.app.processEvents()
                    image = control.grab().toImage()
                    colours.append(image.pixelColor(image.width() // 6, image.height() // 2).name())
                self.assertNotEqual(colours[0], colours[1], f'{control.text()} in {name}')
        window.theme_combo.setCurrentText('Teal')
        window.deleteLater()

    def test_themes_apply_to_the_whole_application(self):
        from gmxtransplant.gui.window import MainWindow, THEMES
        window = MainWindow()
        self.assertEqual(len(THEMES), 3)
        for name, theme in THEMES.items():
            window.theme_combo.setCurrentText(name)
            self.assertIn(theme['accent'], self.app.styleSheet())
            self.assertEqual(self.app.palette().window().color().name(), theme['bg'])
        window.theme_combo.setCurrentText('Teal')
        window.deleteLater()

if __name__ == '__main__':
    unittest.main()
