"""PyMOL and VMD are found where each system installs them, including macOS app bundles."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import viewers
from viewers import find_viewer


def executable(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n")
    path.chmod(0o755)
    return path


class ViewerLookupTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.env = patch.dict(os.environ, {"PATH": "/nonexistent"})
        self.env.start()
        for name in viewers.OVERRIDE.values():
            os.environ.pop(name, None)
        self.python = patch.object(sys, "executable", str(self.root / "env/bin/python"))
        self.python.start()

    def tearDown(self):
        self.python.stop()
        self.env.stop()
        self.tmp.cleanup()

    def test_macos_application_bundles(self):
        apps = self.root / "Applications"
        pymol = executable(apps / "PyMOL.app/Contents/MacOS/PyMOL")
        vmd = executable(apps / "VMD 1.9.4a57-arm64-Rev12.app/Contents/vmd/vmd_MACOSXARM64")
        with patch.object(sys, "platform", "darwin"), patch.object(viewers, "MAC_APPLICATION_FOLDERS", [apps]):
            self.assertEqual(find_viewer("pymol").executable, str(pymol))
            found = find_viewer("vmd")
            self.assertEqual(found.executable, str(vmd))
            self.assertEqual(found.env["VMDDIR"], str(vmd.parent))
            self.assertEqual(found.command("-e", "view.vmd"), [str(vmd), "-e", "view.vmd"])

    def test_bundles_are_only_searched_on_macos(self):
        apps = self.root / "Applications"
        executable(apps / "PyMOL.app/Contents/MacOS/PyMOL")
        with patch.object(sys, "platform", "linux"), patch.object(viewers, "MAC_APPLICATION_FOLDERS", [apps]):
            self.assertIsNone(find_viewer("pymol"))

    def test_environment_folder_and_explicit_override(self):
        sibling = executable(self.root / "env/bin/pymol")
        self.assertEqual(find_viewer("pymol").executable, str(sibling))
        custom = executable(self.root / "elsewhere/my-vmd")
        with patch.dict(os.environ, {"GMXTRANSPLANT_VMD": str(custom)}):
            self.assertEqual(find_viewer("vmd").executable, str(custom))
        with patch.dict(os.environ, {"GMXTRANSPLANT_VMD": str(self.root / "missing")}):
            self.assertIsNone(find_viewer("vmd"))


class LaunchTests(unittest.TestCase):
    def test_vmd_gets_an_open_input_pipe_and_pymol_does_not(self):
        import os
        import subprocess
        from unittest import mock
        from viewers import Viewer, launch
        calls = []

        class Process:
            pid = 4321

        def popen(command, **kwargs):
            calls.append((command, kwargs))
            return Process()
        with mock.patch('subprocess.Popen', side_effect=popen):
            launch(Viewer('pymol', '/x/pymol'), ['view.pml'], '/tmp')
            self.assertEqual(calls[-1][1]['stdin'], subprocess.DEVNULL)
            calls.clear()
            launch(Viewer('vmd', '/x/vmd', {'VMDDIR': '/x'}), ['-e', 'view.vmd'], '/tmp')
        if os.name != 'posix':
            return
        (vmd, vmd_kwargs), (watcher, watcher_kwargs) = calls
        self.assertEqual(vmd, ['/x/vmd', '-e', 'view.vmd'])
        # VMD quits at end of input: it reads a pipe that the watcher holds open while VMD runs.
        self.assertIsInstance(vmd_kwargs['stdin'], int)
        self.assertIsInstance(watcher_kwargs['stdout'], int)
        self.assertIn('kill -0 4321', watcher[-1])
        self.assertTrue(vmd_kwargs['start_new_session'] and watcher_kwargs['start_new_session'])
        self.assertEqual(vmd_kwargs['env']['VMDDIR'], '/x')


if __name__ == "__main__":
    unittest.main()
