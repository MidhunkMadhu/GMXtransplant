"""Open Babel installed by pip is found even when its folder is not on PATH."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from charmm36_cholesterol import PipelineError, check_obabel


class ObabelLookupTests(unittest.TestCase):
    def test_found_next_to_python_when_not_on_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake_python = Path(tmp) / "python"
            fake_python.write_text("")
            obabel = Path(tmp) / "obabel"
            obabel.write_text("#!/bin/sh\n")
            obabel.chmod(0o755)
            with patch.dict(os.environ, {"PATH": "/nonexistent"}), patch.object(sys, "executable", str(fake_python)):
                self.assertEqual(check_obabel("obabel"), str(obabel))

    def test_missing_obabel_explains_how_to_install(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {"PATH": "/nonexistent"}), \
                patch.object(sys, "executable", str(Path(tmp) / "python")):
            with self.assertRaisesRegex(PipelineError, "openbabel-wheel.*brew install open-babel"):
                check_obabel("obabel")


if __name__ == "__main__":
    unittest.main()
