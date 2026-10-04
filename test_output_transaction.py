import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from output import OutputError, _promote_staged_outputs


class OutputTransactionTests(unittest.TestCase):
    def test_second_promotion_failure_restores_both_old_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "final.pdb"
            second = root / "final.gro"
            staged_first = root / ".final.pdb.staged"
            staged_second = root / ".final.gro.staged"
            first.write_text("old pdb\n")
            second.write_text("old gro\n")
            staged_first.write_text("new pdb\n")
            staged_second.write_text("new gro\n")

            real_replace = os.replace

            def fail_second(source, destination):
                if str(source) == str(staged_second):
                    raise OSError("simulated second rename failure")
                return real_replace(source, destination)

            with patch("output.os.replace", side_effect=fail_second):
                with self.assertRaisesRegex(OutputError, "prior outputs were restored"):
                    _promote_staged_outputs(
                        [(str(staged_first), str(first)), (str(staged_second), str(second))]
                    )

            self.assertEqual(first.read_text(), "old pdb\n")
            self.assertEqual(second.read_text(), "old gro\n")
            self.assertFalse(any(root.glob("*.backup.*")))


if __name__ == "__main__":
    unittest.main()
