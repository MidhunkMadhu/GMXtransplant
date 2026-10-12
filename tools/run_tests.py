#!/usr/bin/env python3
"""Run the GMXtransplant test suite from the checkout.

    python3 tools/run_tests.py            # everything (about 10 minutes)
    python3 tools/run_tests.py --fast     # skip the @slow tests that build whole systems (about 1 minute)
    python3 tools/run_tests.py test_salt  # a pattern, as unittest -k
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--fast", action="store_true", help="skip tests marked @slow (GMXT_FAST=1)")
parser.add_argument("patterns", nargs="*", help="only tests matching these (unittest -k)")
args = parser.parse_args()
env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(root / "src"), str(root)]),
           QT_QPA_PLATFORM=os.environ.get("QT_QPA_PLATFORM", "offscreen"))
if args.fast:
    env["GMXT_FAST"] = "1"
command = [sys.executable, "-B", "-m", "unittest", "discover", "-s", "tests", "-t", "."]
for pattern in args.patterns:
    command += ["-k", pattern]
sys.exit(subprocess.call(command, cwd=root, env=env))
