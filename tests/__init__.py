"""GMXtransplant tests. Run all: python -m unittest discover -s tests (from the checkout, PYTHONPATH=src).

Tests that build whole example systems are marked @slow; GMXT_FAST=1 skips them
(tools/run_tests.py --fast).
"""
import os
import unittest

FAST = os.environ.get("GMXT_FAST", "").strip().lower() in {"1", "true", "yes"}


def slow(test):
    """Skip a test (or test class) that builds a full system when GMXT_FAST=1."""
    return unittest.skipIf(FAST, "slow test: unset GMXT_FAST to run it")(test)
