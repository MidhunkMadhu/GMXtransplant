"""Output publication and report regression checks for cholesterol mode."""
import io
import unittest
from contextlib import ExitStack, redirect_stdout, redirect_stderr
from types import SimpleNamespace
from unittest.mock import patch

from config import Config
from run_pipeline import _run_cholesterol
from topology import TopologyError


class CholesterolOutputTests(unittest.TestCase):
    def run_workflow(self, topology_error=None):
        cfg = Config(mode="chl")
        cfg.topology.enabled = topology_error is not None
        cfg.ndx.enabled = False
        result = SimpleNamespace(
            n_experimental_cholesterols=1, n_experimental_cholesterol_atoms=74,
            heavy_atom_audit=[], alignment=SimpleNamespace(rmsd_after=0.0),
            clash_result=SimpleNamespace(counts_by_class={}),
            composition=SimpleNamespace(
                lipid_residues_removed=[], ions_removed=[],
                estimated_final_concentration_molar=None, warnings=[]),
            diagnostic_paths=[],
        )
        with ExitStack() as stack:
            stack.enter_context(redirect_stdout(io.StringIO()))
            stack.enter_context(redirect_stderr(io.StringIO()))
            stack.enter_context(patch("run_pipeline.restore_cholesterols",
                                     return_value=(object(), object(), result)))
            stack.enter_context(patch("run_pipeline.assemble_topology",
                                     side_effect=topology_error))
            writer = stack.enter_context(patch("run_pipeline.write_outputs",
                                               return_value={"checked": "both formats"}))
            builder = stack.enter_context(patch("run_pipeline.build_cholesterol_report"))
            report = stack.enter_context(patch("run_pipeline.write_cholesterol_report",
                                               return_value=("report.txt", "report.json")))
            index = stack.enter_context(patch("run_pipeline.write_ndx"))
            status = _run_cholesterol(cfg)
        return status, writer, builder, report, index

    def test_topology_failure_reports_error_without_publishing_coordinates(self):
        status, writer, builder, report, index = self.run_workflow(
            TopologyError("coordinate/topology atom count mismatch"))
        self.assertEqual(status, 6)
        writer.assert_not_called()
        index.assert_not_called()
        report.assert_called_once()
        self.assertIn("atom count mismatch", builder.call_args.kwargs["topology_error"])

    def test_success_retains_coordinate_validation_in_report(self):
        status, writer, builder, report, _ = self.run_workflow()
        self.assertEqual(status, 0)
        writer.assert_called_once()
        report.assert_called_once()
        self.assertEqual(builder.call_args.kwargs["output_validation"],
                         {"checked": "both formats"})
        self.assertTrue(builder.call_args.kwargs["coordinates_written"])


if __name__ == "__main__":
    unittest.main()
