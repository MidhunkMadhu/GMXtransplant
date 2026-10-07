"""Allow ``python -m gmxtransplant`` to run the installed command."""

from run_pipeline import cli


raise SystemExit(cli())
