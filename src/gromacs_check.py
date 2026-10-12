"""Check a built system with GROMACS itself: `gmx grompp` with zero warnings allowed.

grompp reads the coordinates, the whole topology and every force-field parameter,
so it catches what the Python checks cannot: a bonded term without parameters, an
atom-name or count mismatch, a charge-group or cut-off problem. The check uses a
plain minimization mdp (no restraints), writes nothing next to the system, and is
skipped quietly when no GROMACS executable is found (`--gromacs-check auto`).
"""
from __future__ import annotations

from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile

MDP = """integrator              = steep
nsteps                  = 10
emtol                   = 1000.0
cutoff-scheme           = Verlet
nstlist                 = 10
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
constraints             = h-bonds
constraint_algorithm    = LINCS
"""


def find_gmx(explicit=None):
    """The GROMACS command: explicit, $GMXTRANSPLANT_GMX, $GMX, then gmx / gmx_mpi / gmx_d on PATH."""
    for candidate in (explicit, os.environ.get("GMXTRANSPLANT_GMX"), os.environ.get("GMX")):
        if candidate:
            return shlex.split(candidate)
    for name in ("gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d"):
        if shutil.which(name):
            return [name]
    return None


def check_folder(folder, gmx=None, timeout=900, top=None, gro=None, index=None):
    """Run grompp on folder/topol.top + folder/step5_input.gro (or the files given); returns a result dict."""
    folder = Path(folder).resolve()
    top = Path(top) if top else folder / "topol.top"
    gro = Path(gro) if gro else folder / "step5_input.gro"
    index = Path(index) if index else folder / "index.ndx"
    result = {"folder": str(folder), "status": "skipped"}
    if not top.is_file() or not gro.is_file():
        result["reason"] = "no topol.top and step5_input.gro"
        return result
    command = find_gmx(gmx)
    if command is None:
        result["reason"] = "no GROMACS executable found (set GMXTRANSPLANT_GMX)"
        return result
    with tempfile.TemporaryDirectory(prefix="gmxtransplant-grompp-") as tmp:
        mdp = Path(tmp) / "check.mdp"
        mdp.write_text(MDP)
        args = command + ["grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top), "-o", str(Path(tmp) / "check.tpr"),
                          "-po", str(Path(tmp) / "mdout.mdp"), "-maxwarn", "0"]
        if index.is_file():
            args += ["-n", str(index)]
        try:
            run = subprocess.run(args, cwd=tmp, capture_output=True, text=True, timeout=timeout)
        except (OSError, subprocess.TimeoutExpired) as exc:
            result.update(status="error", reason=str(exc))
            return result
    output = run.stdout + run.stderr
    notes = len(re.findall(r"^NOTE \d+", output, re.M))
    warnings_ = len(re.findall(r"^WARNING \d+", output, re.M))
    errors = re.findall(r"(?ms)^(?:ERROR \d+|Fatal error:).*?(?=^\s*$)", output)
    charge = re.search(r"System has non-zero total charge: ([-\d.e+]+)", output)
    result.update(status="passed" if run.returncode == 0 else "failed", command=" ".join(command + ["grompp"]),
                  notes=notes, warnings=warnings_, returncode=run.returncode,
                  errors=[e.strip()[:800] for e in errors][:5])
    if charge:
        result["net_charge_reported"] = float(charge[1])
    return result


def describe(result):
    if result["status"] == "skipped":
        return f"[gromacs] check skipped for {result['folder']}: {result['reason']}"
    if result["status"] == "passed":
        return (f"[gromacs] grompp passed for {result['folder']} ({result['warnings']} warnings, "
                f"{result['notes']} notes)")
    detail = "; ".join(result.get("errors") or [result.get("reason", "")])
    return f"[gromacs] grompp FAILED for {result['folder']}: {detail}"


def main(argv=None):
    import argparse
    import json
    import sys
    parser = argparse.ArgumentParser(prog="gmxtransplant check-gromacs",
                                     description="Run gmx grompp (no warnings allowed) on built system folders.")
    parser.add_argument("folders", nargs="+", help="folders with topol.top and step5_input.gro")
    parser.add_argument("--gmx", help="GROMACS command, e.g. 'gmx_mpi' or 'srun -n 1 gmx_mpi'")
    parser.add_argument("--json", action="store_true", help="print the results as JSON")
    args = parser.parse_args(argv)
    results = [check_folder(folder, args.gmx) for folder in args.folders]
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        for result in results:
            print(describe(result))
    return 1 if any(r["status"] in ("failed", "error") for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
