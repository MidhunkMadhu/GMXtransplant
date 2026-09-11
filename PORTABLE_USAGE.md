# Portable usage

The project is now installable; the full installation, input requirements, YAML
guide, and examples are maintained in [README.md](README.md).

In short:

```bash
cd /path/to/GMXtransplant
python3 -m pip install .
gmxtransplant --help
```

Install OpenMM support when energy minimization is wanted:

```bash
python3 -m pip install '.[minimize]'
python3 -m openmm.testInstallation
```

Write one of the bundled YAML files to a separate job directory, edit the
grouped `paths:` block and system-specific selections, validate it, then run it:

```bash
gmxtransplant --show-example receptor > receptor_replace.yaml
gmxtransplant --mode receptor -i receptor_replace.yaml --dry-run
gmxtransplant --mode receptor -i receptor_replace.yaml
gmxtransplant --mode receptor -i receptor_replace.yaml --minimize
```

Use mode `lig` or `chl` with the matching example. NumPy, SciPy, PyYAML, and
MDAnalysis are installed from package metadata. Open Babel is an external
requirement for CHARMM36 cholesterol conversion; PyMOL is optional for the two
PyMOL alignment methods. OpenMM is optional and used only by `--minimize`.
