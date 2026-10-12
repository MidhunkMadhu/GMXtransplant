"""Configuration basics: the error type, the strict YAML loader and the path rules."""

from __future__ import annotations

import os
import math
import re
import sys
import sysconfig
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import yaml
from yaml.nodes import MappingNode

from mode_registry import EDIT_MODES, path_kind as edit_path_kind


class ConfigError(ValueError):
    pass


class _UniqueKeySafeLoader(yaml.SafeLoader):
    """Safe YAML loader that refuses ambiguous duplicate mapping keys."""


def _construct_unique_mapping(loader, node, deep=False):
    if not isinstance(node, MappingNode):
        raise ConfigError("expected a YAML mapping")
    loader.flatten_mapping(node)
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in mapping
        except TypeError as exc:
            raise ConfigError(
                f"YAML mapping key at line {key_node.start_mark.line + 1} "
                "is not a scalar value"
            ) from exc
        if duplicate:
            raise ConfigError(
                f"Duplicate YAML key {key!r} at line "
                f"{key_node.start_mark.line + 1}; duplicate keys are unsafe because "
                "they make the active option ambiguous"
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeySafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def _resolve_path_references(raw):
    """Resolve whole-value ${name} references from the literal paths registry."""
    paths = _require_mapping(raw.get("paths"), "paths")
    for key, value in paths.items():
        if not isinstance(key, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", key):
            raise ConfigError("paths keys must use lower-case snake_case names")
        if value is not None and not isinstance(value, str):
            raise ConfigError(f"paths.{key} must be a path string or null")
        if isinstance(value, str) and "${" in value:
            raise ConfigError(f"paths.{key} must be a literal path, not a reference")

    def resolve(value):
        if isinstance(value, dict):
            return {key: resolve(item) for key, item in value.items()}
        if isinstance(value, list):
            return [resolve(item) for item in value]
        if isinstance(value, str) and "${" in value:
            match = re.fullmatch(r"\$\{([a-z][a-z0-9_]*)\}", value)
            if not match:
                raise ConfigError(f"Invalid path reference {value!r}; use a whole value like '${{name}}'")
            name = match.group(1)
            if name not in paths:
                raise ConfigError(f"Unknown path reference {value!r}: paths.{name} is missing")
            return paths[name]
        return value

    return {key: value if key == "paths" else resolve(value) for key, value in raw.items()}


# ---------------------------------------------------------------------------
# Explicit paths
#
# No path in this program defaults to a location on the host system. Every
# input path must be spelled out in full: nothing is resolved against the
# working directory, and there is no shared input root to assemble paths from.
# Generated files are the single exception -- a bare output filename lands in
# the run's output folder, which is the --output directory when one is given
# and otherwise the directory the command was started from.
#
# Configurations shipped inside the installed package are the other exception:
# their relative paths resolve against the example folder that ships with the
# code, so the bundled examples run without referring to the user's system.
# ---------------------------------------------------------------------------

# Bundled templates use this prefix where a real path belongs, so the GUI can
# present an empty field instead of a path that does not exist.
PLACEHOLDER_PATH_PREFIX = "/path/to/"

# Cholesterol-mode conversion diagnostics; generated, not supplied.
CHOLESTEROL_DIAGNOSTIC_PATHS = frozenset({
    "raw_protonated_pdb_path", "converted_pdb_path", "aligned_reference_pdb_path",
    "reference_placed_pdb_path", "merged_pdb_path",
})

# Names of files or folders written inside an output folder that is itself anchored.
_UNANCHORED_OUTPUTS = {("minimization", "output_gro_path"), ("minimization", "report_path"),
                       ("cholesterol", "repack", "output_dir"), ("popc", "equilibration", "output_dir"),
                       *((key, "equilibration", "output_dir") for key in EDIT_MODES)}


def path_kind(path):
    """Classify a configuration field as a file, directory or executable path."""
    key = path[-1] if path else ""
    if isinstance(key, int):
        key = path[-2] if len(path) > 1 else ""
    if not isinstance(key, str):
        return None
    if path[:1] == ("charmprot",) and key in {"reference", "transplant"}:
        return "directory"
    if path[:1] == ("charmprot",) and key in {"reference_gro", "transplant_gro"}:
        return "file"
    if path[:1] in (("addbinder",), ("popc",)) and key == "host":
        return "directory"
    if path[:1] == ("popc",) and key in {"host_gro", "reference_gro"}:
        return "file"
    if path[:1] == ("addbinder",) and key in {
            "host_gro", "binder_coordinates", "binder_itp", "binder_forcefield"}:
        return "file"
    if path and path[0] in EDIT_MODES:
        return edit_path_kind(path[0], key)
    if key == "obabel_command":
        return "executable"
    if key.endswith(("_dir", "_toppar")) or key in {"include_dir", "system_toppar"}:
        return "directory"
    if (path and path[0] == "paths") or key.endswith(("_path", "_paths")) or key in {
        "path", "refgro", "reference_topol", "protein_template_top", "topology_path"}:
        return "file"
    return None


def is_output_path(path):
    """True when this field names something the pipeline writes."""
    if not path:
        return False
    if path[0] in EDIT_MODES:
        return path[1:] == ("output_dir",)
    return (path[:2] in (("charmprot", "output_dir"), ("addbinder", "output_dir"), ("popc", "output_dir"))
            or path[0] == "output" or path[:2] in [
        ("topology", "output_dir"), ("ndx", "output_path"),
        ("minimization", "output_dir"), ("minimization", "output_gro_path"),
        ("minimization", "report_path")]
        or (path[0] == "cholesterol" and path[-1] in CHOLESTEROL_DIAGNOSTIC_PATHS))


def _package_roots():
    """Directories in which the installed package keeps its own example inputs."""
    return [
        Path(__file__).resolve().parent.parent,
        Path(sys.prefix) / "share" / "gmxtransplant",
        Path(sysconfig.get_path("data")) / "share" / "gmxtransplant",
    ]


def bundled_config_root(config_path):
    """Directory a *bundled* configuration resolves its relative inputs against.

    Returns None for every configuration outside the installed package, which
    must therefore give each input path in full.
    """
    if not config_path:
        return None
    try:
        resolved = Path(config_path).resolve()
    except OSError:
        return None
    for root in _package_roots():
        try:
            root = root.resolve()
        except OSError:
            continue
        if root.is_dir() and root in resolved.parents:
            return resolved.parent
    return None


def resolve_output_root(output_root=None) -> Path:
    """The one folder generated files land in: --output, else the launch directory."""
    root = Path(output_root).expanduser() if output_root else Path.cwd()
    return root.absolute()


def anchor_config_paths(raw, config_path=None, output_root=None, input_root=None):
    """Return `raw` with every path made explicit, rejecting relative inputs.

    `input_root` is only ever a folder that travels with the configuration --
    a bundled example, or a prepared minimization folder -- never a location
    chosen from the host system.
    """
    example_root = input_root or bundled_config_root(config_path)
    out_root = resolve_output_root(output_root)
    relative_inputs = []

    def visit(value, path):
        if isinstance(value, dict):
            return {key: visit(item, path + (key,)) for key, item in value.items()}
        if isinstance(value, list):
            return [visit(item, path + (index,)) for index, item in enumerate(value)]
        if not isinstance(value, str) or not value.strip():
            return value
        kind = path_kind(path)
        if kind is None or path in _UNANCHORED_OUTPUTS:
            return value
        # A bare command name is looked up on PATH, not treated as a path.
        if kind == "executable" and os.sep not in value and "/" not in value:
            return value
        candidate = Path(value).expanduser()
        if candidate.is_absolute():
            return str(candidate)
        if is_output_path(path):
            return str(out_root / candidate)
        if example_root is not None:
            return str((example_root / candidate).resolve())
        relative_inputs.append(f"{'.'.join(map(str, path))} = {value!r}")
        return value

    resolved = {key: (value if key == "paths" else visit(value, (key,)))
                for key, value in raw.items()}
    if relative_inputs:
        raise ConfigError(
            "Every input path must be given in full. GMXtransplant does not resolve an input "
            "against the working directory or any shared input folder, so these relative "
            "paths cannot be used:\n  - " + "\n  - ".join(sorted(relative_inputs))
        )
    return resolved


def anchor_output_defaults(cfg, output_root=None):
    """Place output names that came from dataclass defaults in the run's output folder.

    Without this a default such as ``step5_input.gro`` would silently fall back
    to the working directory rather than the folder the run writes to.
    """
    root = resolve_output_root(output_root)

    def place(owner, attribute):
        value = getattr(owner, attribute, None)
        if isinstance(value, str) and value and not Path(value).expanduser().is_absolute():
            setattr(owner, attribute, str(root / Path(value).expanduser()))

    for attribute in ("pdb_path", "gro_path", "report_path", "inspection_pdb_path"):
        place(cfg.output, attribute)
    place(cfg.ndx, "output_path")
    place(cfg.topology, "output_dir")
    place(cfg.minimization, "output_dir")
    for attribute in sorted(CHOLESTEROL_DIAGNOSTIC_PATHS):
        place(cfg.cholesterol, attribute)
    return cfg


def _is_finite_number(value) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
    )

def _require_mapping(value, label: str) -> dict:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a mapping")
    return value
