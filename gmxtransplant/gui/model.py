"""Configuration/path handling shared by the GUI and its isolated worker (no Qt)."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import uuid
import webbrowser

import yaml
import config as pipeline_config
from config import path_kind, is_output_path as is_output, PLACEHOLDER_PATH_PREFIX

MODES = {'charmprot': 'CHARMM-GUI protein transplant', 'protein': 'Protein insertion', 'lig': 'Ligand replacement',
         'chl': 'Cholesterol restoration', 'addbinder': 'Add binder outside the membrane protein',
         'minimize': 'Minimization preparation'}
FOLDERS = {'charmprot': 'charmprot', 'protein': 'protein', 'lig': 'ligand', 'chl': 'cholesterol',
           'addbinder': 'addbinder', 'minimize': 'minimization'}
EXAMPLES = {'charmprot': ('charmprot', 'charmprot.yaml'), 'protein': ('protein_insertion', 'protein_replace.yaml'),
            'lig': ('ligand_replacement', 'ligand_replace.yaml'),
            'chl': ('cholesterol_restoration', 'cholesterol_restore.yaml'),
            'addbinder': ('addbinder', 'addbinder.yaml')}
# Further examples of a mode: key -> (mode, example directory, YAML, output folder).
EXTRA_EXAMPLES = {'addbinder_gprotein': ('addbinder', 'addbinder', 'addbinder_gprotein.yaml', 'addbinder_gprotein')}
# Every example card, in display order.
EXAMPLE_KEYS = [*EXAMPLES, *EXTRA_EXAMPLES]
# A folder every complete example ships with, proving its inputs are installed.
EXAMPLE_INPUT_FOLDER = {'charmprot': 'reference', 'addbinder': 'host'}
# Sections that are a mode's own self-contained configuration.
OWN_SECTION_MODES = ('charmprot', 'addbinder')
TEMPLATES = {**{k: v[1] for k, v in EXAMPLES.items()}, 'minimize': 'minimization.yaml'}
DIAGNOSTICS = {'raw_protonated_pdb_path': 'cholesterol_protonated_raw.pdb',
               'converted_pdb_path': 'cholesterol_charmm36.pdb',
               'aligned_reference_pdb_path': 'cholesterol_reference_aligned.pdb',
               'reference_placed_pdb_path': 'cholesterol_reference_placed.pdb',
               'merged_pdb_path': 'merged_system.pdb'}


def data_roots():
    return [Path(__file__).resolve().parents[2],
            Path(sys.prefix) / 'share/gmxtransplant',
            Path(sysconfig.get_path('data')) / 'share/gmxtransplant']


def template_path(mode):
    for root in data_roots():
        for path in [root / TEMPLATES[mode], root / 'examples' / TEMPLATES[mode]]:
            if path.is_file():
                return path
    raise FileNotFoundError('Configuration templates are missing. Reinstall GMXtransplant.')


def example_mode(key):
    """The mode an example card runs in."""
    return EXTRA_EXAMPLES[key][0] if key in EXTRA_EXAMPLES else key


def example_folder(key):
    """The folder under <output>/examples/ an example card writes to."""
    return EXTRA_EXAMPLES[key][3] if key in EXTRA_EXAMPLES else FOLDERS[key]


def example_path(key):
    mode = example_mode(key)
    directory, filename = EXTRA_EXAMPLES[key][1:3] if key in EXTRA_EXAMPLES else EXAMPLES[key]
    for root in data_roots():
        path = root / 'examples' / directory / filename
        if path.is_file() and (path.parent / EXAMPLE_INPUT_FOLDER.get(mode, 'environment')).is_dir():
            return path
    raise FileNotFoundError('Example inputs are missing. Reinstall the current GMXtransplant package.')


DOCUMENTATION_PDF = 'GMXtransplant.pdf'


def documentation_path():
    """The documentation PDF shipped with the installed package."""
    for root in data_roots():
        for candidate in (root / 'docs' / DOCUMENTATION_PDF, root / DOCUMENTATION_PDF):
            if candidate.is_file():
                return candidate
    raise FileNotFoundError(
        f'The documentation PDF ({DOCUMENTATION_PDF}) was not found in this GMXtransplant '
        'installation. Reinstall the package to restore it.')


def running_in_wsl():
    return 'WSL_DISTRO_NAME' in os.environ or 'microsoft' in os.uname().release.lower()


def _viewer_command(path):
    """The command that opens a file with this desktop's default application."""
    if sys.platform == 'darwin':
        return ['open', str(path)]
    if running_in_wsl():
        # Under WSL the Windows default viewer is the one the user expects;
        # a Linux xdg-open association may be an image editor or missing.
        if shutil.which('wslview'):
            return ['wslview', str(path)]
        if shutil.which('explorer.exe') and shutil.which('wslpath'):
            windows = subprocess.run(['wslpath', '-w', str(path)], capture_output=True,
                                     text=True, timeout=5).stdout.strip()
            if windows:
                return ['explorer.exe', windows]
    return ['xdg-open', str(path)] if shutil.which('xdg-open') else None


def open_with_default_app(path):
    """Open a file or folder with the default application, without waiting (never blocks the GUI)."""
    if os.name == 'nt':
        os.startfile(str(path))
        return True
    command = _viewer_command(path)
    if command is None:
        return False
    subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL, start_new_session=True)
    return True


_system_viewer = open_with_default_app


def open_documentation(viewer=None, browser=None):
    """Open the documentation PDF: system viewer, then web browser, then report why not."""
    path = documentation_path()
    reasons = []
    for label, action in (('the system PDF viewer', viewer or _system_viewer),
                          ('the default web browser', browser or (lambda p: webbrowser.open(p.as_uri())))):
        try:
            if action(path):
                return path
            reasons.append(f'{label} is not available')
        except Exception as exc:
            reasons.append(f'{label} failed: {exc}')
    raise RuntimeError(f'Could not open the documentation at {path}: ' + '; '.join(reasons) + '.')


def read_yaml(text):
    raw = yaml.load(text, Loader=pipeline_config._UniqueKeySafeLoader)
    if not isinstance(raw, dict) or not raw:
        raise ValueError('The configuration must be a non-empty YAML mapping.')
    if 'reference' in raw and 'transplant' in raw and 'charmprot' not in raw:
        common = {k: v for k, v in raw.items() if k in {'paths', 'minimization'}}
        raw = {**common, 'charmprot': {k: v for k, v in raw.items() if k not in common}}
    return raw


def load_yaml(path):
    return read_yaml(Path(path).read_text(encoding='utf-8'))


def dump_yaml(raw):
    return yaml.safe_dump(raw, sort_keys=False, allow_unicode=True)


def defaults(mode):
    if mode == 'charmprot':
        from charmprot import CharmProtSpec
        return {'charmprot': asdict(CharmProtSpec(reference='', transplant='')),
                'minimization': asdict(pipeline_config.MinimizationSpec())}
    if mode == 'addbinder':
        from addbinder import AddBinderSpec
        return {'addbinder': asdict(AddBinderSpec(host='', binder_coordinates='', binder_itp=''))}
    if mode == 'minimize':
        return {'paths': {}, 'minimization': asdict(pipeline_config.MinimizationSpec())}
    values = asdict(pipeline_config.Config())
    selected = {k: v for k, v in values.items() if k in pipeline_config._MODE_SECTIONS[mode]}
    for section, value in selected.items():
        if isinstance(value, dict) and section in pipeline_config._SECTION_KEYS:
            selected[section] = {k: v for k, v in value.items()
                                 if k in pipeline_config._SECTION_KEYS[section]}
    return {'paths': {}, **selected}


def resolve_inputs(raw, example_base=None):
    """Resolve ${...} references, then require every input path to be explicit.

    `example_base` is only ever the folder of a configuration shipped with the
    package. A document the user wrote must name each input in full: no path is
    resolved against a shared input folder or the working directory.
    """
    resolved = pipeline_config._resolve_path_references(deepcopy(raw))
    resolved.pop('paths', None)
    inputs, relative = [], []

    def visit(value, path=()):
        if isinstance(value, dict):
            return {k: visit(v, path + (k,)) for k, v in value.items()}
        if isinstance(value, list):
            return [visit(v, path + (i,)) for i, v in enumerate(value)]
        if isinstance(value, str) and value and path_kind(path) and not is_output(path):
            if path_kind(path) == 'executable' and '/' not in value:
                return value
            p = Path(value).expanduser()
            if not p.is_absolute():
                if example_base is None:
                    relative.append('.'.join(map(str, path)) + f' = {value!r}')
                    return value
                p = Path(example_base) / p
            p = p.resolve()
            inputs.append(p)
            return str(p)
        return value

    result = visit(resolved)
    if relative:
        raise ValueError(
            'Every input must be a complete path. Use the Browse buttons to choose these '
            'inputs; nothing is resolved against a shared input folder:\n  '
            + '\n  '.join(sorted(relative)))
    return result, inputs


def example_document(mode):
    """A bundled example with its inputs rewritten as the absolute paths they are."""
    path = example_path(mode)
    raw = deepcopy(load_yaml(path))
    registry = raw.get('paths')
    if isinstance(registry, dict):
        kinds = _registry_kinds(raw)
        for key, value in list(registry.items()):
            if (isinstance(value, str) and value and kinds.get(key) == 'input'
                    and not Path(value).expanduser().is_absolute()):
                registry[key] = str((path.parent / value).resolve())
    resolved, _ = resolve_inputs(raw, example_base=path.parent)
    for section, value in resolved.items():
        raw[section] = value
    return raw


def _registry_kinds(raw):
    """Classify each ${...} registry key by how the document actually uses it."""
    kinds = {}

    def visit(value, path=()):
        if isinstance(value, dict):
            for key, item in value.items():
                visit(item, path + (key,))
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, path + (index,))
        elif isinstance(value, str) and value.startswith('${') and value.endswith('}'):
            if path_kind(path):
                name = value[2:-1]
                kind = 'output' if is_output(path) else 'input'
                kinds[name] = 'output' if kinds.get(name) == 'output' else kind

    visit({k: v for k, v in raw.items() if k != 'paths'})
    return kinds


def template_document(mode):
    """A blank starting document: relative example inputs become empty fields."""
    raw = deepcopy(load_yaml(template_path(mode)))
    kinds = _registry_kinds(raw)
    registry = raw.get('paths')
    if isinstance(registry, dict):
        for key, value in list(registry.items()):
            if isinstance(value, str) and value and kinds.get(key) != 'output' and (
                    value.startswith(PLACEHOLDER_PATH_PREFIX)
                    or not Path(value).expanduser().is_absolute()):
                registry[key] = ''

    def blank(value, path=()):
        if isinstance(value, dict):
            return {k: blank(v, path + (k,)) for k, v in value.items()}
        if isinstance(value, list):
            return [blank(v, path + (i,)) for i, v in enumerate(value)]
        if (isinstance(value, str) and value and not value.startswith('${')
                and path_kind(path) and not is_output(path)
                and path_kind(path) != 'executable'
                and (value.startswith(PLACEHOLDER_PATH_PREFIX)
                     or not Path(value).expanduser().is_absolute())):
            return ''
        return value

    return {key: (value if key == 'paths' else blank(value, (key,)))
            for key, value in raw.items()}


def destination(output_root, mode, example=False):
    """example: False, True (the mode's example folder) or the name of an example's own folder."""
    if not str(output_root).strip():
        raise ValueError('Choose an output folder first.')
    root = Path(output_root).expanduser().absolute()
    folder = example if isinstance(example, str) else FOLDERS[mode]
    result = root / 'examples' / folder if example else root / FOLDERS[mode]
    # A symlink must never redirect a mode folder outside the selected root.
    if result.resolve() != root.resolve().joinpath(*result.relative_to(root).parts):
        raise ValueError('The output mode folder must not be a symbolic link.')
    return result.resolve()


def _filename(value, fallback):
    name = Path(str(value or fallback)).name
    if name in {'', '.', '..'}:
        raise ValueError('Output filenames must have a filename, not just a directory.')
    return name


def build_job(raw, output_root, mode, example_base=None, prepare=False, example_folder=None):
    """Build a CLI-ready snapshot, enforcing the GUI's flat, contained output layout."""
    out = destination(output_root, mode, example_folder or example_base is not None)
    if mode == 'charmprot' and 'charmprot' not in raw:
        raw = read_yaml(dump_yaml(raw))
    values, inputs = resolve_inputs(raw, example_base)
    for p in inputs:
        if p == out or out in p.parents or (p.is_dir() and p in out.parents):
            raise ValueError(f'Input and output locations overlap: {p}. Choose a separate output folder.')
    if mode == 'addbinder':
        section = values.get('addbinder')
        if not isinstance(section, dict):
            raise ValueError('addbinder must be a mapping of settings.')
        section['output_dir'] = str(out)
        poses = section.get('poses') or [{}]
        names = [(p.get('name') if isinstance(p, dict) else None) or f'binderpose{i + 1}'
                 for i, p in enumerate(poses)]
        count = section.get('random_poses', 3)
        names += [f'binderpose{len(names) + i + 1}' for i in range(count if isinstance(count, int) else 0)]
        for name in names:
            if not isinstance(name, str) or Path(name).name != name or name in {'.', '..'}:
                raise ValueError(f'Pose names must be simple folder names: {name!r}')
        return {'id': uuid.uuid4().hex, 'mode': mode, 'directory': str(out), 'config': values,
                'inputs': [str(p) for p in inputs],
                'artifacts': ['run.yaml', 'run.log', 'run-status.json', 'summary.txt', 'summary.json',
                              'visualization', 'view.pml', 'view.vmd', 'view.pse', *dict.fromkeys(names)]}
    if mode == 'charmprot':
        from charmprot import GENERATED
        values['charmprot']['output_dir'] = str(out)
        if prepare or values.get('minimization'):
            mini = values.setdefault('minimization', {})
            mini['output_dir'] = str(out / 'openmm_minimization')
            mini['output_gro_path'] = _filename(mini.get('output_gro_path'), 'minimized.gro')
            mini['report_path'] = _filename(mini.get('report_path'), 'minimization_report')
            if prepare:
                mini['enabled'] = True
        return {'id': uuid.uuid4().hex, 'mode': mode, 'directory': str(out), 'config': values,
                'inputs': [str(p) for p in inputs],
                'artifacts': ['run.yaml', 'run.log', 'run-status.json', 'openmm_minimization', *sorted(GENERATED)]}
    for section in (['minimization'] if mode == 'minimize' else
                    ['output', 'ndx', 'topology', 'minimization'] + (['cholesterol'] if mode == 'chl' else [])):
        if values.get(section) is None:
            values[section] = {}
        if not isinstance(values[section], dict):
            raise ValueError(f'{section} must be a mapping of settings.')
    artifacts = ['run.yaml', 'run.log', 'run-status.json', 'visualization', 'view.pml', 'view.vmd', 'view.pse']
    reserved = set(artifacts + ['.gui-manifest.json', '.gui.lock', 'toppar', 'openmm_minimization'])

    def output_file(section, key, fallback, report=False):
        mapping = values.setdefault(section, {})
        name = _filename(mapping.get(key), fallback)
        names = [name + '.txt', name + '.json'] if report else [name]
        if any(n in reserved for n in names):
            raise ValueError(f'Output filename is reserved or used twice: {name}')
        reserved.update(names)
        artifacts.extend(names)
        mapping[key] = str(out / name)

    if mode != 'minimize':
        for key, name in [('pdb_path', 'step5_input.pdb'), ('gro_path', 'step5_input.gro'),
                          ('inspection_pdb_path', 'aligned_replacement_inspection.pdb')]:
            output_file('output', key, name)
        output_file('output', 'report_path', FOLDERS[mode] + '_report', report=True)
        output_file('ndx', 'output_path', 'index.ndx')
        if 'topol.top' in reserved:
            raise ValueError('topol.top is reserved for the generated topology.')
        reserved.add('topol.top')
        values.setdefault('topology', {})['output_dir'] = str(out)
        artifacts.extend(['topol.top', 'toppar'])
        if mode == 'chl':
            for key, name in DIAGNOSTICS.items():
                output_file('cholesterol', key, name)
            repack = values['cholesterol'].get('repack') or {}
            artifacts.append(repack.get('output_dir') or 'repack')
    minimization = values.setdefault('minimization', {})
    minimization['output_dir'] = str(out / 'openmm_minimization')
    # The exported runner's products must stay in its bundle.
    minimization['output_gro_path'] = _filename(minimization.get('output_gro_path'), 'minimized.gro')
    minimization['report_path'] = _filename(minimization.get('report_path'), 'minimization_report')
    if prepare:
        minimization['enabled'] = True
    artifacts.append('openmm_minimization')
    return {'id': uuid.uuid4().hex, 'mode': mode, 'directory': str(out), 'config': values,
            'inputs': [str(p) for p in inputs], 'artifacts': artifacts}


def validate_job(job, check_paths=False):
    with tempfile.TemporaryDirectory(prefix='gmxtransplant-check-') as tmp:
        path = Path(tmp) / 'config.yaml'
        path.write_text(dump_yaml(job['config']), encoding='utf-8')
        if job['mode'] == 'charmprot':
            from charmprot import load_charmprot_config
            cfg = load_charmprot_config(path, check_paths=check_paths)
            if cfg.minimization.enabled:
                from minimization_bundle import validate_preparation
                validate_preparation(cfg.minimization)
            return cfg
        if job['mode'] == 'addbinder':
            from addbinder import load_addbinder_config
            return load_addbinder_config(path, check_paths=check_paths)
        if job['mode'] == 'minimize':
            cfg = pipeline_config.load_minimization_config(str(path), check_paths=check_paths)
            from minimization_bundle import validate_preparation
            validate_preparation(cfg)
        else:
            cfg = pipeline_config.load_config(str(path), mode=job['mode'], check_paths=check_paths)
            if cfg.minimization.enabled:
                from run_pipeline import _post_pipeline_minimization_spec
                from minimization_bundle import validate_preparation
                validate_preparation(_post_pipeline_minimization_spec(cfg))
        return cfg
