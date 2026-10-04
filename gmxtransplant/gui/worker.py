"""An isolated, cancellable CLI worker; no Qt dependency."""
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import fcntl
import json
import os
import shutil
import signal
import sys
import traceback

from .model import dump_yaml, validate_job


class Tee:
    def __init__(self, stream, log):
        self.stream, self.log = stream, log

    def write(self, text):
        self.stream.write(text)
        self.log.write(text)
        self.flush()
        return len(text)

    def flush(self):
        self.stream.flush()
        self.log.flush()


def cleanup(job):
    """Replace only known generated products; unrelated files are left alone."""
    root = Path(job['directory'])
    names = set(job['artifacts'])
    manifest = root / '.gui-manifest.json'
    if manifest.is_symlink():
        raise ValueError('The output manifest must not be a symbolic link.')
    if manifest.is_file():
        old = json.loads(manifest.read_text())
        names.update(old['artifacts'])
    paths = []
    for name in names:
        if not isinstance(name, str) or Path(name).name != name or name in {'', '.', '..', '.gui.lock'}:
            raise ValueError('Invalid generated-output manifest entry.')
        p = root / name
        if p.is_symlink():
            raise ValueError(f'Output is a symbolic link: {p}')
        for source in map(Path, job['inputs']):
            if source == p or p in source.parents:
                raise ValueError(f'Output replacement would affect input {source}')
        paths.append(p)
    # Complete all safety checks before removing anything.
    for p in paths:
        if p.is_dir():
            shutil.rmtree(p)
        elif p.exists():
            p.unlink()
    manifest.write_text(json.dumps({'artifacts': job['artifacts']}, indent=2))


def execute(job, action='run'):
    print('Checking configuration' + (' and input paths…' if action != 'schema' else '…'), flush=True)
    validate_job(job, check_paths=action != 'schema')
    if action != 'run':
        print('Configuration check passed.' if action == 'schema' else 'Configuration and input-path checks passed.')
        print('These checks do not validate assembled coordinates or simulate the system.')
        return 0
    root = Path(job['directory'])
    if root.resolve() != root:
        raise ValueError('The output directory changed or contains a symbolic link.')
    root.mkdir(parents=True, exist_ok=True)
    # Do not follow an existing lock symlink. The worker holds the lock until exit.
    lock_fd = os.open(root / '.gui.lock', os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    with os.fdopen(lock_fd, 'w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another run is already using this output folder.')
        cleanup(job)
        snapshot = root / 'run.yaml'
        snapshot.write_text(dump_yaml(job['config']), encoding='utf-8')
        status_path = root / 'run-status.json'
        status_path.write_text(json.dumps({'status': 'running', 'mode': job['mode'], 'id': job.get('id')}))
        previous_cwd = Path.cwd()
        code = 1
        status = 'failed'
        try:
            os.chdir(root)
            with (root / 'run.log').open('w', encoding='utf-8') as log:
                with redirect_stdout(Tee(sys.stdout, log)), redirect_stderr(Tee(sys.stderr, log)):
                    from run_pipeline import cli
                    args = (['prepare-minimization'] if job['mode'] == 'minimize'
                            else ['--mode', job['mode']])
                    try:
                        if job['mode'] == 'minimize':
                            print('Preparing portable minimization inputs…')
                        code = cli(args + ['-i', str(snapshot)])
                        status = 'completed' if code == 0 else 'failed'
                    except KeyboardInterrupt:
                        code, status = 130, 'cancelled'
                        print('\nRun cancelled. Partial products may remain; rerunning replaces them.')
                    except Exception:
                        traceback.print_exc()
        finally:
            os.chdir(previous_cwd)
            status_path.write_text(json.dumps({'status': status, 'mode': job['mode'], 'exit_code': code, 'id': job.get('id')}))
        return code


def main():
    # Each worker and its external tools get their own process group for cancellation.
    os.setsid()
    try:
        job = json.loads(Path(sys.argv[1]).read_text())
        return execute(job, sys.argv[2])
    except KeyboardInterrupt:
        print('Cancelled.', flush=True)
        return 130
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr, flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
