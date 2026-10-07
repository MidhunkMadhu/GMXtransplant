"""Translate complete worker log lines into concise, factual GUI progress."""
import re
import ast
from pathlib import Path


class RunProgress:
    def __init__(self, mode='protein'):
        self.mode = mode
        self.pending = ''
        self.stage = 'Ready to start'
        self.facts = {}
        self.notices = []
        self.error = ''

    def feed(self, text, final=False):
        self.pending += text
        lines = self.pending.split('\n')
        self.pending = lines.pop()
        if final and self.pending:
            lines.append(self.pending)
            self.pending = ''
        for line in lines:
            self.consume(line.strip())

    def consume(self, line):
        if line.startswith('Checking configuration'):
            self.stage = 'Checking your settings and input files'
        elif line == 'Periodic-box preflight':
            self.stage = 'Checking the simulation box and periodic geometry'
        elif line.startswith('Classification:'):
            if '(repaired)' in line:
                self.facts['Simulation box'] = 'Orientation repaired; geometry check passed.'
            elif '(passed)' in line:
                self.facts['Simulation box'] = 'Geometry check passed.'
        elif re.match(r'^\[\d+(?:-\d+)?\]', line):
            self.stage = self.stage_description(line)
        elif line.startswith('STEP 1: DETECT'):
            self.stage = 'Preparing cholesterol hydrogen atoms'
        elif line.startswith('STEP 2: MAP'):
            self.stage = 'Matching cholesterol to CHARMM36 parameters'
        elif line.startswith('STEP 3: INDEPENDENT'):
            self.stage = 'Verifying the converted cholesterol molecules'
        elif line == 'Preparing portable minimization inputs…':
            self.stage = 'Preparing the portable minimization inputs'
        elif line.startswith('Minimization inputs prepared:'):
            self.facts['Minimization'] = 'Portable inputs prepared; no simulation was run.'
        elif line.startswith('Configuration check passed.'):
            self.facts['Validation'] = 'Configuration is valid. Molecular input files were not opened.'
        elif line.startswith('Configuration and input-path checks passed.'):
            self.facts['Validation'] = 'Settings and input paths passed checks. Full assembly has not been validated.'

        charm_stages = {'[reference]': 'Reading and classifying the reference system',
                        '[transplant]': 'Reading and classifying the transplant system',
                        '[alignment]': 'Proteins aligned; preparing replacement',
                        '[clashes]': 'Checking overlaps with the retained environment',
                        '[charge]': 'Checking charge and adjusting counterions',
                        '[topology]': 'Building and validating topology and coordinates',
                        '[VISUALIZATION]': 'Preparing molecular comparison views',
                        '[host]': 'Reading and classifying the host system',
                        '[frame]': 'Locating the membrane and the protein tip',
                        '[binder]': 'Reading the binder and checking its parameters'}
        for prefix, description in charm_stages.items():
            if line.startswith(prefix):
                self.stage = description
        if line.startswith('[alignment]'):
            fit = re.search(r'RMSD ([\d.]+) A', line)
            if fit:
                self.facts['Alignment'] = f'Protein alignment RMSD: {fit[1]} Å.'
        rmsd = re.search(r'(?:RMSD.*?:|Alignment RMSD:)\s*([\d.]+)\s*A\b', line)
        if rmsd:
            self.facts['Alignment'] = f'Fit RMSD: {rmsd[1]} Å.'
        tip = re.match(r'\[frame\] .*tip (\S+) ([\d.]+) A beyond the (\w+) plane', line)
        if tip:
            self.facts['Protein tip'] = f'{tip[1]}, {tip[2]} Å beyond the {tip[3]} headgroup plane.'
        # "[binderpose1] flat at 20 A: accepted: host ..." (older runs: "[pose01] accepted ...")
        pose = re.match(r'\[([^\]]+)\] (?:.*: )?(accepted|rejected|failed)(?::|$)', line)
        if pose and pose[1] not in {'host', 'frame', 'binder', 'complete'}:
            self.stage = f'Building pose {pose[1]}'
            self.facts.setdefault('Poses', '')
            self.facts['Poses'] = (self.facts['Poses'] + '; ' if self.facts['Poses'] else '') + f'{pose[1]} {pose[2]}'
        if pose and pose[1] not in {'host', 'frame', 'binder', 'complete'}:
            self._last_pose = pose[1]
        salt = re.match(r'salt ([\d.]+) M \(host ([\d.]+) M\): (\d+) ion\(s\) removed, (\d+) added', line)
        if salt:
            entry = f'{getattr(self, "_last_pose", "pose")} {salt[3]} removed / {salt[4]} added'
            previous = self.facts.get('Salt', '').split(': ', 1)
            entries = previous[1].rstrip('.').split('; ') if len(previous) == 2 else []
            self.facts['Salt'] = f'{salt[1]} M (host {salt[2]} M), ions per pose: ' + '; '.join(entries + [entry]) + '.'
        built = re.match(r'\[complete\] (\d+)/(\d+) pose\(s\) built', line)
        if built:
            self.facts['Result'] = f'{built[1]} of {built[2]} poses built; see summary.txt.'
        repack = re.search(r'Wrote (.*?)/: repacking equilibration \(([\d.]+) ns, (\d+) restored', line)
        if repack:
            self.facts['Repacking'] = (f'{Path(repack[1]).name}/ holds a {repack[2]} ns repacking equilibration '
                                       f'with the {repack[3]} restored cholesterols held: topol.top, toppar/, '
                                       'GRO, index and mdp files; a sample run script is also saved.')
        restored = re.match(r'Restored cholesterols: (\d+), final residues (.+?) \(cholesterol (.+?) of (\d+) in topol\.top\)', line)
        if restored:
            # These three lines replace the generic clash/composition facts with simpler ones.
            for key in ('Clash removal', 'Composition', 'Membrane'):
                self.facts.pop(key, None)
            self.facts['Restored cholesterol'] = (f'{restored[1]} restored: final residues {restored[2]} '
                                                  f'(cholesterol {restored[3]} of {restored[4]}).')
        nearby = re.match(r'Removed next to them: (.+?)\.?$', line)
        if nearby:
            self.facts['Removed nearby'] = nearby[1][0].upper() + nearby[1][1:] + '.'
        far = re.match(r'Removed far away: (.+?)\.?$', line)
        if far:
            self.facts['Removed far away'] = far[1][0].upper() + far[1][1:] + '.'
        inserted = re.match(r'Inserted (\d+) cholesterol residue', line)
        if inserted:
            self.facts['Cholesterol'] = f'{inserted[1]} molecules inserted.'
        removed = re.match(r'Removed (\d+) counterion\(s\) -> final net charge ([+\-\d.]+) e', line)
        if removed:
            self.facts['Charge'] = f'{removed[1]} counterion(s) removed; final charge {removed[2]} e.'
        charge = re.search(r'Audited (?:final topology charge:|topology/coordinate.*?net charge)\s*([+\-\d.]+) e', line)
        if charge:
            self.facts['Charge'] = f'Final topology charge verified: {charge[1]} e.'
        composition = re.match(r'Composition removals: (\d+) lipid\(s\), (\d+) ion\(s\)', line)
        if composition:
            self.facts['Composition'] = f'{composition[1]} extra lipids and {composition[2]} ions removed.'
        verification = re.match(r'Summary: (\d+) PASS, (\d+) FAIL, (\d+) checked', line)
        if verification:
            self.facts['Cholesterol conversion'] = f'{verification[1]} passed verification; {verification[2]} failed.'
        classes = re.match(r'Removed by class: (\{.*?\})', line)
        if classes:
            try:
                counts = ast.literal_eval(classes[1])
                if isinstance(counts, dict) and all(isinstance(k, str) and type(v) is int for k, v in counts.items()):
                    self.facts['Clash removal'] = (', '.join(f'{v} {k}' for k, v in counts.items()) + ' molecule(s) removed.'
                                                   if counts else 'No clashing molecules needed removal.')
            except (ValueError, SyntaxError):
                pass
        lipids = re.search(r'Removed by class:.*?\(lipids (\d+)/(\d+)', line)
        if lipids:
            self.facts['Membrane'] = f'{lipids[1]} of {lipids[2]} lipids removed.'
        if line.startswith('WARNING:'):
            notice = line.removeprefix('WARNING:').strip()
            if notice not in self.notices:
                self.notices.append(notice)
        if line.startswith('ERROR:') or re.match(r'^\[[^\]]*(?:ERROR|STOPPED)[^\]]*\]', line):
            self.error = re.sub(r'^(?:ERROR:|\[[^\]]*\])\s*', '', line)
        # Config dumps, package metadata warnings, coordinate arrays, and atom
        # pair listings remain exclusively in the detailed command log.

    def stage_description(self, line):
        content = line.split(']', 1)[1].strip().lower()
        if self.mode == 'chl' and re.match(r'^\[[1-4]\]', line):
            # The CLI announces these four operations together before doing
            # them. Do not present that announcement as four completed stages.
            return 'Restoring cholesterol: preparing, aligning, and adjusting composition'
        if 'repacking' in content:
            return 'Writing the GROMACS repacking equilibration'
        if 'restoring' in content:
            return 'Restoring residue names using matching topology parameters'
        if 'reading' in content:
            return 'Reading your molecular structures'
        if 'aligning' in content:
            return 'Aligning the incoming protein to the target'
        if 'positioning the new ligand' in content:
            return 'Replacing and fitting the incoming ligand'
        if 'inserting' in content:
            return 'Inserting the aligned replacement structure'
        if 'steric clashes' in content:
            return 'Checking overlaps and removing clashing molecules'
        if 'net charge' in content:
            return 'Checking charge and adjusting counterions'
        if 'molecular composition' in content:
            return 'Building and validating the final molecular composition'
        if 'topol' in content:
            return 'Building and checking the GROMACS topology'
        if 'index' in content:
            return 'Writing atom selection groups'
        if 'writing' in content and ('coordinate' in content or 'pdb/gro' in content):
            return 'Writing coordinates and verifying the saved structures'
        if 'reference' in content:
            return 'Comparing the final system with the reference'
        return self.stage

    def finish(self, success, cancelled=False, action='run'):
        self.feed('', final=True)
        if cancelled:
            self.stage = 'Run cancelled. Partial outputs may remain.'
        elif success:
            self.stage = 'Your outputs are ready.' if action == 'run' else 'Checks finished.'
        else:
            self.stage = 'The run could not finish.' if action == 'run' else 'The checks did not pass.'
            if not self.error:
                self.error = 'Open Command Progress for the detailed error.'
