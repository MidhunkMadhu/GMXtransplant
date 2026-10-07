"""Regression coverage for independently parameterized target/replacement systems."""
import tempfile
import unittest
from pathlib import Path

from config import TopologySpec, load_config
from replacement_ligands import resolve_replacement_ligands
from topology import collect_topology_definitions, assemble_topology, TopologyError
from tests.test_topology_ligand_replacement import _write_itp, _final_universe


class SourceSelectionTests(unittest.TestCase):
    def test_source_selection_is_shared_by_charges_and_written_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            incoming = root / 'replacement' / 'toppar'
            target = root / 'environment' / 'toppar'
            incoming.mkdir(parents=True)
            target.mkdir(parents=True)
            for folder, charge in ((incoming, '0.0'), (target, '1.0')):
                _write_itp(folder / 'PROT.itp', 'PROT', [('CT', 'ALA', 'CA', charge)])
                _write_itp(folder / 'LI1.itp', 'LI1', [('N', 'LI1', 'N1', charge), ('C', 'LI1', 'C1', '0')])
                _write_itp(folder / 'SOD.itp', 'SOD', [('NA', 'SOD', 'SOD', charge)])
            (incoming / 'forcefield.itp').write_text('; parameters\n')
            (incoming.parent / 'topol.top').write_text('#include "toppar/forcefield.itp"\n[ molecules ]\nPROT 1\nLI1 1\n')
            cfg = TopologySpec(enabled=True, protein_toppar_dir=str(incoming),
                               environment_toppar_dir=str(target), output_dir=str(root / 'out'),
                               ligand_itp_paths=[str(incoming / 'LI1.itp')])
            definitions, files = collect_topology_definitions(cfg, cfg.ligand_itp_paths)
            self.assertEqual(definitions['PROT'].charge, 0)
            self.assertEqual(definitions['SOD'].charge, 1)
            self.assertEqual(files['LI1'], str(incoming / 'LI1.itp'))
            result = assemble_topology(_final_universe(), cfg, cfg.ligand_itp_paths)
            for name in ('PROT', 'LI1'):
                self.assertEqual((Path(result.toppar_dir) / (name + '.itp')).read_bytes(),
                                 (incoming / (name + '.itp')).read_bytes())
            # Ambiguity within one source remains an error.
            _write_itp(incoming / 'duplicate.itp', 'PROT', [('CT', 'ALA', 'CA', '2')])
            with self.assertRaisesRegex(TopologyError, 'conflicting definitions'):
                collect_topology_definitions(cfg, cfg.ligand_itp_paths)

    def test_minimal_metadata_and_automatic_environment_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ligand = root / 'ligand.itp'
            _write_itp(ligand, 'DistinctType', [('N', 'LIG', 'N1', '1')])
            config = root / 'job.yaml'
            config.write_text(f'''target_box:
  path: {root / 'target.pdb'}
  protein_mask: ':1'
replacement_structure:
  path: {root / 'incoming.gro'}
  protein_mask: ':1'
replacement_ligands:
  - resname: LIG
    itp_path: {ligand}
topology:
  environment_toppar_dir: {root / 'environment/toppar'}
name_restoration:
  enabled: true
  method: itp_atom_count
  pdb_to_full_resname:
    TIP: [TIP3]
''')
            cfg = load_config(str(config), check_paths=False)
            self.assertEqual(cfg.name_restoration.reference_topol,
                             str(root / 'environment/topol.top'))
            resolved = resolve_replacement_ligands(cfg)
            self.assertEqual(resolved[0].moleculetype, 'DistinctType')
            self.assertEqual(resolved[0].net_charge, 1)
