"""Typed form editor for GMXtransplant configurations.

One control per setting. A setting with a default is shown pre-filled and is
written to YAML only when it was in the loaded document or has been changed.
A setting without a default (auto-detected or inferred) has a single checkbox:
unticked means the pipeline decides, ticked reveals the value to enter.
Sections that the pipeline switches on and off have one checkbox, never a
second "enabled" field. Output locations follow the output folder selected at
the top of the window, so they are not shown at all.
"""
from copy import deepcopy
from dataclasses import is_dataclass
from pathlib import Path
import re
from typing import get_args, get_origin

from PySide6.QtCore import Qt, QLocale
from PySide6.QtGui import QDoubleValidator, QIntValidator
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QComboBox, QCheckBox, QScrollArea, QListWidget, QStackedWidget,
    QSplitter, QFileDialog, QFrame)

from .model import defaults, path_kind, is_output
from .fields import annotation_for, concrete, initial_value, residue_mapping

ENUMS = {
    ('charmprot', 'ligands'): ['auto', 'ignore'],
    ('charmprot', 'environment_overrides', '*'): ['lipid', 'sterol', 'water', 'ion', 'detergent', 'solvent'],
    ('addbinder', 'side'): ['upper', 'lower'],
    ('cholesterol', 'repack', 'stages', '*', 'integrator'): ['steep', 'nvt', 'npt'],
    ('addbinder', 'lateral_reference'): ['tip', 'host_center'],
    ('addbinder', 'ion_side'): ['any', 'binder_side'],
    ('addbinder', 'poses', '*', 'distance_to'): ['nearest_atom', 'centroid'],
    ('addbinder', 'poses', '*', 'orientation'): ['auto', 'flat', 'end_on', 'edge', 'as_is', 'euler'],
    ('addbinder', 'environment_overrides', '*'): ['lipid', 'sterol', 'water', 'ion', 'detergent', 'solvent'],
    ('alignment', 'method'): ['mask_fit', 'pymol_align', 'pymol_cealign'],
    ('ligand_replace', 'fit', 'method'): ['autofit', 'pairfit', 'mcsfit', 'nofit'],
    ('name_restoration', 'method'): ['itp_atom_count', 'atom_signature'],
    ('name_restoration', 'reference'): ['replacement_structure', 'target_box'],
    ('name_restoration', 'apply_to'): ['target_box', 'replacement_structure', 'both'],
    ('ndx', 'style'): ['both', 'generic', 'charmm_gui'],
    ('box_validation', 'mode'): ['strict', 'repair', 'off'],
    ('minimization', 'platform'): ['auto', 'CPU', 'CUDA', 'HIP', 'OpenCL', 'Reference'],
    ('minimization', 'restraint_residue_classes', '*'): ['lipid', 'ligand', 'protein', 'water', 'ion'],
}
TITLES = {
    'charmprot': 'Reference and transplant', 'paths': 'Input files', 'target_box': 'Target environment',
    'replacement_structure': 'Incoming structure', 'ndx': 'Index groups', 'refgro': 'Reference GRO',
    'name_restoration': 'Residue name restoration', 'pdb_to_full_resname': 'PDB residue name replacements',
    'ligand_replace': 'Ligand replacement', 'replacement_ligands': 'Incoming ligands',
    'reference': 'Reference system folder', 'transplant': 'Transplant system folder',
    'reference_proteins': 'Reference protein molecules', 'transplant_proteins': 'Transplant protein molecules',
    'reference_ligands': 'Reference ligand molecules', 'transplant_ligands': 'Transplant ligand molecules',
    'ligands': 'Ligand handling', 'reference_gro': 'Reference coordinates (GRO)',
    'transplant_gro': 'Transplant coordinates (GRO)', 'clash_distance': 'Clash distance (Å)',
    'keep_lipids': 'Keep lipids', 'target_net_charge': 'Target net charge (e)',
    'ion_exclusion_distance': 'Ion exclusion distance from protein (Å)',
    'lipid_ion_exclusion_distance': 'Ion exclusion distance from lipids (Å)',
    'environment_overrides': 'Environment overrides',
    'restraint_force_constant_kj_mol_nm2': 'Restraint strength (kJ/mol/nm²)',
    'max_iterations': 'Maximum iterations', 'platform': 'OpenMM platform',
    'restraint_residue_classes': 'Custom residue classes', 'coordinates_path': 'Assembled coordinates (GRO)',
    'topology_path': 'Assembled topology (topol.top)', 'include_dir': 'Topology include folder',
    'box_dimensions': 'Box dimensions', 'protein_mask': 'Protein mask', 'itp_path': 'ITP file',
    'resname': 'Residue name', 'moleculetype': 'Molecule type name',
    'upper': 'Upper leaflet', 'lower': 'Lower leaflet',
    'forcefield_path': 'Force field (forcefield.itp)', 'coord_path': 'Coordinates',
    'addbinder': 'Host and binder',
}
# Label of the one checkbox that switches a section on.
TOGGLES = {
    'topology': 'Write the GROMACS topology (topol.top and toppar/)',
    'name_restoration': 'Restore truncated residue names',
    'minimization': 'Prepare a restrained minimization folder',
}
HELP = {
    'charmprot': 'Choose two prepared CHARMM-GUI system folders, each containing topol.top, toppar/ and '
                 'coordinates. All reference proteins and ligands are replaced from the transplant; the '
                 'membrane, water and ions of the reference are kept.',
    'addbinder': 'Choose the prepared membrane system (host folder) and the binder: its coordinates, ITP and, '
                 'for a CHARMM-GUI ligand, the forcefield.itp generated with it. Each pose places the binder '
                 'the given distance beyond the protein\'s outermost heavy atom on the chosen side, checks it '
                 'against the protein, the membrane and every periodic image, and writes one system per pose. '
                 'Distances are in Å.',
    'paths': 'Give each input as a complete path, or browse for it. Other sections refer to these files.',
    'alignment': 'Fit matching protein regions. Masks use sequential residue positions, not printed residue IDs. Distances are in Å.',
    'output': 'Filenames of the generated files.',
    'replacement_ligands': 'Bound ligands moving with the incoming protein, each with its matching ITP file.',
    'minimization': 'Prepares a portable OpenMM folder; it is not run here. Protein and ligand heavy atoms '
                    'and lipid head groups are restrained, so only lipid tails, water and ions relax.',
    'charge': 'Charge is measured in elementary charges (e); exclusion distances are in Å.',
    'cholesterol': 'Experimental cholesterol is fitted into the prepared target. CHARMM36 conversion requires Open Babel.',
    'topology': 'Select matching parameter folders. A protein template topology can override the sibling topol.top.',
    'box_validation': 'Validate periodic geometry and hard clashes; distances are in Å.',
    'ndx': 'An index.ndx with atom groups for GROMACS is always written with the system.',
}
PATH_TITLES = {
    ('name_restoration', 'reference'): 'Name source',
    ('name_restoration', 'min_jaccard'): 'Minimum name similarity',
    ('name_restoration', 'reference_topol'): 'Reference topology (topol.top)',
    ('output', 'pdb_path'): 'Final PDB filename',
    ('output', 'gro_path'): 'Final GRO filename',
    ('output', 'report_path'): 'Report filename prefix',
    ('output', 'inspection_pdb_path'): 'Alignment inspection PDB filename',
    ('ndx', 'output_path'): 'Index filename',
    ('ndx', 'style'): 'Group naming style',
    ('ndx', 'extra_groups'): 'Extra groups',
    ('topology', 'protein_toppar_dir'): 'Incoming protein toppar folder',
    ('topology', 'environment_toppar_dir'): 'Environment toppar folder',
    ('topology', 'protein_template_top'): 'Incoming topol.top (override)',
    ('topology', 'ligand_itp_paths'): 'Extra ligand ITP files',
    ('topology', 'moleculetype_overrides'): 'Molecule type overrides',
    ('charge', 'tolerance'): 'Charge tolerance (e)',
    ('charge', 'exclusion_distance_from_protein'): 'Ion exclusion distance from protein (Å)',
    ('charge', 'exclusion_distance_from_ligand'): 'Ion exclusion distance from ligands (Å)',
    ('charge', 'exclusion_distance_from_lipid'): 'Ion exclusion distance from lipids (Å)',
    ('charge', 'exclude_membrane_interior'): 'Keep ions near lipids',
    ('charge', 'charge_table_overrides'): 'Residue charge overrides',
    ('clash_detection', 'threshold'): 'Clash distance (Å)',
    ('clash_detection', 'thresholds'): 'Clash distance per class (Å)',
    ('clash_detection', 'use_pbc'): 'Use periodic boundaries',
    ('clash_detection', 'heavy_atoms_only'): 'Heavy atoms only',
    ('clash_detection', 'flag_lipid_removal_fraction'): 'Lipid removal warning fraction',
    ('clash_detection', 'keep_classes'): 'Never remove these classes',
    ('box_validation', 'mode'): 'Box check mode',
    ('box_validation', 'hard_clash_cutoff'): 'Hard clash distance (Å)',
    ('box_validation', 'molecule_protrusion_allowance'): 'Molecule protrusion allowance (Å)',
    ('box_validation', 'wrap_after_validation'): 'Wrap molecules into the box',
    ('box_validation', 'search_sample_size'): 'Repair search: sampled atoms',
    ('box_validation', 'search_starts'): 'Repair search: starting orientations',
    ('box_validation', 'search_max_iterations'): 'Repair search: iterations per start',
    ('box_validation', 'max_candidate_verifications'): 'Repair search: candidates verified',
    ('box_validation', 'max_search_seconds'): 'Repair search: time limit (s)',
    ('addbinder', 'host'): 'Host system folder',
    ('addbinder', 'host_gro'): 'Host coordinates (GRO)',
    ('addbinder', 'binder_coordinates'): 'Binder coordinates (PDB, GRO or MOL2)',
    ('addbinder', 'binder_itp'): 'Binder ITP',
    ('addbinder', 'binder_forcefield'): 'Binder force field (forcefield.itp)',
    ('addbinder', 'side'): 'Membrane side',
    ('addbinder', 'distance'): 'Distance from the protein tip (Å)',
    ('addbinder', 'min_image_gap'): 'Minimum gap to periodic images (Å)',
    ('addbinder', 'min_membrane_gap'): 'Minimum gap to lipids (Å)',
    ('addbinder', 'lateral_reference'): 'Lateral position measured from',
    ('addbinder', 'water_clash_distance'): 'Water/ion removal distance (Å)',
    ('addbinder', 'ion_side'): 'Remove or add ions in',
    ('addbinder', 'salt_concentration'): 'Salt concentration (host, or mol/L)',
    ('addbinder', 'poses'): 'Poses',
    ('addbinder', 'poses', '*', 'name'): 'Pose name (output folder)',
    ('addbinder', 'poses', '*', 'distance'): 'Own distance from the tip (Å)',
    ('addbinder', 'poses', '*', 'lateral_offset'): 'Lateral offset x y (Å)',
    ('addbinder', 'poses', '*', 'angles'): 'Euler angles z y z (degrees)',
    ('addbinder', 'poses', '*', 'spin'): 'Spin about the membrane normal (degrees)',
    ('addbinder', 'poses', '*', 'flip'): 'Flip upside down',
    ('addbinder', 'poses', '*', 'approach'): 'Approach direction: tilt, azimuth (degrees)',
    ('addbinder', 'poses', '*', 'centroid'): 'Binder centre x y z (Å, host GRO frame)',
    ('addbinder', 'poses', '*', 'from_input_position'): 'Start from the input position',
    ('addbinder', 'poses', '*', 'distance_to'): 'Distance measured to',
    ('addbinder', 'poses', '*', 'reduce_distance_by'): 'If it does not fit, reduce the distance by (Å)',
    ('addbinder', 'poses', '*', 'min_distance'): 'Shortest distance to try (Å)',
    ('addbinder', 'random_rotated'): 'Random poses that are also rotated',
    ('addbinder', 'random_from_input_position'): 'Random poses: start from the input position',
    ('addbinder', 'random_poses'): 'Random poses',
    ('addbinder', 'random_distance'): 'Random distance range min max (Å)',
    ('addbinder', 'random_max_tilt'): 'Random poses: largest tilt from the normal (degrees)',
    ('addbinder', 'random_min_angle'): 'Random poses: smallest angle between directions (degrees)',
    ('addbinder', 'random_min_separation'): 'Random poses: smallest distance between binder centres (Å)',
    ('addbinder', 'random_seed'): 'Random seed',
    ('addbinder', 'same_orientation'): 'Same orientation for every pose',
    ('addbinder', 'copy_run_inputs'): 'Copy the host\'s mdp files and README into each pose',
    ('alignment', 'method'): 'Alignment method',
    ('alignment', 'region_mask_original'): 'Fit region in the target',
    ('alignment', 'region_mask_replacement'): 'Fit region in the incoming structure',
    ('alignment', 'cycles'): 'PyMOL outlier-rejection cycles',
    ('alignment', 'cutoff'): 'PyMOL outlier cutoff (Å)',
    ('ligand_replace', 'structure_path'): 'Prepared system (with the old ligand)',
    ('ligand_replace', 'fit', 'method'): 'Placement method',
    ('ligand_replace', 'fit', 'old_ligand_fit_atoms'): 'Fit atoms in the old ligand',
    ('ligand_replace', 'fit', 'new_ligand_fit_atoms'): 'Fit atoms in the new ligand',
    ('cholesterol', 'experimental_structure_path'): 'Experimental structure (with cholesterol)',
    ('cholesterol', 'target_system_path'): 'Prepared target system',
    ('cholesterol', 'target_protein_mask'): 'Protein mask in the target',
    ('cholesterol', 'experimental_protein_mask'): 'Protein mask in the experimental structure',
    ('cholesterol', 'cholesterol_resnames'): 'Cholesterol residue names to look for',
    ('cholesterol', 'charmm_resname'): 'CHARMM cholesterol residue name',
    ('cholesterol', 'convert_to_charmm36'): 'Convert to CHARMM36',
    ('cholesterol', 'max_heavy_atom_fit_rmsd'): 'Maximum reconstruction fit RMSD (Å)',
    ('cholesterol', 'charmm36_reference_path'): 'CHARMM36 cholesterol reference (optional)',
    ('cholesterol', 'obabel_command'): 'Open Babel command',
    ('cholesterol', 'verify_tolerance'): 'Verification tolerance (Å)',
    ('cholesterol', 'target_system_format'): 'Format',
    ('cholesterol', 'composition', 'reference_system_path'): 'Composition reference system',
    ('cholesterol', 'composition', 'distance_from_protein'): 'Protected distance from protein (Å)',
    ('cholesterol', 'composition', 'target_concentration'): 'Salt concentration (mol/L)',
    ('cholesterol', 'composition', 'concentration_tolerance_fraction'): 'Concentration tolerance (fraction)',
    ('cholesterol', 'composition', 'lipid_targets'): 'Lipids per leaflet',
    ('cholesterol', 'repack'): 'Repacking equilibration (GROMACS)',
    ('cholesterol', 'repack', 'enabled'): 'Write the repacking equilibration',
    ('cholesterol', 'repack', 'output_dir'): 'Folder name',
    ('cholesterol', 'repack', 'moleculetype'): 'Name for the restrained cholesterol type',
    ('cholesterol', 'repack', 'lipid_restraint'): 'Other lipids: restraint (kJ/mol/nm²)',
    ('cholesterol', 'repack', 'temperature'): 'Temperature (K)',
    ('cholesterol', 'repack', 'stages'): 'Stages',
    ('cholesterol', 'repack', 'stages', '*', 'name'): 'Stage (mdp file name)',
    ('cholesterol', 'repack', 'stages', '*', 'integrator'): 'Type',
    ('cholesterol', 'repack', 'stages', '*', 'nsteps'): 'Steps',
    ('cholesterol', 'repack', 'stages', '*', 'dt'): 'Time step (ps)',
    ('cholesterol', 'repack', 'stages', '*', 'restraint'): 'Protein, ligands, restored cholesterols (kJ/mol/nm²)',
    ('cholesterol', 'raw_protonated_pdb_path'): 'Diagnostic: protonated cholesterol PDB',
    ('cholesterol', 'converted_pdb_path'): 'Diagnostic: converted cholesterol PDB',
    ('cholesterol', 'aligned_reference_pdb_path'): 'Diagnostic: aligned reference PDB',
    ('cholesterol', 'reference_placed_pdb_path'): 'Diagnostic: placed reference PDB',
    ('cholesterol', 'merged_pdb_path'): 'Diagnostic: merged system PDB',
}
# One line under each setting: what it does and what you can do with it.
DESCRIPTIONS = {
    ('addbinder', 'host'): 'The prepared membrane system the binder joins: a folder with topol.top, toppar/ and step5_input.gro, e.g. a CHARMM-GUI build. Only its water and ions change.',
    ('addbinder', 'host_gro'): "Tick to use a coordinate file other than the host folder's own step5 GRO, e.g. an MD frame.",
    ('addbinder', 'binder_coordinates'): 'One binder molecule with hydrogens. Atom names must match the ITP. Its position does not matter.',
    ('addbinder', 'binder_itp'): 'The ITP with the binder\'s one molecule type. A name already used by a different host molecule (e.g. PROA) is renamed automatically. Binders of several molecules (one ITP each) are set in the YAML as a list.',
    ('addbinder', 'binder_forcefield'): 'The force field holding the binder\'s parameters, e.g. the forcefield.itp CHARMM-GUI generated with the ligand. Untick only when the host force field already covers the binder.',
    ('addbinder', 'side'): 'upper = +z side of the membrane, lower = −z side. Which side is extracellular depends on how your system was built.',
    ('addbinder', 'distance'): 'Gap along the membrane normal from the protein\'s outermost heavy atom to the binder\'s nearest heavy atom. Used by every pose without its own distance.',
    ('addbinder', 'min_image_gap'): 'The binder must stay at least this far from every periodic image of the protein and of itself. A pose that fails is rejected with the box height it would need.',
    ('addbinder', 'min_membrane_gap'): 'The binder must stay at least this far from lipid heavy atoms.',
    ('addbinder', 'lateral_reference'): 'tip: centre the binder over the tip atom. host_center: over the protein\'s centre in x and y.',
    ('addbinder', 'water_clash_distance'): 'Water and ions with a heavy atom closer than this to any binder atom are removed. Protein and lipids are never removed.',
    ('addbinder', 'neutralize'): 'Add or remove counterions so the system reaches the target charge. Untick to keep the current ion difference.',
    ('addbinder', 'salt_concentration'): 'The binder displaces water and ions. host keeps the host\'s salt concentration (pairs × 55.51 / waters) on the remaining water; or give a value in mol/L. Extra ions are removed and missing ones replace bulk waters.',
    ('addbinder', 'same_orientation'): 'Poses with orientation auto each get a different orientation, whatever their distance. Tick to give them all the same (flat) one.',
    ('addbinder', 'copy_run_inputs'): 'Copies the host folder\'s *.mdp files and CHARMM-GUI README run script, so each pose folder is ready for grompp. Its index.ndx has the CHARMM-GUI groups SOLU, MEMB, SOLV and SOLU_MEMB.',
    ('addbinder', 'poses', '*', 'flip'): 'Turn the binder upside down (180° about x). Needs an orientation other than auto.',
    ('addbinder', 'poses', '*', 'distance_to'): 'nearest_atom: the gap between the nearest heavy atoms. centroid: from the tip (or host centre) to the binder\'s centre; with an approach direction this places the binder in spherical coordinates (distance, tilt, azimuth).',
    ('addbinder', 'poses', '*', 'reduce_distance_by'): 'When the pose does not fit the box (gap to periodic images, membrane), try again this much closer, down to the shortest distance. 0 tries only the one distance.',
    ('addbinder', 'random_rotated'): 'Untick for all random poses to get a random rotation. Tick and give a number to rotate only that many (the last ones); the others keep the input orientation, e.g. a G protein in its bound arrangement.',
    ('addbinder', 'random_from_input_position'): 'Tick to move each random pose out from the input position, so only the side and tilt are random.',
    ('addbinder', 'poses', '*', 'from_input_position'): 'For a binder whose coordinates are already in the host\'s frame, e.g. a G protein in its receptor-bound arrangement: start there and move it straight out along the membrane normal until it is the distance from the protein.',
    ('addbinder', 'poses', '*', 'approach'): 'Tick to come in tilted from the membrane normal by tilt, towards azimuth (0 = +x, 90 = +y), instead of straight down. The distance is then to the nearest protein heavy atom.',
    ('addbinder', 'poses', '*', 'centroid'): 'Tick to put the binder\'s heavy-atom centre exactly at x y z, in Å in the host GRO\'s coordinates (GRO nm × 10). Distance, approach and lateral offset are then not used; the pose is rejected only if it overlaps the host (under 3 Å) or fails the box checks.',
    ('addbinder', 'random_poses'): 'Extra poses, numbered on after the listed ones (binderpose2, binderpose3 ...), each from a random direction, distance and orientation. A draw that fails the box checks, or comes too close in angle or position to another pose, is drawn again. 0 for none.',
    ('addbinder', 'random_distance'): 'Each random pose takes a distance in this range from the nearest protein heavy atom. Untick for distance to distance + 5 Å.',
    ('addbinder', 'random_seed'): 'The same seed gives the same random poses. Each pose report records the values drawn, so you can copy one into the poses list.',
    ('addbinder', 'ion_side'): 'any: bulk water anywhere in the box (the two water layers are connected through the periodic boundary, so ions redistribute during equilibration). binder_side: the binder\'s side of the membrane first.',
    ('addbinder', 'poses'): 'Each pose is built and checked on its own and written to its own folder, ready to run; a rejected pose keeps only its report. No poses: one flat pose.',
    ('addbinder', 'poses', '*', 'orientation'): 'auto: a different orientation for each pose (flat, end-on, end-on flipped, edge ...), or the same one when "Same orientation for every pose" is ticked. flat: thinnest axis along the membrane normal (uses the least box height). end_on: longest axis along the normal, its far-reaching end towards the protein. edge: middle axis along the normal. as_is: keep the input orientation. euler: rotate by the angles below.',
    ('addbinder', 'poses', '*', 'distance'): 'Tick to override the distance set above for this pose.',
    ('addbinder', 'poses', '*', 'angles'): 'Used only with orientation euler.',
    ('charmprot', 'reference'): 'The prepared system whose membrane, water and ions are kept. Choose its CHARMM-GUI folder (topol.top, toppar/, coordinates).',
    ('charmprot', 'transplant'): 'The system the new protein and ligands come from. Choose its CHARMM-GUI folder.',
    ('charmprot', 'reference_proteins'): 'Normally detected automatically. Tick to name the protein molecule types to replace, e.g. PROA, PROB.',
    ('charmprot', 'transplant_proteins'): 'Normally detected automatically. Tick to name the incoming protein molecule types, e.g. PROA.',
    ('charmprot', 'reference_ligands'): 'Normally detected automatically. Tick to name exactly which reference ligands are removed.',
    ('charmprot', 'transplant_ligands'): 'Normally detected automatically. Tick to name exactly which transplant ligands are inserted.',
    ('charmprot', 'reference_gro'): "Tick to use a coordinate file other than the reference folder's own step5 GRO.",
    ('charmprot', 'transplant_gro'): "Tick to use a coordinate file other than the transplant folder's own step5 GRO.",
    ('charmprot', 'clash_distance'): 'Environment molecules with a heavy atom closer than this to the inserted protein are removed. Lower keeps more.',
    ('charmprot', 'neutralize'): 'Remove counterions to bring the system to the target charge. Untick to only report the charge.',
    ('charmprot', 'target_net_charge'): 'The net charge the system is brought to; 0 for a neutral system.',
    ('charmprot', 'ion_exclusion_distance'): 'Ions closer than this to the protein or ligands are never removed during neutralization.',
    ('charmprot', 'lipid_ion_exclusion_distance'): 'Ions closer than this to lipids are never removed during neutralization.',
    ('minimization', 'restraint_force_constant_kj_mol_nm2'): 'How strongly restrained atoms are held in place. Higher keeps them closer; 1000 is a common choice.',
    ('minimization', 'max_iterations'): 'Upper limit on minimizer steps. Increase if the report says it did not converge.',
    ('minimization', 'platform'): 'Hardware used when the folder is run later. auto picks the fastest available: CUDA, HIP, OpenCL, then CPU.',
    ('minimization', 'coordinates_path'): 'The assembled system to minimize, e.g. step5_input.gro from a previous run.',
    ('minimization', 'topology_path'): 'The topol.top that matches those coordinates.',
    ('target_box', 'path'): 'The prepared system the new protein goes into (GRO recommended: it keeps full names and the box).',
    ('target_box', 'protein_mask'): 'The old protein (and bound ligands) to remove, e.g. :1-963. Positions count residues in file order.',
    ('target_box', 'format'): 'auto detects from the file extension; choose gro or pdb to force it.',
    ('target_box', 'box_dimensions'): 'Only needed if the file has no box. Tick and give a b c α β γ in Å and degrees.',
    ('replacement_structure', 'path'): 'The structure containing the new protein and any ligands bound to it.',
    ('replacement_structure', 'protein_mask'): 'The block to insert, including bound ligands, e.g. :1-963.',
    ('replacement_structure', 'format'): 'auto detects from the file extension; choose gro or pdb to force it.',
    ('replacement_structure', 'box_dimensions'): 'Only needed if the file has no box. Tick and give a b c α β γ in Å and degrees.',
    ('replacement_ligands', '*', 'resname'): 'Residue name of this ligand in the incoming structure. Use NAME:N to pick one of several copies.',
    ('replacement_ligands', '*', 'itp_path'): 'The ligand ITP (atoms, charges, bonds). Browse to the file from CHARMM-GUI.',
    ('replacement_ligands', '*', 'charge'): 'Read from ITP sums the ITP charges; choose Enter charge only to override it.',
    ('replacement_ligands', '*', 'moleculetype'): 'Leave empty: found from the residue name. Only needed if the ITP defines several molecule types with that residue name.',
    ('refgro',): 'Optional check: a GRO whose non-protein molecule types should all appear in the output. Missing ones are noted in the report.',
    ('alignment', 'method'): 'mask_fit fits the two regions below atom for atom (no PyMOL needed). pymol_align / pymol_cealign need PyMOL.',
    ('alignment', 'region_mask_original'): 'Atoms of the old protein used for the fit, e.g. :1-961@CA. Leave unticked to use the whole protein mask.',
    ('alignment', 'region_mask_replacement'): 'Matching atoms of the new protein, same count and order as the target region.',
    ('alignment', 'cycles'): 'PyMOL methods only: rounds of outlier rejection.',
    ('alignment', 'cutoff'): 'PyMOL methods only: atoms further apart than this after fitting are rejected as outliers.',
    ('name_restoration', 'method'): 'atom_signature compares atom names with a full-name structure; itp_atom_count uses atom counts from the ITPs.',
    ('name_restoration', 'reference'): 'Which structure supplies the full residue names (atom_signature).',
    ('name_restoration', 'apply_to'): 'Which structure gets its truncated names restored.',
    ('name_restoration', 'min_jaccard'): 'How similar atom-name sets must be to accept a name (0-1). Lower accepts looser matches.',
    ('name_restoration', 'reference_topol'): 'itp_atom_count only: the topol.top listing the full molecule names.',
    ('name_restoration', 'pdb_to_full_resname'): 'Truncated PDB name → allowed full names, e.g. POP → POPC, POPE. Add one row per truncated name.',
    ('box_validation', 'mode'): 'strict only checks the box; repair may rotate a wrongly oriented frame; off skips the check.',
    ('box_validation', 'hard_clash_cutoff'): 'Atom pairs closer than this in the final system are reported as hard clashes.',
    ('box_validation', 'molecule_protrusion_allowance'): 'Report molecules sticking out of the box by more than this; it never rejects the system.',
    ('box_validation', 'wrap_after_validation'): 'Move whole molecules back inside the box after checking. Off keeps coordinates as they are.',
    ('box_validation', 'random_seed'): 'Seed for the repair search, so repeated runs give the same result.',
    ('box_validation', 'search_sample_size'): 'Repair mode: atoms sampled when testing orientations. Rarely needs changing.',
    ('box_validation', 'search_starts'): 'Repair mode: number of starting orientations tried. Rarely needs changing.',
    ('box_validation', 'search_max_iterations'): 'Repair mode: refinement steps per starting orientation. Rarely needs changing.',
    ('box_validation', 'max_candidate_verifications'): 'Repair mode: best orientations checked in full. Rarely needs changing.',
    ('box_validation', 'max_search_seconds'): 'Repair mode: time limit for the orientation search.',
    ('clash_detection', 'threshold'): 'Molecules closer than this to the inserted block are removed (whole molecules). Lower keeps more.',
    ('clash_detection', 'thresholds'): 'Tick to use a different clash distance for a class, e.g. water → 0.6.',
    ('clash_detection', 'use_pbc'): 'Measure distances across periodic boundaries. Only turn off together with box check mode off.',
    ('clash_detection', 'heavy_atoms_only'): 'Check only non-hydrogen atoms. Untick to include hydrogens (stricter).',
    ('clash_detection', 'flag_lipid_removal_fraction'): 'Warn when more than this fraction of lipids is removed; it is not a limit.',
    ('clash_detection', 'keep_classes'): 'Classes that are reported but never removed, e.g. lipid to keep the membrane intact.',
    ('charge', 'target_net_charge'): 'The net charge the system is brought to; 0 for a neutral system.',
    ('charge', 'neutralize'): 'Remove counterions to reach the target charge. Untick to only report the charge.',
    ('charge', 'tolerance'): 'How far from the target the final charge may be.',
    ('charge', 'exclusion_distance_from_protein'): 'Ions closer than this to the protein are never removed.',
    ('charge', 'exclusion_distance_from_ligand'): 'Ions closer than this to ligands are never removed.',
    ('charge', 'exclusion_distance_from_lipid'): 'Ions closer than this to lipids are never removed (with the option below).',
    ('charge', 'exclude_membrane_interior'): 'Protect ions within the lipid distance above from removal. Untick to allow removing them.',
    ('charge', 'random_seed'): 'Seed for choosing which eligible ions are removed, so repeated runs give the same result.',
    ('charge', 'charge_table_overrides'): 'Net charge of residues the pipeline does not know, e.g. MYLIG → -1. Add one row per residue.',
    ('topology', 'protein_toppar_dir'): "The incoming protein's CHARMM-GUI toppar/ folder; its topol.top must sit next to it.",
    ('topology', 'environment_toppar_dir'): 'A toppar/ folder with the lipid, water and ion ITPs of the kept environment.',
    ('topology', 'protein_template_top'): 'Only if the incoming topol.top is not next to its toppar/ folder: choose it here.',
    ('topology', 'ligand_itp_paths'): 'Tick to add ITPs of cofactors that are not already in the toppar folders.',
    ('topology', 'moleculetype_overrides'): 'For residues whose molecule type has a different name, e.g. HEM → HEME. Add one row per residue.',
    ('ndx', 'output_path'): 'Filename of the index file written with the system.',
    ('ndx', 'style'): 'generic (Protein, Lipid, ...), charmm_gui (SOLU, MEMB, SOLV) or both.',
    ('ndx', 'extra_groups'): 'Your own groups: a name and a selection mask, e.g. Pocket → :120-130.',
    ('output', 'pdb_path'): 'Filename of the final PDB, written in the output folder.',
    ('output', 'gro_path'): 'Filename of the final GRO, written in the output folder.',
    ('output', 'report_path'): 'Reports are written as this name with .txt and .json.',
    ('output', 'inspection_pdb_path'): 'PDB overlaying the old and new placement, for checking the fit.',
    ('ligand_replace', 'structure_path'): 'The complete prepared system; everything is kept except the old ligand.',
    ('ligand_replace', 'format'): 'auto detects from the file extension; choose gro or pdb to force it.',
    ('ligand_replace', 'box_dimensions'): 'Only needed if the file has no box. Tick and give a b c α β γ in Å and degrees.',
    ('ligand_replace', 'original_ligand', 'resname'): 'Residue name of the ligand to remove. Use NAME:N to pick one of several copies.',
    ('ligand_replace', 'new_ligand', 'coord_path'): 'Coordinates of the incoming ligand (MOL2, PDB, ...); atom order must match its ITP.',
    ('ligand_replace', 'new_ligand', 'resname'): 'Residue name the new ligand gets in the output.',
    ('ligand_replace', 'new_ligand', 'itp_path'): 'The new ligand ITP (atoms, charges, bonds). It is added to the topology automatically.',
    ('ligand_replace', 'new_ligand', 'format'): 'auto detects from the file extension; choose one to force it.',
    ('ligand_replace', 'fit', 'old_ligand_fit_atoms'): 'pairfit only: atom names in the old ligand, e.g. N4, C10, C9. Ignored by other methods.',
    ('ligand_replace', 'fit', 'new_ligand_fit_atoms'): 'pairfit only: the matching atom names in the new ligand, in the same order.',
    ('cholesterol', 'experimental_structure_path'): 'The experimental structure (e.g. a PDB entry) containing the cholesterol to restore.',
    ('cholesterol', 'target_system_path'): 'The prepared membrane system the cholesterol is added to.',
    ('cholesterol', 'target_system_format'): 'auto detects from the file extension; choose gro or pdb to force it.',
    ('cholesterol', 'box_dimensions'): 'Only needed if the file has no box. Tick and give a b c α β γ in Å and degrees.',
    ('cholesterol', 'target_protein_mask'): 'Protein atoms used to align the two structures, in the target, e.g. :1-328.',
    ('cholesterol', 'experimental_protein_mask'): 'The corresponding protein atoms in the experimental structure.',
    ('cholesterol', 'cholesterol_resnames'): 'Residue names that count as cholesterol in the experimental structure.',
    ('cholesterol', 'charmm_resname'): 'Residue name the restored cholesterol gets (CHL1 in CHARMM36).',
    ('cholesterol', 'convert_to_charmm36'): 'Add hydrogens and CHARMM36 names with Open Babel. Untick if the input is already CHARMM36 CHL1.',
    ('cholesterol', 'reconstruct_missing_heavy_atoms'): 'Rebuild a few missing heavy atoms from the CHARMM36 template. Untick to reject incomplete molecules.',
    ('cholesterol', 'max_missing_heavy_atoms'): 'The most heavy atoms that may be rebuilt per cholesterol.',
    ('cholesterol', 'max_heavy_atom_fit_rmsd'): 'Rebuilding is refused if the template fits worse than this.',
    ('cholesterol', 'write_diagnostics'): 'Also keep the intermediate PDB files listed below, for inspection.',
    ('cholesterol', 'charmm36_reference_path'): 'Leave empty to use the bundled CHARMM36 cholesterol; choose a file to use your own.',
    ('cholesterol', 'obabel_command'): 'The Open Babel program. Leave obabel to use the one installed with GMXtransplant, or browse to another.',
    ('cholesterol', 'verify_tolerance'): 'Converted atoms must stay within this distance of the experimental positions.',
    ('cholesterol', 'composition', 'reference_system_path'): 'Tick to count lipids and salt from another system; otherwise the target is used.',
    ('cholesterol', 'composition', 'distance_from_protein'): 'Lipids closer than this to the protein are never removed when adjusting composition.',
    ('cholesterol', 'composition', 'target_concentration'): 'Salt concentration of the reference system, kept in the output.',
    ('cholesterol', 'composition', 'concentration_tolerance_fraction'): 'Allowed relative deviation from that concentration, e.g. 0.05 for 5%.',
    ('cholesterol', 'composition', 'salt_cation'): 'Residue name of the salt cation, e.g. SOD.',
    ('cholesterol', 'composition', 'salt_anion'): 'Residue name of the salt anion, e.g. CLA.',
    ('cholesterol', 'composition', 'random_seed'): 'Seed for choosing which lipids and ions are removed, so repeated runs give the same result.',
    ('cholesterol', 'composition', 'lipid_targets'): 'Wanted number of each lipid in the upper and lower leaflet. Extra lipids far from the protein are removed.',
    ('cholesterol', 'repack'): 'Why: the lipids that overlapped the restored cholesterols were removed, so the cholesterols start with few neighbours and drift off their poses in free MD. What is written: repack/ with its own topol.top, toppar/, GRO, index and mdp files for a short equilibration (about 1.4 ns) that holds every heavy atom of the protein, ligands and restored cholesterols while the other lipids fill the gaps. A sample run script is also saved. Afterwards continue from its last GRO with the main topol.top.',
    ('cholesterol', 'repack', 'lipid_restraint'): 'Sets POSRES_FC_LIPID and DIHRES_FC for every other lipid. 0 leaves them free to repack.',
    ('cholesterol', 'repack', 'moleculetype'): 'The restored cholesterols get a copy of the cholesterol ITP under this name, with every heavy atom restrained in x, y and z. The main topol.top is unchanged.',
    ('cholesterol', 'repack', 'stages', '*', 'integrator'): 'steep = minimization, nvt = constant volume, npt = constant pressure (semi-isotropic).',
    ('cholesterol', 'raw_protonated_pdb_path'): 'Written only when diagnostics are on.',
    ('cholesterol', 'converted_pdb_path'): 'Written only when diagnostics are on.',
    ('cholesterol', 'aligned_reference_pdb_path'): 'Written only when diagnostics are on.',
    ('cholesterol', 'reference_placed_pdb_path'): 'Written only when diagnostics are on.',
    ('cholesterol', 'merged_pdb_path'): 'Written only when diagnostics are on.',
    ('charmprot', 'keep_lipids'): 'Keep every lipid and sterol, even where it overlaps the incoming protein. Only clashing water, ions and other molecules are removed. Minimize before simulating.',
    ('charmprot', 'environment_overrides'): 'Classify molecule types the automatic detector does not recognise, for example a custom lipid MYLIP → lipid. Unrecognised types are treated as ligands.',
    ('charmprot', 'ligands'): 'auto: reference proteins and ligands are removed; transplant proteins and ligands are inserted. ignore: reference proteins and ligands are removed; only transplant proteins are inserted (its ligands are left out).',
    ('minimization', 'restraint_residue_classes'): 'Only needed for residues that are not recognised automatically, for example a custom lipid MYLIP → lipid.',
    ('ligand_replace', 'new_ligand', 'forcefield_path'): 'The forcefield.itp CHARMM-GUI generated with this ligand ITP; it holds the bond, angle and dihedral parameters. Leave empty to use a forcefield.itp next to the ITP.',
    ('replacement_ligands', '*', 'forcefield_path'): 'The forcefield.itp CHARMM-GUI generated with this ligand ITP. Leave empty to use a forcefield.itp next to the ITP.',
    ('ligand_replace', 'fit', 'method'): 'autofit: pair atoms with the same name. pairfit: use the atom lists below. mcsfit: find the largest shared substructure automatically and report the pairs it chose. nofit: keep the incoming coordinates as they are.',
    ('minimization', 'include_dir'): 'Folder that the topology #include lines are relative to. Leave empty to use the topology folder.',
}
# The only minimization settings shown; every other one keeps its default.
MINIMIZATION_FIELDS = {'restraint_force_constant_kj_mol_nm2', 'max_iterations', 'platform',
                       'restraint_residue_classes'}
STANDALONE_MINIMIZATION_FIELDS = {'coordinates_path', 'topology_path', 'include_dir'}
PLACEHOLDERS = {
    'box_dimensions': 'a b c α β γ   (Å and degrees)',
    'membrane_z_override': 'lower upper   (Å)',
}


def title(key):
    text = str(key).replace('_', ' ')
    return TITLES.get(str(key), text[:1].upper() + text[1:])


def lookup(table, path):
    """Exact path match, then with list indices, then the last key, as '*' wildcards."""
    if not path:
        return None
    path = tuple(path)
    indexless = tuple('*' if isinstance(k, int) else k for k in path)
    for candidate in (path, indexless, path[:-1] + ('*',)):
        if candidate in table:
            return table[candidate]
    return None


def referenced_path_types(raw):
    """How each ${name} registry entry is used: file, directory, or a GUI-managed output."""
    roles, usages = {}, {}

    def visit(value, path=()):
        if isinstance(value, dict):
            for key, child in value.items():
                visit(child, path + (key,))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, path + (index,))
        elif isinstance(value, str) and value.startswith('${') and value.endswith('}'):
            role = path_kind(path)
            if role:
                name = value[2:-1]
                roles[name] = role
                usages.setdefault(name, []).append(is_output(path))
    visit(raw)
    for name, outputs in usages.items():
        if all(outputs):
            roles[name] = 'managed_directory' if roles[name] == 'directory' else 'output'
    return roles


def choose_path(parent, kind, initial, callback, save=False):
    """A non-modal Qt file chooser.

    On WSLg a modal dialog can open behind the main window and leave it
    unresponsive; a non-modal chooser never blocks the rest of the window.
    """
    dialog = QFileDialog(parent.window(), 'Choose folder' if kind == 'directory' else 'Choose file', initial)
    dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)
    dialog.setWindowModality(Qt.WindowModality.NonModal)
    if kind == 'directory':
        dialog.setFileMode(QFileDialog.FileMode.Directory)
        dialog.setOption(QFileDialog.Option.ShowDirsOnly, True)
    elif save:
        dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
    else:
        dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        if kind == 'yaml':
            dialog.setNameFilters(['YAML (*.yaml *.yml)', 'All files (*)'])
        elif kind != 'executable':
            dialog.setNameFilters(['Molecular inputs (*.gro *.pdb *.mol2 *.sdf *.itp *.top)', 'All files (*)'])
    dialog.fileSelected.connect(callback)
    dialog.finished.connect(dialog.deleteLater)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return dialog


def initial_folder(text):
    """The chooser opens on the value already in the field, never an assumed folder."""
    if text and not text.startswith('${'):
        p = Path(text).expanduser()
        if p.is_absolute():
            return str(p if p.is_dir() else p.parent)
    return ''


def clear_layout(layout):
    while layout.count():
        item = layout.takeAt(0)
        if item.widget():
            item.widget().deleteLater()
        elif item.layout():
            clear_layout(item.layout())
            item.layout().deleteLater()


def _list_element(annotation):
    annotation = concrete(annotation)
    if get_origin(annotation) is list:
        args = get_args(annotation)
        return concrete(args[0]) if args else None
    return None


class Context:
    def __init__(self, mode=None, path_types=None):
        self.mode, self.path_types = mode, path_types or {}


# Settings shown first in a mode's own section; the rest go under Advanced settings.
PRIMARY = {
    ('charmprot',): ['reference', 'transplant'],
    ('addbinder',): ['host', 'binder_coordinates', 'binder_itp', 'binder_forcefield', 'side',
                     'distance', 'min_image_gap', 'same_orientation', 'poses',
                     'random_poses', 'random_distance'],
}
# Inputs that are always written, never behind a checkbox.
REQUIRED = {('charmprot', 'reference'), ('charmprot', 'transplant'),
            ('addbinder', 'host'), ('addbinder', 'binder_coordinates'), ('addbinder', 'binder_itp')}
NOUNS = {('replacement_ligands',): 'ligand', ('addbinder', 'poses'): 'pose'}


# Sections whose `enabled` flag is not a user choice: index.ndx is always written.
ALWAYS_ENABLED = {('ndx',)}


def is_toggle_section(path, default, ctx):
    return (isinstance(default, dict) and isinstance(default.get('enabled'), bool)
            and tuple(path) not in ALWAYS_ENABLED
            and not (ctx.mode == 'minimize' and path == ('minimization',)))


# Deprecated settings: rejected or ignored by the pipeline, so never offered.
# Values in a loaded YAML are still passed through unchanged.
DEPRECATED = {('charge', 'membrane_z_margin'), ('charge', 'membrane_z_override'), ('alignment', 'atoms')}


def is_hidden(path, ctx, toggle_parent=False):
    key = path[-1]
    if tuple(path) in DEPRECATED or (tuple(path[:-1]) in ALWAYS_ENABLED and key == 'enabled'):
        return True
    if toggle_parent and key == 'enabled':
        return True
    if path in {('charmprot', 'output_dir'), ('addbinder', 'output_dir')}:
        return True
    if path[:1] == ('minimization',) and len(path) == 2:
        shown = MINIMIZATION_FIELDS | (STANDALONE_MINIMIZATION_FIELDS if ctx.mode == 'minimize' else set())
        return key not in shown
    if path[:1] == ('paths',) and len(path) == 2:
        return ctx.path_types.get(key) == 'managed_directory'
    return is_output(path) and path_kind(path) == 'directory'


class ValueEditor(QWidget):
    """The control for one value: a group of settings, a list, a switch, a number or text."""

    def __init__(self, value, default=None, path=(), parent=None, path_types=None, mode=None,
                 label=None, context=None, present=True):
        super().__init__(parent)
        self.ctx = context or Context(mode, path_types)
        # False when the value is only a default: none of its keys came from the document.
        self.present = present
        self.path, self.default = tuple(path), default
        self.original_value = deepcopy(value)
        self.annotation = annotation_for(self.path)
        self.label = label or (title(self.path[-1]) if self.path else '')
        self.rows, self.scalar, self.charge_source = [], None, None
        self.outer = QVBoxLayout(self)
        self.outer.setContentsMargins(0, 0, 0, 0)
        self.outer.setSpacing(6)
        if self.path[:1] == ('replacement_ligands',) and self.path[-1:] == ('charge',):
            self.charge_source = QComboBox()
            self.charge_source.addItems(['Read charge from ITP', 'Enter charge (e)'])
            self.charge_source.setCurrentIndex(0 if value == 'from_itp' else 1)
            self.outer.addWidget(self.charge_source)
        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 0, 0, 0)
        self.body_layout.setSpacing(6)
        self.outer.addWidget(self.body)
        if self.charge_source:
            self.render(0.0 if value == 'from_itp' else value)
            self.body.setVisible(self.charge_source.currentIndex() == 1)
            self.charge_source.currentIndexChanged.connect(lambda i: self.body.setVisible(i == 1))
        else:
            self.render(initial_value(self.annotation, self.path) if value is None else value)

    # -- construction -----------------------------------------------------
    def classify(self, value):
        hint = concrete(self.annotation)
        element = _list_element(self.annotation)
        if self.charge_source:
            return 'Number'
        if isinstance(value, dict):
            if is_dataclass(hint) or self.path in {('paths',), ('charmprot',), ('addbinder',)} or self.lipid_target():
                return 'Group'
            return 'Entries'
        # Residue replacement entries accept one name or a list of names.
        if (self.path[:2] == ('name_restoration', 'pdb_to_full_resname') and len(self.path) == 3
                and isinstance(value, (str, list))):
            return 'Names'
        if isinstance(value, list):
            if is_dataclass(element) or (value and all(isinstance(v, dict) for v in value)):
                return 'Items'
            if path_kind(self.path) and not is_output(self.path):
                return 'Paths'
            if element in (int, float) or (value and all(type(v) in (int, float) for v in value)):
                return 'Numbers'
            return 'Names'
        if hint is bool or isinstance(value, bool):
            return 'Boolean'
        if hint is int or (hint is not float and type(value) is int):
            return 'Integer'
        if hint is float or isinstance(value, float):
            return 'Number'
        return 'Text'

    def lipid_target(self):
        return len(self.path) == 4 and self.path[:3] == ('cholesterol', 'composition', 'lipid_targets')

    def render(self, value):
        clear_layout(self.body_layout)
        self.rows, self.scalar = [], None
        self.kind = self.classify(value)
        if self.kind == 'Group':
            self.render_group(value)
        elif self.kind == 'Entries':
            self.render_entries(value)
        elif self.kind == 'Items':
            self.render_items(value)
        elif self.kind == 'Paths':
            self.render_paths(value)
        else:
            self.render_scalar(value)

    def render_group(self, value):
        default = self.default if isinstance(self.default, dict) else {}
        keys = list(dict.fromkeys([*value, *default]))
        toggle_parent = is_toggle_section(self.path, default or value, self.ctx)
        self.passthrough = {}
        layout = self.body_layout
        primary = PRIMARY.get(self.path)
        if primary:
            keys = primary + [k for k in keys if k not in primary]
        for key in keys:
            child_path = self.path + (key,)
            if is_hidden(child_path, self.ctx, toggle_parent):
                if self.present and key in value:
                    # An always-on section keeps a loaded `enabled` key but never as false.
                    self.passthrough[key] = (True if self.path in ALWAYS_ENABLED and key == 'enabled'
                                             else deepcopy(value[key]))
                continue
            if primary and key not in primary and layout is self.body_layout:
                layout = self.add_disclosure('Advanced settings')
            child_default = default.get(key)
            if child_default is None and is_dataclass(concrete(annotation_for(child_path))):
                child_default = initial_value(annotation_for(child_path), child_path)
            required = (self.path == ('paths',) or self.lipid_target()
                        or child_path in REQUIRED)
            field = Field(key, value.get(key, deepcopy(child_default)), child_default,
                          self.present and key in value, self.path, self.ctx,
                          kind='required' if required else None, known=key in default or required)
            self.rows.append(field)
            layout.addWidget(field)

    def add_disclosure(self, text):
        toggle = QPushButton('▸  ' + text)
        toggle.setObjectName('disclosure')
        toggle.setCheckable(True)
        panel = QWidget()
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(12, 0, 0, 0)
        panel.hide()

        def expand(shown):
            panel.setVisible(shown)
            toggle.setText(('▾  ' if shown else '▸  ') + text)
        toggle.toggled.connect(expand)
        self.body_layout.addWidget(toggle, 0, Qt.AlignmentFlag.AlignLeft)
        self.body_layout.addWidget(panel)
        self.advanced_toggle = toggle
        return panel_layout

    def render_entries(self, value):
        """Named entries (residue maps, overrides) with an inline add row."""
        self.items_layout = QVBoxLayout()
        self.items_layout.setSpacing(4)
        self.body_layout.addLayout(self.items_layout)
        for key, entry in value.items():
            self.add_entry(key, entry)
        row = QHBoxLayout()
        self.new_key = QLineEdit()
        self.new_key.setPlaceholderText('Residue name, e.g. POPC' if residue_mapping(self.path)
                                        else 'Molecule type, e.g. MYLIP' if self.path in {
                                            ('charmprot', 'environment_overrides'), ('addbinder', 'environment_overrides')}
                                        else 'Name')
        self.new_key.setMaximumWidth(260)
        self.add_button = QPushButton('Add')
        self.add_button.setEnabled(False)
        self.new_key.textChanged.connect(lambda text: self.add_button.setEnabled(bool(text.strip())))
        self.new_key.returnPressed.connect(lambda: self.add_item(self.new_key.text()))
        self.add_button.clicked.connect(lambda: self.add_item(self.new_key.text()))
        self.add_note = QLabel('')
        self.add_note.setObjectName('muted')
        row.addWidget(self.new_key)
        row.addWidget(self.add_button)
        row.addWidget(self.add_note)
        row.addStretch()
        self.body_layout.addLayout(row)

    def add_entry(self, key, value, present=True):
        field = Field(key, value, None, present, self.path, self.ctx, kind='entry')
        field.remove_button.clicked.connect(lambda: self.remove_row(field))
        self.rows.append(field)
        self.items_layout.addWidget(field)
        return field

    def add_item(self, key=None):
        if self.kind == 'Entries':
            key = (key or '').strip()
            if residue_mapping(self.path):
                key = key.upper()
            if not key:
                return None
            if key in [r.key for r in self.rows]:
                self.add_note.setText(f'{key} is already listed.')
                return None
            self.add_note.setText('')
            self.new_key.clear()
            path = self.path + (key,)
            choices = lookup(ENUMS, path)
            return self.add_entry(key, choices[0] if choices else initial_value(annotation_for(path), path), False)
        if self.kind == 'Items':
            return self.add_list_item(initial_value(annotation_for(self.path + (0,)), self.path + (0,)), False)
        if self.kind == 'Paths':
            return self.add_path(key or '')
        return None

    def render_items(self, value):
        self.items_layout = QVBoxLayout()
        self.body_layout.addLayout(self.items_layout)
        for entry in value:
            self.add_list_item(entry)
        noun = NOUNS.get(self.path, 'entry')
        add = QPushButton(f'Add {noun}')
        add.clicked.connect(lambda: self.add_item())
        self.body_layout.addWidget(add, 0, Qt.AlignmentFlag.AlignLeft)

    def add_list_item(self, entry, present=True):
        index = len(self.rows)
        noun = NOUNS.get(self.path, 'entry').capitalize()
        field = Field(index, entry, initial_value(annotation_for(self.path + (index,)), self.path + (index,)),
                      present, self.path, self.ctx, kind='entry', label=f'{noun} {index + 1}')
        field.remove_button.clicked.connect(lambda: self.remove_row(field))
        self.rows.append(field)
        self.items_layout.addWidget(field)
        return field

    def render_paths(self, value):
        self.items_layout = QVBoxLayout()
        self.body_layout.addLayout(self.items_layout)
        for entry in value:
            self.add_path(entry)
        add = QPushButton('Add file…')
        add.clicked.connect(lambda: choose_path(self, 'file', '', lambda p: self.add_path(p)))
        self.body_layout.addWidget(add, 0, Qt.AlignmentFlag.AlignLeft)

    def add_path(self, text):
        field = Field(len(self.rows), text, None, True, self.path, self.ctx, kind='entry', label='')
        field.remove_button.clicked.connect(lambda: self.remove_row(field))
        self.rows.append(field)
        self.items_layout.addWidget(field)
        return field

    def remove_row(self, field):
        if field in self.rows:
            self.rows.remove(field)
            field.deleteLater()

    def render_scalar(self, value):
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        self.body_layout.addLayout(row)
        kind = self.kind
        if kind == 'Boolean':
            self.scalar = QCheckBox(self.label)
            self.scalar.setChecked(bool(value))
        elif kind in {'Names', 'Numbers'}:
            items = value if isinstance(value, list) else [value]
            self.scalar = QLineEdit(', '.join(str(v) for v in items))
            self.scalar.setPlaceholderText(PLACEHOLDERS.get(self.path[-1] if self.path else '',
                                                            'Comma-separated, e.g. PROA, PROB'
                                                            if kind == 'Names' else 'Numbers separated by spaces'))
        else:
            choices = lookup(ENUMS, self.path)
            if self.path[-1:] in [('format',), ('target_system_format',)]:
                choices = ['auto', 'gro', 'pdb'] + (['mol2', 'sdf'] if self.path == (
                    'ligand_replace', 'new_ligand', 'format') else [])
            if choices and kind == 'Text':
                self.scalar = QComboBox()
                self.scalar.addItems(choices)
                if str(value) not in choices:
                    self.scalar.addItem(str(value))
                self.scalar.setCurrentText(str(value))
            else:
                self.scalar = QLineEdit(str(value))
                if kind == 'Integer':
                    self.scalar.setValidator(QIntValidator(self.scalar))
                elif kind == 'Number':
                    validator = QDoubleValidator(self.scalar)
                    validator.setLocale(QLocale.c())
                    self.scalar.setValidator(validator)
            self.scalar.setMinimumWidth(180)
        row.addWidget(self.scalar, 1)
        role = (self.ctx.path_types.get(self.path[-1], 'file') if self.path[:1] == ('paths',) and len(self.path) == 2
                else path_kind(self.path))
        if kind == 'Text' and role in {'file', 'directory', 'executable'} and not is_output(self.path):
            browse = QPushButton('Folder…' if role == 'directory' else 'Browse…')
            browse.clicked.connect(lambda: self.browse(role == 'directory'))
            row.addWidget(browse)
        self.setToolTip('.'.join(map(str, self.path)))

    # -- interaction ------------------------------------------------------
    def set_path(self, path):
        if path and isinstance(self.scalar, QLineEdit):
            self.scalar.setText(path)

    def browse(self, directory=False):
        kind = 'directory' if directory else (path_kind(self.path) or 'file')
        return choose_path(self, kind, initial_folder(self.scalar.text()), self.set_path)

    # -- value ------------------------------------------------------------
    def value(self):
        if self.charge_source and self.charge_source.currentIndex() == 0:
            return 'from_itp'
        kind = self.kind
        if kind == 'Group':
            result = {}
            order = list(self.original_value) if isinstance(self.original_value, dict) else []
            values = {f.key: f.value() for f in self.rows if f.include()}
            values.update(self.passthrough)
            for key in order + [k for k in values if k not in order]:
                if key in values:
                    result[key] = values[key]
            return result
        if kind == 'Entries':
            return {f.key: f.value() for f in self.rows}
        if kind in {'Items', 'Paths'}:
            return [f.value() for f in self.rows]
        if kind == 'Boolean':
            return self.scalar.isChecked()
        text = self.scalar.currentText() if isinstance(self.scalar, QComboBox) else self.scalar.text()
        name = '.'.join(map(str, self.path)) or 'value'
        if kind in {'Names', 'Numbers'}:
            parts = [p for p in re.split(r'[,\s]+', text.strip()) if p]
            if kind == 'Numbers':
                try:
                    parts = [int(p) if type(self.original_value) is list and p.lstrip('-').isdigit()
                             and all(type(v) is int for v in self.original_value) else float(p) for p in parts]
                except ValueError as exc:
                    raise ValueError(f'{name}: enter numbers separated by spaces.') from exc
            if isinstance(self.original_value, str) and len(parts) == 1:
                return parts[0]
            return parts
        try:
            if kind in {'Number', 'Integer'} and type(self.original_value) in (int, float) and text == str(self.original_value):
                return self.original_value
            return int(text) if kind == 'Integer' else float(text) if kind == 'Number' else text
        except ValueError as exc:
            raise ValueError(f'{name}: enter a valid {kind.lower()}.') from exc


class Field(QFrame):
    """One labelled setting and the rule for whether it is written to YAML.

    kinds: required (always written), entry (a user-added, removable entry),
    toggle (a section switched on by its one checkbox), optional (no default;
    written only when ticked), boolean and fixed (written when loaded or changed).
    """

    def __init__(self, key, value, default, present, parent_path, ctx, kind=None, known=True, label=None):
        super().__init__()
        self.setObjectName('fieldRow')
        self.key, self.present, self.default = key, present, deepcopy(default)
        self.path = tuple(parent_path) + (key,)
        self.original = deepcopy(value)
        self.was_null = present and value is None
        annotation = concrete(annotation_for(self.path))
        if kind is None:
            if is_toggle_section(self.path, default if isinstance(default, dict) else value, ctx):
                kind = 'toggle'
            elif annotation is bool or isinstance(default, bool) or (not known and isinstance(value, bool)):
                kind = 'boolean'
            elif default is None and known:
                kind = 'optional'
            elif (default == {} and not is_dataclass(annotation)) or (
                    default == [] and path_kind(self.path) and self.path != ('replacement_ligands',)):
                kind = 'optional'
            else:
                kind = 'fixed'
        self.kind = kind
        text = label if label is not None else (
            str(key) if residue_mapping(parent_path) else lookup(PATH_TITLES, self.path) or title(key))
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 4, 0, 4)
        layout.setSpacing(4)
        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(header)
        self.toggle = None
        self.remove_button = None
        if kind in {'toggle', 'optional'}:
            self.toggle = QCheckBox(TOGGLES.get(str(key), text) if kind == 'toggle' else text)
            self.toggle.setObjectName('fieldToggle')
            self.toggle.setChecked(bool(value.get('enabled')) if kind == 'toggle' and isinstance(value, dict)
                                   else present and value is not None)
            header.addWidget(self.toggle)
        elif kind != 'boolean' and text:
            caption = QLabel(text)
            caption.setObjectName('fieldLabel')
            header.addWidget(caption)
        header.addStretch()
        if kind == 'entry':
            self.remove_button = QPushButton('Remove')
            self.remove_button.setObjectName('quiet')
            header.addWidget(self.remove_button)
        child_default = default
        if isinstance(value, dict) and child_default is None and is_dataclass(annotation):
            child_default = initial_value(annotation, self.path)
        self.editor = ValueEditor(value, child_default, self.path, context=ctx, label=text, present=present)
        if kind == 'boolean':
            self.editor.scalar.setObjectName('fieldToggle')
        if kind == 'entry' and self.editor.kind not in {'Group', 'Entries', 'Items', 'Paths'}:
            # A named entry fits on one line: name, value, Remove.
            header.insertWidget(1 if text else 0, self.editor, 1)
        else:
            body = QHBoxLayout()
            indent = 26 if kind in {'toggle', 'optional'} else 0
            body.setContentsMargins(indent, 0, 0, 0)
            body.addWidget(self.editor)
            layout.addLayout(body)
        description = lookup(DESCRIPTIONS, self.path)
        if description:
            note = QLabel(description)
            note.setWordWrap(True)
            note.setObjectName('muted')
            note.setContentsMargins(26 if kind in {'toggle', 'optional', 'boolean'} else 0, 0, 0, 0)
            layout.insertWidget(1 if kind in {'toggle', 'optional'} else layout.count(), note)
        if self.toggle is not None:
            self.editor.setVisible(self.toggle.isChecked())
            self.toggle.toggled.connect(self.editor.setVisible)

    def include(self):
        if self.kind in {'required', 'entry'}:
            return True
        if self.kind == 'toggle':
            return self.toggle.isChecked() or self.present
        if self.kind == 'optional':
            return self.toggle.isChecked() or self.was_null
        if self.editor.kind == 'Group':
            return self.present or bool(self.editor.value())
        return self.present or self.editor.value() != self.default

    def value(self):
        if self.kind == 'optional' and not self.toggle.isChecked():
            return None
        result = self.editor.value()
        if self.kind == 'toggle':
            if self.toggle.isChecked() or (isinstance(self.original, dict) and 'enabled' in self.original):
                result['enabled'] = self.toggle.isChecked()
        return result


class ConfigEditor(QWidget):
    def __init__(self, raw, mode):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter()
        self.sections = QListWidget()
        self.sections.setMinimumWidth(200)
        self.pages = QStackedWidget()
        splitter.addWidget(self.sections)
        splitter.addWidget(self.pages)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter)
        self.fields = {}
        schema = defaults(mode)
        ctx = Context(mode, referenced_path_types(raw))
        # Nested dataclass fields omitted by the YAML are still accessible.
        for key in dict.fromkeys([*raw, *schema]):
            if key == 'paths' and not raw.get('paths'):
                continue  # No file registry in this document: nothing to show.
            page = QWidget()
            page.setObjectName('page')
            box = QVBoxLayout(page)
            box.setContentsMargins(18, 14, 18, 14)
            heading = QLabel(title(key))
            heading.setObjectName('sectionTitle')
            box.addWidget(heading)
            help_label = QLabel(HELP.get(key, 'Settings are pre-filled with their defaults.'))
            help_label.setWordWrap(True)
            help_label.setObjectName('muted')
            box.addWidget(help_label)
            # The file registry and the CHARMM-GUI folders are the inputs themselves:
            # always shown and always written, never behind a checkbox.
            required = key in {'charmprot', 'addbinder', 'paths'} or (mode == 'minimize' and key == 'minimization')
            field = Field(key, raw.get(key, deepcopy(schema.get(key))), schema.get(key), key in raw, (),
                          ctx, kind='required' if required else None, known=key in schema,
                          label=title(key) if key not in schema or schema.get(key) is None else '')
            box.addWidget(field)
            box.addStretch()
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(page)
            self.pages.addWidget(scroll)
            self.sections.addItem(title(key))
            self.fields[key] = field
        self.sections.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.sections.setCurrentRow(0)
        splitter.setSizes([215, 700])

    def value(self):
        return {key: field.value() for key, field in self.fields.items() if field.include()}
