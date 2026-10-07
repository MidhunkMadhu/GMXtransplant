import unittest
import MDAnalysis as mda
from clashes import detect_and_remove_clashes
from config import ClashDetectionSpec


class RoundingClashTests(unittest.TestCase):
    def run_contact(self, resname, keep_classes=(), offset=0):
        u = mda.Universe.empty(4, n_residues=2, atom_resindex=[0, 1, 1, 1], trajectory=True)
        u.add_TopologyAttr('names', ['O', 'H2', 'OH2', 'H1'])
        u.add_TopologyAttr('resnames', ['ALA', resname])
        u.add_TopologyAttr('resids', [168, 28087])
        u.atoms.positions = [[24.7305393219, 78.4682617188, 40.4241104126],
                            [24.872625351, 79.66217804 + offset, 40.3577461243],
                            [26, 81, 42], [27, 81, 42]]
        u.dimensions = [105.828, 105.828, 192.017, 90, 90, 90]
        before = u.atoms.positions.copy()
        result = detect_and_remove_clashes(u.atoms[1:], u.atoms[:1], u.dimensions,
            ClashDetectionSpec(threshold=1.2, thresholds={'lipid': 0.8},
                               keep_classes=list(keep_classes), heavy_atoms_only=False))
        self.assertTrue((before == u.atoms.positions).all())
        return result

    def test_gro_rounding_removes_complete_water(self):
        result = self.run_contact('TIP3')
        self.assertEqual(len(result.removed_ag), 3)
        self.assertEqual(result.removed_molecules[0].contact_context, 'inserted_block:gro')
        self.assertLess(result.removed_molecules[0].min_distance, 1.2)

    def test_clear_water_stays(self):
        self.assertEqual(len(self.run_contact('TIP3', offset=0.05).removed_ag), 0)

    def test_kept_water_stays_and_is_reported(self):
        result = self.run_contact('TIP3', ['water'])
        self.assertEqual(len(result.removed_ag), 0)
        self.assertEqual(len(result.flagged_but_kept), 1)

    def test_lipid_uses_its_lower_cutoff(self):
        self.assertEqual(len(self.run_contact('POPC').removed_ag), 0)


if __name__ == '__main__':
    unittest.main()
