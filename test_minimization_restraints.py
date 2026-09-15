from types import SimpleNamespace
import unittest

import numpy as np

from config import MinimizationSpec
from minimization_restraints import add_heavy_atom_restraints, molecule_atom_roles


class ClassificationTests(unittest.TestCase):
    def test_protein_caps_custom_lipids_and_ligands(self):
        self.assertEqual(molecule_atom_roles(["ACE", "ALA", "NME"], {}), ["protein"] * 3)
        self.assertEqual(molecule_atom_roles(["LIG"], {}), ["ligand"])
        for name in ["POPC", "CHL1", "CLR", "TIP3", "SOD"]:
            self.assertNotIn(molecule_atom_roles([name], {})[0], ["protein", "ligand"])
        self.assertEqual(molecule_atom_roles(["CUSTOM", "CAP"], {"CUSTOM": "lipid"}), ["lipid"] * 2)
        self.assertEqual(molecule_atom_roles(["ALA", "MOD"], {"MOD": "ligand"}), ["protein", "ligand"])


class RestraintPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import openmm
        except ImportError:
            raise unittest.SkipTest("OpenMM is not installed")

    def test_clash_relaxes_while_protein_and_ligand_heavy_atoms_are_restrained(self):
        import openmm
        from openmm import app, unit
        topology = app.Topology()
        chain = topology.addChain()
        definitions, molecules = {}, []
        # Includes a heavy-mass hydrogen to ensure element-based selection.
        for rn, names, elements in [("ALA", ["CA", "HG"], [app.element.carbon, app.element.hydrogen]),
                                    ("LIG", ["C1"], [app.element.carbon]),
                                    ("POPC", ["C2"], [app.element.carbon]),
                                    ("TIP3", ["OH2"], [app.element.oxygen]),
                                    ("SOD", ["NA"], [app.element.sodium])]:
            residue = topology.addResidue(rn, chain)
            for name, element in zip(names, elements):
                topology.addAtom(name, element, residue)
            definitions[rn] = SimpleNamespace(atoms=[[None, None, None, rn] for _ in names])
            molecules.append((rn, 1))
        top = SimpleNamespace(topology=topology, _molecules=molecules, _moleculeTypes=definitions)
        system = openmm.System()
        nonbonded = openmm.NonbondedForce()
        for index in range(6):
            system.addParticle(4 if index == 1 else 12)
            nonbonded.addParticle(0, 0.3, 1 if index in (0, 3) else 0)
        system.addForce(nonbonded)
        points = np.array([[0, 0, 0], [0, 1, 0], [2, 0, 0], [.23, 0, 0], [2, 2, 0], [3, 3, 0]])
        reference = points * unit.nanometer
        report, selected = add_heavy_atom_restraints(openmm, unit, top, system, reference, MinimizationSpec())
        self.assertEqual(selected, [0, 2])
        self.assertEqual(report["protein_heavy_atoms"], 1)
        self.assertEqual(report["ligand_heavy_atoms"], 1)
        integrator = openmm.VerletIntegrator(.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
        context.setPositions(reference)
        initial = context.getState(getEnergy=True).getPotentialEnergy()
        openmm.LocalEnergyMinimizer.minimize(context, 0.001, 5000)
        state = context.getState(getEnergy=True, getPositions=True)
        final = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        self.assertLess(state.getPotentialEnergy(), initial)
        self.assertGreater(np.linalg.norm(final[3] - final[0]), .33)
        self.assertLess(np.linalg.norm(final[0] - points[0]), .03)
        self.assertLess(np.linalg.norm(final[0] - points[0]),
                        np.linalg.norm(final[3] - points[3]) / 5)
        self.assertGreater(np.linalg.norm(final[3] - points[3]), .09)
        np.testing.assert_allclose(final[2], points[2], atol=1e-6)

    def test_periodic_restraint_uses_nearest_image(self):
        import openmm
        from openmm import app, unit
        topology = app.Topology()
        residue = topology.addResidue("LIG", topology.addChain())
        topology.addAtom("C1", app.element.carbon, residue)
        top = SimpleNamespace(topology=topology, _molecules=[("LIG", 1)],
                              _moleculeTypes={"LIG": SimpleNamespace(atoms=[[None, None, None, "LIG"]])})
        system = openmm.System()
        system.addParticle(12)
        system.setDefaultPeriodicBoxVectors(openmm.Vec3(4, 0, 0), openmm.Vec3(0, 4, 0), openmm.Vec3(0, 0, 4))
        nb = openmm.NonbondedForce()
        nb.setNonbondedMethod(nb.CutoffPeriodic)
        nb.addParticle(0, .3, 0)
        system.addForce(nb)
        add_heavy_atom_restraints(openmm, unit, top, system, np.array([[.1, 0, 0]]) * unit.nanometer, MinimizationSpec())
        integrator = openmm.VerletIntegrator(.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
        context.setPositions([[4.2, 0, 0]])
        state = context.getState(getEnergy=True, getForces=True)
        self.assertAlmostEqual(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 5, places=6)
        self.assertAlmostEqual(state.getForces()[0][0].value_in_unit(unit.kilojoule_per_mole / unit.nanometer), -100, places=6)
