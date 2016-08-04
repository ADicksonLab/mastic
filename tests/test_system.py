import doctest
import unittest

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.tests.data as mastdata

import mast.config.molecule as mastmolconfig
import mast.config.system as mastsysconfig

class TestSystemType(unittest.TestCase):

    def setUp(self):
        self.mock_atom1_attrs = {}
        self.mock_atom1_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_atom1_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock1AtomType = mastmol.AtomType.factory("Mock1AtomType", **self.mock_atom1_attrs)
        self.mock_atom2_attrs = {}
        self.mock_atom2_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute_2"
        self.mock_atom2_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock2AtomType = mastmol.AtomType.factory("Mock2AtomType", **self.mock_atom2_attrs)

        self.atom_types = (self.Mock1AtomType, self.Mock2AtomType)
        self.mock_bond_attrs = {}
        self.mock_bond_attrs[mastmolconfig.BOND_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_bond_attrs['undefined_attribute'] = "undefined_mock_attribute"

        self.MockBondType = mastmol.BondType.factory("MockBondType",
                                                     atom_types=self.atom_types,
                                                     **self.mock_bond_attrs)

        self.bond_types = [self.MockBondType]
        self.mock_mol_attrs = {}
        self.bond_map = {0:(0,1)}
        self.mock_mol_attrs[mastmolconfig.MOLECULE_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_mol_attrs['undefined_attribute'] = "undefined_mock_attribute"

        self.MockMoleculeType = mastmol.MoleculeType.factory("MockMoleculeType",
                                                        atom_types=self.atom_types,
                                                        bond_types=self.bond_types,
                                                        bond_map=self.bond_map,
                                                        **self.mock_mol_attrs)

        self.system_attributes = {'name' : 'mock_system'}
        self.molecule_types = [self.MockMoleculeType, self.MockMoleculeType]
        self.atom_types = [self.Mock1AtomType, self.Mock2AtomType]
        self.member_types = self.atom_types + self.molecule_types
        self.system_attributes = {'name' : 'mock_system'}

    def tearDown(self):
        pass

    def test_factory_molecules_only(self):
        MockSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.molecule_types,
                                                    **self.system_attributes)
        for attribute in mastsysconfig.SYSTEM_ATTRIBUTES:
            self.assertIn(attribute, MockSystemType.attributes)

        self.assertEqual(len(MockSystemType.member_types), 2)
        for member_type in MockSystemType.member_types:
            self.assertIs(member_type, self.MockMoleculeType)

        for molecule_type in MockSystemType.molecule_types():
            self.assertIs(molecule_type, self.MockMoleculeType)

        self.assertFalse(MockSystemType.atom_types())
        self.assertFalse(MockSystemType.association_types)
        self.assertEqual(MockSystemType.member_type_library, {self.MockMoleculeType})
        self.assertEqual(MockSystemType.name, self.system_attributes['name'])

    def test_factory_atoms_only(self):
        MockSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.atom_types,
                                                    **self.system_attributes)

        self.assertEqual(len(MockSystemType.member_types), 2)

        for i, member_type in enumerate(MockSystemType.member_types):
            self.assertIs(member_type, self.atom_types[i])

        for i, atom_type in enumerate(MockSystemType.atom_types()):
            self.assertIs(atom_type, self.atom_types[i])

        self.assertFalse(MockSystemType.molecule_types())
        self.assertFalse(MockSystemType.association_types)
        self.assertEqual(MockSystemType.member_type_library, set(self.atom_types))
        self.assertEqual(MockSystemType.name, self.system_attributes['name'])

    def test_factory_atoms_and_molecules(self):
        MockSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.member_types,
                                                    **self.system_attributes)

        self.assertEqual(len(MockSystemType.member_types), 4)

        for i, member_type in enumerate(MockSystemType.member_types):
            self.assertIs(member_type, self.member_types[i])

        for i, atom_type in enumerate(MockSystemType.atom_types()):
            self.assertIs(atom_type, self.atom_types[i])

        self.assertFalse(MockSystemType.association_types)
        self.assertEqual(MockSystemType.member_type_library, set(self.member_types))
        self.assertEqual(MockSystemType.name, self.system_attributes['name'])

class TestSystem(unittest.TestCase):

    def setUp(self):
        self.mock_atom1_attrs = {}
        self.mock_atom1_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_atom1_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock1AtomType = mastmol.AtomType.factory("Mock1AtomType", **self.mock_atom1_attrs)
        self.mock_atom2_attrs = {}
        self.mock_atom2_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute_2"
        self.mock_atom2_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock2AtomType = mastmol.AtomType.factory("Mock2AtomType", **self.mock_atom2_attrs)

        self.atom_types = (self.Mock1AtomType, self.Mock2AtomType)
        self.mock_bond_attrs = {}
        self.mock_bond_attrs[mastmolconfig.BOND_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_bond_attrs['undefined_attribute'] = "undefined_mock_attribute"

        self.MockBondType = mastmol.BondType.factory("MockBondType",
                                                     atom_types=self.atom_types,
                                                     **self.mock_bond_attrs)

        self.bond_types = [self.MockBondType]
        self.mock_mol_attrs = {}
        self.bond_map = {0:(0,1)}
        self.mock_mol_attrs[mastmolconfig.MOLECULE_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_mol_attrs['undefined_attribute'] = "undefined_mock_attribute"

        self.MockMoleculeType = mastmol.MoleculeType.factory("MockMoleculeType",
                                                        atom_types=self.atom_types,
                                                        bond_types=self.bond_types,
                                                        bond_map=self.bond_map,
                                                        **self.mock_mol_attrs)

        self.system_attributes = {'name' : 'mock_system'}
        self.molecule_types = [self.MockMoleculeType, self.MockMoleculeType]
        self.atom_types = [self.Mock1AtomType, self.Mock2AtomType]
        self.member_types = self.atom_types + self.molecule_types

        self.nonoverlapping_coords = [np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
                                      np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 0.0]])]

        self.overlapping_coords = [np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
                                   np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])]

        self.nonoverlapping_molecules = []
        for coords in self.nonoverlapping_coords:
            self.molecules.append(self.MockMoleculeType.to_molecule(coords))

        self.overlapping_molecules = []
        for coords in self.overlapping_coords:
            self.molecules.append(self.MockMoleculeType.to_molecule(coords))

        self.MockAtomsSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.atom_types,
                                                    **self.system_attributes)
        self.MockMoleculesSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.molecule_types,
                                                    **self.system_attributes)
        self.MockSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.member_types,
                                                    **self.system_attributes)

        def tearDown(self):
            pass

        def test_nonoverlapping_molecules_system(self):
            system = self.MockMoleculesSystemType.to_system(self.nonoverlapping_coords)
            self.assertTrue(isinstance(system, mastsys.System))
            self.assertIs(system.system_type, self.MockMoleculesSystemType)
            # things that should be in the system
            self.assertTrue(len(system.members) == 2)
            self.assertTrue(len(system.molecules) == 2)
            # things that shouldn't be in the system
            self.assertFalse(system.associations)
            self.assertFalse(system.association_types)
            self.assertFalse(system.atoms)
            self.assertFalse(system.atom_types)
            # methods
            # selection of molecules
            mol_sel = system.molecules_sel()
            self.assertTrue(isinstance(mol_sel, mastsel.IndexedSelection))
            for idx, selected_mol in mol_sel.items():
                self.assertIs(selected_mol, system.molecules[idx])
            # flags
            self.assertIn('system', system.flags)
            for member in system.members:
                self.assertIn('system', member.flags)

        def test_overlapping_system(self):
            with self.assertRaises(AssertionError):
                self.MockMoleculesSystemType.to_system(self.overlapping_coords)

        def test_molecules_and_atom_system(self):
            system = self.MockSystemType.to_system(self.nonoverlapping_coords)
            self.assertTrue(isinstance(system, mastsys.System))
            self.assertIs(system.system_type, self.MockMoleculesSystemType)
            # things that should be in the system
            self.assertTrue(len(system.members) == 4)
            self.assertTrue(len(system.molecules) == 2)
            self.assertTrue(len(system.atoms) == 2)
            # things that shouldn't be in the system
            self.assertFalse(system.associations)
            self.assertFalse(system.association_types)
            # methods
            # selection of molecules
            mol_sel = system.molecules_sel()
            self.assertTrue(isinstance(mol_sel, mastsel.IndexedSelection))
            for idx, selected_mol in mol_sel.items():
                self.assertIs(selected_mol, system.molecules[idx])
            # selection of atoms
            atom_sel = system.atoms_sel()
            self.assertTrue(isinstance(atom_sel, mastsel.IndexedSelection))
            for idx, selected_atom in atom_sel.items():
                self.assertIs(selected_atom, system.atoms[idx])
            # flags
            self.assertIn('system', system.flags)
            for member in system.members:
                self.assertIn('system', member.flags)

        def test_substantiated_member_system_uniqueness(self):
            """Test to make sure you cannot add a substantiated member to two
            different systems."""

            # this shouldn't be a big deal if you only ever use the
            # Type.to_xxx method, (Type.substantiate in the future)
            # will wait to write these tests then
            self.assertTrue(False)


class TestSystem(unittest.TestCase):
    def setUp(self):
        pass

if __name__ == "__main__":

    from mast import system

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(system, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
