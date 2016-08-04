import doctest
import unittest

import mast.molecule as mastmol
import mast.system as mastsys
import mast.tests.data as mastdata

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
        self.mock_attrs = {}
        self.bond_map = {0:(0,1)}
        self.mock_attrs[mastmolconfig.MOLECULE_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_attrs['undefined_attribute'] = "undefined_mock_attribute"

        self.MockMoleculeType = mastmol.MoleculeType.factory("MockMoleculeType",
                                                        atom_types=self.atom_types,
                                                        bond_types=self.bond_types,
                                                        bond_map=self.bond_map,
                                                        **self.mock_attrs)

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

        self.system_attributes = {'name' : 'mock_system'}
        self.member_types = [self.MockMoleculeType, self.MockMoleculeType]

    def tearDown(self):
        pass

    def test_factory_molecules_only(self):
        MockSystemType = mastsys.SystemType.factory("MockSystemType",
                                                    member_types=self.member_types,
                                                    **system_attrs)
        for attribute in mastsysconfig.SYSTEM_ATTRIBUTES:
            self.assertIn(MockSystemType.attributes)

        self.assertEqual(len(MockSystemType.member_types), 2)
        for member_type in MockSystemType.member_types:
            self.assertIs(member_type, self.MockMoleculeType)

        for molecule_type in MockSystemType.molecule_types():
            self.assertIs(molecule_type, self.MockMoleculeType)

        self.assertFalse(MockSystemType.atom_types())
        self.assertFalse(MockSystemType.association_types)
        self.assertIs(MockSystemType.member_type_library, {self.MockMoleculeType})
        self.assertEqual(MockSystemType.name, self.system_attributes['name'])

    def test_factory_atoms_only(self):
        pass

    def test_factory_atoms_and_molecules(self):
        pass


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
