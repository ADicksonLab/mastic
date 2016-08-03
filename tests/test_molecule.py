import doctest
import unittest

import itertools

import numpy as np
import numpy.testing as npt

from mast import molecule

import mast.molecule as mastmol
import mast.selection as mastsel
import mast.tests.config.molecule as mastmolconfig
import mast.tests.data as mastdata

class TestAtomType(unittest.TestCase):
    def setUp(self):
        self.mock_attrs = {}
        self.mock_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_attrs['undefined_attribute'] = "undefined_mock_attribute"

    def tearDown(self):
        pass

    def test_factory(self):
        MockAtomType = mastmol.AtomType.factory("MockAtomType", **self.mock_attrs)
        self.assertIsInstance(MockAtomType, type)

class TestBondType(unittest.TestCase):
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

    def tearDown(self):
        pass

    def test_factory(self):
        MockBondType = mastmol.BondType.factory("MockBondType",
                                                atom_types=self.atom_types,
                                                **self.mock_bond_attrs)

        self.assertIsInstance(MockBondType, type)

        # test the non-domain specific attributes work
        self.assertEqual(MockBondType.atom_types, (self.Mock1AtomType, self.Mock2AtomType))

class TestFakeMoleculeType(unittest.TestCase):
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

    def tearDown(self):
        pass

    def test_factory(self):
        MockMoleculeType = mastmol.MoleculeType.factory("MockMoleculeType",
                                                        atom_types=self.atom_types,
                                                        bond_types=self.bond_types,
                                                        bond_map=self.bond_map,
                                                        **self.mock_attrs)
        self.assertIsInstance(MockMoleculeType, type)

        # test that the non-domain specific attributes and functions
        # work
        self.assertEqual(MockMoleculeType.atom_types, self.atom_types)
        self.assertEqual(MockMoleculeType.atom_type_library,
                          set(MockMoleculeType.atom_types))
        self.assertEqual(MockMoleculeType.bond_types, self.bond_types)
        self.assertEqual(MockMoleculeType.bond_type_library,
                          set(MockMoleculeType.bond_types))
        self.assertEqual(MockMoleculeType.bond_map, self.bond_map)
        # make sure we get the correct AtomTypes from the bond map
        begin_atom_type = MockMoleculeType.atom_types[MockMoleculeType.bond_map[0][0]]
        end_atom_type = MockMoleculeType.atom_types[MockMoleculeType.bond_map[0][1]]
        self.assertEqual(begin_atom_type, self.Mock1AtomType)
        self.assertEqual(end_atom_type, self.Mock2AtomType)
        # we didn't set these so make sure they are empty
        self.assertFalse(MockMoleculeType.features)
        self.assertFalse(MockMoleculeType.feature_families())
        self.assertFalse(MockMoleculeType.feature_families_map())
        self.assertFalse(MockMoleculeType.feature_types())
        self.assertFalse(MockMoleculeType.feature_types_map())

        # test the domain specific stuff is the same as in the mock
        # config files
        for attr in mastmolconfig.MOLECULE_ATTRIBUTES:
            self.assertIn(attr, MockMoleculeType.attributes)


class testAtom(unittest.TestCase):

    def setUp(self):
        # make a MoleculeType
        self.mock_attrs = {}
        self.mock_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.MockAtomType = mastmol.AtomType.factory("MockAtomType", **self.mock_attrs)
        self.coords = np.array((0.0, 0.0, 0.0))

    def tearDown(self):
        pass

    def test_constructors(self):
        atom1 = mastmol.Atom(coords=self.coords, atom_type=self.MockAtomType)
        atom2 = self.MockAtomType.to_atom(self.coords)
        npt.assert_array_almost_equal(atom1.coords, atom2.coords)
        self.assertEqual(atom1.atom_type, atom2.atom_type)

class testBond(unittest.TestCase):

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

        self.coords = (np.array((0.0, 0.0, 0.0)), np.array((0.0, 0.0, 1.0)))
        self.atoms = (self.Mock1AtomType.to_atom(self.coords[0]),
                      self.Mock2AtomType.to_atom(self.coords[1]))
    def tearDown(self):
        pass
    def test_constructors(self):
        bonds = []
        bonds.append(mastmol.Bond(atom_container=self.atoms, bond_type=self.MockBondType))
        bonds.append(mastmol.Bond(atom_container=self.atoms, atom_ids=(0,1),
                                  bond_type=self.MockBondType))
        bonds.append(self.MockBondType.to_bond(*self.coords))

        for bond_a, bond_b in itertools.combinations(bonds, 2):
            npt.assert_array_almost_equal(bond_a.coords, bond_b.coords)
            self.assertEqual(bond_a.bond_type, bond_b.bond_type)

if __name__ == "__main__":

    from mast import molecule

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(molecule, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()