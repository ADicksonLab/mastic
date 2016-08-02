import doctest
import unittest

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
        atom1_attrs = {'name' : "FAKE1"}
        self.Atom1Type = mastmol.AtomType.factory("Atom1Type", **atom1_attrs)

        atom2_attrs = {'name' : "FAKE2"}
        self.Atom2Type = mastmol.AtomType.factory("Atom2Type", **atom2_attrs)

        self.atom_types = [self.Atom1Type, self.Atom2Type]
        self.bond_map = {0 : (0,1)}
        self.bond1_attrs = {'bond_type' : "TRIPLE"}

    def tearDown(self):
        pass

    def test_factory(self):
        MockBondType = mastmol.BondType.factory("MockBondType",
                                                atom_types=tuple(self.atom_types),
                                                bond_map=self.bond_map,
                                                **self.bond1_attrs)
        self.assertIsInstance(MockBondType, type)

        # test the non-domain specific attributes work
        self.assertEqual(MockBondType.atom_types, (self.Atom1Type, self.Atom2Type))

class TestFakeMoleculeType(unittest.TestCase):
    def setUp(self):

        atom1_attrs = {'name' : "FAKE1"}
        self.Atom1Type = mastmol.AtomType.factory("Atom1Type", **atom1_attrs)

        atom2_attrs = {'name' : "FAKE2"}
        self.Atom2Type = mastmol.AtomType.factory("Atom2Type", **atom2_attrs)

        self.atom_types = [self.Atom1Type, self.Atom2Type]
        self.bond_map = {0 : (0,1)}
        self.bond1_attrs = {'bond_type' : "TRIPLE"}
        MockBondType = mastmol.BondType.factory("MockBondType",
                                                atom_types=tuple(self.atom_types),
                                                bond_map=self.bond_map,
                                                **self.bond1_attrs)
        self.bond_types = [MockBondType]
        self.mock_attrs = {}


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
        self.assertEqual(begin_atom_type, self.Atom1Type)
        self.assertEqual(end_atom_type, self.Atom2Type)
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



if __name__ == "__main__":

    from mast import molecule

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(molecule, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
