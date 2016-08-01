import doctest
import unittest

import numpy as np
import numpy.testing as npt

from mast import molecule

import mast.molecule as mastmol
import mast.selection as mastsel
import mast.test.config.molecule as mastmolconfig

class TestAtomType(unittest.TestCase):
    def setUp(self):
        mock_attrs = {}
        mock_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        mock_attrs['undefined_attribute'] = "undefined_mock_attribute"

    def tearDown(self):
        pass

    def test_factory(self):
        MockAtomType = mastmol.AtomType.factory("MockAtomType", **mock_attrs)
        self.assertIsInstance(MockAtom, type)

class TestBondType(unittest.TestCase):
    def setUp(self):
        mock_attrs = {}
        mock_attrs[mastmolconfig.BOND_ATTRIBUTES[0]] = "mock_attribute"
        mock_attrs['undefined_attribute'] = "undefined_mock_attribute"

    def tearDown(self):
        pass

    def test_factory(self):
        MockBondType = mastmol.BondType.factory("MockBondType", **mock_attrs)
        self.assertIsInstance(MockBond, type)

class TestMoleculeType(unittest.TestCase):
    def setUp(self):
        mock_attrs = {}
        mock_attrs[mastmolconfig.BOND_ATTRIBUTES[0]] = "mock_attribute"
        mock_attrs['undefined_attribute'] = "undefined_mock_attribute"

    def tearDown(self):
        pass

    def test_factory(self):
        MockMoleculeType = mastmol.MoleculeType.factory("MockMoleculeType", **mock_attrs)
        self.assertIsInstance(MockMolecule, type)



if __name__ == "__main__":

    from mast import molecule

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(molecule, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
