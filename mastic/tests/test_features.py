import doctest
import unittest

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.features as mastfeat
import mast.interactions as mastinx

import mast.tests.data as mastdata

import mast.config.molecule as mastmolconfig
import mast.config.system as mastsysconfig
import mast.config.features as mastfeatconfig
import mast.config.interactions as mastinxconfig


class TestFeatureType(unittest.TestCase):

    def setUp(self):

        self.mock_atom1_attrs = {}
        self.mock_atom1_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_atom1_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock1AtomType = mastmol.AtomType("Mock1AtomType", **self.mock_atom1_attrs)
        self.mock_atom2_attrs = {}
        self.mock_atom2_attrs[mastmolconfig.ATOM_ATTRIBUTES[0]] = "mock_attribute_2"
        self.mock_atom2_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.Mock2AtomType = mastmol.AtomType("Mock2AtomType", **self.mock_atom2_attrs)

        self.atom_types = (self.Mock1AtomType, self.Mock2AtomType)
        self.mock_bond_attrs = {}
        self.mock_bond_attrs[mastmolconfig.BOND_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_bond_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.MockBondType = mastmol.BondType("MockBondType",
                                                atom_types=self.atom_types,
                                                **self.mock_bond_attrs)

        self.bond_types = [self.MockBondType]
        self.mock_attrs = {}
        self.bond_map = {0:(0,1)}
        self.mock_attrs[mastmolconfig.MOLECULE_ATTRIBUTES[0]] = "mock_attribute"
        self.mock_attrs['undefined_attribute'] = "undefined_mock_attribute"
        self.MockMoleculeType = mastmol.MoleculeType("MockMoleculeType",
                                                             atom_types=self.atom_types,
                                                             bond_types=self.bond_types,
                                                             bond_map=self.bond_map,
                                                             **self.mock_attrs)


    def tearDown(self):
        pass

    def test_factory(self):
        atom_idxs = [0]
        feature_attrs = {"mock_attribute" : 35}
        MockFeatureType = mastfeat.FeatureType("MockFeatureType",
                                               molecule_type=self.MockMoleculeType,
                                               atom_idxs=atom_idxs)

    def test_find_features_make_molecule_type(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from mast.interfaces.rdkit import RDKitMoleculeWrapper

        # load a string from the data (small-molecule :: sml) as a
        # file-like object for reading
        sml_path = mastdata.BEN_path
        sml_rdkit = Chem.MolFromPDBFile(sml_path, removeHs=False)
        wrapper = RDKitMoleculeWrapper(sml_rdkit)

        BENType = wrapper.make_molecule_type(find_features=True)
        for feature_id, feature in BENType.feature_types.items():
            self.assertTrue(feature.attributes == mastfeatconfig.FEATURE_ATTRIBUTES)
            for atom_type in feature.atom_types:
                self.assertIn(atom_type, BENType.atom_types)
            for bond_type in feature.bond_types:
                self.assertIn(bond_type, BENType.bond_types)


if __name__ == "__main__":

    from mast import features

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(features, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
