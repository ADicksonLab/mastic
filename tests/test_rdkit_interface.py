import unittest
import doctest

import mast.molecule as mastmol
from mast.interfaces import RDKitMoleculeWrapper

import mast.tests.config.molecule as mastmolconfig
import mast.tests.data as mastdata

from rdkit import Chem

class TestRDKitWrapper(unittest.TestCase):

    def setUp(self):
        # load a string from the data (small-molecule :: sml) as a
        # file-like object for reading
        self.sml_path = mastdata.BEN_path
        self.sml_rdkit = Chem.MolFromPDBFile(self.sml_path, removeHs=False)
        self.wrapper = RDKitMoleculeWrapper(self.sml_rdkit)
        self.atom = self.wrapper
    def tearDown(self):
        pass

    def test_wrapper_conformers(self):
        pass

    def test_bond_map(self):
        pass

    def test_atom_data(self):
        pass

    def test_bond_data(self):
        pass


    # test only the attributes we need for extracting useful information
    def test_molecule_data(self):
        ring_info = self.wrapper.GetRingInfo()
        ring_info.NumRings()
        self.wrapper.rdkit_molecule.GetNumAtoms()
        self.wrapper.rdkit_molecule.GetNumBonds()
        self.wrapper.rdkit_molecule.GetNumHeavyAtoms()

    def test_find_features(self):
        pass


if __name__ == "__main__":

    from mast.interfaces import rdkit

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(rdkit, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
