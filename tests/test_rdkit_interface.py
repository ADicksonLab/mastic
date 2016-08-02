import unittest
import doctest

import numpy as np

import mast.molecule as mastmol
from mast.interfaces.rdkit import RDKitMoleculeWrapper

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
        self.atom_idx = 0
        self.bond_idx = 0
        self.atom = self.wrapper.atoms[self.atom_idx]
        self.bond = self.wrapper.bonds[self.bond_idx]
    def tearDown(self):
        pass

    def test_wrapper_conformers(self):
        self.assertIsInstance(self.wrapper.get_conformer_coords(0), np.ndarray)


    def test_atom_data(self):
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['atomic_num'],
                          self.atom.GetAtomicNum())

        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['bond_degree_no_Hs'],
                          self.atom.GetDegree())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['bond_degree'],
                          self.atom.GetDegree())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['bond_degree_with_Hs'],
                          self.atom.GetTotalDegree())
        # same but want a convenience attribute
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['total_bond_degree'],
                          self.atom.GetTotalDegree())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['explicit_valence'],
                          self.atom.GetExplicitValence())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['implicit_valence'],
                          self.atom.GetImplicitValence())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['total_valence'],
                          self.atom.GetTotalValence())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['formal_charge'],
                          self.atom.GetFormalCharge())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['hybridization'],
                          self.atom.GetHybridization())

        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['is_aromatic'],
                          self.atom.GetIsAromatic())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['in_ring'],
                          self.atom.IsInRing())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['isotope'],
                          self.atom.GetIsotope())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['mass'],
                          self.atom.GetMass())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['num_radical_electrons'],
                          self.atom.GetNumRadicalElectrons())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['element'],
                          self.atom.GetSymbol())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['num_Hs'],
                          self.atom.GetTotalNumHs())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['rdkit_mol_idx'],
                          self.atom.GetIdx())
        monomer_info = self.atom.GetMonomerInfo()
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['monomer_type'],
                          monomer_info.GetMonomerType())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_name'],
                          monomer_info.GetName().strip())
        # self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_chain_id'],
        # monomer_info.GetChainID())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_insertion_code'],
                          monomer_info.GetInsertionCode())
        # self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_heteroatom'],
        # monomer_info.IsHeteroAtom())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_occupancy'],
                          monomer_info.GetOccupancy())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_residue_name'],
                          monomer_info.GetResidueName())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_residue_number'],
                          monomer_info.GetResidueNumber())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_serial_number'],
                          monomer_info.GetSerialNumber())
        # self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_segment_number'],
        # monomer_info.GetSegmentNumber())
        self.assertEqual(self.wrapper.atom_data(self.atom_idx)['pdb_temp_factor'],
                          monomer_info.GetTempFactor())

    def test_bond_data(self):
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['bond_type'],
                         str(self.bond.GetBondType()))
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['bond_type_number'],
                         str(self.bond.GetBondTypeAsDouble()))
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['is_aromatic'],
                         self.bond.GetIsAromatic())
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['in_ring'],
                         self.bond.IsInRing())
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['stereo'],
                         str(self.bond.GetStereo()))
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['is_conjugated'],
                         self.bond.GetIsConjugated())
        atom1_idx = self.bond.GetBeginAtomIdx()
        atom2_idx = self.bond.GetEndAtomIdx()
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['rdkit_atom_idxs'],
                         (atom1_idx, atom2_idx))
        self.assertEqual(self.wrapper.bond_data(self.bond_idx)['rdkit_mol_idx'],
                         self.bond.GetIdx())


    def test_bonds_map(self):
        for bond_idx, atom_tuple in self.wrapper.bonds_map().items():
            atom1_idx = atom_tuple[0]
            atom2_idx = atom_tuple[1]
            bond = self.wrapper.bonds[bond_idx]
            self.assertEqual(bond.GetBeginAtomIdx(), atom1_idx)
            self.assertEqual(bond.GetEndAtomIdx(), atom2_idx)

    # test only the attributes we need for extracting useful information
    def test_molecule_data(self):
        ring_info = self.wrapper.rdkit_molecule.GetRingInfo()
        self.assertEqual(self.wrapper.molecule_data()['num_rings'],
                          ring_info.NumRings())
        self.assertEqual(self.wrapper.molecule_data()['num_atoms'],
                    self.wrapper.rdkit_molecule.GetNumAtoms())
        self.assertEqual(self.wrapper.molecule_data()['num_bonds'],
                          self.wrapper.rdkit_molecule.GetNumBonds())
        self.assertEqual(self.wrapper.molecule_data()['num_heavy_atoms'],
                          self.wrapper.rdkit_molecule.GetNumHeavyAtoms())

    def test_find_features(self):
        features = self.wrapper.find_features()
        atom_idxs = range(len(self.wrapper.atoms))
        for feature_dict in features.values():
            for atom_id in feature_dict['atom_ids']:
                self.assertIn(atom_id, atom_idxs)
            self.assertIsInstance(feature_dict['position'], tuple)
            self.assertTrue(len(feature_dict['position']), 3)
            for dim in feature_dict['position']:
                self.assertIsInstance(dim, float)

    def test_make_molecule_type(self):
        test_mol_type = self.wrapper.make_molecule_type()
        

if __name__ == "__main__":

    from mast.interfaces import rdkit

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(rdkit, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
