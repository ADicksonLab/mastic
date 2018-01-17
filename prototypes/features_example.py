from rdkit import Chem
from rdkit.Chem import AllChem
import os.path as osp

import mastic.selection as masticsel
import mastic.molecule as masticmol
import mastic.features as masticfeat

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

import mastic.config.interactions as masticinxconfig

from rdkit import Chem

tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
pka_rdkit = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
pka_rdkit_wrapper = RDKitMoleculeWrapper(pka_rdkit, mol_name="PKA")

PKAType = pka_rdkit_wrapper.make_molecule_type(find_features=True)
pka_features = pka_rdkit_wrapper.find_features()

feature1 = pka_features[1]
atom_idxs = feature1['atom_ids']
feature1_attrs = {}
feature1_attrs['rdkit_family'] = feature1['family']
feature1_attrs['rdkit_type'] = feature1['type']
feature1_attrs['rdkit_position'] = feature1['position']
MyPKAFeature1Type = masticfeat.FeatureType("MyPKAFeature1Type",
                                               molecule_type=PKAType,
                                               atom_idxs=atom_idxs,
                                               **feature1_attrs)


PKAType.add_feature_type('mine', MyPKAFeature1Type)

pka_coords = pka_rdkit_wrapper.get_conformer_coords(0)
pka_mol = PKAType.to_molecule(pka_coords)

pka_feature_1 = MyPKAFeature1Type.to_feature(molecule=pka_mol)


feature_type = PKAType.feature_types[1]
