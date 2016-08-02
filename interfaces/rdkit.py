import os.path as osp
import collections as col

import numpy as np

import mast.molecule as mastmol

from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures


class RDKitMoleculeWrapper(object):
    def __init__(self, rdkit_molecule, mol_name=None):
        self.rdkit_molecule = rdkit_molecule
        self.mol_name = mol_name

    @property
    def atoms(self):
        return [atom for atom in self.rdkit_molecule.GetAtoms()]

    @property
    def bonds(self):
        return [bond for bond in self.rdkit_molecule.GetBonds()]

    @property
    def num_atoms(self):
        return self.rdkit_molecule.GetNumAtoms()

    @property
    def num_bonds(self):
        return self.rdkit_molecule.GetNumBonds()

    @property
    def num_conformers(self):
        return self.rdkit_molecule.GetNumConformers()

    def atoms_data(self):
        atoms_data = []
        for atom in self.rdkit_molecule.GetAtoms():
            atom_data = self.atom_data(atom.GetIdx())
            atoms_data.append(atom_data)

        return atoms_data

    def atom_data(self, atom_idx):
        """Extracts useful rdkit information about an atom and returns it as a
        dictionary.

        """
        atom = self.atoms[atom_idx]
        atom_dict = {}
        atom_dict['atomic_num'] = atom.GetAtomicNum()
        atom_dict['bond_degree_no_Hs'] = atom.GetDegree()
        # same but want a convenience attribute
        atom_dict['bond_degree'] = atom.GetDegree()
        atom_dict['bond_degree_with_Hs'] = atom.GetTotalDegree()
        # same but want a convenience attribute
        atom_dict['total_bond_degree'] = atom.GetTotalDegree()
        atom_dict['explicit_valence'] = atom.GetExplicitValence()
        atom_dict['implicit_valence'] = atom.GetImplicitValence()
        atom_dict['total_valence'] = atom.GetTotalValence()
        atom_dict['formal_charge'] = atom.GetFormalCharge()
        atom_dict['hybridization'] = atom.GetHybridization()

        atom_dict['is_aromatic'] = atom.GetIsAromatic()
        atom_dict['in_ring'] = atom.IsInRing()
        atom_dict['isotope'] = atom.GetIsotope()
        atom_dict['mass'] = atom.GetMass()
        atom_dict['num_radical_electrons'] = atom.GetNumRadicalElectrons()
        atom_dict['element'] = atom.GetSymbol()
        atom_dict['num_Hs'] = atom.GetTotalNumHs()
        monomer_info = atom.GetMonomerInfo()
        if monomer_info:
            atom_dict['monomer_type'] = monomer_info.GetMonomerType()
            atom_dict['pdb_name'] = monomer_info.GetName().strip()
            # atom_dict['pdb_chain_id'] = monomer_info.GetChainID()
            atom_dict['pdb_insertion_code'] = monomer_info.GetInsertionCode()
            # atom_dict['pdb_heteroatom'] = monomer_info.IsHeteroAtom()
            atom_dict['pdb_occupancy'] = monomer_info.GetOccupancy()
            atom_dict['pdb_residue_name'] = monomer_info.GetResidueName()
            atom_dict['pdb_residue_number'] = monomer_info.GetResidueNumber()
            atom_dict['pdb_serial_number'] = monomer_info.GetSerialNumber()
            # atom_dict['pdb_segment_number'] = monomer_info.GetSegmentNumber()
            atom_dict['pdb_temp_factor'] = monomer_info.GetTempFactor()

        atom_dict['rdkit_mol_idx'] = atom.GetIdx()
        atom_dict['name'] = atom_dict['rdkit_mol_idx']

        return atom_dict

    def bonds_data(self):
        bonds_data = []
        for bond_idx, bond in enumerate(self.bonds):
            bond_data = self.bond_data(bond_idx)
            bonds_data.append(bond_data)

        return bonds_data

    def bond_data(self, bond_idx):
        """Extracts useful rdkit information about an atom and returns it as a
        dictionary.

        """
        bond = self.bonds[bond_idx]
        bond_dict = {}
        bond_dict['bond_type'] = str(bond.GetBondType())
        bond_dict['bond_type_number'] = str(bond.GetBondTypeAsDouble())
        bond_dict['is_aromatic'] = bond.GetIsAromatic()
        bond_dict['in_ring'] = bond.IsInRing()
        bond_dict['stereo'] = str(bond.GetStereo())
        bond_dict['is_conjugated'] = bond.GetIsConjugated()
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        bond_dict['rdkit_atom_idxs'] = (atom1_idx, atom2_idx)
        bond_dict['rdkit_mol_idx'] = bond.GetIdx()
        bond_dict['name'] = bond_dict['rdkit_mol_idx']

        return bond_dict

    def bonds_map(self):
        bond_map_dict = {}
        for bond_idx, bond in enumerate(self.bonds):
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            bond_map_dict[bond_idx] = (atom1_idx, atom2_idx)

        return bond_map_dict

    def molecule_data(self):
        """Extracts useful rdkit information about an atom and returns it as a
        dictionary.

        """
        molecule_dict = {}
        molecule_dict['name'] = self.mol_name
        ring_info = self.rdkit_molecule.GetRingInfo()
        molecule_dict['num_rings'] = ring_info.NumRings()
        molecule_dict['num_atoms'] = self.rdkit_molecule.GetNumAtoms()
        molecule_dict['num_bonds'] = self.rdkit_molecule.GetNumBonds()
        molecule_dict['num_heavy_atoms'] = self.rdkit_molecule.GetNumHeavyAtoms()

        return molecule_dict

    def make_atom_type(self, atom_idx, atom_type_name):
        atom_data = self.atom_data(atom_idx)
        return mastmol.AtomType.factory(atom_type_name, **atom_data)

    def make_atom_types(self):
        atoms_data = self.atoms_data()
        atom_types = []
        for atom_data in atoms_data:
            atom_type_name = "{1}Atom{0}Type".format(atom_data['name'], self.mol_name)
            atom_type = mastmol.AtomType.factory(atom_type_name, **atom_data)
            atom_types.append(atom_type)
        return atom_types

    def make_bond_type(self, bond_idx, bond_type_name):
        bond_data = self.bond_data(bond_idx)
        return mastmol.BondType.factory(bond_type_name,
                                        atom_types=bond_atom_types,
                                        **bond_data)

    def make_bond_types(self):
        # get atom types
        atom_types = self.make_atom_types()
        bonds_data = self.bonds_data()
        bond_types = []
        for bond_data in bonds_data:
            bond_type_name = "{1}Bond{0}Type".format(bond_data['name'], self.mol_name)
            bond_atom_types = (atom_types[bond_data['rdkit_atom_idxs'][0]],
                          atom_types[bond_data['rdkit_atom_idxs'][1]])
            bond_type = mastmol.BondType.factory(bond_type_name,
                                                 atom_types=bond_atom_types,
                                                 **bond_data)
            bond_types.append(bond_type)
        return bond_types

    def make_molecule_type(self, find_features=False):
        # get relevant data
        bond_map = self.bonds_map()
        molecule_data = self.molecule_data()
        if find_features:
            features = self.find_features()
        else:
            features = {}

        # AtomTypes
        atom_types = self.make_atom_types()
        # BondTypes
        bond_types = self.make_bond_types()
        # MoleculeType
        molecule_data.update({"name" : self.mol_name})

        molecule_type_name = "{}Type".format(self.mol_name)
        molecule_type = mastmol.MoleculeType.factory(molecule_type_name,
                                             atom_types=atom_types,
                                             bond_types=bond_types, bond_map=bond_map,
                                             features=features,
                                             **molecule_data)

        return molecule_type

    def find_features(self, fdef="BaseFeatures.fdef"):
        """Uses a feature definition (fdef) database to to find features in
        the molecule.

        Returns a nested dictionary mapping the feature indices to the
        feature dict.

        Examples
        --------

        """

        assert isinstance(fdef, str)
        fdef_path = osp.join(RDConfig.RDDataDir, fdef)
        feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
        factory_features = feature_factory.GetFeaturesForMol(self.rdkit_molecule)
        features = {}
        for feature in factory_features:
            # unpack the coordinates from the rdkit object
            pos_obj = feature.GetPos()
            feature_info = {'family' : feature.GetFamily(),
                            'type' : feature.GetType(),
                            'atom_ids' : feature.GetAtomIds(),
                            'position' : (pos_obj.x, pos_obj.y, pos_obj.z)}
            features[feature.GetId()] = feature_info


        families = col.defaultdict(list)
        types = col.defaultdict(list)
        for idx, info in features.items():
            families[info['family']].append(idx)
            types[info['type']].append(idx)

        return features

    def get_conformer_coords(self, conf_idx):
        assert self.rdkit_molecule.GetNumConformers() > 0, \
            "{0} has no conformers".format(self)

        conformer = self.rdkit_molecule.GetConformer(conf_idx)
        atom_idxs = range(self.rdkit_molecule.GetNumAtoms())
        # make the CoordArray
        coords = []
        for atom_idx in atom_idxs:
            coord = conformer.GetAtomPosition(atom_idx)
            coord = np.array([coord.x, coord.y, coord.z])
            coords.append(coord)
        coords = np.array(coords)
        return coords
