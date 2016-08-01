import numpy as np

from rdkit import Chem
import os.path as osp

from mast.selection import CoordArray, IndexedSelection, Selection, \
    GenericSelection

from mast.molecule import Bond, Molecule, Atom, \
    MoleculeType, AtomType, BondType, RDKitMoleculeWrapper

import mast.config.molecule as mastmolconfig


### Making AtomTypes, BondTypes, and MoleculeTypes
print("making AtomTypes")
atom1_attrs = {'pdb_name' : "FAKE1"}
Atom1Type = AtomType.factory("Atom1Type", **atom1_attrs)

atom2_attrs = {'pdb_name' : "FAKE2"}
Atom2Type = AtomType.factory("Atom2Type", **atom2_attrs)

atom3_attrs = {'pdb_name' : "FAKE3"}
Atom3Type = AtomType.factory("Atom3Type", **atom3_attrs)

print("making BondType")
bond1_attrs = {'bond_type' : "TRIPLE"}
Bond1Type = BondType.factory("Bond1Type", atom_types=(Atom1Type, Atom2Type), **bond1_attrs)

print("making MoleculeType")
molecule1_attrs = {'name' : "FAKE"}
Molecule1Type = MoleculeType.factory("Molecule1Type",
                                     atom_types=[Atom1Type, Atom2Type],
                                     bond_types=[Bond1Type], bond_map = {0 : (0,1)},
                                     **molecule1_attrs)

### From an external representation
print("Read in an external representation")
tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
pka_rdkit = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
# cast external representation to a wrapper
pka_rdkit_wrapper = RDKitMoleculeWrapper(pka_rdkit)
print(pka_rdkit_wrapper)

# extract data from it
# atom type data
pka_atom_data = pka_rdkit_wrapper.atoms_data()
print(pka_atom_data)
# features
pka_features = pka_rdkit_wrapper.find_features()
print(pka_features)
# bond types
pka_bond_data = pka_rdkit_wrapper.bonds_data()
print(pka_bond_data)
# bond map
pka_bond_map = pka_rdkit_wrapper.bonds_map()

# molecule data
pka_molecule_data = pka_rdkit_wrapper.molecule_data()
print(pka_molecule_data)
# get the coordinates from a conformer
pka_coords = pka_rdkit_wrapper.get_conformer_coords(0)

# create types from data sources
# AtomTypes
pka_atom_types = []
for atom_data in pka_atom_data:
    atom_type_name = "PKAAtom{0}Type".format(atom_data['name'])
    atom_type = AtomType.factory(atom_type_name, **atom_data)
    pka_atom_types.append(atom_type)

# BondTypes
pka_bond_types = []
for bond_data in pka_bond_data:
    bond_type_name = "PKABond{0}Type".format(bond_data['name'])
    atom_types = (pka_atom_types[bond_data['rdkit_atom_idxs'][0]],
                  pka_atom_types[bond_data['rdkit_atom_idxs'][1]])
    bond_type = BondType.factory(bond_type_name, atom_types=atom_types, **bond_data)
    pka_bond_types.append(bond_type)

# MoleculeType
pka_molecule_data.update({"name" : "PKA"})

PKA1Type = MoleculeType.factory("PKAType", atom_types=pka_atom_types,
                               bond_types=pka_bond_types, bond_map=pka_bond_map,
                               **pka_molecule_data)
# using the rdkit wrapper converter
PKA2Type = pka_rdkit_wrapper.molecule_type()

# find features
PKA3Type = pka_rdkit_wrapper.molecule_type(find_features=True)


### making actual Atoms, Bonds, and Molecules

print("making an CoordArray for atoms")
array = np.array([[0,0,0], [0,0,1], [1,0,0]])
atom_array = CoordArray(array)
print(atom_array)

print("making atoms")
atom1 = Atom(np.array([5,5,5]))
print(atom1)
print(atom1.coords)

atom2 = Atom(np.array([6,6,6]), atom_array=atom_array, atom_type=atom2_type)
print(atom2)
print(atom2.coords)
print("testing overlap of two atoms")
print(atom2.overlaps(atom1))
atom3 = Atom(atom_array=atom_array, array_idx=0)
print(atom3)
print(atom3.coords)
atoms = [atom1, atom2, Atom(np.array([0,1,0]))]
# # make a selection of those atoms
atomsel = IndexedSelection(atoms, [0,1])
print(atomsel)

# make a selection of atoms for bonds, and angle

# making a bond between each atom and the next one in the index,
# just to test
print("making up bonds")
idx_a = range(len(atoms))[:-1]
idx_b = [a+1 for a in idx_a]
idx = zip(idx_a, idx_b)
bonds = [Bond(atoms, bond_idx) for bond_idx in idx]

print("accessing bonds from an atom")
print("Is an atom in a bond?")
print(atoms[0]._in_bond)
print("how many bonds is it in")
print("first atom", len(atoms[0].bonds))
print("second atom", len(atoms[1].bonds))
print("get the bonds themselves")
print(atoms[0].bonds)
print("get the other atom in the bond")
bond = atoms[0].bonds[0]
other_atom = next((a for a in bond.atoms if a is not atoms[0]))

print("same thing using the method")
print(atoms[0].adjacent_atoms)
# angles still stubbed out for now
angles = [IndexedSelection(atoms, [0,1,2])]
print("making a MoleculeType")
moltype = MoleculeType.factory(mol_type_name="TestMolType", name='test_molecule')
print("making a molecule")
mol = Molecule(atoms, bonds, angles, mol_type=moltype)
print(mol)
print("atom_types in mol")
print(mol.atom_types)


print("Making a mast.Molecule from the RDKitMoleculeWrapper data")
pka_mol = PKAType().to_molecule(pka_coords)
# pka_mol = Molecule(mol_type=pka_type, coords=pka_coords)
print(pka_mol)
print(pka_mol.molecule_type)

print("testing overlap of two molecules")
print(pka_mol.overlaps(mol))


