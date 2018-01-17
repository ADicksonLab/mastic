import numpy as np

from rdkit import Chem
import os.path as osp

from mastic.selection import CoordArray, IndexedSelection, Selection, \
    GenericSelection

from mastic.molecule import Bond, Molecule, Atom, \
    MoleculeType, AtomType, BondType
from mastic.interfaces.rdkit import RDKitMoleculeWrapper
import mastic.config.molecule as masticmolconfig


### Making AtomTypes, BondTypes, and MoleculeTypes
print("making AtomTypes")
atom1_attrs = {'pdb_name' : "FAKE1"}
Atom1Type = AtomType("Atom1Type", **atom1_attrs)

atom2_attrs = {'pdb_name' : "FAKE2"}
Atom2Type = AtomType("Atom2Type", **atom2_attrs)

atom3_attrs = {'pdb_name' : "FAKE3"}
Atom3Type = AtomType("Atom3Type", **atom3_attrs)

print("making BondType")
bond1_attrs = {'bond_type' : "TRIPLE"}
Bond1Type = BondType("Bond1Type", atom_types=(Atom1Type, Atom2Type), **bond1_attrs)

print("making MoleculeType")
molecule1_attrs = {'name' : "FAKE"}
Molecule1Type = MoleculeType("Molecule1Type",
                                     atom_types=[Atom1Type, Atom2Type],
                                     bond_types=[Bond1Type], bond_map = {0 : (0,1)},
                                     **molecule1_attrs)

### From an external representation
print("Read in an external representation")
tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
pka_rdkit = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
# cast external representation to a wrapper
pka_rdkit_wrapper = RDKitMoleculeWrapper(pka_rdkit, mol_name="PKA")
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
    atom_type = AtomType(atom_type_name, **atom_data)
    pka_atom_types.append(atom_type)

# BondTypes
pka_bond_types = []
for bond_data in pka_bond_data:
    bond_type_name = "PKABond{0}Type".format(bond_data['name'])
    atom_types = (pka_atom_types[bond_data['rdkit_atom_idxs'][0]],
                  pka_atom_types[bond_data['rdkit_atom_idxs'][1]])
    bond_type = BondType(bond_type_name, atom_types=atom_types, **bond_data)
    pka_bond_types.append(bond_type)

# MoleculeType
PKA1Type = MoleculeType("PKAType", atom_types=pka_atom_types,
                               bond_types=pka_bond_types, bond_map=pka_bond_map,
                               **pka_molecule_data)
# using the rdkit wrapper converter
PKA2Type = pka_rdkit_wrapper.make_molecule_type()

# find features
PKA3Type = pka_rdkit_wrapper.make_molecule_type(find_features=True)


### making actual Atoms, Bonds, and Molecules

# the CoordArray
array = np.array([[0,0,0], [0,0,1], [1,0,0]])
atom_array = CoordArray(array)
print(atom_array)

# the Atoms
atom1 = Atom(np.array([5,5,5]), atom_type=Atom1Type)
print(atom1)
print(atom1.coords)

atom2 = Atom(np.array([6,6,6]), atom_array=atom_array, atom_type=Atom2Type)
print(atom2)
print(atom2.coords)
print("testing overlap of two atoms")
print(atom2.overlaps(atom1))
atom3 = Atom(atom_array=atom_array, array_idx=0, atom_type=Atom3Type)
print(atom3)
print(atom3.coords)
atoms = [atom1, atom2, Atom(np.array([0,1,0]), atom_type=Atom2Type)]

# selection of atoms is fun
atomsel = IndexedSelection(atoms, [0,1])

# Bonds
bond = Bond(atoms, atom_ids=(0,1), bond_type=Bond1Type)
idx_a = range(len(atoms))[:-1]
idx_b = [a+1 for a in idx_a]
idx = zip(idx_a, idx_b)
bonds = [Bond(atoms, bond_idx, bond_type=Bond1Type) for bond_idx in idx]

# getting the atoms out of bonds
print("accessing bonds from an atom")
print("Is an atom in a bond?")
print(atoms[0].isin_bond)
print("how many bonds is it in")
print("first atom", len(atoms[0].bonds))
print("second atom", len(atoms[1].bonds))
print("get the bonds themselves")
print(atoms[0].bonds)
print("get the other atom in the bond")
bond = atoms[0].bonds[0]
other_atom = next((a for a in bond.atoms if a is not atoms[0]))

# using the class's method
print(atoms[0].adjacent_atoms)

# Molecule
mol = Molecule(atoms, bonds, mol_type=Molecule1Type)
print(mol)
print("atom_types in mol")
print(mol.atom_types)


print("Making a mastic.Molecule from the RDKitMoleculeWrapper data")
pka_mol = PKA3Type.to_molecule(pka_coords)
# pka_mol = Molecule(mol_type=pka_type, coords=pka_coords)
print(pka_mol)
print(pka_mol.molecule_type)

print("testing overlap of two molecules")
print(pka_mol.overlaps(mol))

pka_atom_type = PKA3Type.atom_types[0]
pka_atom = pka_mol.atoms[0]
pka_bond = pka_mol.bonds[0]
pka_feature = pka_mol.features[1]
