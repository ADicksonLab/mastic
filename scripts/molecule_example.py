import numpy as np

from mast.selection import CoordArray, IndexedSelection, Selection, \
    GenericSelection

from mast.molecule import Bond, Molecule, Atom, \
    MoleculeType, AtomType, RDKitMoleculeType, AtomTypeLibrary, MoleculeTypeLibrary

print("making an CoordArray for atoms")
array = np.array([[0,0,0], [0,0,1], [1,0,0]])
atom_array = CoordArray(array)
print(atom_array)

print("making atoms")
atom1 = Atom(np.array([5,5,5]))
print(atom1)
print(atom1.coords)
print("making AtomType")
atom2_type = AtomType({'element':'C'})

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


from rdkit import Chem
import os.path as osp
tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
pka = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
print(pka)


# PKA = RDKitMoleculeType.create_molecule_type(pka, mol_type='PKA')

# PKA = RDKitMoleculeType.create_molecule_type(
pka_type = RDKitMoleculeType(pka, mol_name="PKA")
print(pka_type)

print("getting positions objects from rdchem.Mol for atoms")
conf = list(pka_type.molecule.GetConformers()).pop()
positions = {idx : conf.GetAtomPosition(idx) for idx in range(conf.GetNumAtoms())}
print(positions[0])
print("Making an AtomTypeLibrary and Atom list for pka")
atom_types = []
atom_names = {}
atoms = []
pka_atom_type_library = AtomTypeLibrary()
for atom_idx in range(len(pka_type.atoms)):

    atom_type = AtomType(pka_type.atom_data(atom_idx))
    atom_types.append(atom_type)

    atom_name = atom_type.pdb_name
    if atom_name in atom_names.keys() and \
       not pka_atom_type_library.attributes_match(atom_type):
        atom_names[atom_name] += 1
    elif atom_name not in atom_names.keys():
        atom_names[atom_name] = 0

    if atom_names[atom_name] > 0:
        pka_atom_type_library.add_type(atom_type, atom_name + str(atom_names[atom_name]) )
    else:
        pka_atom_type_library.add_type(atom_type, atom_name)

    coord = np.array([positions[atom_idx].x, positions[atom_idx].y, positions[atom_idx].z])
    atoms.append(Atom(coords=coord, atom_type=atom_type))

print(pka_atom_type_library)
print(atoms)

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
moltype = MoleculeType({'name':'test_molecule'})
print("making a molecule")
mol = Molecule(atoms, bonds, angles, mol_type=moltype)
print(mol)
print("atom_types in mol")
print(mol.atom_types)

print("Making a mast.Molecule from the RDKitMoleculeType")
pka_mol = pka_type.to_molecule(0)
# pka_mol = Molecule(mol_type=pka_type, coords=pka_coords)
print(pka_mol)
print(pka_mol.molecule_type)

print("testing overlap of two molecules")
print(pka_mol.overlaps(mol))

print("finding features using mast.molecule method")
pka_type.find_features()
pka_mol.make_feature_selections()
print(pka_mol.features)

print("using type_constructor")
pka_mol2 = pka_type.to_molecule_from_coords(pka_mol.atom_coords)
