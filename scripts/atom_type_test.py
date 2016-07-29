from copy import copy

atom_attrs = ['name',
              'atomic_num',
              'bond_degree_no_Hs',
              'bond_degree_with_Hs',
              'total_bond_degree',
              'explicit_valence',
              'implicit_valence',
              'total_valence',
              'formal_charge',
              'hybridization',
              'is_aromatic',
              'in_ring',
              'isotope',
              'mass',
              'num_radical_electrons',
              'element',
              'num_Hs',
              'monomer_type',
              'pdb_name',
              'pdb_insertion_code',
              'pdb_occupancy'
              'pdb_residue_name',
              'pdb_residue_number',
              'pdb_serial_number',
              'pdb_temp_factor',]

none_atom_attr_dict = {attr : None for attr in atom_attrs}

AtomType = type('AtomType', (object,), none_atom_attr_dict)
AtomType.attributes = atom_attrs

def atom_type_factory(atom_attrs, atom_type_name):
    attributes = {attr : None for attr in AtomType.attributes}

    # simply keep track of which attributes the input did not provide
    for attr in AtomType.attributes:
        try:
            assert attr in atom_attrs.keys()
        except AssertionError:
            # logging
            print("Attribute {0} not found in atom input.".format(attr))

    # add the attributes that it has into the class
    for attr, value in atom_attrs.items():
        # make sure AtomType has an attribute that matches
        try:
            assert attr in AtomType.attributes
        # if it doesn't then report it
        except AssertionError:
            # logging
            print("Input attribute {0} not in AtomType attributes, ignoring.".format(attr))
        # if no error then add it
        else:
            attributes[attr] = value

    return type(atom_type_name, (AtomType,), attributes)

O_attrs = {'atomic_num' : 16, 'element' : 'O'}
OAtomType = atom_type_factory(O_attrs, 'OAtomType')

H_attrs = {'atomic_num' : 1, 'element' : 'H'}
HAtomType = atom_type_factory(H_attrs, 'HAtomType')

C_attrs = {'atomic_num' : 12, 'element' : 'C'}
CAtomType = atom_type_factory(C_attrs, 'CAtomType')

class MoleculeType(type):

    # attributes determined during factory creation
    _atom_type_library = set()
    _features = None
    _feature_families = None
    _feature_types = None

    _bonds = None

    name = None
    atom_types = None

    @classmethod
    def atom_type_library(self):
        return list(self._atom_type_library)

    @classmethod
    def features(self):
        return self._features

    @classmethod
    def feature_families(self):
        return self._feature_families

    @classmethod
    def feature_types(self):
        return self._feature_types

    @classmethod
    def atom_types(self):
        return self._atom_type_sequence

    @classmethod
    def bonds(self):
        return self._bonds

    @classmethod
    def to_molecule(self, coords=None):
        """ Construct a Molecule using input coordinates with mapped indices"""

        coord_array = CoordArray(coords)

        # Make atoms out of the coord array
        self.make_atom_type_library()
        atom_idxs = range(self.molecule.GetNumAtoms())
        atoms = []
        for atom_idx in atom_idxs:
            atom_type = self.atom_type_library[atom_idx]
            atom = Atom(atom_array=coord_array, array_idx=atom_idx, atom_type=atom_type)
            atoms.append(atom)

        # handle bonds
        bonds = []
        for bond in self.molecule.GetBonds():
            bonds.append(Bond(atoms, (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))

        # TODO handle and create angles
        angles = None

        return Molecule(atoms, bonds, angles, mol_type=self)

def molecule_type_factory(mol_type_name, name=None, atom_types=None):
    mol_type = type(mol_type_name, (MoleculeType,), {})
    mol_type.name = name
    mol_type.atom_types = atom_types
    mol_type._atom_type_library = set(mol_type.atom_types)
    return mol_type


water_attrs = {'atom_types' : [HAtomType, OAtomType, HAtomType],
               'name' : 'water'}

methanol_attrs = {'atom_types' : [HAtomType, OAtomType, CAtomType,
                                  HAtomType, HAtomType, HAtomType],
                  'name' : 'methanol'}

WaterType = molecule_type_factory('WaterType', **water_attrs)
MethanolType = molecule_type_factory('MethanolType', **methanol_attrs)
