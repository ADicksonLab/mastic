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


def _atom_type_factory(atom_attrs, atom_type_name):
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

AtomType = type('AtomType', (object,), {attr : None for attr in atom_attrs})
AtomType.attributes = atom_attrs
AtomType.factory = _atom_type_factory


O_attrs = {'atomic_num' : 16, 'element' : 'O'}
O_ATOM_TYPE = AtomType.factory(O_attrs, 'OAtomType')

H_attrs = {'atomic_num' : 1, 'element' : 'H'}
H_ATOM_TYPE = AtomType.factory(H_attrs, 'HAtomType')

C_attrs = {'atomic_num' : 12, 'element' : 'C'}
C_ATOM_TYPE = AtomType.factory(C_attrs, 'CAtomType')

class MoleculeType(object):

    def __init__(self, name=None, atom_types=None, ):
        self._atom_type_library = set(atom_types)
        self._features = None
        self._feature_families = None
        self._feature_types = None
        self._bonds = None

        self.name = name
        self._atom_type_sequence = atom_types


    @property
    def atom_type_library(self):
        return list(self._atom_type_library)

    @property
    def features(self):
        return self._features

    @property
    def feature_families(self):
        return self._feature_families

    @property
    def feature_types(self):
        return self._feature_types

    @property
    def atom_types(self):
        return self._atom_type_sequence

    @property
    def bonds(self):
        return self._bonds

    @property
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

    @classmethod
    def factory(cls, mol_type_name, name=None, atom_types=None):
        mol_class = type(mol_type_name, (cls,), {})
        mol_type = mol_class(name=name, atom_types=atom_types)

        return mol_type


water_attrs = {'atom_types' : [H_ATOM_TYPE, O_ATOM_TYPE, H_ATOM_TYPE],
               'name' : 'water'}

methanol_attrs = {'atom_types' : [H_ATOM_TYPE, O_ATOM_TYPE, C_ATOM_TYPE,
                                  H_ATOM_TYPE, H_ATOM_TYPE, H_ATOM_TYPE],
                  'name' : 'methanol'}

WATER_TYPE = MoleculeType.factory('WaterType', **water_attrs)
METHANOL_TYPE = MoleculeType.factory('MethanolType', **methanol_attrs)
