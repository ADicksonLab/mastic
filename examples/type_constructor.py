import mast.selection as mastsel

def meta_type_factory(attr_list, class_name, bases):
    attr_dict = {attr : None for attr in attr_list}
    attr_dict['attributes'] = attr_list
    return type(class_name, bases, attr_dict)

mol_attrs = ['pdb_name', 'atom_types', ]
mol_class_name = "MoleculeType"

MoleculeType = meta_type_factory(mol_attrs, mol_class_name, (mastsel.SelectionMember,))

def molecule_type_factory(attributes, mol_type_name):

    for attr in MoleculeType.attributes:
        assert attr in attributes.keys()
    return type(mol_type_name, (MoleculeType,), attributes)

water_attrs = {'atom_types' : ['H', 'O', 'H'],
               'pdb_name' : 'water'}

WaterType = molecule_type_factory(water_attrs, 'WaterType')

class Molecule(object):
    def __init__(self, molecule_type, coords):
        self.molecule_type = molecule_type
        self.coords = coords

water_mol = Molecule(WaterType, [[0,0,1],[0,1,0],[1,0,0]])

