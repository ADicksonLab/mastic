INTERACTION_ATTRIBUTES = ['name', ]

# families and types from the RDKit feature definition database
# (.fdef)

# Feature dictionaries have as their keys the grouping attribute.  The
# list has as many elements as the degree of the interaction.  Each
# element is a list of keywords that can be used to identify the
# feature from rdkit.
# the grouping attributes assign which key to use in the INX_FEATURES
# dictionary when finding appropriate features for each member of the
# interaction

# HydrogenBondType : hydrogen bonds with explicit hydrogens
HYDROGEN_BOND_DEGREE = 2
HYDROGEN_BOND_COMMUTATIVITY = False
HYDROGEN_BOND_PARAM_KEYS = ['distance', 'angle']
HYDROGEN_BOND_FEATURE_KEYS = ['donor', 'acceptor']
HYDROGEN_BOND_RDKIT_FEATURES = {HYDROGEN_BOND_FEATURE_KEYS[0] : ['Donor'],
                        HYDROGEN_BOND_FEATURE_KEYS[1] : ['Acceptor']}
HYDROGEN_BOND_FEATURES = HYDROGEN_BOND_RDKIT_FEATURES


# HalogenBondType : strong interactions with halogen acceptors and
# explicit hydrogens
HALOGEN_BOND_DEGREE = 2
HALOGEN_BOND_COMMUTATIVITY = False
HALOGEN_BOND_PARAM_KEYS = ['distance', 'donor_angle', 'acceptor_angle']
HALOGEN_BOND_FEATURE_KEYS = ['donor', 'halogen_acceptor']
HALOGEN_BOND_GROUPING_ATTRIBUTES = ['rdkit_family', 'rdkit_type']
HALOGEN_BOND_FEATURES = {'rdkit_family' : [['Donor'], []],
                         'rdkit_type' : [[], ['HalogenAcceptor']]}

# Halogen bond -Y-O ~ X-C-
# where Y = C, P, S and X = F, Cl, Br, I
# parameters taken from (Halogen bonds in biological molecules., Auffinger et al., 2004)
# Maximum distance between oxygen and halogen +0.5
HALOGEN_BOND_DIST_MAX = 4.0 # /AA
# Optimal donor angle
HALOGEN_BOND_DON_ANGLE = 165 # the X-C angle in /degrees
# Optimal acceptor angle
HALOGEN_BOND_ACC_ANGLE = 120 # the Y-O angle in /degrees
# Maximum deviation from optimal angle
HALOGEN_BOND_ANGLE_DEV = 30 # /degrees

# PiStackingType
PISTACKING_DEGREE = 2
PISTACKING_COMMUTATIVITY = True
PISTACKING_PARAM_KEYS = ['centroid_distance', 'ring_normal_angle',
                         'centroid_offset', 'stacking_type']
PISTACKING_FEATURE_KEYS = ['arom_1', 'arom_2']
PISTACKING_RDKIT_FEATURES = {PISTACKING_FEATURE_KEYS[0] : ['Aromatic', 'Arom6'],
                             PISTACKING_FEATURE_KEYS[1] : ['Aromatic', 'Arom6']}
PISTACKING_FEATURES = PISTACKING_RDKIT_FEATURES

# Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACKING_CENTROID_DIST_MAX = 7.5 # /AA
# Max. Deviation from parallel or perpendicular orientation
PISTACKING_ANGLE_DEVIATION = 30 # /degrees
# Maximum offset of the two rings (corresponds to the radius of
# benzene + 0.5 A)
PISTACKING_OFFSET_MAX = 2.0 # /AA

# Maximum distance the closest atom must be in a T-stack
PISTACKING_T_DIST = 5.0 #/AA


HYDROPH_DEGREE = 2
HYDROPH_GROUPING_ATTRIBUTES = ['rdkit_family', 'rdkit_family']
HYDROPH_FEATURES = {'rdkit_family' : [['Hydrophobe', 'LumpedHydrophobe'],
                                      ['Hydrophobe', 'LumpedHydrophobe']],
                    'rdkit_type' : [[], []]}

PICATION_DEGREE = 2
PICATION_GROUPING_ATTRIBUTES = ['rdkit_family', 'rdkit_family']
PICATION_FEATURES = {'rdkit_family' : [['Aromatic'], ['PosIonizable']],
                     'rdkit_type' : [[], []]}

SALTBRIDGE_DEGREE = 2
SALTBRIDGE_GROUPING_ATTRIBUTES = ['rdkit_family', 'rdkit_family']
SALTBRIDGE_FEATURES = {'rdkit_family' : [['PosIonizable'], ['NegIonizable']],
                       'rdkit_type' : [[], []]}


WATER_BRIDGE_FEATURE_FAMILIES = []
WATER_BRIDGE_FEATURE_TYPES = []

METAL_FEATURE_FAMILIES = []
METAL_FEATURE_TYPES = []

# Hydrogen Bond -D-H ~ A-
# parameters taken from (Hubbard & Haider, 2001)
# Max. distance between hydrogen bond donor and acceptor  + 0.6 A
HYDROGEN_BOND_DIST_MAX = 4.1 # /AA
# Min. angle at the hydrogen bond donor + 10
HYDROGEN_BOND_DON_ANGLE_MIN = 100 # /AA




# Determines allowed deviation from planarity in aromatic rings
AROMATIC_PLANARITY = 5.0 # /AA


# Some distance thresholds were extended (max. 1.0A) if too
# restrictive too account for low-quality structures Distance cutoff
# for detection of hydrophobic contacts
HYDROPH_DIST_MAX = 4.0 # /AA

# Max. distance between charged atom and aromatic ring center
# (Gallivan and Dougherty, 1999)
PICATION_DIST_MAX = 6.0 # /AA

# Max. distance between centers of charge for salt bridges (Barlow and
# Thornton, 1983) + 1.5
SALTBRIDGE_DIST_MAX = 5.5 # /AA


# Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
WATER_BRIDGE_MINDIST = 2.5 # /AA
# Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.4
WATER_BRIDGE_MAXDIST = 4.0 # /AA
# Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005) - 5
WATER_BRIDGE_OMEGA_MIN = 75 # /degrees
# Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
WATER_BRIDGE_OMEGA_MAX = 140 # /degrees
# Min. angle between water oxygen, donor hydrogen and donor atom (Jiang et al., 2005)
WATER_BRIDGE_THETA_MIN = 100 # /degrees

# Max. distance between metal ion and interacting atom (Harding, 2001)
METAL_DIST_MAX = 3.0 # /AA
