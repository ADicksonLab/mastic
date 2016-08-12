INTERACTION_ATTRIBUTES = ['name']

# families and types from the RDKit feature definition database
# (.fdef)
HBOND_FEATURES = {'rdkit_family' : ['Donor', 'Acceptor'],
                  'rdkit_type' : []}

PISTACK_FEATURE_FAMILIES = ['Aromatic']
PISTACK_FEATURE_TYPES = []

HYDROPH_FEATURE_FAMILIES = ['Hydrophobe', 'LumpedHydrophobe']
HYDROPH_FEATURE_TYPES = []

PICATION_FEATURE_FAMILIES = ['Aromatic', 'PosIonizable']
PICATION_FEATURE_TYPES = []

SALTBRIDGE_FEATURE_FAMILIES = ['PosIonizable', 'NegIonizable']
SALTBRIDGE_FEATURE_TYPES = []

HALOGEN_FEATURE_FAMILIES = []
HALOGEN_FEATURE_TYPES = ['HalogenAcceptor']

WATER_BRIDGE_FEATURE_FAMILIES = []
WATER_BRIDGE_FEATURE_TYPES = []

METAL_FEATURE_FAMILIES = []
METAL_FEATURE_TYPES = []

# Determines allowed deviation from planarity in aromatic rings
AROMATIC_PLANARITY = 5.0 # /AA

# Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
HBOND_DIST_MAX = 4.1 # /AA
# Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
HBOND_DON_ANGLE_MIN = 100 # /AA


# Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACK_DIST_MAX = 7.5 # /AA
# Max. Deviation from parallel or perpendicular orientation (in degrees)
PISTACK_ANG_DEV = 30 # /degrees
# Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
PISTACK_OFFSET_MAX = 2.0 # /AA

# Some distance thresholds were extended (max. 1.0A) if too restrictive too account for low-quality structures
# Distance cutoff for detection of hydrophobic contacts
HYDROPH_DIST_MAX = 4.0 # /AA

# Max. distance between charged atom and aromatic ring center (Gallivan and Dougherty, 1999)
PICATION_DIST_MAX = 6.0 # /AA

# Max. distance between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.5
SALTBRIDGE_DIST_MAX = 5.5 # /AA

# Max. distance between oxy. and halogen (Halogen bonds in biological molecules., Auffinger)+0.5
HALOGEN_DIST_MAX = 4.0 # /AA
# Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_ACC_ANGLE = 120 # /degrees
# Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_DON_ANGLE = 165 # /degrees
# Max. deviation from optimal angle
HALOGEN_ANGLE_DEV = 30 # /degrees

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
