from mast.tests.data import Trypsin_Benzamidine_System_cryst, Trypsin_Benzamidine_SystemType
from mast.interactions.hydrogen_bond import HydrogenBondType

# use the association polynomial function of the system
assoc_terms = Trypsin_Benzamidine_SystemType.association_polynomial(
    # input the degree of the interaction
    interaction_degree=HydrogenBondType.degree,
    # return the indices of the system members instead of the members
    # themselves
    return_idxs=True,
    # whether or not the interaction is symmetric or not
    commutative=False)

hbond_inx_classes = Trypsin_Benzamidine_SystemType.interaction_space(
    assoc_terms, HydrogenBondType)

tryp_ben_prof_results = new_system.associations[0].profile_interactions([HydrogenBondType],
                                                      interaction_classes=hbond_inx_classes)
