for assoc_term, assoc_profile in system_profile.association_profiles.items():
    for hit_idx in assoc_profile.hits:
        assoc_idx = sEH_TPPU_SystemType.assoc_member_idxs.index(assoc_term)
        inx_class = sEH_TPPU_SystemType.association_types[assoc_idx].interaction_subspace[hit_idx]
        print("interaction class index {0}".format(hit_idx))
        for feature_type in inx_class.feature_types:
            pdb_serials = [atom_type.pdb_serial_number for atom_type in
                           feature_type.atom_types]
            print("    Feature {0}, {1}: {2}".format(feature_type.name,
                                                 feature_type.rdkit_family,
                                                 " ".join([str(i) for i in pdb_serials])))
