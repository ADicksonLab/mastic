import os.path as osp
import pickle
import itertools as it
import mdtraj as mdj

from mast.interactions.hydrogen_bond import HydrogenBondType
from mast.interactions.pi_stacking import PiStackingType

def profile_frame_inxs(frame_xyz, lig_idxs, receptor_idxs, system_type, inx_type, inx_classes):
    # get the coords for each member and convert to Angstroms
    lig_coords = frame_xyz[lig_idxs,:]
    receptor_coords = frame_xyz[receptor_idxs,:]
    member_coords = [lig_coords, receptor_coords]

    # substantiate the system with these coordinates
    system = system_type.to_system(member_coords)

    # profile the interactions
    inxs = []
    rec_lig_inxs = system.associations[0].profile_interactions(
        [inx_type],
        interaction_classes=inx_classes)[inx_type]
    inxs.extend(rec_lig_inxs)

    print(len(inxs))
    return inxs


if __name__ == "__main__":
    # load the system type pickle in
    system_pkl_path = osp.join(".", "Trypsin_Benzamidine_SystemType.pkl")
    with open(system_pkl_path, 'rb') as rf:
        Trypsin_Benzamidine_SystemType = pickle.load(rf)

    # use the association polynomial function of the system
    assoc_terms = Trypsin_Benzamidine_SystemType.association_polynomial(
        # input the degree of the interaction
        interaction_degree=HydrogenBondType.degree,
        # return the indices of the system members instead of the members
        # themselves
        return_idxs=True,
        # whether or not the interaction is symmetric or not
        commutative=False)

    # for pi-stacking we only want ligand-protein interactions, and it is
    # commutative so we only need one inter- term
    pistack_inx_classes = Trypsin_Benzamidine_SystemType.interaction_space(
        [(0,1)], PiStackingType)[(0,1)]

    # paths for the trajectory
    # get the coordinates from the trajectory
    traj_path = osp.join(".", "frames.dcd")
    top_path = osp.join(".", "topology.pdb")
    traj =mdj.load_dcd(traj_path, top=top_path)
    n_frames = traj.n_frames

    # indices for the ligand and receptor
    top_df = traj.top.to_dataframe()[0]
    lig_idxs = list(top_df[top_df['segmentID'] == 'SML'].index)
    receptor_idxs = list(top_df[top_df['segmentID'] == 'PROA'].index)

    # for input to mapping function
    # inputs to the mapping function
    # separate frames individually
    traj_xyzs = [frame.xyz[0,:,:] for frame in traj]
    lig_idxs_inputs = [lig_idxs for i in range(n_frames)]
    receptor_idxs_inputs = [receptor_idxs for i in range(n_frames)]
    system_type_inputs = [Trypsin_Benzamidine_SystemType for i in range(n_frames)]
    inx_type_inputs = [PiStackingType for i in range(n_frames)]
    inx_classes_inputs = [pistack_inx_classes for i in range(n_frames)]

    inputs = [traj_xyzs, lig_idxs_inputs, receptor_idxs_inputs,
              system_type_inputs, inx_type_inputs, inx_classes_inputs]

    print("Done with set up")
    # this is what actually runs the function on the inputs
    cluster_inxs = list(futures.map(profile_frame_inxs, *inputs))

    # now persist to disk
    inx_pickle_path = osp.join(".", "pi_cation_inxs.pkl")
    with open(inx_pickle_path, 'rb') as wf:
        pickle.dump(cluster_inxs, wf)
