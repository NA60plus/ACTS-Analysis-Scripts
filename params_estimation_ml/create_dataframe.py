import uproot
import pandas as pd
import math
import numpy as np


# Function to find the index of a sublist by id
def find_particle_index(particles, target_id):
    index_list = [] 
    for index, particle in enumerate(particles):
        if particle[0] == target_id:
            index_list.append(index)
    return index_list  # Return -1 if the id is not found


def create_training_set():
    # Open the ROOT file
    file = uproot.open("../output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6_realTarget_beam0.5_twosteps_rej0.114_perigeeZ400_suffix/hits_ms.root")

    # Access the TTree (adjust the tree name as needed)
    tree = file["hits"]

    # Convert the tree to a pandas DataFrame
    # You can specify the branches you want to extract or leave it empty to get all branches
    df = tree.arrays(library="pd")


    cols_to_drop = ['event_id', 'geometry_id', 'tt', 'tpx',
        'tpy', 'tpz', 'te', 'deltapx', 'deltapy', 'deltapz', 'deltae', 'index',
        'volume_id', 'boundary_id', 'layer_id', 'approach_id', 'sensitive_id']

    df = df.drop(columns=cols_to_drop)


    grouped = df.groupby("particle_id")

    combinations = [
                        [0,1,2],
                        [0,1,3],
                        [0,1,4],
                        [0,1,5],
                        [0,2,3],
                        [0,2,4],
                        [0,2,5],
                        [0,3,4],
                        [0,3,5],
                        [0,4,5],
                        [1,2,3],
                        [1,2,4],
                        [1,2,5],
                        [1,3,4],
                        [1,3,5],
                        [1,4,5],
                        [2,3,4],
                        [2,3,5],
                        [2,4,5],
                        [3,4,5]
                    ]

    sp_triplets = []
    for particle_id, group in grouped:
        # Loop over each row within the group
        trip = []
        trip.append(particle_id)
        for index, row in group.iterrows():
            trip.append(row["tx"])
            trip.append(row["ty"])
            trip.append(row["tz"])
        for comb in combinations:
            if 3+3*comb[2] > len(trip) -1 or 3+3*comb[1] > len(trip) -1:
                continue
            sp_triplets.append([trip[0],
                                
                                trip[1+3*comb[0]],
                                trip[2+3*comb[0]],
                                trip[3+3*comb[0]],

                                trip[1+3*comb[1]],
                                trip[2+3*comb[1]],
                                trip[3+3*comb[1]],

                                trip[1+3*comb[2]],
                                trip[2+3*comb[2]],
                                trip[3+3*comb[2]]])    


    """
    PARTICLE ID
    """


    # Open the ROOT file
    file = uproot.open("../output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6_realTarget_beam0.5_twosteps_rej0.114_perigeeZ400_suffix/particles_simulation_ms.root")

    # Access the TTree (adjust the tree name as needed)
    tree = file["particles"]

    # Convert the tree to a pandas DataFrame
    # You can specify the branches you want to extract or leave it empty to get all branches
    df_part = tree.arrays(library="pd")

    cols_to_drop = ['event_id', 'particle_type', 'process', 'pt',
        'vertex_primary', 'vertex_secondary', 'particle', 'generation',
        'sub_particle', 'e_loss', 'total_x0', 'total_l0', 'number_of_hits',
        'outcome']
    df_part = df_part.drop(columns=cols_to_drop)

    # Loop over each row in the DataFrame
    # Loop over the events

    for index, row in df_part.iterrows():
        num_entries = len(row['particle_id'])  # Assuming each column has the same number of entries

        for i in range(num_entries):
            d0 = math.sqrt(row["vx"][i]**2+row["vx"][i]**2)
            z0 = row["vz"][i]
            theta = np.arctan2(math.sqrt(row["px"][i]**2 + row["py"][i]**2), row["pz"][i])
            phi = row["phi"][i]
            qop = row["q"][i]/row["p"][i]
            indexes = find_particle_index(sp_triplets,row["particle_id"][i])
            if indexes:
                for index in indexes:
                    sp_triplets[index].append(d0)
                    sp_triplets[index].append(z0)
                    sp_triplets[index].append(theta)
                    sp_triplets[index].append(phi)
                    sp_triplets[index].append(qop)

    return sp_triplets


# Python's entry point: if the script is run (not imported), main() will be executed
if __name__ == "__main__":
    list = create_training_set()