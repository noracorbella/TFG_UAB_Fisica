#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 11:17:25 2025

@author: nora
"""

import MDAnalysis as mda
from MDAnalysis.analysis import msd
import numpy as np
from scipy.spatial.distance import pdist, squareform

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/DB8_traj_unwrapped.xtc"


topology_pdb = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/DB8_frame0.pdb"


output_file_DB = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/MDAnalysis_MSD_DB8_filtered.dat"
output_file_DB_original = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/MDAnalysis_MSD_DB8.dat"

# Distance threshold for considering beads as stuck together (in pixels)
DISTANCE_THRESHOLD = 20.0

# Load the universe
u = mda.Universe(topology_pdb, trajectory_file)

# Select Dynabeads
DBs = u.select_atoms('name cell')

def compute_and_save_msd(atomgroup, output_file, label):
    msd_calc = msd.EinsteinMSD(atomgroup, msd_type='xy', fft=True, apply_pbc=True)
    msd_calc.run()

    nframes = msd_calc.n_frames
    timestep = 100  # One frame every 100 MCS
    lagtimes = np.arange(nframes) * timestep
    msd_values = msd_calc.results.timeseries

    with open(output_file, 'w') as f:
        f.write(f"# MSD of {label}\n")
        f.write("LagTime(MCS)\tMSD\n")
        for t, msd_val in zip(lagtimes, msd_values):
            f.write(f"{t}\t{msd_val}\n")

def find_non_stuck_beads(atomgroup, distance_threshold):
    """
    Find beads that are not stuck to any other bead (distance > threshold).
    Returns a boolean mask for atoms to keep.
    """
    # Get positions for all beads
    positions = atomgroup.positions[:, :2]  # Use only x,y coordinates (for 2D analysis)
    
    # Calculate pairwise distances between all beads
    distances = squareform(pdist(positions))
    
    # Set diagonal to a large value so we don't consider self-distances
    np.fill_diagonal(distances, 999999)
    
    # Find beads that are too close to any other bead
    # For each bead, check if any distance to another bead is less than threshold
    mask = np.all(distances > distance_threshold, axis=1)
    
    return mask

# Calculate MSD for all beads (original approach)
compute_and_save_msd(DBs, output_file_DB_original, label="All Dynabeads (H)")

# Now calculate MSD only for beads that are not stuck to others
print(f"Total number of Dynabeads: {len(DBs)}")

# Create a list to track how many beads were excluded in each frame
excluded_count = []
filtered_indices = []

# Loop through each frame to identify non-stuck beads
for ts in u.trajectory:
    # Find beads that are not stuck to others
    keep_mask = find_non_stuck_beads(DBs, DISTANCE_THRESHOLD)
    
    excluded_count.append(np.sum(~keep_mask))
    print(f"Frame {ts}: Keeping {len(filtered_indices)} beads, excluding {np.sum(~keep_mask)} beads")
    # If this is the first frame, save the indices of beads to include
    if ts.frame == 0:
        filtered_indices = np.where(keep_mask)[0]
        print(f"Frame 0: Keeping {len(filtered_indices)} beads, excluding {np.sum(~keep_mask)} beads")

# Select only the beads that were not stuck at the beginning
# Note: MDAnalysis typically requires consistent atom selections across frames
if len(filtered_indices) > 0:
    filtered_DBs = DBs[filtered_indices]
    print(f"Selected {len(filtered_DBs)} non-stuck beads for MSD analysis")
    
    # Calculate MSD for non-stuck beads
    compute_and_save_msd(filtered_DBs, output_file_DB, label="Non-stuck Dynabeads (H)")
else:
    print("Warning: All beads appear to be too close to each other in the first frame.")
    
print(f"Average number of excluded beads per frame: {np.mean(excluded_count)}")


'''

msd_calc = msd.EinsteinMSD(DBs, msd_type='xy', fft = True, apply_pbc=True)
msd_calc.run()



nframes = msd_calc.n_frames
timestep = 100 
lagtimes = np.arange(nframes)*timestep
print(lagtimes)

MSDs = msd_calc.results.timeseries


with open(output_file, 'w') as f1:
    for t, msd_value in zip(lagtimes, MSDs):
        f1.write(f"{t}\t{msd_value}\n")

'''
