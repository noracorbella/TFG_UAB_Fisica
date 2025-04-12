#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 11:17:25 2025

@author: nora
"""


import MDAnalysis as mda
from MDAnalysis.analysis import msd
import numpy as np

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_4/DB4_traj_unwrapped.xtc"


topology_pdb = r"mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_4/DB4_frame0.pdb"


output_file_DB = r"mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_4/MDAnalysis_MSD_DB4.dat"

output_file_TC = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker_2/TcellsDynabeadSystem_DiffusionTracker_def/MDAnalysis_MSD_TC.dat"



u = mda.Universe(topology_pdb, trajectory_file)



DBs = u.select_atoms('name H')
TCs = u.select_atoms('name O')

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

# Run for each type
compute_and_save_msd(DBs, output_file_DB, label="Dynabeads (H)")
compute_and_save_msd(TCs, output_file_TC, label="T cells (O)")

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
