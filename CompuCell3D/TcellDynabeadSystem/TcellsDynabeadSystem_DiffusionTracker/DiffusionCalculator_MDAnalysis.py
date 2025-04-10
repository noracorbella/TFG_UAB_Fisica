#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 11:17:25 2025

@author: nora
"""


import MDAnalysis as mda
from MDAnalysis.analysis import msd
from MDAnalysis.coordinates.XTC import XTCReader
import MDAnalysis.transformations as trans
import numpy as np

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DBTracker/DB_trajectories.xtc"

topology_pdb = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DBTracker/DB_frame0.pdb"


output_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DBTracker/MDAnalysis_MSD_DB.dat"



#We only have a trajectory file so we make the topology file from frame 0 of the trajectory file

with XTCReader(trajectory_file) as trajectory:
    n_atoms = trajectory.n_atoms 
u = mda.Universe.empty(n_atoms, trajectory=True, atom_resindex=[0]*n_atoms, residue_segindex=[0])    

u.add_TopologyAttr("name", ["H"]*n_atoms)
u.add_TopologyAttr("type", ["H"]*n_atoms)
u.add_TopologyAttr("resnames", ["MET"])
u.add_TopologyAttr("resids", [1])
u.add_TopologyAttr("segids", ["SYSTEM"]) 


u.load_new(trajectory_file)
u.trajectory[0]
u.atoms.write(topology_pdb)

print("Created topology.pdb with first frame coordinates!")


u = mda.Universe(topology_pdb, trajectory_file)
u.trajectory.add_transformations(trans.wrap(u.atoms)) # Unwrap trajectory to correct for periodic boundary effects

dt = 1.0 

times = [i * dt for i in range(len(u.trajectory))]  
print(times)
DBs = u.select_atoms('name H')

msd_calc = msd.EinsteinMSD(DBs, msd_type='xyz', fft = True)
msd_calc.run()




MSDs = msd_calc.results.timeseries


with open(output_file, 'w') as f1:
    for t, msd_value in zip(times, MSDs):
        f1.write(f"{t}\t{msd_value}\n")
