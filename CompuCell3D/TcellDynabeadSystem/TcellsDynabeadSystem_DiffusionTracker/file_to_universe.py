#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:31:55 2025

@author: nora
"""



positions_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.txt"

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/trajectory_vmd.xtc"

frame0_xyz = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/frame0_xyz.xyz"



trajectory = open(trajectory_file, "w")
frame0_xyz = open(frame0_xyz, "w")

n_TC = 24
n_DB = 24

frame0_xyz.write(f"{n_TC+n_DB}\n")
frame0_xyz.write("TCells Dynabeads frame 0 \n")



MCS = []


with open(positions_file,'r') as f:
    next(f)
    current_mcs = -1
    for l in f:
        mcs, cell_id, cell_type, x, y = map(float, l.split())
        cell_id = int(cell_id)
        if int(cell_type) == 1:
            name = 'O'
        if int(cell_type) == 2:
            name = 'H'
        

        
        if mcs == 0:
            frame0_xyz.write(f"{name} \t {x} \t {y} \t 0.0 \n")
        
        if mcs != current_mcs:
            if current_mcs != -1:
                trajectory.write("\n")  # End previous frame

            current_mcs = mcs
            MCS.append(int(mcs))
            trajectory.write("48\n")
            trajectory.write(f"Frame {int(mcs/100)}\n")
        
        # Write the atom line
        trajectory.write(f"{name} \t {x} \t {y} \t 0.0 \n")

# Ensure the last frame ends properly
trajectory.write("\n")
trajectory.close()


frame0_xyz.close()

frame0_xyz = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/frame0_xyz.xyz"

frame0_gro = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/frame0_gro.gro"




with open(frame0_xyz, "r") as xyz:
    next(xyz)  # Skip the first line
    next(xyz)  # Skip the second line
    atoms = []
    for line in xyz:
        parts = line.split()
        atom_name = parts[0]

        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        atoms.append((atom_name, x * 0.01, y * 0.01, z * 0.01))
        
n_atoms = len(atoms)
      
box_size = (5.000, 5.000, 0.001)

    
with open(frame0_gro, "w") as gro:
    gro.write("MD of 16 waters, t=0.0\n")
    gro.write(f"    {n_atoms}\n")
    atom_id = 1
    residue_name = '1MET'
    for atom_name, x, y, z in atoms:
        gro.write(f"   {residue_name:>5}{atom_name:>4}{atom_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        atom_id += 1
        
    gro.write(f"{box_size[0]:8.5f}{box_size[1]:8.5f}{box_size[2]:8.5f}\n")

gro.close()





























            
                    



