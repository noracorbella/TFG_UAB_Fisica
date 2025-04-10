#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:31:55 2025

@author: nora
"""

import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import XYZ



trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/trajectory_vmd.xtc"
trajectory_gro = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/trajectory_gro.gro"


frames = []

box_size = (5.000, 5.000, 0.001)

with open(trajectory_file,'r') as f:
    
    lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("Frame") or line == 48:
            i += 1
            continue
        
        num_atoms = int(line)
        i += 1
        
        frame_data = []
        while len(frame_data) < num_atoms:
            line = lines[i].strip()
            parts = line.split()
            if len(parts) >= 4: 
                atom_type = parts[0]
                try:
                    x, y, z = map(float, parts[1:])
                    frame_data.append((atom_type, x/100, y/100, z/100))
                except ValueError:
                    print(f"Skipping invalid line: {line}")
            i += 1
        frames.append(frame_data)
        
        i += 1
        
f.close()



with open(trajectory_gro, "w") as gro:
    for frame_idx, frame_data in enumerate(frames):
        
        t = (frame_idx) / 10
        
        gro.write(f"MD of 16 waters, t={t:.1f}\n")        
        gro.write(f"{num_atoms}\n")

        for atom_idx, (atom_type, x, y, z) in enumerate(frame_data):
            atom_name = atom_type.strip()[:5]  
            residue_name = "1MET"  
            residue_id = 1  
            atom_id = atom_idx + 1  
            gro.write(f"   {residue_name:>5}{atom_name:>4}{atom_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
  
        gro.write(f"{box_size[0]:8.5f}{box_size[1]:8.5f}{box_size[2]:8.5f}\n")

gro.close()















            
                    



