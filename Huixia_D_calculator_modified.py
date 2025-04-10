#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 15:44:19 2025

@author: nora
"""

#-------------------------------------------------------------------------------------------
# calculate diffusion of "OW", get "initial_frames" every "delta_frame" from DCD file
#-------------------------------------------------------------------------------------------

import MDAnalysis as mda
import numpy as np

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/TrajectoryFiles/trajectories_mol2.mol2"


topology_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/frame0_gro.gro"
output_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/MSD_DB.dat"



u = mda.Universe(topology_file, trajectory_file)

print(u.trajectory)


G = u.select_atoms('name H') #Select DB (accounted as H in the .gro file)
num_G = len(G)

print(num_G)


#Parameters for the MSD calculation
dt = u.trajectory.dt # Time step

print(dt)

steps=6000  # Number of time steps to calculate MSD
initial_frames=9  # Number of initial frames to average over
delta_frame=100 # Frame interval for sampling

Lx, Ly = 5.00000, 5.00000
#Lz = 0.00100

#Periodic boundary condition correction with minimum image convention
def pbc(dr, box):
    return dr - box * np.round(dr / box)

# Coordinate lists
x, y = [], []
#z = []

# Initial positions
u.trajectory[0]
x0 = G.positions[:,0].copy()
y0 = G.positions[:,1].copy()
#z0 = G.positions[:,2].copy()

x.append(x0)
y.append(y0)
#z.append(z0)

#Unwrap trajectory
for i in range(1, len(u.trajectory)):
    
    u.trajectory[i]

    
    #Current positions
    xn = G.positions[:,0]
    yn = G.positions[:,1]
    #zn = G.positions[:,2]
    
    # Calculate displacements with PBC corrections
    dx = xn - x[-1]
    dx = pbc(dx, Lx)
    
    dy = yn - y[-1]
    dy = pbc(dy, Ly)
    
    #dz = zn - z[-1]
    #dz = pbc(dz, Lz)

	
    # Update Coordinates
    x.append(x[-1] + dx)
    y.append(y[-1] + dy)
    #z.append(z[-1] + dz)
				

# Initialize avg MSD values
ave = np.zeros(steps)
count = 0

# Calculate MSD
for r in range(initial_frames):
	
    first_frame = r * delta_frame 
    end_frame = first_frame + steps
    
    
    if end_frame > len(u.trajectory):
        end_frame = len(u.trajectory)
        print("end_frame > len(u.trajectory)")
        print(len(u.trajectory))
        steps = end_frame - first_frame #Adjust steps

    if steps <= 0:
        print("steps <= 0")
        continue

    # Reference positions
    xr = x[first_frame]
    yr = y[first_frame]
    #zr = z[first_frame]

    # Temporary MSD 
    result = np.zeros(steps)

    # MSD for each time window
    for i in range(first_frame + 1, end_frame):
        dx = x[i] - xr
        dy = y[i] - yr
        #zi = z[i]
        result[i-first_frame] = np.mean(dx**2 + dy**2)

    ave[:steps] += result
    count += 1


ave /= count

with open(output_file, 'w') as f1:
    for i in range(steps):
        f1.write(f"{i*dt}\t{ave[i]}\n")






'''
    # Average over all DB
    for i in range(steps):
        result[i] /= num_G
        
# Average over all initial frames
for i in range(steps):
    ave[i] /= initial_frames


# Write MSD
for i in range(steps):
    f1.write(str(i * dt) + '\t' + str(ave[i]) + '\n')


f1.close()
sys.exit()
'''
