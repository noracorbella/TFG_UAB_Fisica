# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import MDAnalysis as mda
import matplotlib.pyplot as plt
import math
import random

trajectory_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/DB8_traj_unwrapped.xtc"


topology_pdb = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/DB8_frame0.pdb"

DB_positions_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_8/DB8_traj_plot.png"


u = mda.Universe(topology_pdb, trajectory_file)

DB_atoms = u.select_atoms("name cell")

DB_positions = {atom.index: ([], []) for atom in DB_atoms}

for ts in u.trajectory:
    for atom in DB_atoms:
        DB_positions[atom.index][0].append(atom.position[0])  
        DB_positions[atom.index][1].append(atom.position[1])  


plt.figure(figsize=(10, 8))

def plot_traj(data, marker, linewidth, label_prefix="Dynabead"):
    num_beads = len(data)
    colors = [
        (random.random(), random.random(), random.random()) for _ in range(num_beads)
    ]  # Generate random colors
    bead_ids = sorted(data.keys())  # Ensure consistent order

    for i, bead_id in enumerate(bead_ids):
        x_coords, y_coords = data[bead_id]
        x_plot, y_plot = [], []
        color = colors[i]  # Get the color for this bead


        for j in range(len(x_coords)):
            if j == 0 or math.sqrt(
                (x_coords[j] - x_coords[j - 1]) ** 2
                + (y_coords[j] - y_coords[j - 1]) ** 2
            ) <= 100:
                x_plot.append(x_coords[j])
                y_plot.append(y_coords[j])
            else:
                if x_plot:
                    plt.plot(
                        x_plot,
                        y_plot,
                        color=color,
                        marker=marker,
                        linewidth=linewidth,
                    )  # Plot with specific color and label
                x_plot, y_plot = [x_coords[j]], [y_coords[j]]
        if x_plot:
            plt.plot(
                x_plot,
                y_plot,
                color=color,
                marker=marker,
                linewidth=linewidth,
            )  # Plot the last segment

plot_traj(DB_positions, None, 1)  # Plot all beads with different colors

plt.grid(True)
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("Trajectories of Dynabeads")
plt.savefig(DB_positions_plot, format="png", dpi=300)
plt.show()


            
'''


plt.figure(figsize=(10, 8))

def plot_traj(data, color, marker, linewidth, label):
    for cell_id, (x_coord, y_coord) in data.items():
        x_plot, y_plot = [], []
        for i in range(len(x_coord)):
            if i == 0 or math.sqrt((x_coord[i] - x_coord[i-1])**2 + (y_coord[i] - y_coord[i-1])**2) <= 100:
                x_plot.append(x_coord[i])
                y_plot.append(y_coord[i])
            else:
                if x_plot:
                    plt.plot(x_plot, y_plot, color=color, marker=marker, linewidth=linewidth, label=label if cell_id == list(data.keys())[0] else "")
                x_plot, y_plot = [x_coord[i]], [y_coord[i]]
        if x_plot:
            plt.plot(x_plot, y_plot, color=color, marker=marker, linewidth=linewidth, label=label if cell_id == list(data.keys())[0] else "")

plot_traj(DB_positions, 'green', None, 1, 'Dynabeads')

plt.grid(True)
plt.savefig(DB_positions_plot, format='png', dpi=300)
plt.show()


'''