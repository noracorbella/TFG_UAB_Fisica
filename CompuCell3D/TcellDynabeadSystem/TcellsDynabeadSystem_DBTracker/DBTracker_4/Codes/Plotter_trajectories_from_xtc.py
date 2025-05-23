
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import math

traj_file = r"DB4_traj_unwrapped.xtc"
top_file = r"DB4_frame0.pdb"
plot = r"DB4_traj_plot.png"

u = mda.Universe(top_file, traj_file)

print(f"Number of frames: {len(u.trajectory)}")
print(f"Number of atoms: {len(u.atoms)}")
print(f"Atom info: {u.atoms[0] if len(u.atoms) > 0 else 'No atoms'}")

atom = u.atoms[0]

def plot_traj(data, color, marker, linewidth):
    for cell_id, (x_coord, y_coord) in data.items():
        x_plot, y_plot = [], []
        for i in range(len(x_coord)):
            if i == 0 or math.sqrt((x_coord[i] - x_coord[i-1])**2 + (y_coord[i] - y_coord[i-1])**2) <= 100:
                x_plot.append(x_coord[i])
                y_plot.append(y_coord[i])
            else:
                if x_plot:
                    plt.plot(x_plot, y_plot, color=color, marker=marker, linewidth=linewidth)
                x_plot, y_plot = [x_coord[i]], [y_coord[i]]
        if x_plot:
            plt.plot(x_plot, y_plot, color=color, marker=marker, linewidth=linewidth)

positions = []
for ts in u.trajectory:
    positions.append(atom.position[:2].copy()) 

positions = np.array(positions)
data = {0: (positions[:, 0], positions[:, 1])}

plot_traj(data, 'green', None, 1)

plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('DB Trajectory')
plt.xlim(0, 500)
plt.ylim(0, 500)
plt.savefig(plot, format='png', dpi=300)
plt.show()




