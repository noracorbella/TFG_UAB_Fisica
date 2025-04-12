"""
The .xtc file is required to have the coordinates in the unwrapped convention.
This can be done with GROMACS software, with the following command line:
gmx trjconv -f DB_traj.xtc -s DB_frame0.pdb -pbc nojump -o DB_traj_unwrapped.xtc
where DB_traj.xtc is the wrapped trajectory file, DB_frame0.pdb is the topology file and DB_traj_unwrapped.xtc is the converted trajectory file with the unwrapped convention.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import msd
import numpy as np

trajectory_file = r"DB4_traj_unwrapped.xtc"
topology_pdb = r"DB4_frame0.pdb"
output_file_DB = r"MDAnalysis_MSD_DB4.dat"

u = mda.Universe(topology_pdb, trajectory_file)

DBs = u.select_atoms('all')

print(f"Number of DB: {len(DBs)}") #it should be 1

msd_calc = msd.EinsteinMSD(DBs, msd_type='xy', fft = True, apply_pbc=True)
msd_calc.run()

nframes = msd_calc.n_frames
timestep = 100 
lagtimes = np.arange(nframes)*timestep
print(lagtimes)

MSDs = msd_calc.results.timeseries

with open(output_file_DB, 'w') as f1:
    for t, msd_value in zip(lagtimes, MSDs):
        f1.write(f"{t}\t{msd_value}\n")
