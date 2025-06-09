from cc3d.core.PySteppables import *
import numpy as np
import random
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.coordinates.XTC import XTCWriter


class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=10):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
        Nt = 24 #number of Tcells
        Nd = Nt #number of dynabeads
        Dynabead_size = 8 #diameter
        

        for _ in range(Nd):
            x, y = self.get_random_position(Dynabead_size)
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = np.pi * (Dynabead_size/2) ** 2 
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.sqrt(cell.targetVolume * np.pi )
            cell.lambdaSurface = 5.0
    
    
    def get_random_position(self, cell_size):
        x = random.randint(0, self.dim.x - cell_size)
        y = random.randint(0, self.dim.y - cell_size)
        return x, y

     
        
        

class TrajectoryTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.trajectory_writer = None
        self.universe = None

    def start(self):
        # Define the initial positions of cells in your simulation space
        self.num_atoms = len(self.cell_list)  # Assuming each cell is an atom
        
        # Create an empty Universe for MDAnalysis with the same number of cells
        self.universe = MDAnalysis.Universe.empty(self.num_atoms, trajectory=True)  # No bonds, just positions
        
        # Create the xtc file writer
        self.trajectory_writer = XTCWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\TcellsDynabeadSystem_DBTracker\DB_trajectories.xtc", self.num_atoms)
        
        # Prepare the atoms and set the coordinates as we don't have topology data
        self.universe.add_TopologyAttr("name")
        for i, cell in enumerate(self.cell_list):
            self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]  # Assuming 2D positions, z=0 for simplicity
            self.universe.atoms[i].name = f"cell_{cell.id}"

    def step(self, mcs):
        if mcs % 100 == 0:
            # Update the positions of cells and write to the trajectory file
            for i, cell in enumerate(self.cell_list):
                self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]  # Update 2D positions
            # Write the current frame to the XTC file
            self.trajectory_writer.write(self.universe)

    def finish(self):
        # Ensure the trajectory writer finalizes the xtc file correctly
        self.trajectory_writer.close() 
