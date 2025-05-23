from cc3d.core.PySteppables import *
import numpy as np
import random
import MDAnalysis
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.coordinates.PDB import PDBWriter



class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=10):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
        Nt = 24 #number of Tcells
        Nd = 24 #number of dynabeads
        Tcell_size = 10 #diameter
        Dynabead_size = 8 #diameter
        
        for _ in range(Nt):
            x, y = self.get_random_position(Tcell_size)
            cell = self.new_cell(self.TCELL)
            self.cell_field[x:x + Tcell_size, y:y + Tcell_size, 0] = cell
            cell.targetVolume = np.pi * (Tcell_size/2) ** 2 
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.sqrt(np.pi * cell.targetVolume)
            cell.lambdaSurface = 5.0

        for _ in range(Nd):
            x, y = self.get_random_position(Dynabead_size)
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = np.pi * (Dynabead_size/2) ** 2 
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
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
        self.num_atoms = len(self.cell_list)
        self.universe = MDAnalysis.Universe.empty(self.num_atoms, trajectory=True)
        
        self.trajectory_writer = XTCWriter(
            r"cell_trajectories.xtc",self.num_atoms)
                
        self.universe.add_TopologyAttr("name")
        self.universe.add_TopologyAttr("resname")

        for i, cell in enumerate(self.cell_list):
            self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]

            
            if cell.type == self.TCELL:
                self.universe.atoms[i].name = "O"
            if cell.type == self.DYNABEAD:
                self.universe.atoms[i].name = "H"            
    
        with PDBWriter(r"frame0.pdb") as pdb_writer:
            pdb_writer.write(self.universe)




    
    def step(self, mcs):
        if mcs % 100 == 0:
            for i, cell in enumerate(self.cell_list):
                self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]
            
            self.trajectory_writer.write(self.universe)

    def finish(self):
        if self.trajectory_writer:
            self.trajectory_writer.close()