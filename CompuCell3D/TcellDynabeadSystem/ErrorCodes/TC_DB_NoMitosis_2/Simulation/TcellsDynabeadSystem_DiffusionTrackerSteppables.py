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
        Nd = 24 #number of dynabeads
        Dynabead_size = 8 #diameter
        Tcell_size = 10

        for _ in range(Nd):
            while True:
                x, y = self.get_random_position(Dynabead_size)
                if self.is_area_free(x, y, Dynabead_size):
                    break
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = np.pi * (Dynabead_size / 2) ** 2
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
            cell.lambdaSurface = 5.0
    
        for _ in range(Nt):
            while True:
                x, y = self.get_random_position(Tcell_size)
                if self.is_area_free(x, y, Tcell_size):
                    break
            cell = self.new_cell(self.TCELL)
            self.cell_field[x:x + Tcell_size, y:y + Tcell_size, 0] = cell
            cell.targetVolume = np.pi * (Tcell_size / 2) ** 2
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
            cell.lambdaSurface = 5.0
        
    
    def get_random_position(self, cell_size):
        x = random.randint(0, self.dim.x - cell_size)
        y = random.randint(0, self.dim.y - cell_size)
        return x, y

    def is_area_free(self, x, y, size):
        for i in range(x, x + size):
            for j in range(y, y + size):
                if self.cell_field[i, j, 0] is not None:
                    return False
        return True    
   
   
        
class TrajectoryTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.trajectory_writer = None
        self.universe = None

    def start(self):
        self.num_atoms = len(self.cell_list)
        atom_resindex = np.arange(self.num_atoms)
        
        self.universe = MDAnalysis.Universe.empty(n_atoms=self.num_atoms, n_residues=self.num_atoms, n_segments=1,atom_resindex=atom_resindex, trajectory=True)
        
        
        
        for i, cell in enumerate(self.cell_list):
            atom = self.universe.atoms[i]
            residue = atom.residue

            atom.position = [cell.xCOM, cell.yCOM, 0.0]
            
            if cell.type == self.TCELL:
                atom.name = "T"
                residue.resname = "TCE"
            elif cell.type == self.DYNABEAD:
                atom.name = "D"
                residue.resname = "DYN"
        
        self.trajectory_writer = XTCWriter("cell_trajectories.xtc", self.num_atoms)


    def step(self, mcs):
        if mcs % 100 == 0:
            for i, cell in enumerate(self.cell_list):
                self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]
            self.trajectory_writer.write(self.universe)


            
            
