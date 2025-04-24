from cc3d.core.PySteppables import *
import numpy as np
import random
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.coordinates.XTC import XTCWriter


class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=10):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
        Nt = 1 #number of Tcells
        Nd = Nt #number of dynabeads
        Dynabead_size = 8 #diameter
        

        for _ in range(Nd):
            x, y = self.get_random_position(Dynabead_size)
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = np.pi * (Dynabead_size/2) ** 2 
            cell.lambdaVolume = 20.0
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
        self.universe.add_TopologyAttr("name")
        
        for i, cell in enumerate(self.cell_list):
            self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]  
            self.universe.atoms[i].name = f"cell_{cell.id}"
            
        with PDBWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTracker\DBTracker_6\DB6_frame0.pdb") as pdb_writer:
            pdb_writer.write(self.universe.atoms)

        self.trajectory_writer = XTCWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTracker\DBTracker_6\DB6_trajectories.xtc", self.num_atoms)

        self.plot_win = self.add_new_plot_window(title = "Volume", x_axis_title = "MCS", y_axis_title = "V(pixels^2)", x_scale_type = "linear", y_scale_type = "linear", grid = True)
    
        self.plot_win.add_plot("Volume", style = "Donts", color = "red", size = 5)
        
        
        
    def step(self, mcs):
        if mcs % 100 == 0:
            for i, cell in enumerate(self.cell_list):
                self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]  
                
            self.trajectory_writer.write(self.universe)
            
        Volume = 0
        L = len(self.cell_list)
        for cell in self.cell_list:
            Volume += cell.volume/L
            
        self.plot_win.add_data_point("Volume", mcs, Volume)

