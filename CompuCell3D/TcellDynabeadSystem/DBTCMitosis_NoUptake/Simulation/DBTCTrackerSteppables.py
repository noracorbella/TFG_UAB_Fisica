from cc3d.core.PySteppables import *
import numpy as np
import random
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.coordinates.XTC import XTCWriter


class CellInitialiser(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        Nt = 30 #number of Tcells
        Nd = 30 #number of dynabeads
        Tcell_size = 10.7
        Dynabead_size = 8
        
        for _ in range(Nt):
            x, y = self.get_random_position(Tcell_size)
            cell = self.new_cell(self.TCELL)
            self.cell_field[x:x + Tcell_size, y:y + Tcell_size, 0] = cell
            cell.targetVolume = np.pi * (Tcell_size/2) ** 2 
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
            cell.lambdaSurface = 5.0

        for _ in range(Nd):
            x, y = self.get_random_position(Dynabead_size)
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = np.pi * (Dynabead_size/2) ** 2 
            cell.lambdaVolume = 5.0
            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
            cell.lambdaSurface = 5.0
        
        
    def get_random_position(self, cell_size, max_tries = 1000):
        cell_size = int(cell_size)
        for _ in range(max_tries):   
            x = random.randint(0, self.dim.x - cell_size)
            y = random.randint(0, self.dim.y - cell_size)
            
            if all(self.cell_field[xi, yi, 0] is None 
                for xi in range(x, x + cell_size) 
                for yi in range(y, y + cell_size)
            ):
                
                return x, y
                

        raise RuntimeError("Could not find a non-overlapping position after many tries.")



class TrajectoryTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.trajectory_writer = None
        self.universe = None

    def start(self):
        self.num_atoms = len(self.cell_list)  
        self.universe = MDAnalysis.Universe.empty(self.num_atoms, n_residues=3, trajectory=True)
        
        self.universe.add_TopologyAttr("name")
        self.universe.add_TopologyAttr("type")
        self.universe.add_TopologyAttr("resname")

        
        for i, cell in enumerate(self.cell_list):
            atom = self.universe.atoms[i]
            atom.position = [cell.xCOM, cell.yCOM, 0.0]
            if cell.type == self.TCELL:
                atom.name = "C"
                atom.type = "C"
                atom.residue.resname = "TCL"

            elif cell.type == self.DYNABEAD:
                atom.name = "H"
                atom.type = "H"
                atom.residue.resname = "DYN"

            else:
                atom.name = "C"
                atom.type = "C"
                atom.residue.resname = "ATC"

        with PDBWriter(r"DBTC_InitialFrame.pdb") as pdb_writer:
            pdb_writer.write(self.universe.atoms)

        self.trajectory_writer = XTCWriter(r"DBTC_trajectories.xtc", self.num_atoms)

    def step(self, mcs):
        if mcs % 100 == 0:
            for i, cell in enumerate(self.cell_list):
                if i < len(self.universe.atoms):
                    self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0] 
            self.trajectory_writer.write(self.universe)



class TcellGrowthSteppable(MitosisSteppableBase):
    def __init__(self, frequency = 1):
        MitosisSteppableBase.__init__(self, frequency)
        self.Tcell_maxsize = 21.3 #max size a Tcell can achieve 
    
    def step(self, mcs):
        growth_rate = (1.0 / 100.0 ) # pixel / MCS
        max_volume = np.pi * (self.Tcell_maxsize/2) ** 2
        
        for cell in self.cell_list:
            if cell.type == self.TCELL or cell.type == self.ACTIVATEDTCELL:
                if cell.volume <= max_volume:
                    cell.targetVolume += growth_rate
                    cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)    
                
class TCellMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self, frequency)  
        self.Tcell_minsize = 17.8 #minimum size for division

    def start(self):
        self.TC_count_file = open(r"TC_count.txt", "w")
        self.TC_count_file.write("MCS\tNormalTCells\tActivatedTCells\tTotalTCells\n")

        self.TC_count_file.flush() 



    def step(self, mcs):
        regular_count = len(self.cell_list_by_type(self.TCELL))
        activated_count = len(self.cell_list_by_type(self.ACTIVATEDTCELL))
        total_count = regular_count + activated_count
        
        self.TC_count_file.write(f"{mcs}\t{regular_count}\t{activated_count}\t{total_count}\n")
        self.TC_count_file.flush()


        for cell in self.cell_list_by_type(self.TCELL, self.ACTIVATEDTCELL):
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            if cell.volume >= min_volume and self.dynabead_neighbor(cell):
                self.divide_cell_random_orientation(cell)
        


        
    def dynabead_neighbor(self, cell):
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor and neighbor.type == self.DYNABEAD and common_surface_area > 0:
                return True
        return False

    
    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0
        self.parent_cell.targetSurface = 2 * np.pi * np.sqrt(self.parent_cell.targetVolume)
        
        self.clone_parent_2_child()
        
        self.child_cell.type = self.parent_cell.type
        
    def finish(self):
        if hasattr(self, 'TC_count_file'):
            self.TC_count_file.close()



class TCellActivationSteppable(SteppableBasePy):
    def __init__(self, frequency = 1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        self.activation_count_file = open("activated_cells.txt", "w")
        self.activation_count_file.write("MCS\tNormalTCells\tActivatedTCells\n")

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.TCELL):
            if self.is_touching_dynabead(cell):
                
                original_target_volume = cell.targetVolume
                original_target_surface = cell.targetSurface
                original_lambda_volume = cell.lambdaVolume
                original_lambda_surface = cell.lambdaSurface
                
                cell.type = self.ACTIVATEDTCELL
                
                cell.targetVolume = original_target_volume
                cell.targetSurface = original_target_surface
                cell.lambdaVolume = original_lambda_volume
                cell.lambdaSurface = original_lambda_surface

        if mcs % 100 == 0:
            normal_count = len(self.cell_list_by_type(self.TCELL))
            activated_count = len(self.cell_list_by_type(self.ACTIVATEDTCELL))
            self.activation_count_file.write(f"{mcs}\t{normal_count}\t{activated_count}\n")
            self.activation_count_file.flush()


    def is_touching_dynabead(self, cell):
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor and neighbor.type == self.DYNABEAD and common_surface_area > 0:
                return True
        return False

    def finish(self):
        if hasattr(self, 'activation_count_file'):
            self.activation_count_file.close()



