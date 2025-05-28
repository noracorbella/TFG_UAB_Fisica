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

            elif cell.type == self.ACTIVATEDTCELL:
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

class NutrientFieldSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.nutrient_field = None
        self.base_uptake_rate = 0.1  # Base uptake rate for regular T cells
        self.activated_uptake_multiplier = 5.0  # Activated T cells consume more nutrients
        self.nutrient_replenishment_rate = 0.05  # Add nutrient replenishment
        self.max_nutrient_concentration = 1.0  # Maximum nutrient concentration
    
    def start(self):
        self.nutrient_field = self.field.Nutrient
        
        # Create a file to track average nutrient levels
        self.nutrient_level_file = open(r"nutrient_levels.txt", "w")
        self.nutrient_level_file.write("MCS\tAvgNutrient\n")
        self.nutrient_level_file.flush()

    def step(self, mcs):
        
        # Replenish nutrients throughout the field
        for x in range(self.dim.x):
            for y in range(self.dim.y):
                current_conc = self.nutrient_field[x, y, 0]
                if current_conc < self.max_nutrient_concentration:
                    self.nutrient_field[x, y, 0] = min(
                        self.max_nutrient_concentration,
                        current_conc + self.nutrient_replenishment_rate
                    )
                    
        for cell in self.cell_list:
            if cell.type == self.TCELL or cell.type == self.ACTIVATEDTCELL:
                raw_x, raw_y = cell.xCOM, cell.yCOM #center of mass coordinates
                
                # PBC by wrapping coordinates
                x = int(raw_x) % self.dim.x
                y = int(raw_y) % self.dim.y
                
                # Uptake rate based on cell type
                if cell.type == self.TCELL:
                    uptake_rate = self.base_uptake_rate
                else:  # Activated T cell
                    uptake_rate = self.base_uptake_rate * self.activated_uptake_multiplier                
                
                # Current nutrient concentration at cell location
                current_conc = self.nutrient_field[x, y, 0]
                
                # Actual uptake (can't consume more than what's available)
                actual_uptake = min(uptake_rate, current_conc)
                
                # Update nutrient field
                self.nutrient_field[x, y, 0] = current_conc - actual_uptake

                # Track nutrients consumed by the cell
                if not hasattr(cell, 'dict') or cell.dict is None:
                    cell.dict = {}
                if 'nutrients_consumed' not in cell.dict:
                    cell.dict['nutrients_consumed'] = 0
                
                cell.dict['nutrients_consumed'] += actual_uptake

        # Record average nutrient level periodically
        if mcs % 100 == 0:
            avg_nutrient = self.calculate_average_nutrient()
            self.nutrient_level_file.write(f"{mcs}\t{avg_nutrient}\n")
            self.nutrient_level_file.flush()

    def calculate_average_nutrient(self):
        total = 0.0
        count = 0
        
        for x in range(self.dim.x):
            for y in range(self.dim.y):
                total += self.nutrient_field[x, y, 0]
                count += 1
        
        return total / count if count > 0 else 0

    def finish(self):
        if hasattr(self, 'nutrient_level_file'):
            self.nutrient_level_file.close()


class TcellGrowthSteppable(MitosisSteppableBase):
    def __init__(self, frequency = 1):
        MitosisSteppableBase.__init__(self, frequency)
        self.Tcell_maxsize = 21.3 #max size a Tcell can achieve 
        
        # Nutrient-dependent growth parameters
        self.nutrient_consumption_threshold = 0.01  # Minimum nutrients needed for growth
        self.max_growth_nutrients = 0.1  # Nutrients needed for maximum growth rate
        self.base_growth_rate = 1.0 / 20.0  # Base growth rate (pixels/MCS) - slower than original
        self.max_growth_rate = 1.0 / 10.0  # Maximum growth rate with optimal nutrients
    
        self.activated_growth_multiplier = 1.0  # Activated cells grow 50% faster

    
    def step(self, mcs):
        max_volume = np.pi * (self.Tcell_maxsize/2) ** 2
        
        for cell in self.cell_list:
            if cell.type == self.TCELL or cell.type == self.ACTIVATEDTCELL:
                if cell.volume <= max_volume:
                    # Get nutrients consumed (default to 0 if not available)
                    nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
                    
                    # Reset nutrients consumed for next step
                    if hasattr(cell, 'dict'):
                        cell.dict['nutrients_consumed'] = 0
                    
                    # Calculate growth based on nutrients
                    if nutrients_consumed >= self.nutrient_consumption_threshold:
                        # Scale growth rate based on nutrient consumption
                        nutrient_factor = min(1.0, nutrients_consumed / self.max_growth_nutrients)
                        
                        # Calculate growth rate interpolating between base and max
                        growth_rate = self.base_growth_rate + (self.max_growth_rate - self.base_growth_rate) * nutrient_factor
                        
                        # Apply activation multiplier if applicable
                        if cell.type == self.ACTIVATEDTCELL:
                            growth_rate *= self.activated_growth_multiplier
                        
                        # Apply growth
                        cell.targetVolume += growth_rate
                        cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
    
class TCellMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self, frequency)  
        self.Tcell_minsize = 17.8 #minimum size for division
        # Add nutrient threshold for division
        self.nutrient_division_threshold = 0.05  # Minimum recent nutrient consumption for division


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


        for cell in self.cell_list_by_type(self.TCELL):
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
            
            if (cell.volume >= min_volume and 
                self.dynabead_neighbor(cell) and 
                nutrients_consumed >= self.nutrient_division_threshold):
                self.divide_cell_random_orientation(cell)
        
        for cell in self.cell_list_by_type(self.TCELL, self.ACTIVATEDTCELL):
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
            
            if (cell.volume >= min_volume and 
                nutrients_consumed >= self.nutrient_division_threshold and
                (cell.type == self.ACTIVATEDTCELL or 
                 (cell.type == self.TCELL and self.dynabead_neighbor(cell)))):
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
        
        self.parent_cell.type = self.TCELL
        self.child_cell.type = self.TCELL
        
    def finish(self):
        if hasattr(self, 'TC_count_file'):
            self.TC_count_file.close()

class TCellActivationSteppable(SteppableBasePy):
    def __init__(self, frequency = 1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        self.activation_count_file = open(r"activated_cells.txt", "w")
        self.activation_count_file.write("MCS\tNormalTCells\tActivatedTCells\n")
        
        # Initialize activation tracking in cell dictionaries
        for cell in self.cell_list_by_type(self.TCELL, self.ACTIVATEDTCELL):
            if not hasattr(cell, 'dict') or cell.dict is None:
                cell.dict = {}
            cell.dict['was_activated'] = cell.type == self.ACTIVATEDTCELL

    def step(self, mcs):
        # Check regular T cells for new activations
        for cell in self.cell_list_by_type(self.TCELL):
            is_touching = self.is_touching_dynabead(cell)
            
            # If in contact with dynabead, activate it
            if is_touching:
                # Preserve cell properties when changing type
                original_target_volume = cell.targetVolume
                original_target_surface = cell.targetSurface
                original_lambda_volume = cell.lambdaVolume
                original_lambda_surface = cell.lambdaSurface
                
                cell.type = self.ACTIVATEDTCELL
                
                # Track that this cell has been activated
                if not hasattr(cell, 'dict') or cell.dict is None:
                    cell.dict = {}
                cell.dict['was_activated'] = True
                
                # Restore properties
                cell.targetVolume = original_target_volume
                cell.targetSurface = original_target_surface
                cell.lambdaVolume = original_lambda_volume
                cell.lambdaSurface = original_lambda_surface
            
            # No code to deactivate cells - they stay activated

        # Log activation stats every 100 MCS
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
