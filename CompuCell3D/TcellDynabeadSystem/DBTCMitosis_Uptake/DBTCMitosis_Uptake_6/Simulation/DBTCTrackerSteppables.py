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
        Nt = 15 #number of Tcells
        Nd = 15 #number of dynabeads
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
        
        tcell_count = len(self.cell_list_by_type(self.TCELL))
        dynabead_count = len(self.cell_list_by_type(self.DYNABEAD))
        activated_count = len(self.cell_list_by_type(self.ACTIVATEDTCELL)) 
        
        self.universe = MDAnalysis.Universe.empty(self.num_atoms, n_residues=3, trajectory=True)
        
        self.universe.add_TopologyAttr("name")
        self.universe.add_TopologyAttr("type")
        self.universe.add_TopologyAttr("resname")
        self.universe.add_TopologyAttr("resid")
        
        residues = self.universe.residues
        residues[0].resname = "TCL"
        residues[1].resname = "DYN"
        residues[2].resname = "ATC"
        
        tcell_atoms = []
        dynabead_atoms = []
        activated_atoms = []
        
        for i, cell in enumerate(self.cell_list):
            atom = self.universe.atoms[i]
            atom.position = [cell.xCOM, cell.yCOM, 0.0]
            
            if cell.type == self.TCELL:
                atom.name = "C"
                atom.type = "C"
                tcell_atoms.append(atom)

            elif cell.type == self.DYNABEAD:
                atom.name = "H"
                atom.type = "H"
                dynabead_atoms.append(atom)

            elif cell.type == self.ACTIVATEDTCELL:
                atom.name = "C"
                atom.type = "C"
                activated_atoms.append(atom)

        for atom in tcell_atoms:
            atom.residue = residues[0]
        for atom in dynabead_atoms:
            atom.residue = residues[1]
        for atom in activated_atoms:
            atom.residue = residues[2]

        with PDBWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\DBTC_frame0.pdb") as pdb_writer:
            pdb_writer.write(self.universe.atoms)

        self.trajectory_writer = XTCWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\DBTC_traj.xtc", self.num_atoms)

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
        self.base_uptake_rate = 0.02  # Base uptake rate for regular T cells
        self.activated_uptake_multiplier = 1.5  # Activated T cells consume more nutrients
        self.nutrient_replenishment_rate = 0.03  # Add nutrient replenishment
        self.max_nutrient_concentration = 1.0 # Maximum nutrient concentration
        # Track cell level scalar attribute for visualization
        self.track_cell_level_scalar_attribute(field_name='uptake', attribute_name='nutrients_consumed')
    
    def start(self):
        self.nutrient_field = self.field.Nutrient
        
        # Create a file to track average nutrient levels
        self.nutrient_uptake_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\nutrient_uptake_data", "w")
        self.nutrient_uptake_file.write("MCS\tAvgNutrientConcentration\tAvgUptake\n")
        self.nutrient_uptake_file.flush()
        

        
        # Add plot windows for visualization
        self.plot_win4 = self.add_new_plot_window(title='Average Uptake',
                                                x_axis_title='MonteCarlo Step (MCS)',
                                                y_axis_title='Uptake', x_scale_type='linear', y_scale_type='linear',
                                                grid=True)
        
        self.plot_win5 = self.add_new_plot_window(title='Average Nutrient Concentration',
                                                x_axis_title='MonteCarlo Step (MCS)',
                                                y_axis_title='Uptake', x_scale_type='linear', y_scale_type='linear',
                                                grid=True) 
                                                
        self.plot_win4.add_plot("AvgUptake", style='Dots', color='lime', size=5)
        self.plot_win5.add_plot("AvgConc", style='Dots', color='orange', size=5)

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
        
        # Track total uptake for averaging
        total_uptake = 0.0
        cell_count = 0
                    
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
                
                cell.dict['nutrients_consumed'] = actual_uptake  # Set to current uptake (not cumulative)
                
                # Add to total for averaging
                total_uptake += actual_uptake
                cell_count += 1

        # Calculate average uptake
        avg_uptake = total_uptake / cell_count if cell_count > 0 else 0
        

        avg_nutrient = self.calculate_average_nutrient()
        
        self.plot_win4.add_data_point("AvgUptake", mcs, avg_uptake)
        self.plot_win5.add_data_point("AvgConc", mcs, avg_nutrient)
  
        # Record average nutrient level and uptake periodically
        if mcs % 100 == 0:  # More frequent updates for better plots            
            # Write to files
            self.nutrient_uptake_file.write(f"{mcs}\t{avg_nutrient}\t{avg_uptake}\n")
            self.nutrient_uptake_file.flush()
            


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
        if hasattr(self, 'uptake_file'):
            self.uptake_file.close()


class TcellGrowthSteppable(MitosisSteppableBase):
    def __init__(self, frequency = 1):
        MitosisSteppableBase.__init__(self, frequency)
        self.Tcell_maxsize = 21.3 #max size a Tcell can achieve 
        
        # Nutrient-dependent growth parameters
        self.nutrient_consumption_threshold = 0.005  # Minimum nutrients needed for growth
        self.max_growth_nutrients = 0.08  # Nutrients needed for maximum growth rate
        self.base_growth_rate = 1.0 / 400.0  # Base growth rate (pixels/MCS) - slower than original
        self.max_growth_rate = 1.0 / 180.0  # Maximum growth rate with optimal nutrients
        
        '''
        self.activated_growth_multiplier = 1.5  # Activated cells grow 50% faster
        '''
        
        self.pressure_threshold = 70.0  # Pressure threshold above which cells won't grow
    
    def step(self, mcs):
        max_volume = np.pi * (self.Tcell_maxsize/2) ** 2
        total_pressure = 0
        num_cells = len(self.cell_list)
        if num_cells == 0:
            return
            
        for cell in self.cell_list:
            if cell.type == self.TCELL or cell.type == self.ACTIVATEDTCELL:
                cell.dict['pressure'] = 2 * cell.lambdaVolume * (cell.targetVolume - cell.volume)
                total_pressure += cell.dict['pressure']        
        
        avg_pressure = total_pressure / num_cells

        for cell in self.cell_list:
            if cell.type == self.TCELL or cell.type == self.ACTIVATEDTCELL:
                if cell.volume <= max_volume:
                    current_pressure = cell.dict['pressure']    
                    if current_pressure < self.pressure_threshold:    
                        # Get nutrients consumed (default to 0 if not available)
                        nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
                        
                        # Calculate growth based on nutrients (using current uptake, not cumulative)
                        if nutrients_consumed >= self.nutrient_consumption_threshold:
                            # Scale growth rate based on nutrient consumption
                            nutrient_factor = min(1.0, nutrients_consumed / self.max_growth_nutrients)
                            
                            # Calculate growth rate interpolating between base and max
                            growth_rate = self.base_growth_rate + (self.max_growth_rate - self.base_growth_rate) * nutrient_factor
                            '''
                            # Apply activation multiplier if applicable
                            if cell.type == self.ACTIVATEDTCELL:
                                growth_rate *= self.activated_growth_multiplier
                            '''
                            # Apply growth
                            cell.targetVolume += growth_rate
                            cell.targetSurface = 2 * np.pi * np.sqrt(cell.targetVolume)
    
class TCellMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self, frequency)  
        self.Tcell_minsize = 17.8 #minimum size for division
        # Add nutrient threshold for division
        self.nutrient_division_threshold = 0.03  # Minimum recent nutrient consumption for division
        self.max_divisions = 10

    def start(self):
        self.TC_count_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\TC_count.txt", "w")
        self.TC_count_file.write("MCS\tNormalTCells\tActivatedTCells\tTotalTCells\n")

        self.TC_count_file.flush() 
        # Initialize division counter for T cells
        for cell in self.cell_list_by_type(self.TCELL, self.ACTIVATEDTCELL):
            if not hasattr(cell, 'dict') or cell.dict is None:
                cell.dict = {}
            if 'division_count' not in cell.dict:
                cell.dict['division_count'] = 0


    def step(self, mcs):
        regular_count = len(self.cell_list_by_type(self.TCELL))
        activated_count = len(self.cell_list_by_type(self.ACTIVATEDTCELL))
        total_count = regular_count + activated_count
        
        self.TC_count_file.write(f"{mcs}\t{regular_count}\t{activated_count}\t{total_count}\n")
        self.TC_count_file.flush()


        for cell in self.cell_list_by_type(self.TCELL):
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
            division_count = cell.dict.get('division_count', 0) if hasattr(cell, 'dict') else 0
            
            if (cell.volume >= min_volume and 
                self.dynabead_neighbor(cell) and 
                nutrients_consumed >= self.nutrient_division_threshold and
                division_count < self.max_divisions):
                self.divide_cell_random_orientation(cell)
        
        for cell in self.cell_list_by_type(self.TCELL, self.ACTIVATEDTCELL):
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            nutrients_consumed = cell.dict.get('nutrients_consumed', 0) if hasattr(cell, 'dict') else 0
            division_count = cell.dict.get('division_count', 0) if hasattr(cell, 'dict') else 0
            
            if (cell.volume >= min_volume and 
                nutrients_consumed >= self.nutrient_division_threshold and
                division_count < self.max_divisions and
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
        
        # Increment division counter for the parent cell
        if not hasattr(self.parent_cell, 'dict') or self.parent_cell.dict is None:
            self.parent_cell.dict = {}
        
        parent_division_count = self.parent_cell.dict.get('division_count', 0)
        self.parent_cell.dict['division_count'] = parent_division_count + 1
        
        # Clone parent attributes to child including division history
        self.clone_parent_2_child()
        
        # Ensure the child cell also has the division count
        if not hasattr(self.child_cell, 'dict') or self.child_cell.dict is None:
            self.child_cell.dict = {}
        self.child_cell.dict['division_count'] = parent_division_count + 1
        
        self.child_cell.type = self.parent_cell.type
        
    def finish(self):
        if hasattr(self, 'TC_count_file'):
            self.TC_count_file.close()
            
            
class TCellActivationSteppable(SteppableBasePy):
    
    def __init__(self, frequency = 1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        
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

    def is_touching_dynabead(self, cell):
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor and neighbor.type == self.DYNABEAD and common_surface_area > 0:
                return True
        return False


            
            
 
class EnergyTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.output_file = None
        self.data_file = None
        
        # Define cell types as constants
        self.MEDIUM_ID = 0
        self.TCELL_ID = 1
        self.DYNABEAD_ID = 2
        self.ACTIVATEDTCELL_ID = 3
        
        # Define contact energy table
        self.contact_energy_table = {
            self.MEDIUM_ID: {
                self.MEDIUM_ID: 0.0,
                self.TCELL_ID: 20.0,
                self.DYNABEAD_ID: 20.0,
                self.ACTIVATEDTCELL_ID: 20.0
            },
            self.TCELL_ID: {
                self.MEDIUM_ID: 20.0,
                self.TCELL_ID: 30.0,
                self.DYNABEAD_ID: 15.0,
                self.ACTIVATEDTCELL_ID: 30.0
            },
            self.DYNABEAD_ID: {
                self.MEDIUM_ID: 20.0,
                self.TCELL_ID: 15.0,
                self.DYNABEAD_ID: 100.0,
                self.ACTIVATEDTCELL_ID: 15.0
            },
            self.ACTIVATEDTCELL_ID: {
                self.MEDIUM_ID: 20.0,
                self.TCELL_ID: 30.0,
                self.DYNABEAD_ID: 15.0,
                self.ACTIVATEDTCELL_ID: 30.0
            }
        }

        self.neighbor_order = 4
        
    def start(self):
        # Create and open file for writing energy components
        self.output_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\energy_file.txt", "w")
        self.output_file.write("MCS\tContactEnergy\tVolumeEnergy\tSurfaceEnergy\tTotalEnergy\n")
        
        # Create and open file for writing cell data (volume and pressure)
        self.data_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTCMitosis_Uptake\DBTCMitosis_Uptake_5\cell_data.txt", "w")
        self.data_file.write("MCS\tNormalTCellAvgVolume\tNormalTCellAvgPressure\tActivatedTCellAvgVolume\tActivatedTCellAvgPressure\n")
        
        # Setup plot windows
        self.plot_win = self.add_new_plot_window(title="Energy Components", x_axis_title="MCS", y_axis_title="Energy", 
                                              x_scale_type="linear", y_scale_type="linear", grid=True)

        self.plot_win2 = self.add_new_plot_window(title="Average Pressure", x_axis_title="MCS", y_axis_title="Pressure", 
                                              x_scale_type="linear", y_scale_type="linear", grid=True)        
        
        self.plot_win3 = self.add_new_plot_window(title="Average Volume", x_axis_title="MCS", y_axis_title="Volume", 
                                              x_scale_type="linear", y_scale_type="linear", grid=True)                                              
        
        # Add plots for energy components
        self.plot_win.add_plot("Contact Energy", style="Lines", color="red", size=1)
        self.plot_win.add_plot("Volume Energy", style="Lines", color="cyan", size=1)
        self.plot_win.add_plot("Surface Energy", style="Lines", color="lime", size=1)
        self.plot_win.add_plot("Total Energy", style="Lines", color="orange", size=2)
        
        # Add plots for cell data
        self.plot_win2.add_plot("Normal T Cell Pressure", style="Lines", color="cyan", size=1)
        self.plot_win2.add_plot("Activated T Cell Pressure", style="Lines", color="red", size=1)
         
        self.plot_win3.add_plot("Normal T Cell Volume", style="Lines", color="cyan", size=1)
        self.plot_win3.add_plot("Activated T Cell Volume", style="Lines", color="red", size=1)
        
    def step(self, mcs):
        if mcs % self.frequency == 0:
            # Calculate energy components
            contact_energy = self.calculate_contact_energy()
            volume_energy = self.calculate_volume_energy()
            surface_energy = self.calculate_surface_energy()
            total_energy = contact_energy + volume_energy + surface_energy
            
            # Write energy components to file
            self.output_file.write(f"{mcs}\t{contact_energy}\t{volume_energy}\t{surface_energy}\t{total_energy}\n")
            
            # Update energy plot
            self.plot_win.add_data_point("Contact Energy", mcs, contact_energy)
            self.plot_win.add_data_point("Volume Energy", mcs, volume_energy)
            self.plot_win.add_data_point("Surface Energy", mcs, surface_energy)
            self.plot_win.add_data_point("Total Energy", mcs, total_energy)
            
            # Calculate and track T cell data
            normal_tcell_volume, normal_tcell_pressure = self.calculate_cell_data(self.TCELL_ID)
            activated_tcell_volume, activated_tcell_pressure = self.calculate_cell_data(self.ACTIVATEDTCELL_ID)
            
            # Write cell data to file
            self.data_file.write(f"{mcs}\t{normal_tcell_volume}\t{normal_tcell_pressure}\t{activated_tcell_volume}\t{activated_tcell_pressure}\n")
            
            # Update cell data plots
            self.plot_win2.add_data_point("Normal T Cell Pressure", mcs, normal_tcell_pressure)
            self.plot_win2.add_data_point("Activated T Cell Pressure", mcs, activated_tcell_pressure)
            
            self.plot_win3.add_data_point("Normal T Cell Volume", mcs, normal_tcell_volume)
            self.plot_win3.add_data_point("Activated T Cell Volume", mcs, activated_tcell_volume)
            
            # Flush files to ensure data is written
            self.output_file.flush()
            self.data_file.flush()
            
            # Print status every 1000 MCS
            if mcs % 1000 == 0:
                print(f"MCS: {mcs}, Total Energy: {total_energy}")
                tcell_count = len(self.cell_list_by_type(self.TCELL_ID))
                activated_count = len(self.cell_list_by_type(self.ACTIVATEDTCELL_ID))
                dynabead_count = len(self.cell_list_by_type(self.DYNABEAD_ID))
                print(f"Cell counts - T cells: {tcell_count}, Activated T cells: {activated_count}, Dynabeads: {dynabead_count}")
                print(f"Normal T Cell - Avg Volume: {normal_tcell_volume:.2f}, Avg Pressure: {normal_tcell_pressure:.2f}")
                print(f"Activated T Cell - Avg Volume: {activated_tcell_volume:.2f}, Avg Pressure: {activated_tcell_pressure:.2f}")
    
    def calculate_cell_data(self, cell_type):
        """Calculate average volume and pressure for cells of a specific type"""
        total_volume = 0
        total_pressure = 0
        cell_count = 0
        
        for cell in self.cell_list_by_type(cell_type):
            total_volume += cell.volume
            
            # Calculate pressure as difference between current and target volume
            lambda_vol = cell.lambdaVolume
            target_vol = cell.targetVolume
            current_vol = cell.volume
            pressure = lambda_vol * 2 * (current_vol - target_vol)  # Derivative of λ(V-V_t)²
            
            total_pressure += abs(pressure)  # Using absolute value for pressure
            cell_count += 1
        
        # Calculate averages
        avg_volume = total_volume / cell_count if cell_count > 0 else 0
        avg_pressure = total_pressure / cell_count if cell_count > 0 else 0
        
        return avg_volume, avg_pressure
    
    def calculate_contact_energy(self):
        '''Calculate contact energy based on cell-cell interfaces'''
        contact_energy = 0
        processed_interfaces = set()
        
        for cell in self.cell_list:
            cell_id = cell.id
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    # Create a unique identifier for this cell-cell interface
                    # Sort IDs to ensure the same interface is identified regardless of order
                    interface_id = tuple(sorted([cell_id, neighbor.id]))
                    
                    # Skip if we've already processed this interface
                    if interface_id in processed_interfaces:
                        continue
                    
                    processed_interfaces.add(interface_id)
                    
                    # Get energy from our table
                    type1 = cell.type
                    type2 = neighbor.type
                    energy = self.contact_energy_table[type1][type2]
                    contact_energy += energy * common_surface_area
                else:
                    # Cell-medium interface
                    type1 = cell.type
                    type2 = self.MEDIUM_ID
                    energy = self.contact_energy_table[type1][type2]
                    contact_energy += energy * common_surface_area
        
        return contact_energy
        
    def calculate_volume_energy(self):
        '''Calculate the volume constraint energy component of the Hamiltonian'''
        volume_energy = 0.0
        
        for cell in self.cell_list:
            lambda_vol = cell.lambdaVolume
            target_vol = cell.targetVolume
            current_vol = cell.volume
            volume_energy += lambda_vol * (current_vol - target_vol)**2
                
        return volume_energy
        
    def calculate_surface_energy(self):
        '''Calculate the surface constraint energy component of the Hamiltonian'''
        surface_energy = 0.0
        
        for cell in self.cell_list:
            lambda_surf = cell.lambdaSurface
            target_surf = cell.targetSurface
            current_surf = cell.surface
            surface_energy += lambda_surf * (current_surf - target_surf)**2
                
        return surface_energy
        
    def finish(self):
        # Close the files when simulation ends
        if self.output_file:
            self.output_file.close()
            print("Energy components data saved to energy_file.txt")
            
        if self.data_file:
            self.data_file.close()
            print("Cell data saved to cell_data.txt")
            
    def on_stop(self):
        self.finish()