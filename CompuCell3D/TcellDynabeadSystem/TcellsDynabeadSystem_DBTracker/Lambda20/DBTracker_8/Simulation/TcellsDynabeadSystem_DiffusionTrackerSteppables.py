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
        
        Nt = 24 #number of Tcells
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

        self.num_atoms = len(self.cell_list)  # Assuming each cell is an atom
        
        # Create an empty Universe for MDAnalysis with the same number of cells
        self.universe = MDAnalysis.Universe.empty(self.num_atoms, trajectory=True)  # No bonds, just positions
        self.universe.add_TopologyAttr("name")

        # Initialize positions for all cells in the Universe
        for i, cell in enumerate(self.cell_list):
            self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]  # Assuming 2D positions, z=0 for simplicity
            self.universe.atoms[i].name = f"cell_{cell.id}"
            
        # Write initial structure to PDB file
        with PDBWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTracker\DBTracker_8\DB8_frame0.pdb") as pdb_writer:
            pdb_writer.write(self.universe.atoms)            
            
        # Initialize XTC trajectory writer
        self.trajectory_writer = XTCWriter(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTracker\DBTracker_8\DB8_traj.xtc", self.num_atoms)
            
            
        # Create plot window for volume tracking
        self.plot_win = self.add_new_plot_window(title="Avg Total Volume", x_axis_title="MCS", y_axis_title="V(pixels^2)", 
                                              x_scale_type="linear", y_scale_type="linear", grid=True)
        
        self.plot_win2 = self.add_new_plot_window(title="Avg DYN Volume", x_axis_title="MCS", y_axis_title="V(pixels^2)", 
                                              x_scale_type="linear", y_scale_type="linear", grid=True)
                                              
        self.plot_win.add_plot("Volume", style="Dots", color="blue", size=5)
        
        # Add separate plots for T cells and Dynabeads if you want to track them separately
        self.plot_win2.add_plot("Dynabead Volume", style="Dots", color="green", size=5)
        
    def step(self, mcs):
        # Update trajectory every 100 MCS
        if mcs % 100 == 0:         
            # Update positions for all cells
            for i, cell in enumerate(self.cell_list):
                self.universe.atoms[i].position = [cell.xCOM, cell.yCOM, 0.0]
            
            # Write frame to trajectory
            self.trajectory_writer.write(self.universe)
        
        # Calculate average volumes
        total_volume = 0
        dynabead_volume = 0
        dynabead_count = 0
        
        for cell in self.cell_list:
            total_volume += cell.volume
            
            # Track dynabead volumes separately
            if cell.type == self.DYNABEAD:
                dynabead_volume += cell.volume
                dynabead_count += 1
        
        # Calculate average volumes
        avg_total_volume = total_volume / len(self.cell_list) if len(self.cell_list) > 0 else 0
        avg_dynabead_volume = dynabead_volume / dynabead_count if dynabead_count > 0 else 0
        
        # Add data points to plot
        self.plot_win.add_data_point("Volume", mcs, avg_total_volume)
        self.plot_win2.add_data_point("Dynabead Volume", mcs, avg_dynabead_volume)
    
    def finish(self):
        # Close the trajectory writer when simulation ends
        if self.trajectory_writer:
            self.trajectory_writer.close()
            print("Trajectory saved successfully")




class EnergyTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.output_file = None
        
        # Define cell types as constants
        self.MEDIUM_ID = 0
        self.DYNABEAD_ID = 2
        
        # Define contact energy table based on your XML configuration
        self.contact_energy_table = {
            self.MEDIUM_ID: {
                self.MEDIUM_ID: 0.0,
                self.DYNABEAD_ID: 2.0
            },
            self.DYNABEAD_ID: {
                self.MEDIUM_ID: 2.0,
                self.DYNABEAD_ID: 5.0
            }
        }

        # Store the neighbor order
        self.neighbor_order = 4
        
    def start(self):
        # Create and open file for writing energy components
        self.output_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\DBTracker\DBTracker_8\energy_components.txt", "w")
        # Write header
        self.output_file.write("MCS\tContactEnergy\tVolumeEnergy\tSurfaceEnergy\tTotalEnergy\n")
        
    def step(self, mcs):
        # Calculate energy components using the working approach
        contact_energy = self.calculate_contact_energy_new()
        volume_energy = self.calculate_volume_energy()
        surface_energy = self.calculate_surface_energy()
        total_energy = contact_energy + volume_energy + surface_energy
        
        # Write to file
        self.output_file.write(f"{mcs}\t{contact_energy}\t{volume_energy}\t{surface_energy}\t{total_energy}\n")
        
        # Flush to ensure data is written immediately
        self.output_file.flush()
        
        if mcs % 1000 == 0:  # Print updates occasionally
            print(f"MCS: {mcs}, Total Energy: {total_energy}")
    
    def calculate_contact_energy_new(self):
        '''Calculate contact energy'''
        contact_energy = 0
        
        for cell in self.cell_list:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    type1 = cell.type
                    type2 = neighbor.type
                else:
                    type1 = cell.type
                    type2 = self.MEDIUM_ID
                
                # Get energy from our table
                energy = self.contact_energy_table[type1][type2]
                contact_energy += energy * common_surface_area
        
        # Should the contact energy be divided by 2?
        return contact_energy
        
    def calculate_volume_energy(self):
        """Calculate the volume constraint energy component of the Hamiltonian"""
        volume_energy = 0.0
        
        for cell in self.cell_list:
            if cell:
                lambda_vol = cell.lambdaVolume
                target_vol = cell.targetVolume
                current_vol = cell.volume
                volume_energy += lambda_vol * (current_vol - target_vol)**2
                
        return volume_energy
        
    def calculate_surface_energy(self):
        """Calculate the surface constraint energy component of the Hamiltonian"""
        surface_energy = 0.0
        
        for cell in self.cell_list:
            if cell:
                lambda_surf = cell.lambdaSurface
                target_surf = cell.targetSurface
                current_surf = cell.surface
                surface_energy += lambda_surf * (current_surf - target_surf)**2
                
        return surface_energy
        
    def finish(self):
        # Close the file when simulation ends
        if self.output_file:
            self.output_file.close()
            print("Energy components data saved to energy_components.txt")
            
    def on_stop(self):
        self.finish()