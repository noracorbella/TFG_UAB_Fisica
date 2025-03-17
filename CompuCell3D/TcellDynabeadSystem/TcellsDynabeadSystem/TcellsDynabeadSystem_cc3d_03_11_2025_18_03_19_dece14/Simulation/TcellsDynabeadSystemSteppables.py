from cc3d.core.PySteppables import *
import numpy as np
import random




class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
        Nt = 12 #number of Tcells
        Nd = Nt #number of dynabeads
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


class DiffusionCalculatorSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.cell_positions = {}
        self.time_steps = []
        self.msd_values = []
    def start(self):
        self.diffusion_data = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\Diffusion_data.txt", "w")
        self.diffusion_data.write("Prova")
        self.diffusion_data.flush()
        
    def step(self, mcs):
        for cell in self.cell_list_by_type(self.DYNABEAD):
            cell_id = cell.id
         
            if cell_id not in self.cell_positions:
                self.cell_positions[cell_id] = (cell.xCOM, cell.yCOM)   
            
            initial_position = np.array(self.cell_positions[cell_id])
            current_position = np.array([cell.xCOM, cell.yCOM])
            displacement_squared = np.sum((current_position - initial_position) ** 2)

            self.msd_values.append(displacement_squared)
            self.time_steps.append(mcs)
  
    def Diffusion(self):
        if len(self.time_steps) > 1:
            time_steps = np.array(self.time_steps)
            msds = np.array(self.msd_values)
            
            fit = np.polyfit(time_steps, msds, 1)
            D = fit[0] / 4
            return D
            
        return None
           
            
    def finalise(self):
        D = self.Diffusion()
        if D is not None:
            print(f"Diffusion constant = {D}")
            
        self.diffusion_data.write(f"Diffusion constant = {D}")


























'''
class TCellGrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.Tcell_maxsize = 12 

        
    def step(self, mcs):
    
        for cell in self.cell_list_by_type(self.TCELL):
            
            max_volume = np.pi * (self.Tcell_maxsize/2) ** 2
            
            if cell.volume < max_volume:
                cell.targetVolume += 0.1    

                
                cell.targetSurface =  2 * np.pi * np.sqrt(cell.targetVolume)
                cell.lambdaSurface = 100.0
                
                
                
                
class TCellMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        self.Tcell_minsize = 10 

    def step(self, mcs):
        
        
        for cell in self.cell_list_by_type(self.TCELL):
            
            min_volume = np.pi * (self.Tcell_minsize/2) ** 2
            
            if cell.volume >= min_volume:
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.DYNABEAD:    
                        self.divide_cell_random_orientation(cell)
                        break
        
                    # Other valid options
                    # self.divide_cell_orientation_vector_based(cell,1,1,0)
                    # self.divide_cell_along_major_axis(cell)
                    # self.divide_cell_along_minor_axis(cell)

        
    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0                  
        self.clone_parent_2_child()            


'''