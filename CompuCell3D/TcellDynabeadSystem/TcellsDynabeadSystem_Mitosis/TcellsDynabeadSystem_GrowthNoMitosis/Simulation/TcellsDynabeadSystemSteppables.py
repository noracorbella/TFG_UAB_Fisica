from cc3d.core.PySteppables import *
import numpy as np
import random




class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.prev_com = {} 
        self.Db_vel_lst = [] 
        self.TC_vel_lst = [] 
        self.mcs_lst = [] 
    def start(self):
        
        self.velocity_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\TcellsDynabeadSystem_RealisticData\velocity_file.txt", "w")
        self.velocity_file.write("MCS \t Avg Db vel \t Avg TC vel \n")
        self.velocity_file.flush()
        
        self.plot_win = self.add_new_plot_window(title = "Velocity", x_axis_title = "MCS", y_axis_title = "Average velocity", x_scale_type = "linear", y_scale_type = "linear", grid = True)
        self.plot_win.add_plot("Dynabead Velocity", style = "Lines", color = "green", size = 5)
        self.plot_win.add_plot("TCell Velocity", style = "Lines", color = "red", size = 5)
        
        Nt = 24 #number of Tcells
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



    def step(self, mcs):
        total_Db_vel = 0
        total_TC_vel = 0
        c_Db = 0
        c_TC = 0
        
        for Db in self.cell_list_by_type(self.DYNABEAD):
            cm_Db = np.array([Db.xCOM, Db.yCOM])
            
            if Db.id in self.prev_com:
                Db_vel = np.linalg.norm(cm_Db - self.prev_com[Db.id])
                total_Db_vel += Db_vel
                c_Db += 1
                
            self.prev_com[Db.id] = cm_Db
        
        avg_Db_vel = total_Db_vel / c_Db if c_Db > 0 else 0
        self.Db_vel_lst.append(avg_Db_vel)
            
            
        for TC in self.cell_list_by_type(self.TCELL):
            cm_TC = np.array([TC.xCOM, TC.yCOM])
            
            if TC.id in self.prev_com:
                TC_vel = np.linalg.norm(cm_TC - self.prev_com[TC.id])
                total_TC_vel += TC_vel
                c_TC += 1
                
            self.prev_com[TC.id] = cm_TC
        
        avg_TC_vel = total_TC_vel / c_TC if c_TC > 0 else 0
        self.TC_vel_lst.append(avg_TC_vel)
        self.mcs_lst.append(mcs)         
        
        self.velocity_file.write(f"{mcs} \t {avg_Db_vel} \t {avg_TC_vel}")
        
        self.plot_win.add_data_point("Dynabead Velocity", mcs, avg_Db_vel)
        self.plot_win.add_data_point("TCell Velocity", mcs, avg_TC_vel)
        
        
        
        
        
class DiffusionCalculatorSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.cell_positions = {"DYNABEAD": {}, "TCELL": {}}
        self.msd_values = {"DYNABEAD": [], "TCELL": []}
        self.diffusion_constants = {"DYNABEAD": [], "TCELL": []}
        #self.cell_types_lst = [self.DYNABEAD, self.TCELL]

    def start(self):
        self.diffusion_data = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\TcellsDynabeadSystem_RealisticData\Diffusion_Data.txt", "w")
        self.diffusion_data.write("MCS \t D_Dynabead \t D_TCell \n")
        self.diffusion_data.flush()

    def step(self, mcs):
        for cell_type in [self.DYNABEAD, self.TCELL]:
            type_name = "DYNABEAD" if cell_type == self.DYNABEAD else "TCELL"
            
            for cell in self.cell_list_by_type(cell_type):      
                cell_id = cell.id
         
                if cell_id not in self.cell_positions[type_name]:
                    self.cell_positions[type_name][cell_id] = (cell.xCOM, cell.yCOM) #cell's position is recorded only the first time it appears
                
                initial_position = np.array(self.cell_positions[type_name][cell_id])
                
                current_position = np.array([cell.xCOM, cell.yCOM])
                displacement_squared = np.sum((current_position - initial_position) ** 2)

                self.msd_values[type_name].append(displacement_squared)
        
        D_Dynabead = self.calculate_diffusion(self.msd_values["DYNABEAD"], mcs)
        D_TCell = self.calculate_diffusion(self.msd_values["TCELL"], mcs)  
           
        self.diffusion_constants["DYNABEAD"].append(D_Dynabead)
        self.diffusion_constants["TCELL"].append(D_TCell)
        
        self.diffusion_data.write(f"{mcs}\t{D_Dynabead}\t{D_TCell}\n")
        self.diffusion_data.flush()

        
    def calculate_diffusion(self, msd_values, mcs):
        if len(msd_values) > 1:
            time_steps = np.arange(1, len(msd_values) + 1)
            msds = np.array(msd_values)
              
            fit = np.polyfit(time_steps, msds, 1)
            D = fit[0] / 4
            return D
             
        return None
         
    def finish(self):          
        self.diffusion_data.close()
        
        
class TCellGrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.Tcell_maxsize = 30
        self.TCellVolume = []
    
    def start(self):
        self.TCellVolume_plot = self.add_new_plot_window(
            title="TCell Volume",
            x_axis_title="MCS",
            y_axis_title="Volume",
            x_scale_type="linear",
            y_scale_type="linear",
            grid=True
        )
        
        self.TCellVolume_plot.add_plot("Average Volume", style="Lines", color="blue", size=5)
        
    def step(self, mcs):
        total_volume = 0
        cell_count = 0
    
        for cell in self.cell_list_by_type(self.TCELL):
            
            max_volume = np.pi * (self.Tcell_maxsize/2) ** 2
            
            if cell.volume < max_volume:
                cell.targetVolume += 1
                cell.targetSurface =  2 * np.pi * np.sqrt(cell.targetVolume)
                cell.lambdaVolume += 1
                cell.lambdaSurface += 1
            
            total_volume += cell.volume
            cell_count += 1
        
        avg_volume = total_volume / cell_count if cell_count > 0 else 0
        
        self.TCellVolume_plot.add_data_point("Average Volume", mcs, avg_volume)
        self.TCellVolume.append(avg_volume)
        
        
        
        
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        