from cc3d.core.PySteppables import *
import numpy as np
import random




class CellInitialiser(SteppableBasePy):
    def __init__(self,frequency=10):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
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

     
        
        

class TrajectoryTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.trajectory_data = None
    
    def start(self):
        self.trajectory_data = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\ExperimentalData\TcellsDynabeadSystem\TcellsDynabeadSystem_DiffusionTracker\cell_positions.txt", "w")
        self.trajectory_data.write("MCS \t CellId \t CellType \t Coordx \t Coordy \n")
        
    def step(self, mcs):
        if mcs % 100 == 0:
            for cell in self.cell_list:
                self.trajectory_data.write(f"{mcs} \t {cell.id} \t {cell.type} \t {cell.xCOM} \t {cell.yCOM} \n")
                    
            self.trajectory_data.flush()
            
    def finish(self):
        self.trajectory_data.close()
       