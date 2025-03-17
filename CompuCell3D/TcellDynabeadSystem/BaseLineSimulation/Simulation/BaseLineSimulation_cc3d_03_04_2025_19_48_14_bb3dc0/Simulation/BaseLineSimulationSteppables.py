from cc3d.core.PySteppables import *
import numpy as np
import random
from random import gauss



class CellInitialiser(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        Nt = 12 #number of Tcells
        Nd = Nt #number of dynabeads
        Tcell_size = 8
        Dynabead_size = 5
        
        for _ in range(Nt):
            x, y = self.get_random_position(Tcell_size)
            cell = self.new_cell(self.TCELL)
            self.cell_field[x:x + Tcell_size, y:y + Tcell_size, 0] = cell
            cell.targetVolume = Tcell_size ** 2
            cell.lambdaVolume = 10.0
            cell.targetSurface = 4 * (cell.targetVolume ** (1/2))       
            cell.lambdaSurface = 10.0

        for _ in range(Nd):
            x, y = self.get_random_position(Dynabead_size)
            cell = self.new_cell(self.DYNABEAD)
            self.cell_field[x:x + Dynabead_size, y:y + Dynabead_size, 0] = cell
            cell.targetVolume = Dynabead_size ** 2
            cell.lambdaVolume = 20.0
            cell.targetSurface = 4 * (cell.targetVolume ** (1/2))
            cell.lambdaSurface = 20.0
        
        
    def get_random_position(self, cell_size):
        x = random.randint(0, self.dim.x - cell_size)
        y = random.randint(0, self.dim.y - cell_size)
        return x, y
            



class TCellMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def step(self, mcs):
        for cell in self.cell_list_by_type(self.TCELL):
            cell.lambdaVecX = 10.0 * random.uniform(-0.5,0.5)
            cell.lambdaVecY = 10.0 * random.uniform(-0.5,0.5)

            if cell.volume < 100:
                cell.targetVolume += 0.1
                cell.targetSurface = 4 * (cell.targetVolume ** (1/2))

class DynabeadMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        brownian_strength = 100.0  # Controls random motion strength
        
        for cell in self.cell_list_by_type(2):  # Dynabead Type = 2
            dx = random.uniform(-brownian_strength, brownian_strength)
            dy = random.uniform(-brownian_strength, brownian_strength)

            # Move the Dynabead by applying the random displacement
            cell.lambdaVecX = dx
            cell.lambdaVecY = dy
            
            
class TCellMitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)
        
    def step(self, mcs):
        for cell in self.cell_list_by_type(self.TCELL):
            if cell.volume>=64:
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.DYNABEAD:
                        self.divide_cell_along_major_axis(cell)
                        break

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()            

  
'''        
class DynabeadMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.temperature = 310 #K
        self.viscosity = 0.001 #Pa·s
        self.bead_radius = 1e-6 #m
        self.timestep_seconds = 1  
        
        
    def step(self, mcs):
        D = self.calculate_diffusion_coefficient()
        for cell in self.cell_list_by_type(self.DYNABEAD):
            dx = sqrt(2 * D * self.timestep_seconds) * gauss(0, 1)
            dy = sqrt(2 * D * self.timestep_seconds) * gauss(0, 1)
            displacement_vector = [int(dx), int(dy), 0]
            self.move_cell(cell, displacement_vector)
            
    def calculate_diffusion_coefficient(self):
        k_B = 1.380649e-23
        return k_B * self.temperature / (6 * np.pi * self.viscosity * self.bead_radius)





'''














    
    