from cc3d.core.PySteppables import *
import numpy as np
from math import pi


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for cell in self.cell_list:

            cell.targetVolume = 25
            cell.lambdaVolume = 3.0
            cell.targetSurface = 20
            cell.lambdaSurface = 3.0
#lambdaVolume and lambdaSurface should be of the same magnitud so the importance of the constraints is comparable
#If one is much bigger than the other it will dominate the growth

class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
    
        for cell in self.cell_list:
            cell.targetVolume += 25/100 #the cell is going to grow 25 pixels every 100 MCS        
            cell.targetSurface = 2*pi*sqrt(cell.targetVolume)
            #we assume that the cells occupy the area of a circle and we calculate the surface from the volume
            
            
            
        # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
            # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # you can use here any fcn of concentrationAtCOM
            # cell.targetVolume += 0.01 * concentrationAtCOM       

        