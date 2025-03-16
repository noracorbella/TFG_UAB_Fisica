from cc3d.core.PySteppables import *
import numpy as np

class CellSortingSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        pass

    def step(self, mcs):
        pass
    def finish(self):
        pass
    def on_stop(self):
        pass
        
        
        
 
class ContactEnergyTracker(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        # Define cell types as constants
        self.MEDIUM_ID = 0
        self.TYPE_A = 1
        self.TYPE_B = 2
        
        # Create contact energy table based on XML configuration
        self.contact_energy_table = {
            self.MEDIUM_ID: {
                self.MEDIUM_ID: 0.0,
                self.TYPE_A: 16.0,
                self.TYPE_B: 16.0
            },
            self.TYPE_A: {
                self.MEDIUM_ID: 16.0,
                self.TYPE_A: 2.0,
                self.TYPE_B: 11.0
            },
            self.TYPE_B: {
                self.MEDIUM_ID: 16.0,
                self.TYPE_A: 11.0,
                self.TYPE_B: 16.0
            }
        }

        # Store the neighbor order
        self.neighbor_order = 4
        
    def start(self):
        self.contact_energy_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\CellSorting\Energy\CellSorting_Energy\contact_energy_data.txt", "w")
        self.contact_energy_file.write("MCS \t ContactEnergy \n")
    def step(self, mcs):

        contact_energy = 0
        for cell in self.cell_list:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    type1 = cell.type
                    type2 = neighbor.type
                else:
                    type1 = cell.type
                    type2 = self.MEDIUM_ID
                
                energy = self.contact_energy_table[type1][type2]
                contact_energy += energy * common_surface_area
        
        self.contact_energy_file.write(f"{mcs} \t {contact_energy}\n")
        self.contact_energy_file.flush()  # we flush the file after each write operation

        print(f"MCS {mcs}: Total Contact Energy = {contact_energy}")

    def finish(self):
        if hasattr(self, "contact_energy_file"): 
            '''
            hasattr(), function, returns "True" if the object 'self' has an attribute named "contact_energy_file", returns "False" otherwise
            '''
            self.contact_energy_file.close()


    def on_stop(self):
        self.finish()







