#Aquest codi em dona error. AttributeError: Could not find NeighborTrackerPlugin

from cc3d.core.PySteppables import *
import numpy as np

class cellsortingSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        
        SteppableBasePy.__init__(self,frequency)
        
        #Creem les variables (equivalent a fer a=0 o lst=[])
        self.energy_file = None
        self.contact_energy_file = None
        #self.neighbor_count_file = None
        
    def start(self):
        self.energy_file = open("/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting_0/cellsorting/energy_log.txt", "w")
        self.contact_energy_file = open("/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting_0/cellsorting/contact_energy_log.txt", "w")
        #self.neighbor_count_file = open("/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting_0/cellsorting/neighbor_count_log.txt", "w")
        print("Logging done.")
        
    
    def step(self, mcs):
        
        total_energy = self.simulator.getPotts().getEnergy()
        self.energy_file.write(f"{mcs} \t {total_energy} \n") 
        print(f"MCS {mcs}: Total Energy {total_energy}")
        
        total_contact_energy = 0
        for cell in self.cell_list:
            for neighbor, common_area in self.getCellNeighborDataList(cell):
                if neighbor:
                    cell_type = cell.type
                    neighbor_type = neighbor.type
                    energy = self.contactEnergy(cell_type, neighbor_type) 
                    total_contact_energy += energy * common_area  

        self.contact_energy_file.write(f"{mcs} \t {total_contact_energy} \n")
        print(f"MCS {mcs}: Total Contact Energy {total_contact_energy}")

    def finish(self):
        if self.energy_file:
            self.energy_file.close()
        if self.contact_energy_file:
            self.contact_energy_file.close()
        #if self.neighbor_count_file:
        #    self.neighbor_count_file.close()
        print("Log files closed.")
            
    def on_stop(self):
        if self.energy_file:
            self.energy_file.close()
        if self.contact_energy_file:
            self.contact_energy_file.close()
        #if self.neighbor_count_file:
        #    self.neighbor_count_file.close()
        print("Simulation stopped. Log files closed.")
        
        
        
        
        '''        
        sorted_neigbors = 0
        for cell in self.cell_list:
            for neighbor, area in self.getCellNeighborDataList(cell): #List of all neighbors for the current cell
            #neighbor is the neighboring cell
            #area is the area of the shared boundary between cell and neighbor
            
                if neighbor and cell.type != neighbor.type:
                #neighbor = none vol dir que estem al límit.
                    sorted_neigbors += 1 #si les dues celules son diferents vol dir que no estan sorted per tant sumem al numero de neighbors
        self.neighbor_count_file.write(f"{mcs} \t {sorted_neighbors} \n")
        print(f"MCS {mcs}: Total Sorted Neighbors {sorted_neighbors}")
'''
