from cc3d.core.PySteppables import *
import numpy as np

class CellSortingSteppable(SteppableBasePy):

    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self,frequency)
        self.prev_com = {}  # Store center of mass of cells (the previous one)
        self.vel_lst = []  # Store velocities
        self.mcs_lst = []  # Store MCS values

    def start(self):
        self.contact_areabytype_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\CellSorting\Velocities_0\contact_areabytype_data.txt", "w")
        self.contact_areawithmedium_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\CellSorting\Velocities_0\contact_areabymedium_data.txt", "w")
        self.velocity_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\CellSorting\Velocities_0\velocity_file.txt", "w")
        
        self.contact_areabytype_file.write("MCS \t ACC \t ACN \t ANN \n")s
        self.contact_areawithmedium_file.write("MCS \t ACM \t ANM \n")
        self.velocity_file.write("MCS \t AverageVelocity \n")
        
        self.plot_win = self.add_new_plot_window(title='Contact Area by Type',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Contact Area', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        self.plot_win2 = self.add_new_plot_window(title='Contact Area with Medium',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Contact Area', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)                                                    
        
        self.plot_win.add_plot("CondensingCondensing", style='Lines', color='blue', size=5)
        self.plot_win.add_plot("CondensingNonCondensing", style='Lines', color='red', size=5)
        self.plot_win.add_plot("NonCondensingNonCondensing", style='Lines', color='green', size=5)
        self.plot_win2.add_plot("CondensingMedium", style='Lines', color='blue', size=5)
        self.plot_win2.add_plot("NonCondensingMedium", style='Lines', color='green', size=5)

        self.plot_velocity = self.add_new_plot_window(title = 'Vel vs MCS', x_axis_title = 'MCS', y_axis_title = 'Avg vel', x_scale_type = 'linear', y_scale_type = 'linear', grid = True)
        self.plot_velocity.add_plot('Velocity', style='Lines', color='blue', size=5)

    def step(self, mcs):
        total_vel = 0
        c = 0
        
        for cell in self.cell_list:
            print(f"MCS: {mcs}, Cell ID: {cell.id}, COM: ({cell.xCOM}, {cell.yCOM})")
            cm = np.array([cell.xCOM, cell.yCOM]) #Current CM
        
            if cell.id in self.prev_com:
                velocity = np.linalg.norm(cm - self.prev_com[cell.id])
                total_vel += velocity
                c += 1
            
            self.prev_com[cell.id] = cm
            
        avg_vel = total_vel / c if c > 0 else 0
        self.vel_lst.append(avg_vel)
        self.mcs_lst.append(mcs)
        self.velocity_file.write(f"{mcs} \t {avg_vel} \n")

        ACC=0.0
        ACN=0.0
        ANN=0.0
        ACM=0.0
        ANM=0.0
        for cell in self.cell_list: #Loop over all cells
            if cell: #this test fails if the cell is of type medium
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor: #this test fails if the neighbor is of type medium
                        if cell.type == self.CONDENSING: #Cell types referenced in ALL CAPS
                            if neighbor.type == self.CONDENSING: #in this case both cells are type dark
                                ACC+=common_surface_area
                            else: # in this case neighbor type is light
                                ACN+=common_surface_area
                        else: #in this case cell type is light 
                            if neighbor.type == self.NONCONDENSING: #in this case both cells are type light
                                ANN+=common_surface_area
                            else: # in this case neighbor type is dark
                               ACN+=common_surface_area                   
                    else: #in this case neighbor is medium
                        if cell.type == self.CONDENSING: #Then Dark-Medium Contact
                            ACM+=common_surface_area
                        else: # then Light-Medium Contact
                            ANM+=common_surface_area
        #Correct for Double Counting of cell-cell contact
        ACC/=2.0
        ACN/=2.0
        ANN/=2.0
        #Plot the values of the five contacts
        self.plot_win.add_data_point("CondensingCondensing", mcs, ACC)
        self.plot_win.add_data_point("CondensingNonCondensing", mcs, ACN)
        self.plot_win.add_data_point("NonCondensingNonCondensing", mcs, ANN)
        self.plot_win2.add_data_point("CondensingMedium", mcs, ACM)
        self.plot_win2.add_data_point("NonCondensingMedium", mcs, ANM) 
        self.plot_velocity.add_data_point("Velocity", mcs, avg_vel)
          
        self.contact_areabytype_file.write(f"{mcs} \t {ACC} \t {ACN} \t {ANN} \n")
        self.contact_areawithmedium_file.write(f"{mcs} \t {ACM} \t {ANM} \n")
          
            
    def finish(self):
        self.contact_areabytype_file.close()
        self.contact_areawithmedium_file.close()
        self.velocity_file.close()
   
    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
