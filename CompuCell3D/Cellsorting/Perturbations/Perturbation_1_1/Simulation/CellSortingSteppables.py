from cc3d.core.PySteppables import *
import numpy as np

class CellSortingSteppable(SteppableBasePy):

    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self,frequency)
        self.prev_com = {}  # Store center of mass of cells (the previous one)
        self.vel_lst = []  # Store velocities
        self.mcs_lst = []  # Store MCS values
        self.movie_maker = None
        
    def start(self):     
        
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
            
            neighbor_types = []
            same_type_count = 0
            
            for neighbor, _ in self.get_cell_neighbor_data_list(cell):
                if neighbor: #ignore medium
                    neighbor_types.append(neighbor.type)
                    if neighbor.type == cell.type:
                        same_type_count += 1
                    
            if len(set(neighbor_types)) ==  1 and cell.type not in neighbor_types:
                cell.type = neighbor_types[0]
                
            elif same_type_count == 1:
                new_type = [t for t in neighbor_types if t != cell.type][0]
                cell.type = new_type
            
            
        avg_vel = total_vel / c if c > 0 else 0
        self.vel_lst.append(avg_vel)
        self.mcs_lst.append(mcs)
        
        if not mcs%10 and mcs>1: #Only plot every ten time steps
            """
            Calculate contact areas for each type pair
            """
            # Create placeholders to accumulate areas
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
       
    def finish(self):
        """
        Called after the last MCS to wrap up the simulation
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
