from cc3d.core.PySteppables import *
import  numpy as np

class VolumeParamSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.volume_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\Mitosis\Mitosis_1\volume_file.txt", "w")
        self.volume_file.write("MCS \t Volume \n")
        self.pressure_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\Mitosis\Mitosis_1\pressure_file.txt", "w")
        self.pressure_file.write("MCS \t Pressure \n")
        
        
        cell_0 = self.new_cell(self.CONDENSING)

        self.cell_field[97:102,97:102,0] = cell_0

        for cell in self.cell_list:
            cell.targetVolume = 25
            cell.lambdaVolume = 2.0
            cell.dict['pressure'] = 0.0

        self.plot_win = self.add_new_plot_window(title='Average Volume',
                                                     x_axis_title='MonteCarlo Step (MCS)',
                                                     y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                     grid=True)
            
        

        self.plot_win2 = self.add_new_plot_window(title='Average Pressure',
                                         x_axis_title='MonteCarlo Step (MCS)',
                                         y_axis_title='Pressure', x_scale_type='linear', y_scale_type='linear',
                                         grid=True)
                                            
        
        self.plot_win.add_plot("AVol", style='Dots', color='red', size=5)
        self.plot_win2.add_plot("APre", style='Dots', color='red', size=5)

    def step(self, mcs):
        AVol = 0 
        APre = 0
        L = len(self.cell_list)
        for cell in self.cell_list:
            cell.targetVolume += 1  #Gradually increases the target volume of each cell by 1 unit per Monte Carlo Step (MCS).
            AVol += cell.volume/L
            cell.dict['pressure'] = 2*cell.lambdaVolume*(cell.targetVolume-cell.volume)
            APre += cell.dict['pressure']/L


        self.plot_win.add_data_point("AVol", mcs, AVol)
        self.plot_win2.add_data_point("APre", mcs, APre)
        
        self.volume_file.write(f"{mcs} \t {AVol}\n")
        self.pressure_file.write(f"{mcs} \t {APre}\n")
    
    
    def finish(self):
        if self.volume_file and self.pressure_file:
            self.pressure_file.close()
            self.volume_file.close()


  
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        
        # Optional: Customize where to place the new cell.
        # 0 - parent and child positions will be randomized between mitosis events
        # negative integer - parent appears on the 'left' of the child
        # positive integer - parent appears on the 'right' of the child
        self.set_parent_child_position_flag(-1)
        

    def start(self):
        self.cell_count_file = open(r"C:\Users\norac\OneDrive - UAB\Escritorio\uab\5\TFGJordi\CC3D_Tutorials\Mitosis\Mitosis_1\cell_count.txt", "w")
        self.cell_count_file.write("MCS \t CellCount \n")
        print("Logging done.")


    def step(self, mcs):
             
        
        self.cell_count_file.write(f"{mcs} \t {len(self.cell_list)}\n")
        print(f"MCS {mcs}; Number of cells {self.cell_list} \n")

        cells_to_divide=[]
        for cell in self.cell_list:
            if cell.volume>50: #the critical mass, chosen arbitrarily
                cells_to_divide.append(cell)

        for cell in cells_to_divide:
            '''
            Loop over all cells_to_divide and perform mitosis with a built-in function
            '''
            # To change mitosis mode, leave one of the below lines uncommented

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        '''
        Automatically called BEFORE each division. It controls the data of the parent and child cells.
        '''
        # Reduce parent target volume BEFORE cloning
        
        self.parent_cell.targetVolume /= 2.0                  
        self.child_cell.targetVolume = self.parent_cell.targetVolume
        self.child_cell.lambdaVolume = self.parent_cell.lambdaVolume      

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        # Make the cell type change every time the cell divides
        if self.parent_cell.type == self.CONDENSING:
            self.child_cell.type = self.CONDENSING
        else:
            self.child_cell.type = self.CONDENSING
            
            
    def finish(self):
        if self.cell_count_file:
            self.cell_count_file.close()
        print("Log files closed.")

        