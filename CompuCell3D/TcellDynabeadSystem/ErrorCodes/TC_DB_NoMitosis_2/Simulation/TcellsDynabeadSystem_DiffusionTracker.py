
from cc3d import CompuCellSetup
        


from TcellsDynabeadSystem_DiffusionTrackerSteppables import CellInitialiser
#from TcellsDynabeadSystemSteppables import TCellGrowthSteppable
#from TcellsDynabeadSystemSteppables import TCellMitosisSteppable
#from TcellsDynabeadSystem_DiffusionTrackerSteppables import DiffusionCalculatorSteppable
from TcellsDynabeadSystem_DiffusionTrackerSteppables import TrajectoryTrackerSteppable


CompuCellSetup.register_steppable(steppable=CellInitialiser(frequency=1))
#CompuCellSetup.register_steppable(steppable=TCellGrowthSteppable(frequency=1))
#CompuCellSetup.register_steppable(steppable=TCellMitosisSteppable(frequency=1))
#CompuCellSetup.register_steppable(steppable=DiffusionCalculatorSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TrajectoryTrackerSteppable(frequency=1))


CompuCellSetup.run()
