from cc3d import CompuCellSetup
        


from TcellsDynabeadSystem_DiffusionTrackerSteppables import CellInitialiser
from TcellsDynabeadSystem_DiffusionTrackerSteppables import TrajectoryTrackerSteppable


CompuCellSetup.register_steppable(steppable=CellInitialiser(frequency=1))
CompuCellSetup.register_steppable(steppable=TrajectoryTrackerSteppable(frequency=1))


CompuCellSetup.run()
