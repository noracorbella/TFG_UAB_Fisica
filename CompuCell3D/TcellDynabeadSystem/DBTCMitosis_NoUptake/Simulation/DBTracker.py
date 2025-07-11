
from cc3d import CompuCellSetup
        

from DBTCTrackerSteppables import CellInitialiser
from DBTCTrackerSteppables import TrajectoryTrackerSteppable
from DBTCTrackerSteppables import TcellGrowthSteppable
from DBTCTrackerSteppables import TCellMitosisSteppable
from DBTCTrackerSteppables import TCellActivationSteppable



CompuCellSetup.register_steppable(steppable=CellInitialiser(frequency=1))
CompuCellSetup.register_steppable(steppable=TrajectoryTrackerSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TcellGrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellMitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellActivationSteppable(frequency=1))



CompuCellSetup.run()
