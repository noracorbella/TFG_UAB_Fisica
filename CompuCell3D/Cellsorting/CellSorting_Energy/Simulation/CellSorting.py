
from cc3d import CompuCellSetup
        
from CellSortingSteppables import CellSortingSteppable, ContactEnergyTracker

CompuCellSetup.register_steppable(steppable=CellSortingSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=ContactEnergyTracker(frequency=1))


CompuCellSetup.run()
