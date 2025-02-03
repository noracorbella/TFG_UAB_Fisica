from cc3d import CompuCellSetup
        

from CellSortingSteppables import CellSortingSteppable

CompuCellSetup.register_steppable(steppable=CellSortingSteppable(frequency=1))


CompuCellSetup.run()
