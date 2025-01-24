from cc3d import CompuCellSetup   
from cellsortingSteppables import cellsortingSteppable


CompuCellSetup.register_steppable(steppable=cellsortingSteppable(frequency=1))
CompuCellSetup.run()
