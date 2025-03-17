
from cc3d import CompuCellSetup
        


from TcellsDynabeadSystemSteppables import CellInitialiser
from TcellsDynabeadSystemSteppables import TCellGrowthSteppable
from TcellsDynabeadSystemSteppables import TCellMitosisSteppable
from TcellsDynabeadSystemSteppables import DiffusionCalculatorSteppable

CompuCellSetup.register_steppable(steppable=CellInitialiser(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellGrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellMitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=DiffusionCalculatorSteppable(frequency=1))



CompuCellSetup.run()
