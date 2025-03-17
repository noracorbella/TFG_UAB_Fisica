
from cc3d import CompuCellSetup
from BaseLineSimulationSteppables import CellInitialiser, TCellMotilitySteppable, DynabeadMotilitySteppable, TCellMitosisSteppable

CompuCellSetup.register_steppable(steppable=CellInitialiser(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellMotilitySteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=DynabeadMotilitySteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=TCellMitosisSteppable(frequency=1))
CompuCellSetup.run()
