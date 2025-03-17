
from cc3d import CompuCellSetup
        

from Mitosis_4_1Steppables import VolumeParamSteppable
from Mitosis_4_1Steppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=VolumeParamSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))


CompuCellSetup.run()
