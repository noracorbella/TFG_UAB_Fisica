
from cc3d import CompuCellSetup
        



from Mitosis_3_1Steppables import VolumeParamSteppable
from Mitosis_3_1Steppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=VolumeParamSteppable(frequency=10))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=10))



CompuCellSetup.run()
