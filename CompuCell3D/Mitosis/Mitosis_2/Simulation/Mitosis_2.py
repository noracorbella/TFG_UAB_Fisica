
from cc3d import CompuCellSetup
        


from Mitosis_2Steppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from Mitosis_2Steppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))


CompuCellSetup.run()
