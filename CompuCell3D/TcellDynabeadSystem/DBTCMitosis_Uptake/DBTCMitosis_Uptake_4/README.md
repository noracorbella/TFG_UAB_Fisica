**NutrientFieldSteppable**

	base_uptake_rate	0.05
	activated_uptake_multiplier 	2
	nutrient_replenishment_rate 	0.09
	max_nutrient_concentration 	100
	minimum_nutrient_threshold	0
 
**TcellGrowthSteppable**

	nutrient_consumption_threshold 	0.01
	max_growth_nutrients 	0.1
	activated_growth_multiplier	1
	base_growth_rate 	0.008333333
	max_growth_rate 	0.033333333
 
**TCellMitosisSteppable**
        
	nutrient_division_threshold 	0.04
        max_divisions 	10
        
**DiffusionSolverFE parameters**

	GlobalDiffusionConstant	0.1
	GlobalDecayConstant	0.001
	InitialConcentrationExpression	100
	Secretion	0.01

