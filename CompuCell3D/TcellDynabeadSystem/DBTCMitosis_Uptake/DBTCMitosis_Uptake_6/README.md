**NutrientFieldSteppable**

  base_uptake_rate	0.02
	activated_uptake_multiplier 	1.5
	nutrient_replenishment_rate 	0.04
	max_nutrient_concentration 	100
	minimum_nutrient_threshold	0

**TcellGrowthSteppable**
  
  nutrient_consumption_threshold 	0.005	
	max_growth_nutrients 	0.008
	activated_growth_multiplier	1
	base_growth_rate 	0.0025
	max_growth_rate 	0.005555556

**TCellMitosisSteppable**
  nutrient_division_threshold 	0.03
	max_divisions 	10
 
**DiffusionSolverFE parameters**

  GlobalDiffusionConstant	0.1
	GlobalDecayConstant	0.001
	InitialConcentrationExpression	100
	Secretion	0.01
