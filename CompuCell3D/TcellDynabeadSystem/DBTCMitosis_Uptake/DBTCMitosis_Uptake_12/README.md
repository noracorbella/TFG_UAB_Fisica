**NutrientFieldSteppable**	
  
  	base_uptake_rate	0.001
	activated_uptake_multiplier 	1
	nutrient_replenishment_rate 	0.0001
	max_nutrient_concentration 	10
	minimum_nutrient_threshold	0.2

**TcellGrowthSteppable**	
  
  	nutrient_consumption_threshold 	0
	max_growth_nutrients 	0.05
	activated_growth_multiplier	1
	base_growth_rate 	0.005
	max_growth_rate 	0.05

**TCellMitosisSteppable**	

  	nutrient_division_threshold 	0
	max_divisions 	5

**DiffusionSolverFE parameters**	
  
	GlobalDiffusionConstant	0.1
	GlobalDecayConstant	0.001
	InitialConcentrationExpression	10
	Secretion	0
