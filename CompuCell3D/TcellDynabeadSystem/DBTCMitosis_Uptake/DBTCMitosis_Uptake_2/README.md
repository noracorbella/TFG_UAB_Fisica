First official simulation on the T cell Dynabead system. Trying to find out what's the best way to save the data for posterior analysis.

The simulation parameters are the following

**NutrientFieldSteppable**	
  base_uptake_rate	0.1
	activated_uptake_multiplier 	5
	nutrient_replenishment_rate 	0.06
	max_nutrient_concentration 	100
	minimum_nutrient_threshold	0
**TcellGrowthSteppable**	
  nutrient_consumption_threshold 	0.01
	max_growth_nutrients 	0.1
	activated_growth_multiplier	1
	base_growth_rate 	0.05
	max_growth_rate 	0.10
**TCellMitosisSteppable**
  nutrient_division_threshold 	0.05
	max_divisions 	infty
**.xml parameters**	
  GlobalDiffusionConstant	0.1
	GlobalDecayConstant	0.001
	InitialConcentrationExpression	100
	Secretion	0.01
