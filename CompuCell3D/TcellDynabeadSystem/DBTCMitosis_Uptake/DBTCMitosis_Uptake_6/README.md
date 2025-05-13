        self.base_uptake_rate = 0.02  # Base uptake rate for regular T cells
        self.activated_uptake_multiplier = 1.5  # Activated T cells consume more nutrients
        self.nutrient_replenishment_rate = 0.03  # Add nutrient replenishment
        self.max_nutrient_concentration = 1.0 # Maximum nutrient concentration
        # Track cell level scalar attribute for visualization
        
        # Nutrient-dependent growth parameters
        self.nutrient_consumption_threshold = 0.005  # Minimum nutrients needed for growth
        self.max_growth_nutrients = 0.08  # Nutrients needed for maximum growth rate
        self.base_growth_rate = 1.0 / 400.0  # Base growth rate (pixels/MCS) - slower than original
        self.max_growth_rate = 1.0 / 180.0  # Maximum growth rate with optimal nutrients

        self.nutrient_division_threshold = 0.03  # Minimum recent nutrient consumption for division
        self.max_divisions = 10
