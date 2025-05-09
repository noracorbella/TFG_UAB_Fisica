
        self.nutrient_field = None
        self.base_uptake_rate = 0.1  # Base uptake rate for regular T cells
        self.activated_uptake_multiplier = 5.0  # Activated T cells consume more nutrients
        self.nutrient_replenishment_rate = 0.05  # Add nutrient replenishment
        self.max_nutrient_concentration = 1.0  # Maximum nutrient concentration
        # Track cell level scalar attribute for visualization
        self.track_cell_level_scalar_attribute(field_name='uptake', attribute_name='nutrients_consumed')

        # Nutrient-dependent growth parameters
        self.nutrient_consumption_threshold = 0.01  # Minimum nutrients needed for growth
        self.max_growth_nutrients = 0.1  # Nutrients needed for maximum growth rate
        self.base_growth_rate = 1.0 / 80.0  # Base growth rate (pixels/MCS) - slower than original
        self.max_growth_rate = 1.0 / 40.0  # Maximum growth rate with optimal nutrients
