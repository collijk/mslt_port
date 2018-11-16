class TobaccoFreeGeneration:
    def __init__(self):
        self.exposure = 'tobacco'

    def setup(self, builder):
        self.year = builder.configuration['tobacco_free_generation'].year
        self.clock = builder.time.clock()
        rate_name = '{}_intervention.incidence'.format(self.exposure)
        builder.value.register_value_modifier(rate_name, self.adjust_rate)

    def adjust_rate(self, index, rates):
        this_year = self.clock().year
        if this_year >= self.year:
            return 0.0 * rates
        else:
            return rates


class TobaccoEradication:
    def __init__(self):
        self.exposure = 'tobacco'

    def setup(self, builder):
        self.year = builder.configuration['tobacco_eradication'].year
        self.clock = builder.time.clock()
        rate_name = '{}_intervention.incidence'.format(self.exposure)
        builder.value.register_value_modifier(rate_name, self.adjust_rate)
        self.prevalence_col = '{}_intervention.yes'.format(self.exposure)
        self.remission_col = '{}_intervention.0'.format(self.exposure)
        view_columns = [self.prevalence_col, self.remission_col]
        self.population_view = builder.population.get_view(view_columns)
        # Add a handler to cease all current tobacco use.
        builder.event.register_listener('time_step__prepare',
                                        self.stop_current_use)

    def adjust_rate(self, index, rates):
        this_year = self.clock().year
        if this_year >= self.year:
            return 0.0 * rates
        else:
            return rates

    def stop_current_use(self, event):
        """
        Make all current tobacco users cease their tobacco use.
        """
        this_year = self.clock().year
        if this_year != self.year:
            return
        pop = self.population_view.get(event.index)
        pop.loc[:, self.remission_col] = (pop.loc[:, self.remission_col] +
                                          pop.loc[:, self.prevalence_col])
        pop.loc[:, self.prevalence_col] = 0.0
        self.population_view.update(pop)
