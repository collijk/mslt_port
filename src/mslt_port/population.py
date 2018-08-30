import numpy as np

from mslt_port.data import (get_base_population, get_all_cause_mortality,
                            get_yld_rate, AGE_GROUP_END)


class BasePopulation:

    def setup(self, builder):
        self.pop_data = get_base_population()

        columns = ['age', 'sex', 'population']
        builder.population.initializes_simulants(self.on_initialize_simulants, creates_columns=columns)
        self.population_view = builder.population.get_view(columns + ['tracked'])

        builder.event.register_listener('time_step', self.on_time_step, priority=6)

    def on_initialize_simulants(self, _):
        self.population_view.update(self.pop_data)

    def on_time_step(self, event):
        pop = self.population_view.get(event.index, query='tracked == True')
        pop['age'] += 1
        pop.loc[pop.age > AGE_GROUP_END, 'tracked'] = False
        self.population_view.update(pop)


class Mortality:

    def setup(self, builder):
        mortality_data = builder.lookup.build_table(get_all_cause_mortality())
        self.mortality_rate = builder.value.register_rate_producer('mortality_rate', source=mortality_data)

        builder.event.register_listener('time_step', self.on_time_step)

        self.population_view = builder.population.get_view(['population'])

    def on_time_step(self, event):
        pop = self.population_view.get(event.index)
        probability_of_death = 1 - np.exp(-self.mortality_rate(event.index))
        pop.population *= 1 - probability_of_death
        self.population_view.update(pop)


class Disability:

    def setup(self, builder):
        yld_rate = builder.lookup.build_table(get_yld_rate())
        builder.value.register_rate_producer('yld_rate', source=yld_rate)
