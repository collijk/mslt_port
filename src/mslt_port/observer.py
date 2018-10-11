"""This module provides classes that record various outputs of interest."""

import numpy as np
import pandas as pd
import itertools

from mslt_port.data import AGE_GROUP_END, YEAR_START, YEAR_END

class DiseaseObserver:

    def __init__(self, name, output_path):
        self.name = name
        self.output_path = output_path

    def setup(self, builder):
        columns = [f'{self.name}_S', f'{self.name}_S_previous',
                   f'{self.name}_C', f'{self.name}_C_previous',
                   f'{self.name}_S_intervention', f'{self.name}_S_intervention_previous',
                   f'{self.name}_C_intervention', f'{self.name}_C_intervention_previous']
        self.population_view = builder.population.get_view(columns)

        builder.event.register_listener('collect_metrics', self.on_collect_metrics)
        builder.event.register_listener('simulation_end', self.write_output)

        idx = pd.MultiIndex.from_product([range(AGE_GROUP_END + 1), range(YEAR_START, YEAR_END + 1), ['male', 'female']])
        self.data = pd.DataFrame(columns=['S', 'C', 'S_int', 'C_int'], index=idx)

    def on_collect_metrics(self, event):
        pass

    def write_output(self, event):
        pass


class AdjustedPYandLE:
    """
    This class calculates the adjusted person-years and the adjusted
    life-expectancy for each cohort at each year of the simulation.
    """

    def __init__(self, output_file=None):
        """
        :param output_file: The name of the CSV file in which to record the
            adjusted person-years and adjusted life-expectancy data.
        """
        self.output_file = output_file

    def setup(self, builder):
        """
        """
        self.yld_rate = builder.value.get_value('yld_rate')
        self.acm_rate = builder.value.get_value('mortality_rate')
        view_cols = ['age', 'sex', 'population']
        # TODO: ugly hack, can't extract mortality rate when initialising
        #       simulants unless the diseases have already been initialised.
        extra_cols = ['chd_S']
        # extra_cols = []
        req_cols = view_cols + extra_cols
        builder.population.initializes_simulants(self.on_initialize,
                                                 requires_columns=req_cols)
        self.population_view = builder.population.get_view(view_cols)
        self.clock = builder.time.clock()
        builder.event.register_listener('collect_metrics', self.on_collect_metrics)
        builder.event.register_listener('simulation_end', self.finalise_output)
        self.idx_cols = ['age', 'sex', 'year']

    def create_table(self, min_age, min_year):
        """
        Create an empty data frame to hold all of the adjusted person-year and
        adjusted life-expectancy values that will be produced during a
        simulation.
        """
        ages = range(min_age, AGE_GROUP_END + 1)
        sexes = ['male', 'female']
        years = range(min_year, YEAR_END + 1)
        # NOTE: columns must be in the same order as self.idx_cols.
        rows = list(itertools.product(ages, sexes, years))
        self.data = pd.DataFrame(rows, columns=self.idx_cols)
        self.data.set_index(self.idx_cols)
        self.data['PYadj'] = np.nan
        self.data['population'] = np.nan

    def on_initialize(self, pop_data):
        """
        Calculate adjusted person-years and adjusted life-expectancy before
        the first time-step.
        """
        idx = pop_data.index
        # Record the year at which the simulation started.
        self.year_0 = self.clock().year
        acm_rate_now = self.acm_rate(idx)
        yld_rate_now = self.yld_rate(idx)
        prob_death = 1 - np.exp(- acm_rate_now)
        pop = self.population_view.get(idx)
        pop['year'] = self.year_0
        pop.set_index(self.idx_cols)
        PY = pop['population'] * (1 - 0.5 * prob_death)
        PY_adj = PY * (1 - yld_rate_now)
        self.create_table(pop['age'].min(), self.year_0)
        pop['PYadj'] = PY_adj
        df = self.data.merge(pop, on=self.idx_cols, how='left', suffixes=('', '_new'))
        new_vals = df['PYadj_new'].notna()
        new_popn = df['population_new'].notna()
        if not new_popn.equals(new_vals):
            raise ValueError('New population and PYadj values differ')
        self.data.loc[new_vals, 'PYadj'] = df.loc[new_vals, 'PYadj_new']
        self.data.loc[new_vals, 'population'] = df.loc[new_vals, 'population_new']

    def on_collect_metrics(self, event):
        """
        Calculate adjusted person-years for the current year.
        """
        idx = event.index
        year = self.clock().year
        acm_rate_now = self.acm_rate(idx)
        yld_rate_now = self.yld_rate(idx)
        prob_death = 1 - np.exp(- acm_rate_now)
        pop = self.population_view.get(idx)
        if len(pop.index) == 0:
            # No tracked population remains.
            return
        pop['year'] = year
        # NOTE: PY(a,t_0) equation actually applies to ANY t.
        PY = pop['population'] * (1 - 0.5 * prob_death)
        PY_adj = PY * (1 - yld_rate_now)
        pop['PYadj'] = PY_adj
        df = self.data.merge(pop, on=self.idx_cols, how='left', suffixes=('', '_new'))
        if len(df.index) > 0:
            new_vals = df['PYadj_new'].notna()
            new_popn = df['population_new'].notna()
            if not new_popn.equals(new_vals):
                raise ValueError('New population and PYadj values differ')
            self.data.loc[new_vals, 'PYadj'] = df.loc[new_vals, 'PYadj_new']
            self.data.loc[new_vals, 'population'] = df.loc[new_vals, 'population_new']

    def get_table(self):
        """
        Return a Pandas data frame that contains the adjusted person-years and
        adjusted life-expectancy for each cohort at each year of the
        simulation.
        """
        mask = self.data['PYadj'].notna()
        return self.data.loc[mask].sort_index()

    def to_csv(self, filename):
        """
        Save the adjusted person-years and the adjusted life-expectancy data
        to a CSV file.
        """
        self.get_table().to_csv(filename, index=False)

    def finalise_output(self, event):
        """
        Calculate the adjusted life-expectancy for each cohort at each year of
        the simulation, now that the simulation has finished and the adjusted
        person-years for each cohort at each year have been calculated.

        If an output file name was provided to the constructor, this method
        will also save these data to a CSV file.
        """
        # Identify each generation by their year of birth.
        self.data['year_of_birth'] = self.data['year'] - self.data['age']
        # Sort the table by cohort (i.e., generation and sex), and then by
        # calendar year, so that results are output in the same order as in
        # the spreadsheet models.
        self.data = self.data.sort_values(by=['year_of_birth', 'sex', 'age'],
                                          axis=0)
        self.data = self.data.reset_index(drop=True)
        # Group the adjusted person-years by cohort
        group_cols = ['year_of_birth', 'sex']
        subset_cols = group_cols + ['PYadj']
        grouped = self.data.loc[:, subset_cols].groupby(by=group_cols)['PYadj']
        # Calculate the reverse-cumulative sums of the adjusted person-years
        # (i.e., the present and future person-years) by:
        #   (a) reversing the adjusted person-years values in each cohort;
        #   (b) calculating the cumulative sums in each cohort; and
        #   (c) restoring the original order.
        cumsum = grouped.apply(lambda x: pd.Series(x[::-1].cumsum()).iloc[::-1])
        # Calculate the adjusted life expectancy for each cohort at each year
        # by dividing the remaining person-years by the population size.
        self.data['LEadj'] = cumsum / self.data['population']
        # Re-order the columns to better reflect how the spreadsheet model
        # tables are arranged.
        self.data = self.data[['year_of_birth', 'sex', 'age', 'year',
                               'population', 'PYadj', 'LEadj']]
        if self.output_file is not None:
            self.to_csv(self.output_file)
