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

    def __init__(self, output_file=None, unadjusted=False):
        """
        :param output_file: The name of the CSV file in which to record the
            adjusted person-years and adjusted life-expectancy data.
        :param unadjusted: Whether to also record unadjusted person-years and
            life-expectancy data.
        """
        self.output_file = output_file
        self.unadjusted = unadjusted

    def setup(self, builder):
        """
        """
        self.yld_rate = builder.value.get_value('yld_rate')
        self.acm_rate = builder.value.get_value('mortality_rate')
        view_cols = ['age', 'sex', 'population', 'bau_population']
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
        self.data['PY'] = np.nan
        self.data['population'] = np.nan
        self.data['bau_PYadj'] = np.nan
        self.data['bau_PY'] = np.nan
        self.data['bau_population'] = np.nan

    def record_person_years(self, idx):
        """
        Record the un-adjusted and adjusted person-years for the BAU and the
        intervention.
        """
        pop = self.population_view.get(idx)
        if len(pop.index) == 0:
            # No tracked population remains.
            return

        pop['year'] = self.clock().year

        # Calculate (adjusted) person-years for the intervention.
        int_pr_death = 1 - np.exp(- self.acm_rate(idx))
        pop['PY'] = pop['population'] * (1 - 0.5 * int_pr_death)
        pop['PYadj'] = pop['PY'] * (1 - self.yld_rate(idx))

        # Calculate (adjusted) person-years for the BAU.
        bau_pr_death = 1 - np.exp(- self.acm_rate.source(idx))
        pop['bau_PY'] = pop['bau_population'] * (1 - 0.5 * bau_pr_death)
        pop['bau_PYadj'] = pop['bau_PY'] * (1 - self.yld_rate.source(idx))

        # Determine if any strata statistics need to be updated.
        df = self.data.merge(pop, on=self.idx_cols, how='left',
                             suffixes=('', '_new'))
        if len(df.index) == 0:
            # No strata to update.
            return

        # Determine which strata to update.
        new_vals = df['PYadj_new'].notna()

        # Update (adjusted) person-years and population denominators.
        int_cols = ['PY', 'PYadj', 'population']
        bau_cols = ['bau_{}'.format(c) for c in int_cols]
        for col in int_cols + bau_cols:
            new_col = '{}_new'.format(col)
            self.data.loc[new_vals, col] = df.loc[new_vals, new_col]

    def on_initialize(self, pop_data):
        """
        Calculate adjusted person-years and adjusted life-expectancy before
        the first time-step.
        """
        idx = pop_data.index
        pop = self.population_view.get(idx)
        self.create_table(pop['age'].min(), self.clock().year)
        self.record_person_years(idx)

    def on_collect_metrics(self, event):
        """
        Calculate adjusted person-years for the current year.
        """
        self.record_person_years(event.index)

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

    def calculate_life_expectancy(self, py_col, le_col, denom_col):
        # Group the person-years by cohort.
        group_cols = ['year_of_birth', 'sex']
        subset_cols = group_cols + [py_col]
        grouped = self.data.loc[:, subset_cols].groupby(by=group_cols)[py_col]
        # Calculate the reverse-cumulative sums of the adjusted person-years
        # (i.e., the present and future person-years) by:
        #   (a) reversing the adjusted person-years values in each cohort;
        #   (b) calculating the cumulative sums in each cohort; and
        #   (c) restoring the original order.
        cumsum = grouped.apply(lambda x: pd.Series(x[::-1].cumsum()).iloc[::-1])
        self.data[le_col] = cumsum / self.data[denom_col]

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
        # Calculate (adjusted) life expectancy for BAU and the intervention.
        self.calculate_life_expectancy('PY', 'LE', 'population')
        self.calculate_life_expectancy('PYadj', 'LEadj', 'population')
        self.calculate_life_expectancy('bau_PY', 'bau_LE', 'bau_population')
        self.calculate_life_expectancy('bau_PYadj', 'bau_LEadj', 'bau_population')
        # Calculate differences between the BAU and the intervention.
        self.data['diff_LE'] = self.data['LE'] - self.data['bau_LE']
        self.data['diff_PY'] = self.data['PY'] - self.data['bau_PY']
        self.data['diff_LEadj'] = self.data['LEadj'] - self.data['bau_LEadj']
        self.data['diff_PYadj'] = self.data['PYadj'] - self.data['bau_PYadj']
        # Re-order the columns to better reflect how the spreadsheet model
        # tables are arranged.
        cols = ['year_of_birth', 'sex', 'age', 'year', 'population',
                'PYadj', 'LEadj', 'bau_PYadj', 'bau_LEadj',
                'diff_LEadj', 'diff_PYadj']
        if self.unadjusted:
            cols.extend(['PY', 'LE', 'bau_PY', 'bau_LE', 'diff_PY', 'diff_LE'])
        self.data = self.data[cols]
        if self.output_file is not None:
            self.to_csv(self.output_file)
