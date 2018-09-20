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


class AdjustedPersonYears:
    """
    This class calculates the adjusted person-years for each cohort at each
    year of the simulation. The results can be obtained by calling the
    get_table() method, or saved to a CSV file by calling the to_csv() method.
    """

    def __init__(self, output_py_file=None, output_le_file=None):
        """
        :param output_py_file: The name of the CSV file in which to record the
            adjusted person-years data.
        :param output_le_file: The name of the CSV file in which to record the
            adjusted life expectancy data.
        """
        self.output_py_file = output_py_file
        self.output_le_file = output_le_file

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
        builder.event.register_listener('simulation_end', self.write_output)
        self.idx_cols = ['age', 'sex', 'year']

    def create_PYadj_table(self, min_age, min_year):
        """
        Create an empty data frame to hold all of the adjusted person-year
        values that will be produced during a simulation.
        """
        ages = range(min_age, AGE_GROUP_END + 1)
        sexes = ['male', 'female']
        years = range(min_year, YEAR_END + 1)
        # NOTE: columns must be in the same order as self.idx_cols.
        rows = list(itertools.product(ages, sexes, years))
        self.data_py = pd.DataFrame(rows, columns=self.idx_cols)
        self.data_py.set_index(self.idx_cols)
        self.data_py['PYadj'] = np.nan
        # Also create a table to record the adjusted life expectancy for each
        # cohort.
        le_rows = list(itertools.product(ages, sexes, [min_year]))
        self.data_le = pd.DataFrame(le_rows, columns=self.idx_cols)
        self.data_le.set_index(self.idx_cols)
        self.data_le['LEadj'] = 0.0

    def on_initialize(self, pop_data):
        """
        Calculate adjusted person-years before the first time-step.
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
        self.create_PYadj_table(pop['age'].min(), self.year_0)
        pop['PYadj'] = PY_adj
        df = self.data_py.merge(pop, on=self.idx_cols, how='left', suffixes=('', '_new'))
        new_vals = df['PYadj_new'].notna()
        self.data_py.loc[new_vals, 'PYadj'] = df.loc[new_vals, 'PYadj_new']
        # Calculate the adjusted life expectancy for each cohort.
        pop['LEadj'] = pop['PYadj'] / pop['population']
        df = self.data_le.merge(pop.loc[:, self.idx_cols + ['LEadj']],
                               on=self.idx_cols, how='left',
                               suffixes=('', '_new'))
        new_vals = df['LEadj_new'].notna()
        self.data_le.loc[new_vals, 'LEadj'] = (self.data_le.loc[new_vals, 'LEadj']
                                               + df.loc[new_vals, 'LEadj_new'])

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
        df = self.data_py.merge(pop, on=self.idx_cols, how='left', suffixes=('', '_new'))
        if len(df.index) > 0:
            new_vals = df['PYadj_new'].notna()
            self.data_py.loc[new_vals, 'PYadj'] = df.loc[new_vals, 'PYadj_new']
            # Calculate the adjusted life expectancy for each cohort.
            pop['LEadj'] = pop['PYadj'] / pop['population']
            pop['age'] = pop['age'] - pop['year'] + self.year_0
            pop['year'] = self.year_0
            df = self.data_le.merge(pop.loc[:, self.idx_cols + ['LEadj']],
                                    on=self.idx_cols, how='left',
                                    suffixes=('', '_new'))
            new_vals = df['LEadj_new'].notna()
            self.data_le.loc[new_vals, 'LEadj'] = (self.data_le.loc[new_vals, 'LEadj']
                                                   + df.loc[new_vals, 'LEadj_new'])

    def get_adj_py(self):
        """
        Return a Pandas data frame that contains the adjusted person-years for
        each cohort at each year of the simulation.
        """
        mask = self.data_py['PYadj'].notna()
        return self.data_py.loc[mask].sort_index()

    def get_adj_le(self):
        """
        Return a Pandas data frame that contains the adjusted life expectancy
        for each cohort at the starting year of the simulation.
        """
        return self.data_le

    def adj_py_to_csv(self, filename):
        """Save the adjusted person-years table to a CSV file."""
        self.get_adj_py().to_csv(filename, index=False)

    def adj_le_to_csv(self, filename):
        """Save the adjusted life expectancy table to a CSV file."""
        self.get_adj_le().to_csv(filename, index=False)

    def write_output(self, event):
        """
        Save the adjusted person-years table at the end of the simulation.
        """
        if self.output_py_file is not None:
            self.adj_py_to_csv(self.output_py_file)
        if self.output_le_file is not None:
            self.adj_le_to_csv(self.output_le_file)
