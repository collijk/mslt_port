"""Build population data tables."""

import pandas as pd
import numpy as np
import pathlib


class Population:

    def __init__(self, data_dir, year_start):
        data_file = '{}/base_population.csv'.format(data_dir)
        data_path = str(pathlib.Path(data_file).resolve())

        df = pd.read_csv(data_path)
        df = df.rename(columns={'mortality per 1 rate': 'mortality_rate',
                                'pYLD rate': 'disability_rate',
                                'APC in all-cause mortality': 'mortality_apc',
                                '5-year': 'population'})

        # Use identical populations in the BAU and intervention scenarios.
        df['bau_population'] = df['population'].values

        # Retain only the necessary columns.
        df['year'] = year_start
        df = df[['year', 'age', 'sex', 'population', 'bau_population',
                 'disability_rate', 'mortality_rate', 'mortality_apc']]

        # Remove strata that have already reached the terminal age.
        df = df[~ (df.age == df['age'].max())]

        # Sort the rows.
        df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

        self.year_start = year_start
        self.year_end = year_start + df['age'].max() - df['age'].min()
        self._num_apc_years = 15

        self._data = df

    def years(self):
        """Return an iterator over the simulation period."""
        return range(self.year_start, self.year_end + 1)

    def get_population(self):
        """Return the initial population size for each stratum."""
        cols = ['year', 'age', 'sex', 'population', 'bau_population']
        # Retain only those strata for whom the population size is defined.
        return self._data.loc[self._data['population'].notna(), cols]

    def get_disability_rate(self):
        """Return the disability rate for each stratum."""
        df = self._data[['age', 'sex', 'disability_rate']]
        df = df.rename(columns={'disability_rate': 'rate'})

        tables = []
        for year in self.years():
            df['year'] = year
            tables.append(df.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.reset_index(drop=True)

        return df

    def get_acmr_apc(self):
        """Return the annual percent change (APC) in mortality rate."""
        df = self._data[['year', 'age', 'sex', 'mortality_apc']]
        df = df.rename(columns={'mortality_apc': 'value'})

        tables = []
        for year in self.years():
            df['year'] = year
            tables.append(df.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.reset_index(drop=True)

        return df

    def get_mortality_rate(self):
        """
        Return the mortality rate for each strata.

        :param df_base: The base population data.
        """
        # NOTE: see column IG in ErsatzInput.
        # - Each cohort has a separate APC (column FE)
        # - ACMR = BASE_ACMR * e^(APC * (year - 2011))
        df_apc = self.get_acmr_apc()
        df_acmr = self._data[['year', 'age', 'sex', 'mortality_rate']]
        df_acmr = df_acmr.rename(columns={'mortality_rate': 'rate'})
        base_acmr = df_acmr['rate'].copy()

        tables = []
        df_acmr['year'] = self.year_start - 1
        tables.append(df_acmr.copy())
        for counter, year in enumerate(self.years()):
            if counter <= self._num_apc_years:
                year_apc = df_apc[df_apc.year == year]
                apc = year_apc['value'].values
                scale = np.exp(apc * (year - year_start))
                df_acmr.loc[:, 'rate'] = base_acmr * scale
            else:
                # NOTE: use the same scale for this cohort as per the previous
                # year; shift by 2 because there are male and female cohorts.
                scale[2:] = scale[:-2]
                df_acmr.loc[:, 'rate'] = base_acmr * scale
            df_acmr['year'] = year
            tables.append(df_acmr.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.reset_index(drop=True)

        return df
