"""Build risk-factor data tables."""

import logging
import pandas as pd
import numpy as np
import pathlib

from .uncertainty import (LogNormalRawSD, sample_column, sample_column_long)


def post_cessation_rr(disease, rr_data, rr_col, num_states, gamma):
    # NOTE: this will preserve the condition that draw 0 is the mean.
    for n in range(1, num_states):
        rr_base = rr_data.loc[:, rr_col].values
        col_name = '{}_{}'.format(disease, n)
        rr_data[col_name] = 1 + (rr_base - 1) * np.exp(- gamma.values * n)

    # Add the RR = 1.0 for the final (absorbing) state.
    col_name = '{}_{}'.format(disease, num_states)
    rr_data[col_name] = 1.0

    if np.any(rr_data.isna()):
        raise ValueError('NA values in post-cessation RRs')

    return rr_data


def sample_tobacco_rate(year_start, year_end, data, rate_name, prev_data,
                        apc_data, num_apc_years, prng, rate_dist, n):
    """
    Draw correlated samples for a tobacco rate at each year.

    :param year_start: The year at which the simulation starts.
    :param year_end: The year at which the simulation ends.
    :param data: The data table that contains the rate values.
    :param rate_name: The column name that defines the mean values.
    :param prev_data: The data table that contains the initial prevalence.
    :param apc_data: The data table that contains the annual percent changes
        (set to ``None`` if there are no changes).
    :param num_apc_years: The number of years over which the annual percent
        changes apply (measured from the start of the simulation).
    :param prng: The random number generator (``numpy.random.RandomState``).
    :param rate_dist: The uncertainty distribution for the rate values.
    :param n: The number of samples to draw per row.
    """
    value_col = rate_name

    if value_col == 'incidence':
        # The uptake rate is defined by the initial prevalence.
        data.loc[:, 'incidence'] = 0.0
        initial_rate = prev_data.loc[prev_data['age'] == 20, 'tobacco.yes']
        data.loc[data['age'] == 20, 'incidence'] = initial_rate

    # Sample the initial rate for each cohort.
    df = sample_column_long(data, value_col, prng, rate_dist, n)

    df.insert(0, 'year_start', 0)
    df.insert(1, 'year_end', 0)

    df_index_cols = ['year_start', 'year_end', 'age', 'sex', 'draw']
    apc_index_cols = ['age', 'sex']

    tables = []
    years = range(year_start, year_end + 1)

    if apc_data is not None and value_col in apc_data.columns:
        data_columns = [c for c in df.columns if c not in df_index_cols]
        apc_values = apc_data.loc[:, value_col].values
        base_values = df.loc[:, data_columns].copy().values

        initial_rate = df.loc[:, data_columns].copy().values
        frac = (1 - apc_data.loc[:, value_col].values)

        # Calculate the correlated samples for each cohort at each year.
        for counter, year in enumerate(years):
            df['year_start'] = year
            if counter < num_apc_years:
                df['year_end'] = year + 1
                timespan = year - year_start
                result = initial_rate * (frac ** timespan)[..., np.newaxis]
                df.loc[:, data_columns] = result
                tables.append(df.copy())
            else:
                df['year_end'] = year_end + 1
                tables.append(df.copy())
                break

        df = pd.concat(tables)

    else:
        df['year_start'] = year_start
        df['year_end'] = year_end + 1

    # Replace 'age' with age groups.
    df = df.rename(columns={'age': 'age_group_start'})
    df.insert(df.columns.get_loc('age_group_start') + 1,
              'age_group_end',
              df['age_group_start'] + 1)

    df = df.sort_values(['year_start', 'age_group_start', 'sex'])
    df = df.reset_index(drop=True)

    return df


def sample_tobacco_rate_from(year_start, year_end, data, rate_name, prev_data,
                             apc_data, num_apc_years, rate_dist, samples):
    """
    Draw correlated samples for a tobacco rate at each year.

    :param year_start: The year at which the simulation starts.
    :param year_end: The year at which the simulation ends.
    :param data: The data table that contains the rate values.
    :param rate_name: The column name that defines the mean values.
    :param prev_data: The data table that contains the initial prevalence.
    :param apc_data: The data table that contains the annual percent changes
        (set to ``None`` if there are no changes).
    :param num_apc_years: The number of years over which the annual percent
        changes apply (measured from the start of the simulation).
    :param rate_dist: The uncertainty distribution for the rate values.
    :param samples: Random samples drawn from the half-open interval [0, 1).
    """
    value_col = rate_name

    if value_col == 'incidence':
        # The uptake rate is defined by the initial prevalence.
        data.loc[:, 'incidence'] = 0.0
        initial_rate = prev_data.loc[prev_data['age'] == 20, 'tobacco.yes']
        data.loc[data['age'] == 20, 'incidence'] = initial_rate

    # Sample the initial rate for each cohort.
    df = sample_column_long(data, value_col, rate_dist, samples)

    df.insert(0, 'year_start', 0)
    df.insert(1, 'year_end', 0)

    df_index_cols = ['year_start', 'year_end', 'age', 'sex', 'draw']
    apc_index_cols = ['age', 'sex']

    tables = []
    years = range(year_start, year_end + 1)

    if apc_data is not None and value_col in apc_data.columns:
        data_columns = [c for c in df.columns if c not in df_index_cols]
        apc_df = apc_data.merge(df.loc[:, df_index_cols])
        apc_values = apc_df.loc[:, value_col].values
        base_values = df.loc[:, data_columns].copy().values

        initial_rate = df.loc[:, data_columns].copy().values
        frac = (1 - apc_df.loc[:, value_col].values)

        # Calculate the correlated samples for each cohort at each year.
        for counter, year in enumerate(years):
            df['year_start'] = year
            if counter < num_apc_years:
                df['year_end'] = year + 1
                timespan = year - year_start
                result = initial_rate * (frac ** timespan)[..., np.newaxis]
                df.loc[:, data_columns] = result
                tables.append(df.copy())
            else:
                df['year_end'] = year_end + 1
                tables.append(df.copy())
                break

        df = pd.concat(tables)

    else:
        df['year_start'] = year_start
        df['year_end'] = year_end + 1

    # Replace 'age' with age groups.
    df = df.rename(columns={'age': 'age_group_start'})
    df.insert(df.columns.get_loc('age_group_start') + 1,
              'age_group_end',
              df['age_group_start'] + 1)

    df = df.sort_values(['year_start', 'age_group_start', 'sex', 'draw'])
    df = df.reset_index(drop=True)

    return df


class Tobacco:

    def __init__(self, data_dir, year_start, year_end):
        self._year_start = year_start
        self._year_end = year_end
        self.data_dir = data_dir
        self._initial_rates = self.load_initial_tobacco_rates()
        self._apc = self.load_tobacco_rates_apc()
        self._prev = self.load_initial_tobacco_prevalence()
        self._tax = self.load_tobacco_tax_effects()
        self._mort_rr = self.load_tobacco_mortality_rr()
        self._dis_rr_dict, self._dis_rr_df = self.load_tobacco_diseases_rr()
        self._df_gamma = self.load_tobacco_disease_rr_gamma()
        self._df_dis_rr_sd = self.load_tobacco_disease_rr_sd()

    def sample_price_elasticity_from(self, elast_dist, samples):
        """
        Sample the price elasticity and return the individual draws.

        :param elast_dist: The sampling distribution.
        :param samples: Samples drawn from the half-open interval [0, 1).
        """
        df_elast = sample_column_long(self._df_elast, 'Elasticity',
                                      elast_dist, samples)
        return df_elast

    def scale_price_elasticity_from(self, df_elast, mean_scale, scale_dist,
                                    samples):
        """
        Scale the price elasticity by variable amounts.

        :param df_elast: Baseline samples of the price elasticity.
        :param mean_scale: The average scale to apply (e.g., 1.2).
        :param scale_dist: The sampling distribution for ``mean_scale``.
        :param samples: Samples drawn from the half-open interval [0, 1).
        """
        # Define draw 0 as the expected value.
        samples = np.append([0.5], samples)

        scales = scale_dist.correlated_samples(pd.Series(mean_scale), samples)
        df_scales = pd.DataFrame({'scale': scales.T[0]})
        df_scales['draw'] = list(df_scales.index)

        df_elast = df_elast.copy()
        df_elast = df_elast.merge(df_scales)
        df_elast.loc[:, 'Elasticity'] = (df_elast.loc[:, 'Elasticity'].values
                                         * df_elast.loc[:, 'scale'].values)
        df_elast = df_elast.drop(columns='scale')
        df_elast = df_elast.sort_values(['age', 'sex', 'draw'])
        df_elast = df_elast.reset_index(drop=True)

        return df_elast

    def sample_tax_effects_from(self, elast_dist, samples):
        """
        Sample the price elasticity and return the effects of a tobacco tax on
        uptake and remission for each elasticity draw.

        :param elast_dist: The sampling distribution.
        :param samples: Samples drawn from the half-open interval [0, 1).
        """
        df_elast = self.sample_price_elasticity_from(elast_dist, samples)
        return self.sample_tax_effects_from_elasticity_wide(df_elast)

    def sample_tax_effects_from_elasticity_wide(self, df_elast):
        """
        Calculate the effects of a tobacco tax on uptake and remission, given
        a number of samples for the price elasticity.

        :param df_elast: Samples of the price elasticity.
        """
        df_price = self._df_price
        df_elast = df_elast.copy()

        df_elast.insert(0, 'year_start', self._year_start)
        df_elast.insert(1, 'year_end', self._year_end + 1)

        if np.any(df_elast.isna()):
            raise ValueError('NA values in elast_cols')

        df_tmp = pd.DataFrame(columns=['incidence_effect', 'remission_effect'])
        df_elast = df_elast.join(df_tmp)

        # Loop over the price at each year, and determine the effects of price
        # elasticity on uptake and cessation rates.
        start_price = df_price.loc[0, 'price']
        tables = []
        for i, row in enumerate(df_price.itertuples()):
            df_elast['year_start'] = row.year
            df_elast['year_end'] = row.year + 1
            df_elast['price'] = row.price

            df_elast['incidence_effect'] = np.exp(
                - df_elast['Elasticity'].values
                * np.log(row.price / start_price))

            prev_price = row.price if i == 0 else df_price.loc[i - 1, 'price']
            if row.price > prev_price:
                df_elast['remission_effect'] = np.exp(
                    - df_elast['Elasticity'].values
                    * np.log(row.price / prev_price))
            else:
                df_elast['remission_effect'] = 1.0

            if np.any(df_elast.isna()):
                raise ValueError('NA values found in tobacco tax effects data')

            tables.append(df_elast.copy())

        df = pd.concat(tables)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)
        df = df.loc[:, ['year_start', 'year_end',
                        'age_group_start', 'age_group_end',
                        'sex', 'draw',
                        'incidence_effect', 'remission_effect']]

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco tax effects data')

        return df

    def sample_price_elasticity(self, prng, elast_dist, n):
        """
        Sample the price elasticity and return the individual draws.

        :param prng: The random number generator.
        :param elast_dist: The sampling distribution.
        :param n: The number of samples to draw.
        """
        df_elast = sample_column(self._df_elast, 'Elasticity',
                                 prng, elast_dist, n)
        return df_elast

    def scale_price_elasticity(self, df_elast, mean_scale, prng, scale_dist):
        """
        Scale the price elasticity by variable amounts.

        :param df_elast: Baseline samples of the price elasticity.
        :param mean_scale: The average scale to apply (e.g., 1.2).
        :param prng: The random number generator.
        :param scale_dist: The sampling distribution for ``mean_scale``.
        """
        elast_cols = [c for c in df_elast.columns
                      if c.startswith('Elasticity_draw_')]
        num_samples = len(elast_cols)

        scale_samples = prng.random_sample(size=num_samples)
        scales = scale_dist.correlated_samples(pd.Series(mean_scale),
                                               scale_samples)
        # Remember to define draw 0 as the expected value.
        scales[0] = mean_scale

        df_elast = df_elast.copy()
        df_elast.loc[:, elast_cols] = (df_elast.loc[:, elast_cols].values
                                       * scales.T)

        return df_elast

    def sample_tax_effects_from_elasticity(self, df_elast):
        """
        Calculate the effects of a tobacco tax on uptake and remission, given
        a number of samples for the price elasticity.

        :param df_elast: Samples of the price elasticity.
        """
        df_price = self._df_price
        df_elast = df_elast.copy()
        df_elast.insert(0, 'year', 0)

        elast_cols = [c for c in df_elast.columns
                      if c.startswith('Elasticity_draw_')]
        inc_cols = ['incidence_effect_draw_{}'.format(i)
                    for i in range(len(elast_cols))]
        rem_cols = ['remission_effect_draw_{}'.format(i)
                    for i in range(len(elast_cols))]

        if np.any(df_elast.loc[:, elast_cols].isna()):
            raise ValueError('NA values in elast_cols')

        df_tmp = pd.DataFrame(columns=inc_cols + rem_cols)
        df_elast = df_elast.join(df_tmp)

        # Loop over the price at each year, and determine the effects of price
        # elasticity on uptake and cessation rates.
        start_price = df_price.loc[0, 'price']
        tables = []
        for i, row in enumerate(df_price.itertuples()):
            df_elast['year'] = row.year
            df_elast['price'] = row.price

            df_elast.loc[:, inc_cols] = np.exp(
                - df_elast.loc[:, elast_cols].values
                * np.log(row.price / start_price))

            prev_price = row.price if i == 0 else df_price.loc[i - 1, 'price']
            if row.price > prev_price:
                df_elast.loc[:, rem_cols] = np.exp(
                    - df_elast.loc[:, elast_cols].values
                    * np.log(row.price / prev_price))
            else:
                df_elast.loc[:, rem_cols] = 1.0

            tables.append(df_elast.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.reset_index(drop=True)
        df = df.loc[:, ['year', 'age', 'sex'] + inc_cols + rem_cols]

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco tax effects data')

        return df

    def sample_tax_effects(self, prng, elast_dist, n):
        """
        Sample the price elasticity and return the effects of a tobacco tax on
        uptake and remission for each elasticity draw.

        :param prng: The random number generator.
        :param elast_dist: The sampling distribution.
        :param n: The number of samples to draw.
        """
        df_elast = self.sample_price_elasticity(prng, elast_dist, n)
        return self.sample_tax_effects_from_elasticity(df_elast)

    def sample_disease_rr_from(self, samples_tbl):
        """
        Sample the relative risk of chronic disease incidence for each
        exposure category, in all years of the simulation.

        Note that the sampling distribution for each disease is predefined.

        :param samples_tbl: A dictionary that maps disease names to samples
            from the half-open interval [0, 1).
        """
        logger = logging.getLogger(__name__)
        tables = []

        for key, table in self._dis_rr_dict.items():
            if key not in self._df_gamma.columns:
                logger.info('Ignoring {} (no gamma)'.format(key))
                continue

            if key not in samples_tbl:
                logger.info('Ignoring {} (no samples)'.format(key))
                continue

            # Extract the relative risk for current smokers.
            rr_col = '{}_0'.format(key)
            rr_yes = table[rr_col]

            # Add the standard deviation as a new column.
            df_rr_sd = self._df_dis_rr_sd.loc[
                self._df_dis_rr_sd['disease'] == key]
            if df_rr_sd.empty:
                raise ValueError('No RR distribution for {}'.format(key))

            df_rr_sd = df_rr_sd.rename(columns={'RR': rr_col})
            table = table.merge(df_rr_sd.loc[:, ['sex', rr_col, 'sd']],
                                how='left')
            table['sd'].fillna(0.0, inplace=True)

            rr_dist = LogNormalRawSD(table['sd'])

            # Determine how many post-cessation states there are (not
            # including the 0 years post-cessation state).
            cols = ['age', 'sex', rr_col, 'sd']
            num_states = len([c for c in table.columns if c not in cols])

            # Sample the relative risk for current smokers.
            df_rr = sample_column_long(table.loc[:, cols], rr_col, rr_dist,
                                       samples_tbl[key])
            rr_cols = [col for col in df_rr.columns
                       if col not in ['age', 'sex', 'draw']]

            if len(rr_cols) != 1:
                raise ValueError('Expected one column: {}'.format(rr_cols))

            # Determine the value of gamma for each cohort.
            # It is only indexed by age, so we need to merge it with the index
            # columns of df_rr ('age' and 'sex').
            gamma = self._df_gamma.loc[:, ['age', key]]
            gamma.columns = ['age', 'gamma']
            df_gamma = gamma.merge(df_rr.loc[:, ['age', 'sex']])
            gamma = df_gamma['gamma']

            # Calculate how this RR decays to 1.0, post-cessation.
            df_rr = post_cessation_rr(key, df_rr, rr_col, num_states, gamma)

            # NOTE: add the 'Disease_no' and 'Disease_yes' columns.
            df_rr.insert(df_rr.columns.get_loc(rr_col), '{}_no'.format(key),
                         1.0)
            df_rr.insert(df_rr.columns.get_loc(rr_col), '{}_yes'.format(key),
                         df_rr[rr_col])

            if len(tables) > 0:
                # Remove age and sex columns, so that the RR tables for each
                # disease can be joined (otherwise these columns will be
                # duplicated in each table, and the join will fail).
                df_rr = df_rr.drop(columns=['age', 'sex', 'draw'])

            tables.append(df_rr)

        if len(tables) > 1:
            df = tables[0].join(tables[1:])
        elif len(tables) == 1:
            df = tables[0]
        else:
            raise ValueError('No diseases with RRs')

        df.insert(0, 'year_start', self._year_start)
        df.insert(1, 'year_end', self._year_end + 1)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        return df

    def sample_disease_rr(self, prng, n):
        """
        Sample the relative risk of chronic disease incidence for each
        exposure category, in all years of the simulation.

        Note that the sampling distribution for each disease is predefined.

        :param prng: The random number generator.
        :param n: The number of samples to draw.
        """
        logger = logging.getLogger(__name__)
        tables = []

        for key, table in self._dis_rr_dict.items():
            if key not in self._df_gamma.columns:
                logger.info('Ignoring {}'.format(key))
                continue

            # Extract the relative risk for current smokers.
            rr_col = '{}_0'.format(key)
            rr_yes = table[rr_col]

            # Add the standard deviation as a new column.
            df_rr_sd = self._df_dis_rr_sd.loc[
                self._df_dis_rr_sd['disease'] == key]
            if df_rr_sd.empty:
                raise ValueError('No RR distribution for {}'.format(key))

            df_rr_sd = df_rr_sd.rename(columns={'RR': rr_col})
            table = table.merge(df_rr_sd.loc[:, ['sex', rr_col, 'sd']],
                                how='left')
            table['sd'].fillna(0.0, inplace=True)

            rr_dist = LogNormalRawSD(table['sd'])

            # Determine how many post-cessation states there are (not
            # including the 0 years post-cessation state).
            cols = ['age', 'sex', rr_col, 'sd']
            num_states = len([c for c in table.columns if c not in cols])

            # Sample the relative risk for current smokers.
            df_rr = sample_column(table.loc[:, cols], rr_col,
                                  prng, rr_dist, n)
            rr_cols = [col for col in df_rr.columns
                       if col not in ['age', 'sex']]

            # Determine the value of gamma for each cohort.
            # It is only indexed by age, so we need to merge it with the index
            # columns of df_rr ('age' and 'sex').
            gamma = self._df_gamma.loc[:, ['age', key]]
            gamma.columns = ['age', 'gamma']
            df_gamma = gamma.merge(df_rr.loc[:, ['age', 'sex']])
            gamma = df_gamma['gamma']

            # Calculate how this RR decays to 1.0, post-cessation.
            df_rr = post_cessation_rr(key, df_rr, rr_cols, num_states, gamma)

            if len(tables) > 0:
                # Remove age and sex columns, so that the RR tables for each
                # disease can be joined (otherwise these columns will be
                # duplicated in each table, and the join will fail).
                df_rr = df_rr.drop(columns=['age', 'sex'])

            tables.append(df_rr)

        df = tables[0].join(tables[1:])
        return df

    def get_expected_tax_effects(self):
        """
        Return the effects of a tobacco tax on incidence and remission rates.
        """
        df = self._tax.copy()

        # Replace 'year' with year bins.
        df = df.rename(columns={'year': 'year_start'})
        df.insert(df.columns.get_loc('year_start') + 1,
                  'year_end',
                  df['year_start'] + 1)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco tax effects data')

        return df

    def get_expected_disease_rr(self, disease):
        """
        Return the relative risk of chronic disease incidence for each
        exposure category, in all years of the simulation.

        :param disease: The name of the disease; set to ``None`` to return a
            single table for all diseases.
        """
        if disease is None:
            df = self._dis_rr_df.copy()
        else:
            if disease not in self._dis_rr_dict:
                msg = 'No relative risks for disease {}'.format(disease)
                raise ValueError(msg)

            df = self._dis_rr_dict[disease].copy()
            df.insert(0, 'year', 0)

        # Replace 'year' with year bins.
        df = df.rename(columns={'year': 'year_start'})
        df.insert(df.columns.get_loc('year_start') + 1,
                  'year_end',
                  self._year_end + 1)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco disease RR data')

        return df

    def get_expected_mortality_rr(self):
        """
        Return the relative risk of mortality for each exposure category, in
        all years of the simulation.
        """
        df = self._mort_rr.copy()

        # Copy the relative-risk columns so they apply to the intervention.
        bau_prefix = 'tobacco.'
        int_prefix = 'tobacco_intervention.'
        for col in df.columns:
            if col.startswith(bau_prefix):
                int_col = col.replace(bau_prefix, int_prefix)
                df[int_col] = df[col]

        df = df.drop(columns='year')
        df.insert(0, 'year_start', self._year_start)
        df.insert(1, 'year_end', self._year_end + 1)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco mortality RR data')

        return df

    def get_expected_prevalence(self):
        """Return the initial exposure prevalence."""
        # NOTE: set prevalence to zero at age 20, it will be taken care of by
        # the incidence rate in the first time-step.
        df = self._prev.copy()
        mask = df['age'] == 20
        df.loc[mask, 'tobacco.no'] += df.loc[mask, 'tobacco.yes']
        df.loc[mask, 'tobacco.yes'] = 0.0

        df = df.drop(columns='year')
        df.insert(0, 'year_start', self._year_start)
        df.insert(1, 'year_end', self._year_end + 1)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco mortality RR data')

        return df

    def get_expected_rates(self):
        """
        Return the incidence and remission rates for tobacco use in all years
        of the simulation.
        """
        # The incidence rate is calculated with respect to the initial
        # prevalence of new smokers (i.e., those aged 20).
        initial_prev = self._prev.loc[self._prev['age'] == 20]
        initial_prev = initial_prev.rename(columns={
            'tobacco.yes': 'prevalence'})

        # Set the initial prevalence in 20-year-old cohorts to zero, so that
        # tobacco interventions can have an immediate effect in 2011.
        # Note that this will not affect the 'prevalence' column of df_apc.
        df_prev = self._prev.copy()
        mask = df_prev['age'] == 20
        df_prev.loc[mask, 'tobacco.no'] += df_prev.loc[mask, 'tobacco.yes']
        df_prev.loc[mask, 'tobacco.yes'] = 0.0

        df_apc = self._apc.copy()
        df_apc['prevalence'] = 0.0
        df_apc.update(initial_prev)
        prev = df_apc.loc[:, 'prevalence'].values
        frac = (1 - df_apc.loc[:, 'incidence'].values)

        df = self._initial_rates.copy()
        tables = []
        for year in range(self._year_start, self._year_end + 1):
            df.loc[:, 'incidence'] = prev * (frac ** (year - self._year_start))
            df['year'] = year
            tables.append(df.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco rate data')

        return df

    def sample_i_from(self, rate_dist, samples):
        """Sample the incidence rate."""
        df = sample_tobacco_rate_from(self._year_start, self._year_end,
                                      self._initial_rates, 'incidence',
                                      self._prev, self._apc, 1e3,
                                      rate_dist, samples)
        return df

    def sample_r_from(self, rate_dist, samples):
        """Sample the remission rate."""
        df = sample_tobacco_rate_from(self._year_start, self._year_end,
                                      self._initial_rates, 'remission',
                                      self._prev, None, 0,
                                      rate_dist, samples)
        return df

    def sample_i(self, prng, rate_dist, n):
        """Sample the incidence rate."""
        return sample_tobacco_rate(self._year_start, self._year_end,
                                   self._initial_rates, 'incidence',
                                   self._prev, self._apc, 1e3,
                                   prng, rate_dist, n)

    def sample_r(self, prng, rate_dist, n):
        """Sample the remission rate."""
        return sample_tobacco_rate(self._year_start, self._year_end,
                                   self._initial_rates, 'remission',
                                   self._prev, None, 0,
                                   prng, rate_dist, n)

    def load_initial_tobacco_rates(self):
        data_file = '{}/tobacco_ir_rates.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path)

        df = df.rename(columns={
            'uptake': 'incidence',
            'Cessation': 'remission'})

        df.insert(0, 'year', self._year_start)
        df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in initial tobacco rates data')

        return df

    def load_tobacco_rates_apc(self):
        data_file = '{}/tobacco_uptake_apc.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path)

        apc_col = 'Percentage yearly decrease in uptake in 20 year olds'
        df = df.rename(columns={apc_col: 'incidence'})

        age_min = self._initial_rates['age'].min()
        age_max = self._initial_rates['age'].max()
        apc_tables = []
        for age in range(age_min, age_max + 1):
            df['age'] = age
            apc_tables.append(df.copy())
        df = pd.concat(apc_tables).sort_values(['age', 'sex'])

        # NOTE: only retain non-zero incidence rates for age 20.
        # There is probably a better way to do this.
        df.loc[df['age'] != 20, 'incidence'] = 0.0

        df = df.sort_values(['age', 'sex']).reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco rates APC data')

        return df

    def load_initial_tobacco_prevalence(self):
        data_file = '{}/tobacco_prevalence.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path)

        df = df.fillna(0.0)
        df = df.rename(columns={'never': 'tobacco.no',
                                'current ': 'tobacco.yes',
                                'former': 'tobacco.post'})
        index_cols = ['sex', 'age']
        post_cols = [c for c in df.columns
                     if c not in index_cols and not c.startswith('tobacco.')]

        # Scale each of the post-cessation prevalence columns by the
        # proportion of the population that are former smokers.
        df.loc[:, post_cols] = df.loc[:, post_cols].mul(df['tobacco.post'],
                                                        axis=0)

        rename_to = {c: 'tobacco.{}'.format(str(c).replace('+', '').strip())
                     for c in post_cols}
        df = df.rename(columns=rename_to)

        # Remove the proportion of former smokers, it is no longer required.
        df = df.drop(columns='tobacco.post')

        df.insert(0, 'year', self._year_start)
        df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco prevalence data')

        # Check that each row sums to unity.
        toln = 1e-12
        max_err = (1 - df.iloc[:, 3:].sum(axis=1)).abs().max()
        if max_err > toln:
            raise ValueError('Tobacco prevalence rows do not sum to 1')

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco prevalence data')

        return df

    def load_tobacco_tax_effects(self):
        price_file = '{}/tobacco_tax_price.csv'.format(self.data_dir)
        elast_file = '{}/tobacco_tax_elasticity.csv'.format(self.data_dir)
        price_path = str(pathlib.Path(price_file).resolve())
        elast_path = str(pathlib.Path(elast_file).resolve())
        df_price = pd.read_csv(price_path)
        df_elast = pd.read_csv(elast_path)

        # Retain the tobacco prices and elasticities; they are needed in order
        # to sample the uncertainty about tobacco tax effects.
        self._df_price = df_price
        self._df_elast = df_elast.copy()

        start_price = df_price.loc[0, 'price']
        tables = []
        for i, row in enumerate(df_price.itertuples()):
            df_elast['year'] = row.year
            df_elast['price'] = row.price
            # Tax always has an effect on uptake.
            df_elast['incidence_effect'] = np.exp(- df_elast['Elasticity']
                                                  * np.log(row.price / start_price))
            # Only *tax increases* have an effect on cessation.
            prev_price = row.price if i == 0 else df_price.loc[i - 1, 'price']
            if row.price > prev_price:
                df_elast['remission_effect'] = np.exp(- df_elast['Elasticity']
                                                      * np.log(row.price / prev_price))
            else:
                df_elast['remission_effect'] = 1.0
            tables.append(df_elast.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
        df = df.loc[:, ['year', 'age', 'sex',
                        'incidence_effect', 'remission_effect']]
        df = df.reset_index(drop=True)

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco tax effects data')

        return df

    def load_tobacco_mortality_rr(self):
        data_file = '{}/tobacco_rr_mortality.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path)

        # The first two columns are sex and age.
        num_cols = df.shape[1]
        base_cols = list(df.columns.values[:2])
        post_cols = ['tobacco.{}'.format(n) for n in range(num_cols - 2)]
        df.columns = base_cols + post_cols
        df = df.fillna(1.0)
        final_col = 'tobacco.{}'.format(num_cols - 2)
        df[final_col] = 1.0
        df.insert(0, 'year', self._year_start)
        df.insert(3, 'tobacco.no', 1.0)
        df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

        # NOTE: the relative risk for a current smoker is the same as that of
        # someone who has stopped smoking one year later (i.e., the values in
        # the 'post_0' column, but shifted up by 1. Here, we shift up by two
        # to skip over the strata of the other sex.
        df.insert(4, 'tobacco.yes', df['tobacco.0'].shift(-2))
        df.loc[df['age'] == df['age'].max(), 'tobacco.yes'] = 1.0
        df.loc[df['tobacco.yes'].isna(), 'tobacco.yes'] = 1.0

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco mortality RR data')

        return df

    def load_tobacco_diseases_rr(self):
        data_file = '{}/tobacco_rr_disease.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path, header=None, prefix='C', comment='#')

        # The first two columns are sex and age.
        sex = df.iloc[2:, 0]
        age = df.iloc[2:, 1].astype(int)

        # Create the base disease table.
        out = pd.DataFrame(data={'year': self._year_start, 'age': age, 'sex': sex})

        # Extract the first row, which lists the diseases.
        disease_headers = df.iloc[0, 2:].fillna('').str.strip()
        # Extract the second row, which lists the data columns.
        column_headers = df.iloc[1, 2:]
        # Extract the data values.
        df = df.iloc[2:, 2:].astype(float).fillna(1.0)

        diseases = {}

        # Extract the initial rates for each chronic and acute disease.
        while len(disease_headers) > 0:
            disease = disease_headers[0]
            if not isinstance(disease, str):
                raise ValueError('Invalid disease name: {}'.format(disease))
            # Check where the next disease begins.
            disease_ixs = np.where(disease_headers.str.len() > 0)[0]
            if len(disease_ixs) > 1:
                end_ix = disease_ixs[1]
                dis_df = df.iloc[:, :end_ix]
                dis_cols = column_headers[:end_ix]
                disease_headers = disease_headers[end_ix:]
                column_headers = column_headers[end_ix:]
                df = df.iloc[:, end_ix:]
            else:
                dis_df = df
                dis_cols = column_headers
                disease_headers = []
                column_headers = []

            if disease == 'RR current':
                # NOTE: these are the RRs for current smokers; 'dis_cols' contains
                # the list of diseases.
                for ix, dis_name in enumerate(dis_cols):
                    # NOTE: correct a typographic error in the input file.
                    if dis_name == 'IHD':
                        dis_name = 'CHD'
                    dis_name = dis_name.replace(' ', '').replace('cancer', 'Cancer')
                    no_col = '{}_no'.format(dis_name)
                    yes_col = '{}_yes'.format(dis_name)
                    out[no_col] = 1.0
                    out[yes_col] = dis_df.iloc[:, ix]
                continue

            # NOTE: correct a typographic error in the input file.
            if disease == 'IHD':
                disease = 'CHD'

            disease = disease.replace(' ', '').replace('cancer', 'Cancer')

            if disease in diseases:
                msg = 'Duplicate RRs for disease {}'.format(disease)
                raise ValueError(msg)

            # NOTE: these are the post-cessation disease RRs.
            dis_cols = [c.replace('+', '') if isinstance(c, str) else int(c)
                        for c in dis_cols]
            col_names = ['{}_{}'.format(disease, c) for c in dis_cols]

            dis_df.columns = col_names
            for ix, col_name in enumerate(col_names):
                out[col_name] = dis_df.iloc[:, ix]

            # NOTE: there are values > 1.0 for the RR For younger ages, where
            # the actual RR is set to 1.0 in the full table.
            yes_col = '{}_yes'.format(disease)
            post_0_col = '{}_0'.format(disease)
            out[yes_col] = out[post_0_col]

            dis_df.insert(0, 'sex', sex)
            dis_df.insert(0, 'age', age)
            dis_df = dis_df.sort_values(['age', 'sex']).reset_index(drop=True)

            if np.any(dis_df.isna()):
                raise ValueError('NA values found in tobacco disease RR data')

            diseases[disease] = dis_df

        out = out.sort_values(['age', 'sex']).reset_index(drop=True)

        return (diseases, out)

    def load_tobacco_disease_rr_gamma(self):
        data_file = '{}/tobacco_rr_disease_gamma.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path)

        df.columns = [col.replace(' ', '').replace('cancer', 'Cancer')
                      for col in df.columns]

        df = df.astype(float).fillna(0.0)
        df = df.sort_values(['age']).reset_index(drop=True)

        return df

    def load_tobacco_disease_rr_sd(self):
        data_file = '{}/tobacco_rr_disease_95pcnt_ci.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path, header=[0, 1], index_col=0)

        # Convert sex into its own column.
        df = df.stack(level=0).reset_index()
        df = df.rename(columns={'level_0': 'disease', 'level_1': 'sex'})

        # Rename diseases and sexes as needed.
        df['disease'] = df['disease'].apply(
            lambda dis: dis.replace(' ', '').replace('cancer', 'Cancer'))
        df['sex'] = df['sex'].apply(lambda s: s.lower())

        # Calculate the SD for log-normal distributions.
        df['sd'] = (np.log(df['95%HCI']) - np.log(df['95%LCI'])) / (2 * 1.96)
        df = df.loc[:, ['disease', 'sex', 'RR', 'sd']]

        if np.any(df.isna()):
            raise ValueError('NA values found in tobacco disease RR SD')

        return df
