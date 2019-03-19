"""Build risk-factor data tables."""

import pandas as pd
import numpy as np
import pathlib

from .disease import sample_column


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
    df = sample_column(data, value_col, prng, rate_dist, n)

    df.insert(0, 'year', 0)
    df_index_cols = ['year', 'age', 'sex']
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
            df['year'] = year
            if counter < num_apc_years and year > year_start:
                timespan = year - year_start
                result = initial_rate * (frac ** timespan)[..., np.newaxis]
                df.loc[:, data_columns] = result
            tables.append(df.copy())
    else:
        # Calculate the correlated samples for each cohort at each year.
        for year in years:
            df['year'] = year
            tables.append(df.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
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

    status = """
      [x] 'risk_factor.{}.incidence' and sampling
      [x] 'risk_factor.{}.remission' and sampling
      [x] 'risk_factor.{}.prevalence'
      [~] 'risk_factor.{}.mortality_relative_risk' and sampling
      [~] 'risk_factor.{}.disease_relative_risk' and sampling
      [~] 'risk_factor.{}.tax_effect_incidence' and sampling
      [~] 'risk_factor.{}.tax_effect_remission' and sampling
    """

    def sample_tax_effects(self):
        # TODO
        raise NotImplementedError()

    def sample_disease_rr(self):
        # TODO
        raise NotImplementedError()

    def sample_mortality_rr(self):
        # TODO
        raise NotImplementedError()

    def get_expected_tax_effects(self):
        """
        Return the effects of a tobacco tax on incidence and remission rates.
        """
        df = self._tax.copy()

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

        tables = []
        for year in range(self._year_start, self._year_end + 1):
            df['year'] = year
            tables.append(df.copy())

        df = pd.concat(tables)
        df = df.sort_values(['year', 'age', 'sex']).reset_index(drop=True)

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

        tables = []
        for year in range(self._year_start, self._year_end + 1):
            df['year'] = year
            tables.append(df.copy())

        df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
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
        df = df.sort_values(['year', 'age', 'sex']).reset_index(drop=True)
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

            dis_df.insert(0, 'sex', sex)
            dis_df.insert(0, 'age', age)
            dis_df = dis_df.sort_values(['age', 'sex']).reset_index(drop=True)

            if np.any(dis_df.isna()):
                raise ValueError('NA values found in tobacco disease RR data')

            diseases[disease] = dis_df

        out = out.sort_values(['age', 'sex']).reset_index(drop=True)

        return (diseases, out)
