"""Build disease-specific data tables."""

import pandas as pd
import numpy as np
import pathlib

from .uncertainty import sample_column_long, sample_fixed_rate_from


def sample_disease_rate_from(year_start, year_end, data, rate_name, apc_data,
                             num_apc_years, rate_dist, apc_dist,
                             rate_samples, apc_samples):
    """
    Draw correlated samples for a rate at each year.

    :param year_start: The year at which the simulation starts.
    :param year_end: The year at which the simulation ends.
    :param data: The data table that contains the rate values.
    :param rate_name: The column name that defines the mean values.
    :param apc_data: The data table that contains the annual percent changes
        (set to ``None`` if there are no changes).
    :param num_apc_years: The number of years over which the annual percent
        changes apply (measured from the start of the simulation).
    :param rate_dist: The uncertainty distribution for the rate values.
    :param apc_dist: The uncertainty distribution for the annual percent
        change in this rate.
    :param rate_samples: Samples drawn from the half-open interval [0, 1).
    :param apc_samples: Samples drawn from the half-open interval [0, 1).
    """
    value_col = rate_name

    # Sample the initial rate for each cohort.
    df = sample_column_long(data, value_col, rate_dist, rate_samples)

    df.insert(0, 'year_start', 0)
    df.insert(1, 'year_end', 0)

    df_index_cols = ['year_start', 'year_end', 'age', 'sex', 'draw']
    apc_index_cols = ['age', 'sex']

    tables = []
    years = range(year_start, year_end + 1)

    if apc_data is not None and value_col in apc_data.columns:
        # Sample the annual percent change for each cohort.
        apc = apc_data.loc[:, apc_index_cols + [value_col]]
        # NOTE: this now adds a new column, 'draw'.
        apc = sample_column_long(apc, value_col, apc_dist, apc_samples)

        draw_columns = [c for c in apc.columns if c not in apc_index_cols + ['draw']]
        data_columns = [c for c in df.columns if c not in df_index_cols]
        if set(draw_columns) != set(data_columns):
            raise ValueError('Inconsistent disease parameter draws')

        base_values = df.loc[:, draw_columns].copy().values
        apc_values = apc.loc[:, draw_columns].copy().values

        # Calculate the correlated samples for each cohort at each year.
        for counter, year in enumerate(years):
            df['year_start'] = year
            if counter <= num_apc_years:
                df['year_end'] = year + 1
                timespan = year - year_start
                scale = np.exp(apc_values * timespan)
                df.loc[:, draw_columns] = base_values * scale
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


class Chronic:

    def __init__(self, name, year_start, year_end, data, apc=None):
        self._name = name
        self._year_start = year_start
        self._year_end = year_end
        self._num_apc_years = 15

        data = self._build_data_table(data)
        if apc is not None:
            apc = self._build_apc_table(apc, data['age'].min(),
                                        data['age'].max())

        self._data = data
        self._apc = apc

    def _build_data_table(self, data):
        # NOTE: the name of the remission column was spelt incorrectly for
        # stomach cancer, which prevented it from being renamed to 'r' and
        # ultimately caused it to have a remission rate of zero.
        if 'Remsision' in data.columns:
            data = data.rename(columns={'Remsision': 'Remission'})
        data = data.rename(columns={
            'Incidence': 'i',
            'prevalence': 'prev',
            'Case Fatality': 'f',
            'Remission': 'r',
        })
        if 'r' not in data.columns:
            data['r'] = 0.0 * data['i'].astype(float)

        value_cols = ['i', 'prev', 'f', 'r', 'DR']
        keep_cols = ['age', 'sex'] + value_cols
        data.loc[:, value_cols] = data.loc[:, value_cols].astype(float).fillna(0)
        data = data.loc[:, keep_cols]
        data = data.sort_values(['age', 'sex']).reset_index(drop=True)
        return data

    def _build_apc_table(self, apc, age_min, age_max):
        prefix = '{}_'.format(self._name)
        apc.columns = [c[len(prefix):] if c.startswith(prefix) else c
                       for c in apc.columns]

        tables = []
        for age in range(age_min, age_max + 1):
            apc['age'] = age
            tables.append(apc.copy())

        apc = pd.concat(tables).sort_values(['age', 'sex'])
        apc = apc.reset_index(drop=True)
        return apc


    def get_expected_rates(self):
        # NOTE: need to do this separately for each rate, since some rates
        # may have APCs and other rates will not.
        years = range(self._year_start, self._year_end + 1)
        tables = []
        df_tmp = self._data.copy()

        if self._apc is None:
            df_tmp['year_start'] = self._year_start
            df_tmp['year_end'] = self._year_end + 1
            tables.append(df_tmp.copy())
        else:
            modify_rates = [c for c in self._apc.columns
                            if c in df_tmp.columns and c in ['i', 'r', 'f']]
            base_rates = self._data.loc[:, modify_rates]
            for counter, year in enumerate(years):
                df_tmp['year_start'] = year
                if counter <= self._num_apc_years:
                    df_tmp['year_end'] = year + 1
                    timespan = year - self._year_start
                    scale = np.exp(self._apc.loc[:, modify_rates].values * timespan)
                    df_tmp.loc[:, modify_rates] = base_rates.values * scale
                    tables.append(df_tmp.copy())
                else:
                    df_tmp['year_end'] = self._year_end + 1
                    tables.append(df_tmp.copy())
                    break

        df = pd.concat(tables).sort_values(['year_start', 'age', 'sex'])
        df = df.reset_index(drop=True)

        # Replace 'age' with age groups.
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

        return df

    def sample_i_from(self, rate_dist, apc_dist, rate_samples, apc_samples):
        """Sample the incidence rate."""
        df = sample_disease_rate_from(self._year_start, self._year_end,
                                      self._data, 'i',
                                      self._apc, self._num_apc_years,
                                      rate_dist, apc_dist,
                                      rate_samples, apc_samples)
        return df

    def sample_r_from(self, rate_dist, apc_dist, rate_samples, apc_samples):
        """Sample the remission rate."""
        df = sample_disease_rate_from(self._year_start, self._year_end,
                                      self._data, 'r',
                                      self._apc, self._num_apc_years,
                                      rate_dist, apc_dist,
                                      rate_samples, apc_samples)
        return df

    def sample_f_from(self, rate_dist, apc_dist, rate_samples, apc_samples):
        """Sample the case fatality rate."""
        df = sample_disease_rate_from(self._year_start, self._year_end,
                                      self._data, 'f',
                                      self._apc, self._num_apc_years,
                                      rate_dist, apc_dist,
                                      rate_samples, apc_samples)
        return df

    def sample_yld_from(self, rate_dist, apc_dist, rate_samples, apc_samples):
        """Sample the years lost due to disability rate."""
        df = sample_disease_rate_from(self._year_start, self._year_end,
                                      self._data, 'DR',
                                      self._apc, self._num_apc_years,
                                      rate_dist, apc_dist,
                                      rate_samples, apc_samples)
        return df

    def sample_prevalence_from(self, rate_dist, rate_samples):
        """Sample the initial prevalence of disease."""
        df = sample_column_long(self._data, 'prev', rate_dist, rate_samples)
        df.insert(0, 'year_start', self._year_start)
        df.insert(1, 'year_end', self._year_start + 1)
        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)
        return df


class Acute:

    def __init__(self, name, year_start, year_end, data):
        self._name = name
        self._data = self._build_data_table(data)
        self._year_start = year_start
        self._year_end = year_end

    def _build_data_table(self, data):
        data = data.rename(columns={
            'Mortality': 'excess_mortality',
            'DR': 'disability_rate',
        })

        value_cols = ['excess_mortality', 'disability_rate']
        keep_cols = ['age', 'sex'] + value_cols
        data.loc[:, value_cols] = data.loc[:, value_cols].astype(float)
        data = data.loc[:, keep_cols]
        data = data.sort_values(['age', 'sex']).reset_index(drop=True)
        return data

    def get_expected_rates(self):
        years = range(self._year_start, self._year_end + 1)
        tables = []
        df = self._data.copy()

        df = df.rename(columns={'age': 'age_group_start'})
        df.insert(df.columns.get_loc('age_group_start') + 1,
                      'age_group_end',
                      df['age_group_start'] + 1)

        df.insert(0, 'year_start', self._year_start)
        df.insert(1, 'year_end', self._year_end + 1)

        df = df.sort_values(['year_start', 'age_group_start', 'sex'])
        df = df.reset_index(drop=True)

        return df

    def sample_excess_mortality_from(self, rate_dist, samples):
        """Sample the excess mortality rate."""
        col = 'excess_mortality'
        df = sample_fixed_rate_from(self._year_start, self._year_end,
                                    self._data, col, rate_dist, samples)
        return df

    def sample_disability_from(self, rate_dist, samples):
        """Sample the disability rate."""
        col = 'disability_rate'
        df = sample_fixed_rate_from(self._year_start, self._year_end,
                                    self._data, col, rate_dist, samples)
        return df


class Diseases:

    def __init__(self, data_dir, year_start, year_end):
        self._year_start = year_start
        self._year_end = year_end
        self.data_dir = '{}/diseases/'.format(data_dir)
        self._initial_rates = self.load_initial_disease_rates()
        self._apcs = self.load_chronic_disease_rates_apc()
        self.create_disease_objects()

    def load_initial_disease_rates(self):
        data_file = '{}/disease_rates.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path, header=None, prefix='C', comment='#')
        return df

    def create_disease_objects(self):
        df = self._initial_rates

        # The first two columns are sex and age.
        sex = df.iloc[2:, 0]
        age = df.iloc[2:, 1].astype(int)

        # Extract the first row, which lists the diseases.
        disease_headers = df.iloc[0, 2:].fillna('').str.strip()
        # Extract the second row, which lists the data columns.
        column_headers = df.iloc[1, 2:]
        # Extract the data values.
        df = df.iloc[2:, 2:]

        # Identify the columns that are required for each type of disease.
        chronic_cols = ['Incidence', 'Case Fatality', 'prevalence', 'DR']
        acute_cols = ['Mortality', 'DR']

        self.chronic = {}
        self.acute = {}

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

            if disease == 'Leukaemia - to be completed':
                continue

            disease = disease.replace(' ', '')

            # Extract the relevant disease rates.
            dis_df.columns = dis_cols.values
            dis_df['age'] = age
            dis_df['sex'] = sex

            is_chronic = np.all([c in dis_cols.values for c in chronic_cols])
            is_acute = np.all([c in dis_cols.values for c in acute_cols])

            if is_chronic:
                if disease in self.chronic:
                    raise ValueError('Duplicate disease {}'.format(disease))
                # Check if there are annual percent changes for this disease.
                disease_prefix = disease + '_'
                apc_cols = [col for col in self._apcs.columns
                            if col.startswith(disease_prefix)]
                if apc_cols:
                    dis_df_apc = self._apcs.loc[:, ['sex'] + apc_cols]
                else:
                    dis_df_apc = None
                # Use a consistent name for this disease.
                if disease == 'HeadNeckCancer':
                    disease = 'MouthandoropharynxCancer'
                self.chronic[disease] = Chronic(disease,
                                                self._year_start,
                                                self._year_end,
                                                dis_df,
                                                dis_df_apc)
            elif is_acute:
                if disease in self.acute:
                    raise ValueError('Duplicate disease {}'.format(disease))
                self.acute[disease] = Acute(disease,
                                            self._year_start, self._year_end,
                                            dis_df)
            else:
                msg = 'Invalid columns for disease {}'.format(disease)
                raise ValueError(msg)

    def load_chronic_disease_rates_apc(self):
        data_file = '{}/disease_rates_apc.csv'.format(self.data_dir)
        data_path = str(pathlib.Path(data_file).resolve())
        df = pd.read_csv(data_path, header=None, prefix='C', comment='#')

        rate_suffix = {
            'INCIDENCE TRENDS': '_i',
            'CASE FATALITY TRENDS': '_f',
            'REMISSION TRENDS': '_r',
        }

        # NOTE: column 1 contains the gender, we can ignore its contents.
        # Row 1 contains the rate type (incidence, remission, CFR)
        rate_types = df.iloc[0, 1:].fillna('').str.strip()
        # Row 2 contains the regression slopes for males.
        # NOTE: values for female-only cancers are empty, rather than zero.
        male_slope = (df.iloc[1, 1:].fillna('0')
                      .str.strip()
                      .str.replace('^$', '0')
                      .astype(float))
        # Row 3 contains the regression slopes for females.
        female_slope = df.iloc[2, 1:].astype(float)
        # Row 4 contains the disease names.
        # NOTE: strip all whitespace and newline characters.
        disease_names = df.iloc[3, 1:].str.replace('\s{1,}', '')

        out = pd.DataFrame(data={'sex': ['male', 'female']})

        suffix = rate_suffix[rate_types[0]]
        for ix, disease in enumerate(disease_names):
            if len(rate_types[ix]) > 0:
                suffix = rate_suffix[rate_types[ix]]
            rate_name = disease + suffix
            out[rate_name] = [male_slope[ix], female_slope[ix]]

        return out
