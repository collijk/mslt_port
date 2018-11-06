import itertools
from pathlib import Path

import pandas as pd
import numpy as np
from vivarium.interface import build_simulation_configuration

YEAR_START = 2011
AGE_GROUP_END = 109
YEAR_END = YEAR_START + AGE_GROUP_END

# Path to the data directory, relative to the repository root.
DATA_DIR = './data'

life_table_data_path = str(Path(f'{DATA_DIR}/inputs_take2.csv').resolve())
disease_data_path = str(Path(f'{DATA_DIR}/inputs_selected_diseases.csv').resolve())


def get_default_config():
    config = build_simulation_configuration()
    config.update({
        'population': {
            'population_size': 2 * (AGE_GROUP_END + 1)
        },
        'time': {
            'start': {'year': YEAR_START},
            'end': {'year': YEAR_END},
            'step_size': 365,  # Days
        }
    })

    return config


def get_annual_percent_change():
    return pd.DataFrame(
        list(itertools.product(range(YEAR_START, YEAR_END + 1),
                               range(AGE_GROUP_END + 1),
                               ['male', 'female'],
                               [0])),
        columns=['year', 'age', 'sex', 'value']
    )


def get_base_mortality():
    df = pd.read_csv(life_table_data_path)
    df['year'] = YEAR_START
    df = df.rename(columns={'mortality per 1\nrate': 'rate'})
    acmr = df[['age', 'sex', 'year', 'rate']]
    acmr = acmr[~(acmr.age == AGE_GROUP_END + 1)]
    return acmr.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)


def get_base_population():
    df = pd.read_csv(life_table_data_path)
    pop = df[['age', 'sex', 'pop_avg_5yr']]
    for age_group in range(0, AGE_GROUP_END, 5):
        age_at = age_group + 2
        for age in range(age_group, age_group + 5):
            pop.loc[pop.age == age, 'pop_avg_5yr'] = pop.loc[pop.age == age_at, 'pop_avg_5yr'].values
    pop = pop[~(pop.age == AGE_GROUP_END + 1)]
    pop = pop.rename(columns={'pop_avg_5yr': 'population'})
    pop['bau_population'] = pop['population']
    return pop.sort_values(by=['age', 'sex']).reset_index(drop=True)


def get_yld_rate():
    df = pd.read_csv(life_table_data_path)
    df = df[['age', 'sex', 'pYLD_rate']].rename(columns={'pYLD_rate': 'rate'})
    return df


def get_disease_data(disease):
    key_columns = ['age', 'sex', 'year']
    diseases = ['chd', 'stroke', 'lungC', 'colorectC']
    if disease not in diseases:
        raise ValueError(f'{disease} must be in {diseases}')
    df = pd.read_csv(disease_data_path)

    cols = [c for c in df.columns if (c in key_columns or disease in c)]
    df = df[cols]

    df = df.rename(columns={f'{disease}_i': 'incidence',
                            f'{disease}_r': 'remission',
                            f'{disease}_f': 'excess_mortality',
                            f'{disease}_prev': 'prevalence',
                            f'{disease}_DR': 'disability_rate', }).fillna(0)

    data = []
    for year in range(YEAR_START, YEAR_END + 1):
        df['year'] = year
        data.append(df.copy())
    return pd.concat(data).sort_values(['year', 'age', 'sex']).reset_index()


def get_all_cause_mortality():
    all_apc = get_annual_percent_change()
    current_mortality = get_base_mortality().set_index(['age', 'sex', 'year'])

    data = [current_mortality]
    for year in range(YEAR_START, YEAR_END - 1):
        apc = all_apc[all_apc.year == year].set_index(['age', 'sex', 'year'])
        new_mortality = pd.DataFrame({'rate': current_mortality.rate * (1 + apc.value / 100)}).reset_index()
        new_mortality['year'] = year + 1
        new_mortality = new_mortality.set_index(['age', 'sex', 'year'])
        data.append(new_mortality)
        current_mortality = new_mortality

    return pd.concat(data).reset_index().sort_values(['year', 'age', 'sex']).reset_index(drop=True)


def get_probability_of_death():
    acmr = get_all_cause_mortality().set_index(['age', 'sex', 'year'])
    p = (1 - np.exp(-acmr)).reset_index()
    return p.rename(columns={'rate': 'probability'}).sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)


def get_population():
    current_pop = get_base_population()
    all_p = get_probability_of_death()

    data = [current_pop]
    for year in range(YEAR_START, YEAR_END):
        new_pop = current_pop.copy()
        pop_correct_age = (new_pop.age >= year - YEAR_START) & (new_pop.age < AGE_GROUP_END)

        p = all_p[all_p.year == year].reset_index(drop=True)
        new_pop.loc[pop_correct_age, 'population'] *= 1 - p['probability']

        new_pop['year'] += 1
        new_pop['age'] += 1
        new_pop.loc[~pop_correct_age, 'population'] = 0
        new_pop.loc[new_pop.age == 111, 'age'] = 0

        new_pop = new_pop.sort_values(by=['age']).reset_index(drop=True)
        data.append(new_pop)
        current_pop = new_pop
    return pd.concat(data).sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)


def get_person_years():
    pop = get_population()
    p = get_probability_of_death()
    py = pop.set_index(['year', 'age', 'sex']).population * (1 - p.set_index(['year', 'age', 'sex']).probability)
    py.name = 'person_years'
    return py.reset_index()


def get_adjusted_person_years():
    py = get_person_years()
    ylds = get_yld_rate().set_index(['age'])

    data = []
    for (sex, year), group in py.groupby(['sex', 'year']):
        group.loc[:, 'person_years'] = (group.set_index(['age']).person_years
                                        * (1 - ylds[ylds.sex == sex].rate)).values
        data.append(group)

    return pd.concat(data, ignore_index=True)


def get_adjusted_life_expectancy():
    py_adj = get_adjusted_person_years().set_index(['age', 'year', 'sex'])
    pop = get_population().set_index(['age', 'year', 'sex'])

    ratio = (py_adj.person_years / pop.population).fillna(0)
    ratio.name = 'ratio'
    ratio = ratio.reset_index()

    data = []
    for sex in ['male', 'female']:
        ratio_by_sex = ratio[ratio.sex == sex].drop(columns='sex')
        ratio_by_sex = ratio_by_sex.pivot(index='age', columns='year', values='ratio')
        le = pd.DataFrame({'age': range(AGE_GROUP_END + 1),
                           'life_expectancy': 0})
        for age in range(le.age.max()):
            le.loc[le.age == age, 'life_expectancy'] = np.sum(np.diagonal(ratio_by_sex, -age))

        le['sex'] = sex
        data.append(le)
    return pd.concat(data).reset_index(drop=True)


def get_incidence(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'incidence']]


def get_remission(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'remission']]


def get_excess_mortality(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'excess_mortality']]


def get_prevalence(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'prevalence']]


def get_disability(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'disability_rate']]


def get_acute_disease_data(disease):
    data_path = str(Path('{}/inputs_{}.csv'.format(DATA_DIR, disease)).resolve())
    df = pd.read_csv(data_path)

    key_columns = ['age', 'sex', 'year', 'excess_mortality', 'disability_rate']
    cols = [c for c in df.columns if c in key_columns]
    df = df[cols].fillna(0)

    if 'year' in cols:
        data = df
    else:
        data = []
        for year in range(YEAR_START, YEAR_END + 1):
            df['year'] = year
            data.append(df.copy())
        data = pd.concat(data)

    return data.sort_values(['year', 'age', 'sex']).reset_index()


def get_acute_excess_mortality(disease):
    df = get_acute_disease_data(disease)
    return df[['year', 'age', 'sex', 'excess_mortality']]


def get_acute_disability(disease):
    df = get_acute_disease_data(disease)
    return df[['year', 'age', 'sex', 'disability_rate']]


def get_delayed_prevalence(risk_name):
    """
    Return the initial prevalence of a delayed risk.

    Prevalence is specified as:

    - The fraction of each strata that:

      (a) have never been exposed ('no');
      (b) are currently exposed ('yes'); and
      (c) were previously exposed ('post').

    - The distribution of previously-exposed is specified separately, in
      columns 'post_0', 'post_1', ..., 'post_N', where the final column
      denotes all previous exposures that ended N or more years ago.
    """
    data_path = str(Path('{}/{}_prevalence.csv'.format(DATA_DIR, risk_name))
                    .resolve())
    df = pd.read_csv(data_path)

    # Note: this table may contain duplicate rows for the final age group.
    subset_cols = ['age', 'sex']
    if 'year' in df.columns:
        subset_cols.append('year')
    df = df.drop_duplicates(subset=subset_cols)

    # Scale each of the 'post_X' columns by the proportion of the strata that
    # were previously exposed (df['post']).
    post_columns = [col for col in df.columns if col.startswith('post_')]
    df.loc[:, post_columns] = df.loc[:, post_columns].mul(df['post'], axis=0)
    df = df.drop(columns=['post'])

    # Rename columns to match those expected by the DelayedRisk class.
    rename_to = {c: c.replace('post_', '{}.'.format(risk_name))
                 for c in post_columns}
    rename_to['yes'] = '{}.yes'.format(risk_name)
    rename_to['no'] = '{}.no'.format(risk_name)
    df = df.rename(columns=rename_to)

    if 'year' in df.columns:
        data = df
    else:
        data = []
        for year in range(YEAR_START, YEAR_END + 1):
            df['year'] = year
            data.append(df.copy())
        data = pd.concat(data)

    data = data.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)
    return data


def get_delayed_ir_rates(risk_name):
    """
    Return the incidence (uptake) and remission (cessation) rates for a
    delayed risk.
    """
    data_path = str(Path('{}/{}_ir_rates.csv'.format(DATA_DIR, risk_name))
                    .resolve())
    df = pd.read_csv(data_path)

    # Note: this table may contain duplicate rows for the final age group.
    subset_cols = ['age', 'sex']
    if 'year' in df.columns:
        subset_cols.append('year')
    df = df.drop_duplicates(subset=subset_cols)

    key_columns = ['age', 'sex', 'year', 'incidence', 'remission']
    cols = [c for c in df.columns if c in key_columns]
    df = df[cols].fillna(0)

    if 'year' in cols:
        data = df
    else:
        data = []
        for year in range(YEAR_START, YEAR_END + 1):
            df['year'] = year
            data.append(df.copy())
        data = pd.concat(data)

    data = data.sort_values(['year', 'age', 'sex']).reset_index()
    return data


def get_delayed_incidence(risk_name):
    """
    Return the incidence (uptake) rate for a delayed risk.
    """
    df = get_delayed_ir_rates(risk_name)
    return df[['year', 'age', 'sex', 'incidence']]


def get_delayed_remission(risk_name):
    """
    Return the remission (cessation) rate for a delayed risk.
    """
    df = get_delayed_ir_rates(risk_name)
    return df[['year', 'age', 'sex', 'remission']]


def get_delayed_mortality_rr(risk_name):
    """
    Return the relative risk of mortality associated with each exposure level
    of a delayed risk.
    """
    data_path = str(Path('{}/{}_rr_mortality.csv'.format(DATA_DIR, risk_name))
                    .resolve())
    df = pd.read_csv(data_path)

    # Rename columns to match those expected by the DelayedRisk class.
    post_columns = [col for col in df.columns if col.startswith('post_')]
    rename_to = {c: c.replace('post_', '{}.'.format(risk_name))
                 for c in post_columns}
    rename_to['yes'] = '{}.yes'.format(risk_name)
    rename_to['no'] = '{}.no'.format(risk_name)
    df = df.rename(columns=rename_to)

    # Copy these risk columns so that they also apply to the intervention.
    bau_prefix = '{}.'.format(risk_name)
    int_prefix = '{}_intervention.'.format(risk_name)
    for col in df.columns:
        if col.startswith(bau_prefix):
            int_col = col.replace(bau_prefix, int_prefix)
            df[int_col] = df[col]

    if 'year' in df.columns:
        data = df
    else:
        data = []
        for year in range(YEAR_START, YEAR_END + 1):
            df['year'] = year
            data.append(df.copy())
        data = pd.concat(data)

    data = data.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)
    return data


def get_delayed_disease_rr(risk_name, single_table=False):
    """
    Return the relative risk of disease incidence associated with each
    exposure level of a delayed risk.

    Note that this function returns a dictionary that maps disease names to
    relative risk data frames, rather than returning a single data frame.

    :param risk_name: The name of the risk.
    :param single_table: Whether to return separate tables for each default
        (the default, `False`) or to return a single table for all diseases
        (`True`).
    """
    data_path = str(Path('{}/{}_rr_diseases.csv'.format(DATA_DIR, risk_name))
                    .resolve())
    df = pd.read_csv(data_path)
    key_columns = ['age', 'sex', 'year']

    # Note: this table may contain duplicate rows for the final age group.
    subset_cols = ['age', 'sex']
    if 'year' in df.columns:
        subset_cols.append('year')
    df = df.drop_duplicates(subset=subset_cols)

    if single_table:
        return df

    # Note: RRs are named "{disease}_{no|yes|post_N}"
    tables = {}
    diseases = [c.replace('_no', '') for c in df.columns if c.endswith('_no')]
    for disease in diseases:
        dis_columns = [c for c in df.columns if c.startswith(disease)]
        dis_keys = [c for c in df.columns if c in key_columns]
        dis_df = df[dis_keys + dis_columns]
        dis_prefix = '{}_'.format(disease)
        bau_prefix = '{}.'.format(risk_name)
        int_prefix = '{}_intervention.'.format(risk_name)
        bau_col = {c: c.replace(dis_prefix, bau_prefix).replace('post_', '')
                   for c in dis_columns}
        int_col = {c: c.replace(dis_prefix, int_prefix).replace('post_', '')
                   for c in dis_columns}
        for column in dis_columns:
            dis_df[int_col[column]] = dis_df[column]
        dis_df = dis_df.rename(columns=bau_col)

        if 'year' in dis_df.columns:
            data = dis_df
        else:
            data = []
            for year in range(YEAR_START, YEAR_END + 1):
                dis_df['year'] = year
                data.append(dis_df.copy())
            data = pd.concat(data)

        data = data.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

        tables[disease] = data

    return tables
