import itertools
from pathlib import Path

import pandas as pd
import numpy as np
from vivarium.interface import build_simulation_configuration

YEAR_START = 2011
AGE_GROUP_END = 109
YEAR_END = YEAR_START + AGE_GROUP_END

life_table_data_path = str(Path('../data/inputs_take2.csv').resolve())
disease_data_path = str(Path('../data/inputs_selected_diseases.csv').resolve())


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
    return df[['year', 'age', 'sex', 'incidence']].set_index(['year', 'age', 'sex'])


def get_remission(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'remission']].set_index(['year', 'age', 'sex'])


def get_excess_mortality(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'excess_mortality']].set_index(['year', 'age', 'sex'])


def get_prevalence(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'prevalence']].set_index(['year', 'age', 'sex'])


def get_disability(disease):
    df = get_disease_data(disease)
    return df[['year', 'age', 'sex', 'disability_rate']].set_index(['year', 'age', 'sex'])


def get_intermediary_parameters(disease):
    i = get_incidence(disease)
    r = get_remission(disease)
    f = get_excess_mortality(disease)
    l = i.incidence + r.remission + f.excess_mortality
    l.name = 'l'
    l

    q = np.sqrt(i.incidence ** 2 + r.remission ** 2 + f.excess_mortality ** 2
                + i.incidence * r.remission
                + f.excess_mortality * r.remission
                - i.incidence * f.excess_mortality)
    q.name = 'q'

    w = np.exp(-(l + q) / 2)
    w.name = 'w'
    v = np.exp(-(l - q) / 2)
    v.name = 'v'

    return l, q, w, v


def update_susceptible(S, C, i, r, f, l, q, w, v, t):
    return (2 * (v - w) * (S * (f + r) + C * r) + S * (v * (q - l) + w * (q + l))) / (2 * q)


def update_prevalent(S, C, i, r, f, l, q, w, v, t):
    return (2 * (v - w) * ((f + r) * (S + C) - l * C) - C * q * (v + w)) / (2 * q)
