#!/usr/bin/env python

"""
This script processes the various input data files for the tobacco model and
constructs the data artifacts required for each of the model simulations.
"""

import argparse
import os
import os.path
import pathlib
import sys

import numpy as np
import pandas as pd

from vivarium_public_health.dataset_manager import hdf
from vivarium_public_health.dataset_manager.artifact import Artifact


def parser():
    p = argparse.ArgumentParser()
    p.add_argument('--bin-edges', action='store_true',
                   help='Index data tables by age and year intervals')
    p.add_argument('--data-dir', default='./data', metavar='DIR',
                   help='The directory containing the input data files')
    p.add_argument('--output-dir', default='.', metavar='DIR',
                   help='The destination directory for the output files')
    p.add_argument('--start-year', default=2011, type=int, metavar='YEAR',
                   help='The first year of the simulation')

    return p


def main(args=None):
    """
    Construct the data artifacts requires for the tobacco simulations.
    """

    p = parser()
    args = p.parse_args(args)

    # The location of the input data files.
    data_root_dir = args.data_dir
    # The output directory.
    out_dir = args.output_dir
    # The first year of the simulation.
    year_start = args.start_year
    # Whether to index age and years by single values (False) or by intervals
    # (True); newer versions of Vivarium require the use of intervals.
    bin_edges = args.bin_edges

    # Generate artifacts for the Maori and non-Maori populations.
    populations = ['non-maori', 'maori']
    for population in populations:
        data_dir = '{}/{}'.format(data_root_dir, population)
        artifact_prefix = 'tobacco_data_{}'.format(population)
        if bin_edges:
            artifact_prefix = 'tobacco_bin_data_{}'.format(population)
        build_population_artifacts(data_dir, out_dir, artifact_prefix, year_start, bin_edges)


def build_population_artifacts(data_dir, out_dir, artifact_prefix, year_start, bin_edges):
    """
    Build a set of artifacts for a specific population.

    :param data_dir: The directory containing the input data files.
    :param out_dir: The destination directory for the artifact files.
    :param artifact_prefix: The prefix for artifact file names.
    :param year_start: The first year of the simulation.
    :param bin_edges: Whether to index age and year bins as intervals.
    """
    # Load all of the input data.
    df_base = get_base_population_data(data_dir, year_start)
    df_dis = get_diseases(data_dir, year_start)
    df_tob_prev = get_tobacco_prevalence(data_dir, year_start)
    df_tob = get_tobacco_rates(data_dir, year_start, df_tob_prev)
    df_tob_rr_m = get_tobacco_mortality_rr(data_dir, df_tob)
    df_tob_rr_d = get_tobacco_diseases_rr(data_dir, df_tob)
    df_tob_tax = get_tobacco_tax_effects(data_dir)

    # Build all of the required artifacts.
    artifact_pattern = artifact_prefix + '_{}_{}.hdf'
    bau_label = ['decr']
    delay_label = {True: '0yrs', False: '20yrs'}

    for const_lbl in bau_label:
        for zero_delay, zero_lbl in delay_label.items():
            artifact_file = artifact_pattern.format(const_lbl, zero_lbl)
            artifact_file = os.path.join(out_dir, artifact_file)
            build_artifact(artifact_file, df_base, df_dis,
                           df_tob, df_tob_prev, df_tob_rr_m, df_tob_rr_d,
                           df_tob_tax,
                           zero_delay=zero_delay, bin_edges=bin_edges)

    return 0


def get_base_population_data(data_dir, year_start):
    """
    Load the base population data, which comprises the population sizes,
    disability and mortality rates, and the annual percent change (APC) in
    mortality rate.

    :param data_dir: The directory containing the input data files.
    :param year_start: The first year of the simulation.
    """
    data_path = str(pathlib.Path(f'{data_dir}/base_population.csv').resolve())
    df = pd.read_csv(data_path)

    df = df.rename(columns={'mortality per 1 rate': 'mortality_rate',
                            'pYLD rate': 'disability_rate',
                            'APC in all-cause mortality': 'mortality_apc',
                            '5-year': 'population'})

    # Represent each 5-year cohort as a single stratum.
    df['bau_population'] = df['population'].values

    # Retain only the necessary columns.
    df['year'] = year_start
    df = df[['year', 'age', 'sex', 'population', 'bau_population',
             'disability_rate', 'mortality_rate', 'mortality_apc']]

    # Remove strata that have already reached the terminal age.
    df = df[~ (df.age == df['age'].max())]

    return df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)


def get_population(df_base):
    """
    Return the initial population.

    :param df_base: The base population data.
    """
    # Retain only those strata for whom the population size is defined.
    df = df_base.loc[df_base['population'].notna(),
                     ['year', 'age', 'sex', 'population', 'bau_population']]
    return df


def get_acmr_apc(df_base, year_start):
    """
    Return the annual percent change (APC) in mortality rate.

    :param df_base: The base population data.
    :param year_start: The first year of the simulation.
    """
    year_end = year_start + df_base['age'].max() - df_base['age'].min()
    df = df_base[['year', 'age', 'sex', 'mortality_apc']]
    df = df.rename(columns={'mortality_apc': 'value'})

    tables = []
    for year in range(year_start, year_end + 1):
        df['year'] = year
        tables.append(df.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
    df = df.reset_index(drop=True)

    return df


def get_disability_rate(df_base):
    """
    Return the disability rate for each strata.

    :param df_base: The base population data.
    """
    year_start = df_base['year'].unique()
    if len(year_start) == 1:
        year_start = year_start[0]
    else:
        raise ValueError('Invalid starting year: {}'.format(year_start))
    year_end = year_start + df_base['age'].max() - df_base['age'].min()

    df = df_base[['age', 'sex', 'disability_rate']]
    df = df.rename(columns={'disability_rate': 'rate'})

    tables = []
    for year in range(year_start, year_end + 1):
        df['year'] = year
        tables.append(df.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
    df = df.reset_index(drop=True)

    return df


def get_mortality_rate(df_base, apc_num_years=15):
    """
    Return the mortality rate for each strata.

    :param df_base: The base population data.
    """
    year_start = df_base['year'].unique()
    if len(year_start) == 1:
        year_start = year_start[0]
    else:
        raise ValueError('Invalid starting year: {}'.format(year_start))
    year_end = year_start + df_base['age'].max() - df_base['age'].min()

    # NOTE: see column IG in ErsatzInput.
    # - Each cohort has a separate APC (column FE)
    # - ACMR = BASE_ACMR * e^(APC * (year - 2011))
    df_apc = get_acmr_apc(df_base, year_start)
    df_acmr = df_base[['year', 'age', 'sex', 'mortality_rate']]
    df_acmr = df_acmr.rename(columns={'mortality_rate': 'rate'})
    base_acmr = df_acmr['rate'].copy()

    tables = []
    df_acmr['year'] = year_start - 1
    tables.append(df_acmr.copy())
    for counter, year in enumerate(range(year_start, year_end + 1)):
        if counter <= apc_num_years:
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


def get_base_disease_data(data_dir, year_start):
    """
    Load the base disease data, which comprises the disease rates for the
    first year of the simulation.

    Note that many of these rates will change over time.

    :param data_dir: The directory containing the input data files.
    :param year_start: The first year of the simulation.
    """
    data_path = str(pathlib.Path(f'{data_dir}/disease_rates.csv').resolve())
    df = pd.read_csv(data_path, header=None, prefix='C', comment='#')

    # The first two columns are sex and age.
    sex = df.iloc[2:, 0]
    age = df.iloc[2:, 1].astype(int)

    # Create the base disease table.
    out = pd.DataFrame(data={'year': year_start, 'age': age, 'sex': sex})

    # Extract the first row, which lists the diseases.
    disease_headers = df.iloc[0, 2:].fillna('').str.strip()
    # Extract the second row, which lists the data columns.
    column_headers = df.iloc[1, 2:]
    # Extract the data values.
    df = df.iloc[2:, 2:]

    # Identify the columns that are required for each type of disease.
    chronic_cols = ['Incidence', 'Case Fatality', 'prevalence', 'DR']
    acute_cols = ['Mortality', 'DR']

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

        is_chronic = np.all([c in dis_cols.values for c in chronic_cols])
        is_acute = np.all([c in dis_cols.values for c in acute_cols])

        if is_chronic:
            # Extract the chronic disease rates.
            incidence = dis_df['Incidence'].astype(float).fillna(0)
            prevalence = dis_df['prevalence'].astype(float).fillna(0)
            excess_mortality = dis_df['Case Fatality'].astype(float).fillna(0)
            disability_rate = dis_df['DR'].astype(float).fillna(0)
            if 'Remission' in dis_cols.values:
                remission = dis_df['Remission'].astype(float).fillna(0)
            else:
                remission = 0.0 * incidence
            # Determine the column names for these rates.
            i = '{}_i'.format(disease)
            r = '{}_r'.format(disease)
            f = '{}_f'.format(disease)
            dr = '{}_DR'.format(disease)
            prev = '{}_prev'.format(disease)
            # Add these disease rates to the master table.
            # out = out.assign(i = incidence, r = remission,
            #                  f = excess_mortality, dr = disability_rate,
            #                  prev = prevalence)
            out[i] = incidence
            out[r] = remission
            out[f] = excess_mortality
            out[dr] = disability_rate
            out[prev] = prevalence
        elif is_acute:
            # Extract the acute disease rates.
            excess_mortality = dis_df['Mortality'].astype(float)
            disability_rate = dis_df['DR'].astype(float)
            # Determine the column names for these rates.
            f = '{}_excess_mortality'.format(disease)
            dr = '{}_disability_rate'.format(disease)
            # Add these disease rates to the master table.
            # out = out.assign(f = excess_mortality, dr = disability_rate)
            out[f] = excess_mortality
            out[dr] = disability_rate
        else:
            raise ValueError('Invalid columns for disease {}'.format(disease))

    out = out.sort_values(['year', 'age', 'sex']).reset_index(drop=True)
    if np.any(out.isna()):
        raise ValueError('NA values found in disease data')
    return out


def get_disease_rates_apc(data_dir):
    """
    Return the annual percent change (APC) in disease rates.

    :param data_dir: The directory containing the input data files.
    """
    data_path = str(pathlib.Path(f'{data_dir}/disease_rates_apc.csv').resolve())
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


def get_diseases(data_dir, year_start, apc_num_years=15):
    """
    Return the disease rates for all years of the simulation.

    :param data_dir: The directory containing the input data files.
    :param year_start: The first year of the simulation.
    :param apc_num_years: Apply annual percent changes (APCs) to disease rates
        for the first `apc_num_years` of the simulation, after which they will
        remain constant.
    """
    df_dis = get_base_disease_data(data_dir, year_start)
    df_apc = get_disease_rates_apc(data_dir)

    year_start = df_dis['year'].unique()
    if len(year_start) == 1:
        year_start = year_start[0]
    else:
        raise ValueError('Invalid starting year: {}'.format(year_start))
    year_end = year_start + df_dis['age'].max() - df_dis['age'].min()

    apc_tables = []
    for age in range(df_dis['age'].min(), df_dis['age'].max() + 1):
        df_apc['age'] = age
        apc_tables.append(df_apc.copy())
    df_apc = pd.concat(apc_tables).sort_values(['age', 'sex'])
    df_apc = df_apc.reset_index(drop=True)

    # Determine which rates are affected by an annual percent change, and
    # calculate the scaling factor for each of these rates.
    modify_rates = [c for c in df_apc.columns.values
                    if c in df_dis.columns.values
                    and c not in ['year', 'age', 'sex']]
    base_rates = df_dis.loc[:, modify_rates].copy()

    tables = []
    for counter, year in enumerate(range(year_start, year_end + 1)):
        if counter < apc_num_years:
            scale = np.exp(df_apc.loc[:, modify_rates].values * (year - year_start))
            df_dis.loc[:, modify_rates] = base_rates.values * scale
        df_dis['year'] = year
        tables.append(df_dis.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
    df = df.reset_index(drop=True)

    return df


def get_base_tobacco_rates(data_dir, year_start):
    """
    Return the initial incidence and remission rates for tobacco use.

    :param data_dir: The directory containing the input data files.
    :param year_start: The first year of the simulation.
    """
    data_path = str(pathlib.Path(f'{data_dir}/tobacco_ir_rates.csv').resolve())
    df = pd.read_csv(data_path)
    df = df.rename(columns={'uptake': 'incidence', 'Cessation': 'remission'})
    df['year'] = year_start
    df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)
    return df


def get_tobacco_rates_apc(data_dir):
    """
    Return the annual percent change in tobacco rates.

    :param data_dir: The directory containing the input data files.
    """
    data_path = str(pathlib.Path(f'{data_dir}/tobacco_uptake_apc.csv').resolve())
    df = pd.read_csv(data_path)
    apc_col = 'Percentage yearly decrease in uptake in 20 year olds'
    df = df.rename(columns={apc_col: 'incidence'})
    return df

def get_tobacco_rates(data_dir, year_start, df_tob_prev):
    """
    Return the incidence and remission rates for tobacco use in all years of
    the simulation.

    :param data_dir: The directory containing the input data files.
    :param year_start: The first year of the simulation.
    """
    df = get_base_tobacco_rates(data_dir, year_start)
    df_apc = get_tobacco_rates_apc(data_dir)

    year_end = year_start + df['age'].max() - df['age'].min()

    # The incidence rate is calculated with respect to the initial prevalence
    # of new smokers (i.e., those aged 20).
    initial_prev = df_tob_prev.loc[df_tob_prev['age'] == 20]
    initial_prev = initial_prev.rename(columns={'tobacco.yes': 'prevalence'})
    df_apc = df_apc.merge(initial_prev[['sex', 'prevalence']])

    # Set the initial prevalence in 20-year-old cohorts to zero, so that
    # tobacco interventions can have an immediate effect in 2011.
    # Note that this will not affect the 'prevalence' column of df_apc.
    mask = df_tob_prev['age'] == 20
    df_tob_prev.loc[mask, 'tobacco.no'] += df_tob_prev.loc[mask, 'tobacco.yes']
    df_tob_prev.loc[mask, 'tobacco.yes'] = 0.0

    apc_tables = []
    for age in range(df['age'].min(), df['age'].max() + 1):
        df_apc['age'] = age
        apc_tables.append(df_apc.copy())
    df_apc = pd.concat(apc_tables).sort_values(['age', 'sex'])
    # NOTE: only retain non-zero values for age 20.
    # There is probably a better way to do this.
    df_apc.loc[df_apc['age'] != 20, ['incidence', 'prevalence']] = 0.0
    df_apc = df_apc.reset_index(drop=True)

    prev = df_apc.loc[:, 'prevalence'].values
    frac = (1 - df_apc.loc[:, 'incidence'].values)

    tables = []
    for year in range(year_start, year_end + 1):
        df.loc[:, 'incidence'] = prev * (frac ** (year - year_start))
        df['year'] = year
        tables.append(df.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
    df = df.reset_index(drop=True)

    if np.any(df.isna()):
        raise ValueError('NA values found in disease data')

    return df


def get_tobacco_mortality_rr(data_dir, df_tob):
    """
    Return the relative risk of mortality due to tobacco use.

    :param data_dir: The directory containing the input data files.
    :param df_tob: The tobacco rates for each year of the simulation.
    """
    year_start = df_tob['year'].min()
    year_end = year_start + df_tob['age'].max() - df_tob['age'].min()
    if year_end != df_tob['year'].max():
        year_exp = df_tob['year'].max()
        raise ValueError('Invalid final year, {} != {}'.format(year_end,
                                                               year_exp))

    data_path = str(pathlib.Path(f'{data_dir}/tobacco_rr_mortality.csv')
                    .resolve())
    df = pd.read_csv(data_path)
    # The first two columns are sex and age.
    num_cols = df.shape[1]
    base_cols = list(df.columns.values[:2])
    post_cols = ['tobacco.{}'.format(n) for n in range(num_cols - 2)]
    df.columns = base_cols + post_cols
    df = df.fillna(1.0)
    final_col = 'tobacco.{}'.format(num_cols - 2)
    df[final_col] = 1.0
    df.insert(0, 'year', year_start)
    df.insert(3, 'tobacco.no', 1.0)
    df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

    # NOTE: the relative risk for a current smoker is the same as that of
    # someone who has stopped smoking one year later (i.e., the values in the
    # 'post_0' column, but shifted up by 1. Here, we shift up by two to skip
    # over the strata of the other sex.
    df.insert(4, 'tobacco.yes', df['tobacco.0'].shift(-2))
    df.loc[df['age'] == df['age'].max(), 'tobacco.yes'] = 1.0
    df.loc[df['tobacco.yes'].isna(), 'tobacco.yes'] = 1.0

    # Copy these RR columns so that they also apply to the intervention.
    bau_prefix = 'tobacco.'
    int_prefix = 'tobacco_intervention.'
    for col in df.columns:
        if col.startswith(bau_prefix):
            int_col = col.replace(bau_prefix, int_prefix)
            df[int_col] = df[col]

    tables = []
    for year in range(year_start, year_end + 1):
        df['year'] = year
        tables.append(df.copy())

    df = pd.concat(tables).sort_values(['year', 'age', 'sex'])
    df = df.reset_index(drop=True)

    if np.any(df.isna()):
        raise ValueError('NA values found in tobacco mortality RR data')

    return df


def get_tobacco_diseases_rr(data_dir, df_tob):
    """
    Return the relative risks of disease incidence associated with current and
    past tobacco use.

    :param data_dir: The directory containing the input data files.
    :param df_tob: The tobacco rates for each year of the simulation.
    """
    year_start = df_tob['year'].min()
    year_end = year_start + df_tob['age'].max() - df_tob['age'].min()
    if year_end != df_tob['year'].max():
        year_exp = df_tob['year'].max()
        raise ValueError('Invalid final year, {} != {}'.format(year_end,
                                                               year_exp))

    data_path = str(pathlib.Path(f'{data_dir}/tobacco_rr_disease.csv')
                    .resolve())
    df = pd.read_csv(data_path, header=None, prefix='C', comment='#')

    # The first two columns are sex and age.
    sex = df.iloc[2:, 0]
    age = df.iloc[2:, 1].astype(int)

    year_end = year_start + age.max() - age.min()

    # Create the base disease table.
    out = pd.DataFrame(data={'year': year_start, 'age': age, 'sex': sex})

    # Extract the first row, which lists the diseases.
    disease_headers = df.iloc[0, 2:].fillna('').str.strip()
    # Extract the second row, which lists the data columns.
    column_headers = df.iloc[1, 2:]
    # Extract the data values.
    df = df.iloc[2:, 2:].astype(float).fillna(1.0)

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

        # NOTE: these are the post-cessation disease RRs.
        dis_cols = [c.replace('+', '') if isinstance(c, str) else int(c)
                    for c in dis_cols]
        col_names = ['{}_{}'.format(disease, c) for c in dis_cols]
        for ix, col_name in enumerate(col_names):
            out[col_name] = dis_df.iloc[:, ix]

    tables = []
    for year in range(year_start, year_end + 1):
        out['year'] = year
        tables.append(out.copy())

    out = pd.concat(tables)
    out = out.sort_values(['year', 'age', 'sex']).reset_index(drop=True)
    if np.any(out.isna()):
        raise ValueError('NA values found in tobacco disease RR data')
    return out


def get_tobacco_prevalence(data_dir, year_start):
    """
    Return the initial prevalence of tobacco use.

    :param data_dir: The directory containing the input data files.
    :param df_tob: The tobacco rates for each year of the simulation.
    """
    data_path = str(pathlib.Path(f'{data_dir}/tobacco_prevalence.csv')
                    .resolve())
    df = pd.read_csv(data_path).fillna(0.0)
    df = df.rename(columns={'never': 'tobacco.no',
                            'current ': 'tobacco.yes',
                            'former': 'tobacco.post'})
    index_cols = ['sex', 'age']
    post_cols = [c for c in df.columns
                 if c not in index_cols and not c.startswith('tobacco.')]

    # Scale each of the post-cessation prevalence columns by the proportion of
    # the population that are former smokers.
    df.loc[:, post_cols] = df.loc[:, post_cols].mul(df['tobacco.post'],
                                                    axis=0)

    rename_to = {c: 'tobacco.{}'.format(str(c).replace('+', '').strip())
                 for c in post_cols}
    df = df.rename(columns=rename_to)

    # Remove the proportion of former smokers, it is no longer required.
    df = df.drop(columns='tobacco.post')

    df.insert(0, 'year', year_start)
    df = df.sort_values(by=['year', 'age', 'sex']).reset_index(drop=True)

    # Check that each row sums to unity.
    toln = 1e-12
    max_err = (1 - df.iloc[:, 3:].sum(axis=1)).abs().max()
    if max_err > toln:
        raise ValueError('Tobacco prevalence rows do not sum to 1')

    if np.any(df.isna()):
        raise ValueError('NA values found in tobacco prevalence data')

    return df


def get_tobacco_tax_effects(data_dir):
    """
    Calculate the effects of a tobacco tax on incidence and remission rates.
    """
    price_path = str(pathlib.Path(f'{data_dir}/tobacco_tax_price.csv')
                     .resolve())
    elast_path = str(pathlib.Path(f'{data_dir}/tobacco_tax_elasticity.csv')
                    .resolve())
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
    df = df.loc[:, ['year', 'age', 'sex', 'incidence_effect', 'remission_effect']]
    df = df.reset_index(drop=True)

    return df


def define_bin_edges(df):
    """
    Define the lower (inclusive) and upper (exclusive) bounds for the 'age'
    and 'year' columns.
    """

    df = df.rename(columns={
        'age': 'age_group_start',
        'year': 'year_start',
    })

    if 'age_group_start' in df.columns:
        df.insert(df.columns.get_loc('age_group_start') + 1,
                  'age_group_end',
                  df['age_group_start'] + 1)

    if 'year_start' in df.columns:
        df.insert(df.columns.get_loc('year_start') + 1,
                  'year_end',
                  df['year_start'] + 1)

    return df


def write_table(artifact, path, df, bin_edges=False, verbose=False):
    """
    Add a data table to an existing artifact.
    """
    if bin_edges:
        df = define_bin_edges(df)
    if verbose and False:
        print('{} has columns {}'.format(path, df.columns.values))
        print('Year range: {} to {}'.format(df['year_start'].min(),
                                            df['year_start'].max()))
        print('Age range: {} to {}'.format(df['age_group_start'].min(),
                                           df['age_group_start'].max()))
    artifact.write(path, df)


def build_artifact(artifact_file, df_base, df_dis, df_tob, df_tob_prev,
                   df_tob_rr_m, df_tob_rr_d, df_tob_tax,
                   zero_delay=False, bin_edges=False):
    """
    Build a data artifact.

    :param artifact_file: The filename of the output artifact.
    :param df_base: The base population data.
    :param df_base: The disease rates data.
    :param df_tob: The tobacco rates for each year of the simulation.
    :param df_tob_prev: The initial tobacco-use prevalence.
    :param df_tob_rr_m: The relative risk of mortality associated with tobacco
        use.
    :param df_tob_rr_d: The relative risk of diseases associated with tobacco
        use.
    :param df_tob_tax: The effects of a tobacco tax on uptake and cessation.
    :param zero_delay: Whether past smokers immediately receive the complete
        benefits of no longer using tobacco.
    :param bin_edges: Whether to index age and year bins as intervals.
    """
    if os.path.exists(artifact_file):
        os.remove(artifact_file)

    hdf.touch(artifact_file, append=False)

    art = Artifact(artifact_file)

    # Store the basic population data.
    write_table(art, 'population.structure',
                get_population(df_base), bin_edges=False)
    write_table(art, 'cause.all_causes.disability_rate',
                get_disability_rate(df_base), bin_edges=bin_edges)
    write_table(art, 'cause.all_causes.mortality',
                get_mortality_rate(df_base), bin_edges=bin_edges)

    # Identify the chronic and acute diseases for which we have data.
    chr_suffix = '_f'
    chronic_diseases = [c.replace(chr_suffix, '') for c in df_dis.columns
                        if c.endswith(chr_suffix)]
    acu_suffix = '_excess_mortality'
    acute_diseases = [c.replace(acu_suffix, '') for c in df_dis.columns
                      if c.endswith(acu_suffix)]

    for disease in chronic_diseases:
        inc = '{}_i'.format(disease)
        write_table(art, 'chronic_disease.{}.incidence'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', inc]],
                    bin_edges=bin_edges)
        rem = '{}_r'.format(disease)
        write_table(art, 'chronic_disease.{}.remission'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', rem]],
                    bin_edges=bin_edges)
        cfr = '{}_f'.format(disease)
        write_table(art, 'chronic_disease.{}.mortality'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', cfr]],
                    bin_edges=bin_edges)
        mbd = '{}_DR'.format(disease)
        write_table(art, 'chronic_disease.{}.morbidity'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', mbd]],
                    bin_edges=bin_edges)
        prv = '{}_prev'.format(disease)
        write_table(art, 'chronic_disease.{}.prevalence'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', prv]],
                    bin_edges=bin_edges)

    for disease in acute_diseases:
        emr = '{}_excess_mortality'.format(disease)
        write_table(art, 'acute_disease.{}.mortality'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', emr]],
                    bin_edges=bin_edges)
        mbd = '{}_disability_rate'.format(disease)
        write_table(art, 'acute_disease.{}.morbidity'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', mbd]],
                    bin_edges=bin_edges)

    for exposure in ['tobacco']:
        incidence = df_tob.loc[:, ['year', 'age', 'sex', 'incidence']]
        remission = df_tob.loc[:, ['year', 'age', 'sex', 'remission']]
        write_table(art, 'risk_factor.{}.incidence'.format(exposure),
                    incidence, bin_edges=bin_edges)
        write_table(art, 'risk_factor.{}.remission'.format(exposure),
                    remission, bin_edges=bin_edges)

        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            df_tob_rr_m = update_mortality_rr(df_tob_rr_m, exposure,
                                              years_to_baseline=0)
        write_table(art,
                    'risk_factor.{}.mortality_relative_risk'.format(exposure),
                    df_tob_rr_m, bin_edges=bin_edges)

        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            df_tob_rr_d = update_disease_rr(df_tob_rr_d, years_to_baseline=0)
        write_table(art,
                    'risk_factor.{}.disease_relative_risk'.format(exposure),
                    df_tob_rr_d, bin_edges=bin_edges,
                    verbose=not zero_delay)

        if zero_delay:
            # Collapse all former-smokers into a single bin.
            df_tob_prev = update_prevalence(df_tob_prev, exposure,
                                            years_to_baseline=0)
        write_table(art, 'risk_factor.{}.prevalence'.format(exposure),
                    df_tob_prev, bin_edges=bin_edges)

        tax_inc = df_tob_tax[['year', 'age', 'sex', 'incidence_effect']]
        tax_rem = df_tob_tax[['year', 'age', 'sex', 'remission_effect']]
        write_table(art, 'risk_factor.{}.tax_effect_incidence'.format(exposure),
                    tax_inc, bin_edges=bin_edges)
        write_table(art, 'risk_factor.{}.tax_effect_remission'.format(exposure),
                    tax_rem, bin_edges=bin_edges)

    print(artifact_file)


def always_use_initial_rates(data):
    """
    Use the initial incidence and remission rates in all future years, rather
    than using projected rates into the future.
    """
    year_min = data['year'].min()
    year_max = data['year'].max()

    tables = []
    df = data.loc[data['year'] == year_min].copy()

    for year in range(year_min, year_max + 1):
        df.loc[:, 'year'] = year
        tables.append(df.copy())
    data = pd.concat(tables)

    data = data.sort_values(['year', 'age', 'sex']).reset_index(drop=True)
    return data


def update_prevalence(data, exposure, years_to_baseline=0):
    """
    Reduce the delay between remission and a relative risk of 1.0, collapsing
    the prevalence in the final tunnel states.

    :param data: The prevalence data.
    :param exposure: The name of the exposure.
    :param years_to_baseline: The number of years after which the relative
        risk returns to unity.
    """
    all_columns = data.columns

    prefix = '{}.'.format(exposure)
    suffixes = ['no', 'yes'] + [str(y) for y in range(years_to_baseline + 1)]

    prev_cols = [prefix + suffix for suffix in suffixes]
    idx_cols = [c for c in all_columns if not c.startswith(exposure)]

    want_cols = idx_cols + prev_cols
    keep_cols = [c for c in all_columns if c in want_cols]
    data = data.loc[:, keep_cols]

    # Ensure that each row sums to unity.
    final_col = prefix + str(years_to_baseline)
    other_prev_cols = [c for c in prev_cols if c != final_col]
    data.loc[:, final_col] = 1.0 - data.loc[:, other_prev_cols].sum(axis=1)

    return data


def update_mortality_rr(data, exposure, years_to_baseline):
    """
    Reduce the delay between remission and a relative risk of 1.0. Note that
    this will only **truncate** the relative risks, which is most appropriate
    when `years_to_baseline = 0`.

    :param data: The mortality relative risk data.
    :param exposure: The name of the exposure.
    :param years_to_baseline: The number of years after which the relative
        risk returns to unity.
    """
    all_columns = data.columns

    bau_prefix = '{}.'.format(exposure)
    int_prefix = '{}_intervention.'.format(exposure)
    suffixes = ['no', 'yes'] + [str(y) for y in range(years_to_baseline + 1)]

    bau_cols = [bau_prefix + suffix for suffix in suffixes]
    int_cols = [int_prefix + suffix for suffix in suffixes]
    idx_cols = [c for c in all_columns if not c.startswith(exposure)]

    want_cols = idx_cols + bau_cols + int_cols
    keep_cols = [c for c in all_columns if c in want_cols]
    data = data.loc[:, keep_cols]

    bau_final = bau_prefix + str(years_to_baseline)
    int_final = int_prefix + str(years_to_baseline)
    data.loc[:, bau_final] = 1.0
    data.loc[:, int_final] = 1.0

    return data


def update_disease_rr(data, years_to_baseline):
    """
    Reduce the delay between remission and a relative risk of 1.0. Note that
    this will only **truncate** the relative risks, which is most appropriate
    when `years_to_baseline = 0`.

    :param data: The disease incidence relative risk data.
    :param years_to_baseline: The number of years after which the relative
        risk returns to unity.
    """
    all_columns = data.columns

    diseases = list({c[:-4] for c in all_columns if c.endswith('_yes')})
    suffixes = ['no', 'yes'] + ['post_' + str(y)
                                for y in range(years_to_baseline + 1)]

    for disease in diseases:
        prefix = '{}_'.format(disease)
        dis_cols = [prefix + suffix for suffix in suffixes]
        drop_cols = [c for c in all_columns
                     if c.startswith(prefix) and c not in dis_cols]
        data = data.drop(columns=drop_cols)
        final_col = prefix + 'post_' + str(years_to_baseline)
        data.loc[:, final_col] = 1.0

    return data


if __name__ == "__main__":
    sys.exit(main())
