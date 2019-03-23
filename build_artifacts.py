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

import mslt_port
from mslt_port import Normal, LogNormal, Beta


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
        artifact_prefix = 'test-tobacco_data_{}'.format(population)
        if bin_edges:
            artifact_prefix = 'test-tobacco_bin_data_{}'.format(population)
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

    popn = mslt_port.Population(data_dir, year_start)
    dis = mslt_port.Diseases(data_dir, year_start, popn.year_end)
    tob = mslt_port.Tobacco(data_dir, year_start, popn.year_end)

    for const_lbl in bau_label:
        for zero_delay, zero_lbl in delay_label.items():
            artifact_file = artifact_pattern.format(const_lbl, zero_lbl)
            artifact_file = os.path.join(out_dir, artifact_file)
            build_artifact(artifact_file, df_base, df_dis,
                           df_tob, df_tob_prev, df_tob_rr_m, df_tob_rr_d,
                           df_tob_tax,
                           popn, dis, tob,
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
    # Ignore strata that have already reached the terminal age.
    year_end = year_start + df_dis['age'].max() - df_dis['age'].min() - 1

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

    # Ignore strata that have already reached the terminal age.
    year_end = year_start + df['age'].max() - df['age'].min() - 1

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
    # Ignore strata that have already reached the terminal age.
    year_end = year_start + df_tob['age'].max() - df_tob['age'].min() - 1
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
    # Ignore strata that have already reached the terminal age.
    year_end = year_start + df_tob['age'].max() - df_tob['age'].min() - 1
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
        # NOTE: there are values > 1.0 for the RR For younger ages, where the
        # actual RR is set to 1.0 in the full table.
        yes_col = '{}_yes'.format(disease)
        post_0_col = '{}_0'.format(disease)
        out[yes_col] = out[post_0_col]

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


def print_basic_comparison(df, compare):
    print('df')
    print(df.head().to_string())
    print('compare')
    print(compare.head().to_string())
    merged = df.merge(compare, indicator=True, how='outer')
    left = merged.loc[merged['_merge'] == 'left_only']
    right = merged.loc[merged['_merge'] == 'right_only']
    print(left.head(4).to_string())
    print(left.tail(4).to_string())
    print(right.head(4).to_string())
    print(right.tail(4).to_string())
    print()


def compare_tables(path, df, compare):
    if 'age' in compare.columns and 'year' in compare.columns:
        print('Old comparison: {}'.format(path))
        is_equal = df.equals(compare)
        if not is_equal:
            print('DATA DO NOT MATCH: {}'.format(path))
            print_basic_comparison(df, compare)
            raise ValueError()
    else:
        # NOTE: we need to compare the non-bin-edge tables to draw #0
        # from tables with bin edges. So we need to iterate over each
        # year, pull out the rate(s), and compare them to draw #0 from the
        # correct year bin.
        is_equal = True
        # NOTE: assume that each age group only spans a single year
        # (i.e., that age_group_end == age_group_start + 1)
        compare = compare.rename(columns={
            'age_group_start': 'age',
        })
        compare = compare.drop(columns='age_group_end')
        # Determine which columns we will be comparing.
        cmp_cols = [c for c in compare.columns if c in df.columns]
        if not any(c not in ['age', 'sex'] for c in cmp_cols):
            raise ValueError('No value column to compare')
        # Treat all non-index columns as numeric for the purpose of
        # *approximate* comparison.
        index_cols = ['age', 'sex']
        num_cols = [c for c in cmp_cols if c not in index_cols]
        show_ignored_cols = False
        if show_ignored_cols:
            ign_cols = [c for c in compare.columns if c not in cmp_cols]
            ign_cols2 = [c for c in df.columns if c not in cmp_cols]
            if not set(ign_cols) <= {'year_start', 'year_end', 'draw'}:
                print('{}: {}'.format(path, ign_cols))
            if set(ign_cols2) != {'year'}:
                print('{}: {}'.format(path, ign_cols2))
        # Identify when we're comparing to an expected rate table.
        if 'draw' not in compare.columns:
            print('NO DRAWS FOR {}'.format(path))
        # Compare the data for each year in turn.
        for year in df['year'].unique():
            mask = ((compare['year_start'] <= year)
                    & (compare['year_end'] > year))
            if 'draw' in compare.columns:
                mask = mask & (compare['draw'] == 0)
            # Extract the relevant subset of both tables.
            df_in = df.loc[df['year'] == year, cmp_cols]
            df_in = df_in.reset_index(drop=True)
            df_cmp = compare.loc[mask, cmp_cols].reset_index(drop=True)
            # Check for NA values (such as NaN and None).
            if np.any(df_in.isna()):
                raise ValueError('NAs in input table')
            if np.any(df_cmp.isna()):
                raise ValueError('NAs in comparison table')
            # Ensure the two tables have identical dimensions.
            if df_in.shape != df_cmp.shape:
                print(year)
                print(df_in.head(1).to_string())
                print(df_cmp.head(1).to_string())
                msg ='Incompatible shapes: {} and {}'.format(
                    df_in.shape, df_cmp.shape)
                raise ValueError(msg)
            else:
                # Compare the index columns for exact equality.
                result = df_in.loc[:, index_cols].equals(
                    df_cmp.loc[:, index_cols])
                # Compare the numeric columns for approximate equality.
                result = result and np.allclose(
                    df_in.loc[:, num_cols],
                    df_cmp.loc[:, num_cols])
                if not result:
                    msg = 'Numerical columns differ: {}'.format(num_cols)
                    raise ValueError(msg)
            is_equal = is_equal and result
            if not result:
                print_basic_comparison(df_in, df_cmp)
                raise ValueError('Data for {} do not match'.format(year))


def write_table(artifact, path, df, bin_edges=False, verbose=False, compare=None):
    """
    Add a data table to an existing artifact.
    """
    if compare is not None:
        compare_tables(path, df, compare)
    else:
        print('NOT COMPARING: {}'.format(path))
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
                   popn, dis, tob,
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

    prng = np.random.RandomState(seed=49430)
    num_draws = 5

    def new_samples():
        return prng.random_sample(num_draws)

    yld_dist = LogNormal(sd_pcnt=10)
    yld_smp = new_samples()

    # Store the basic population data.
    write_table(art, 'population.structure',
                get_population(df_base), bin_edges=False,
                compare=popn.get_population())
    write_table(art, 'cause.all_causes.disability_rate',
                get_disability_rate(df_base), bin_edges=bin_edges,
                compare=popn.sample_disability_rate_from(yld_dist, yld_smp))
    write_table(art, 'cause.all_causes.mortality',
                get_mortality_rate(df_base), bin_edges=bin_edges,
                compare=popn.get_mortality_rate())

    # Identify the chronic and acute diseases for which we have data.
    chr_suffix = '_f'
    chronic_diseases = [c.replace(chr_suffix, '') for c in df_dis.columns
                        if c.endswith(chr_suffix)]
    acu_suffix = '_excess_mortality'
    acute_diseases = [c.replace(acu_suffix, '') for c in df_dis.columns
                      if c.endswith(acu_suffix)]

    for disease in chronic_diseases:
        rename_tbl = {
            'i': '{}_i'.format(disease),
            'r': '{}_r'.format(disease),
            'f': '{}_f'.format(disease),
            'DR': '{}_DR'.format(disease),
            'prev': '{}_prev'.format(disease),
        }
        apc_samples = new_samples()
        i_samples = new_samples()
        r_samples = new_samples()
        f_samples = new_samples()
        yld_samples = new_samples()
        prev_samples = new_samples()
        irfprev_dist = Normal(sd_pcnt=5)
        yld_dist = Normal(sd_pcnt=10)
        apc_dist = Normal(sd_pcnt=0.5)

        cmp_cols = ['year_start', 'year_end',
                    'age_group_start', 'age_group_end',
                    'sex', 'draw']
        inc = '{}_i'.format(disease)
        cmp_dis = dis.chronic[disease].sample_i_from(
            irfprev_dist, apc_dist, i_samples, apc_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        write_table(art, 'chronic_disease.{}.incidence'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', inc]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [inc]])
        cmp_dis = dis.chronic[disease].sample_r_from(
            irfprev_dist, apc_dist, r_samples, apc_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        rem = '{}_r'.format(disease)
        write_table(art, 'chronic_disease.{}.remission'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', rem]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [rem]])
        cfr = '{}_f'.format(disease)
        cmp_dis = dis.chronic[disease].sample_f_from(
            irfprev_dist, apc_dist, f_samples, apc_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        write_table(art, 'chronic_disease.{}.mortality'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', cfr]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [cfr]])
        cmp_dis = dis.chronic[disease].sample_yld_from(
            yld_dist, apc_dist, yld_samples, apc_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        mbd = '{}_DR'.format(disease)
        write_table(art, 'chronic_disease.{}.morbidity'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', mbd]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [mbd]])
        cmp_dis = dis.chronic[disease].sample_prevalence_from(
            irfprev_dist, prev_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        prv = '{}_prev'.format(disease)
        # NOTE: only record the prevalence for the first year.
        df_dis = df_dis.loc[df_dis['year'] == df_dis['year'].min()]
        write_table(art, 'chronic_disease.{}.prevalence'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', prv]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [prv]])

    for disease in acute_diseases:
        rename_tbl = {
            'excess_mortality': '{}_excess_mortality'.format(disease),
            'disability_rate': '{}_disability_rate'.format(disease),
        }
        f_dist = Normal(sd_pcnt=10)
        yld_dist = Normal(sd_pcnt=10)
        f_samples = new_samples()
        yld_samples = new_samples()

        cmp_cols = ['year_start', 'year_end',
                    'age_group_start', 'age_group_end',
                    'sex', 'draw']
        cmp_dis = dis.acute[disease].sample_excess_mortality_from(
            f_dist, f_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        emr = '{}_excess_mortality'.format(disease)
        write_table(art, 'acute_disease.{}.mortality'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', emr]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [emr]])
        cmp_dis = dis.acute[disease].sample_disability_from(
            yld_dist, yld_samples)
        cmp_dis = cmp_dis.rename(columns=rename_tbl)
        mbd = '{}_disability_rate'.format(disease)
        write_table(art, 'acute_disease.{}.morbidity'.format(disease),
                    df_dis.loc[:, ['year', 'age', 'sex', mbd]],
                    bin_edges=bin_edges,
                    compare=cmp_dis.loc[:, cmp_cols + [mbd]])

    for exposure in ['tobacco']:
        incidence = df_tob.loc[:, ['year', 'age', 'sex', 'incidence']]
        remission = df_tob.loc[:, ['year', 'age', 'sex', 'remission']]

        ir_dist = Beta(sd_pcnt=20)
        i_samples = new_samples()
        r_samples = new_samples()

        cmp_cols = ['year_start', 'year_end',
                    'age_group_start', 'age_group_end',
                    'sex', 'draw']
        cmp_rates = tob.sample_i_from(ir_dist, i_samples)
        write_table(art, 'risk_factor.{}.incidence'.format(exposure),
                    incidence, bin_edges=bin_edges,
                    compare=cmp_rates[cmp_cols + ['incidence']])
        cmp_rates = tob.sample_r_from(ir_dist, r_samples)
        write_table(art, 'risk_factor.{}.remission'.format(exposure),
                    remission, bin_edges=bin_edges,
                    compare=cmp_rates[cmp_cols + ['remission']])

        cmp_rr_m = tob.get_expected_mortality_rr()
        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            df_tob_rr_m = update_mortality_rr(df_tob_rr_m, exposure,
                                              years_to_baseline=0)
            cmp_rr_m = update_mortality_rr(cmp_rr_m, exposure,
                                           years_to_baseline=0)
        write_table(art,
                    'risk_factor.{}.mortality_relative_risk'.format(exposure),
                    df_tob_rr_m, bin_edges=bin_edges,
                    compare=cmp_rr_m)

        diseases = list(dis.chronic.keys()) + list(dis.acute.keys())
        samples_tbl = {name: new_samples() for name in diseases}
        cmp_rr_d = tob.sample_disease_rr_from(samples_tbl)
        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            df_tob_rr_d = update_disease_rr(df_tob_rr_d, years_to_baseline=0)
            cmp_rr_d = update_disease_rr(cmp_rr_d, years_to_baseline=0)
        write_table(art,
                    'risk_factor.{}.disease_relative_risk'.format(exposure),
                    df_tob_rr_d, bin_edges=bin_edges,
                    compare=cmp_rr_d,
                    verbose=not zero_delay)

        cmp_prev = tob.get_expected_prevalence()
        if zero_delay:
            # Collapse all former-smokers into a single bin.
            df_tob_prev = update_prevalence(df_tob_prev, exposure,
                                            years_to_baseline=0)
            cmp_prev = update_prevalence(cmp_prev, exposure,
                                         years_to_baseline=0)
        write_table(art, 'risk_factor.{}.prevalence'.format(exposure),
                    df_tob_prev, bin_edges=bin_edges,
                    compare=cmp_prev)

        elast_dist = Normal(sd_pcnt=20)
        samples = new_samples()
        cmp_tax = tob.sample_tax_effects_from(elast_dist, samples)
        cmp_cols = ['year_start', 'year_end',
                    'age_group_start', 'age_group_end',
                    'sex', 'draw']
        tax_inc = df_tob_tax[['year', 'age', 'sex', 'incidence_effect']]
        tax_rem = df_tob_tax[['year', 'age', 'sex', 'remission_effect']]
        write_table(art, 'risk_factor.{}.tax_effect_incidence'.format(exposure),
                    tax_inc, bin_edges=bin_edges,
                    compare=cmp_tax[cmp_cols + ['incidence_effect']])
        write_table(art, 'risk_factor.{}.tax_effect_remission'.format(exposure),
                    tax_rem, bin_edges=bin_edges,
                    compare=cmp_tax[cmp_cols + ['remission_effect']])

    print(artifact_file)


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
