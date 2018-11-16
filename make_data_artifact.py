#!/usr/bin/env python
"""This script bundles the input data into a single artifact."""

import os
import os.path
import sys

import pandas as pd

from vivarium_public_health.dataset_manager import hdf
from vivarium_public_health.dataset_manager.artifact import Artifact

import mslt_port.data as data


def main(args=None):
    artifact_pattern = 'input_{}_{}.hdf'
    bau_label = {True: 'const', False: 'decr'}
    delay_label = {True: '0yrs', False: '20yrs'}

    for const_bau in [True, False]:
        for zero_delay in [True, False]:
            artifact_file = artifact_pattern.format(bau_label[const_bau],
                                                    delay_label[zero_delay])

            build_artifact(artifact_file, const_bau=const_bau,
                           zero_delay=zero_delay)

    return 0


def old_build_process():
    artifact_file = 'input_data.hdf'
    if os.path.exists(artifact_file):
        os.remove(artifact_file)

    hdf.touch(artifact_file, append=False)

    art = Artifact(artifact_file)

    # Store the basic population data.
    art.write('population.structure', data.get_base_population())
    # Ensure that the disability rate data are sorted as required.
    yld_data = data.get_yld_rate().sort_values(by=['age', 'sex']).reset_index(drop=True)
    art.write('cause.all_causes.disability_rate', yld_data)
    art.write('cause.all_causes.mortality', data.get_all_cause_mortality())

    # Store data for each chronic disease.
    for disease in ['chd', 'stroke', 'lungC', 'colorectC']:
        art.write('chronic_disease.{}.incidence'.format(disease),
                  data.get_incidence(disease))
        art.write('chronic_disease.{}.remission'.format(disease),
                  data.get_remission(disease))
        art.write('chronic_disease.{}.mortality'.format(disease),
                  data.get_excess_mortality(disease))
        art.write('chronic_disease.{}.morbidity'.format(disease),
                  data.get_disability(disease))
        art.write('chronic_disease.{}.prevalence'.format(disease),
                  data.get_prevalence(disease))

    # Store data for each acute disease/event.
    for disease in ['lrti', 'rta_injury']:
        art.write('acute_disease.{}.mortality'.format(disease),
                  data.get_acute_excess_mortality(disease))
        art.write('acute_disease.{}.morbidity'.format(disease),
                  data.get_acute_disability(disease))

    # Store data for tobacco use.
    for exposure in ['tobacco']:
        art.write('risk_factor.{}.incidence'.format(exposure),
                  data.get_delayed_incidence(exposure))
        art.write('risk_factor.{}.remission'.format(exposure),
                  data.get_delayed_remission(exposure))
        art.write('risk_factor.{}.prevalence'.format(exposure),
                  data.get_delayed_prevalence(exposure))
        art.write('risk_factor.{}.mortality_relative_risk'.format(exposure),
                  data.get_delayed_mortality_rr(exposure))
        disease_rrs = data.get_delayed_disease_rr(exposure, single_table=True)
        art.write('risk_factor.{}.disease_relative_risk'.format(exposure),
                  disease_rrs)

    # TODO: generate and store data for the tobacco tax intervention.
    #
    # * Reduce incidence rate (and increase remission rate?) according to some
    #   lookup table.

    print(art)

    return(0)


def build_artifact(artifact_file, const_bau=False, zero_delay=False):
    if os.path.exists(artifact_file):
        os.remove(artifact_file)

    hdf.touch(artifact_file, append=False)

    art = Artifact(artifact_file)

    # Store the basic population data.
    art.write('population.structure', data.get_base_population())
    # Ensure that the disability rate data are sorted as required.
    yld_data = data.get_yld_rate().sort_values(by=['age', 'sex']).reset_index(drop=True)
    art.write('cause.all_causes.disability_rate', yld_data)
    art.write('cause.all_causes.mortality', data.get_all_cause_mortality())

    # Store data for each chronic disease.
    for disease in ['chd', 'stroke', 'lungC', 'colorectC']:
        art.write('chronic_disease.{}.incidence'.format(disease),
                  data.get_incidence(disease))
        art.write('chronic_disease.{}.remission'.format(disease),
                  data.get_remission(disease))
        art.write('chronic_disease.{}.mortality'.format(disease),
                  data.get_excess_mortality(disease))
        art.write('chronic_disease.{}.morbidity'.format(disease),
                  data.get_disability(disease))
        art.write('chronic_disease.{}.prevalence'.format(disease),
                  data.get_prevalence(disease))

    # Store data for each acute disease/event.
    for disease in ['lrti', 'rta_injury']:
        art.write('acute_disease.{}.mortality'.format(disease),
                  data.get_acute_excess_mortality(disease))
        art.write('acute_disease.{}.morbidity'.format(disease),
                  data.get_acute_disability(disease))

    # Store data for tobacco use.
    for exposure in ['tobacco']:
        incidence = data.get_delayed_incidence(exposure)
        remission = data.get_delayed_remission(exposure)
        if const_bau:
            # Use the initial incidence and remission rates in future years.
            incidence = always_use_initial_rates(incidence)
            remission = always_use_initial_rates(remission)
        art.write('risk_factor.{}.incidence'.format(exposure), incidence)
        art.write('risk_factor.{}.remission'.format(exposure), remission)
        art.write('risk_factor.{}.prevalence'.format(exposure),
                  data.get_delayed_prevalence(exposure))
        table = data.get_delayed_mortality_rr(exposure)
        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            table = update_mortality_rr(table, exposure, years_to_baseline=0)
        art.write('risk_factor.{}.mortality_relative_risk'.format(exposure),
                  table)
        disease_rrs = data.get_delayed_disease_rr(exposure, single_table=True)
        if zero_delay:
            # Instantly revert to a relative risk of 1.0.
            disease_rrs = update_disease_rr(disease_rrs, years_to_baseline=0)
        art.write('risk_factor.{}.disease_relative_risk'.format(exposure),
                  disease_rrs)

    print(artifact_file)


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
    data = data[keep_cols]

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
    suffixes = ['no', 'yes'] + ['post_' + str(y) for y in range(years_to_baseline + 1)]

    for disease in diseases:
        prefix = '{}_'.format(disease)
        dis_cols = [prefix + suffix for suffix in suffixes]
        drop_cols = [c for c in all_columns
                     if c.startswith(prefix) and c not in dis_cols]
        data = data.drop(columns=drop_cols)
        final_col = prefix + 'post_' + str(years_to_baseline)
        data.loc[:, final_col] = 1.0

    return data


def always_use_initial_rates(data):
    """
    Use the initial incidence and remission rates in all future years, rather
    than using projected rates into the future.
    """

    year_min = data['year'].min()
    year_max = data['year'].max()

    print(data.dtypes)
    tables = []
    df = data.loc[data['year'] == year_min]

    for year in range(year_min, year_max + 1):
        df['year'] = year
        tables.append(df.copy())
    data = pd.concat(tables)

    data = data.sort_values(['year', 'age', 'sex']).reset_index(drop=True)
    print(data.dtypes)
    return data


if __name__ == "__main__":
    sys.exit(main())
