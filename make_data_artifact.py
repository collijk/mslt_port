#!/usr/bin/env python
"""This script bundles the input data into a single artifact."""

import os
import os.path
import sys

from vivarium_public_health.dataset_manager import hdf
from vivarium_public_health.dataset_manager.artifact import Artifact

import mslt_port.data as data


def main(args=None):
    artifact_file = 'input_data.hdf'
    if os.path.exists(artifact_file):
        os.remove(artifact_file)

    hdf.touch(artifact_file, append=False)

    art = Artifact(artifact_file)

    # Store the basic population data.
    art.write('population.cohort_size', data.get_base_population())
    art.write('population.morbidity_rate', data.get_yld_rate())
    art.write('population.mortality_rate', data.get_all_cause_mortality())

    # Store data for each chronic disease.
    for disease in ['chd', 'stroke', 'lungC', 'colorectC']:
        art.write('chronic_disease.{}.incidence_rate'.format(disease),
                  data.get_incidence(disease))
        art.write('chronic_disease.{}.remission_rate'.format(disease),
                  data.get_remission(disease))
        art.write('chronic_disease.{}.mortality_rate'.format(disease),
                  data.get_excess_mortality(disease))
        art.write('chronic_disease.{}.morbidity_rate'.format(disease),
                  data.get_disability(disease))
        art.write('chronic_disease.{}.prevalence'.format(disease),
                  data.get_prevalence(disease))

    # Store data for each acute disease/event.
    for disease in ['lrti', 'rta_injury']:
        art.write('acute_disease.{}.mortality_rate'.format(disease),
                  data.get_acute_excess_mortality(disease))
        art.write('acute_disease.{}.morbidity_rate'.format(disease),
                  data.get_acute_disability(disease))

    # Store data for tobacco use.
    for exposure in ['tobacco']:
        art.write('exposure.{}.incidence_rate'.format(exposure),
                  data.get_delayed_incidence(exposure))
        art.write('exposure.{}.remission_rate'.format(exposure),
                  data.get_delayed_remission(exposure))
        art.write('exposure.{}.prevalence'.format(exposure),
                  data.get_delayed_prevalence(exposure))
        art.write('exposure.{}.mortality_relative_risk'.format(exposure),
                  data.get_delayed_mortality_rr(exposure))
        disease_rrs = data.get_delayed_disease_rr(exposure)
        for disease in disease_rrs:
            art.write('exposure.{}.{}_relative_risk'.format(exposure, disease),
                      disease_rrs[disease])

    print(art)

    return(0)


if __name__ == "__main__":
    sys.exit(main())
