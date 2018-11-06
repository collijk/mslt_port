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

    print(art)

    return(0)


if __name__ == "__main__":
    sys.exit(main())
