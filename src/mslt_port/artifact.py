import numpy as np
import os
import os.path

from vivarium_public_health.dataset_manager import hdf
from vivarium_public_health.dataset_manager.artifact import Artifact

from .population import Population
from .disease import Diseases
from .risk_factor import Tobacco
from .uncertainty import Normal, Beta, LogNormal, LogNormalRawSD


def assemble_chd_only(num_draws):
    artifact_file = 'test_reduce_chd.hdf'

    year_start = 2011
    data_dir_maori = './data/maori'
    data_dir_nonmaori = './data/non-maori'
    prng = np.random.RandomState(seed=49430)

    def new_samples():
        return prng.random_sample(num_draws)

    dist_yld = LogNormal(sd_pcnt=10)
    smp_yld = new_samples()

    dist_chd_apc = Normal(sd_pcnt=0.5)
    smp_chd_apc = new_samples()
    dist_chd_i = Normal(sd_pcnt=5)
    smp_chd_i = new_samples()
    dist_chd_r = Normal(sd_pcnt=5)
    smp_chd_r = new_samples()
    dist_chd_f = Normal(sd_pcnt=5)
    smp_chd_f = new_samples()
    dist_chd_yld = Normal(sd_pcnt=10)
    smp_chd_yld = new_samples()
    dist_chd_prev = Normal(sd_pcnt=5)
    smp_chd_prev = new_samples()

    # Instantiate components for the non-Maori population.
    p_nm = Population(data_dir_nonmaori, year_start)
    l_nm = Diseases(data_dir_nonmaori, year_start, p_nm.year_end)

    # Create a new artifact file, replacing any existing file.
    if os.path.exists(artifact_file):
        os.remove(artifact_file)
    hdf.touch(artifact_file, append=False)
    art = Artifact(artifact_file)

    # 'population.structure'
    popn_nm = p_nm.get_population()
    art.write('population.structure',
              popn_nm)
    # 'cause.all_causes.disability_rate'
    popn_yld_nm = p_nm.sample_disability_rate_from(dist_yld, smp_yld)
    art.write('cause.all_causes.disability_rate',
              popn_yld_nm)
    # 'cause.all_causes.mortality'
    popn_acmr_nm = p_nm.get_mortality_rate()
    art.write('cause.all_causes.mortality',
              popn_acmr_nm)

    name = 'CHD'
    disease = l_nm.chronic[name]

    # 'chronic_disease.{}.incidence'
    chd_i = disease.sample_i_from(
        dist_chd_i, dist_chd_apc,
        smp_chd_i, smp_chd_apc)
    art.write('chronic_disease.CHD.incidence',
              check_for_bin_edges(chd_i))
    # 'chronic_disease.{}.remission'
    chd_r = disease.sample_r_from(
        dist_chd_r, dist_chd_apc,
        smp_chd_r, smp_chd_apc)
    art.write('chronic_disease.CHD.remission',
              check_for_bin_edges(chd_r))
    # 'chronic_disease.{}.mortality'
    chd_f = disease.sample_f_from(
        dist_chd_f, dist_chd_apc,
        smp_chd_f, smp_chd_apc)
    art.write('chronic_disease.CHD.mortality',
              check_for_bin_edges(chd_f))
    # 'chronic_disease.{}.morbidity'
    chd_yld = disease.sample_yld_from(
        dist_chd_yld, dist_chd_apc,
        smp_chd_yld, smp_chd_apc)
    art.write('chronic_disease.CHD.morbidity',
              check_for_bin_edges(chd_yld))
    # 'chronic_disease.{}.prevalence'
    chd_prev = disease.sample_prevalence_from(
        dist_chd_prev, smp_chd_prev)
    art.write('chronic_disease.CHD.prevalence',
              check_for_bin_edges(chd_prev))


def check_for_bin_edges(df):
    """
    Check that lower (inclusive) and upper (exclusive) bounds for year and age
    are defined as table columns.
    """

    if 'age_group_start' in df.columns and 'year_start' in df.columns:
        return df
    else:
        raise ValueError('Table does not have bins')
