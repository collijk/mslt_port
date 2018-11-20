#!/usr/bin/env python
"""
Use this script to profile a Vivarium simulation.
"""

import cProfile
import sys

from vivarium.interface import setup_simulation_from_model_specification


def run_sim(sim_file):
    sim = setup_simulation_from_model_specification(sim_file)
    sim.run()
    sim.finalize()


def main(args=None):
    sim_file = 'tobacco_decr_20yrs.yaml'
    command = 'run_sim("{}")'.format(sim_file)
    out_file = sim_file.replace('.yaml', '.stats')
    cProfile.run(command, filename=out_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
