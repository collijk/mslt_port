#!/usr/bin/env python3

from vivarium.interface import setup_simulation

from mslt_port import data
from mslt_port.population import BasePopulation, Mortality, Disability
from mslt_port.disease import Disease
from mslt_port.observer import AdjustedPersonYears


config = data.get_default_config()
py_adj = AdjustedPersonYears('PYadj.csv')
sim = setup_simulation([BasePopulation(), Mortality(), Disability(),
                        Disease('chd'), py_adj], config)

sim.run()
# NOTE: this is required to trigger the 'simulation_end' event.
sim.finalize()
