.. _first_sim:

First simulation
================

Once you have :ref:`installed <installation>` this framework, you can run
multi-state life table (MSLT) simulations by:

1. Defining a simulation file, such as the :ref:`example below <basic-sim-1>`.

2. Running the simulation with the ``simulate`` command provided by
   `Vivarium`_:

   .. code-block:: sh

      simulate run example1.yaml

.. code-block:: yaml
   :caption: An example simulation file, ``example1.yaml``.
   :name: basic-sim-1

   components:
       mslt_port:
           population:
               - BasePopulation()
               - Mortality()
               - Disability()
           observer:
               - AdjustedPYandLE('out_py_le.csv')
               - MorbidityMortality('out_morb_mort.csv')
   configuration:
       population:
           population_size: 220 # Males and females, aged 0 to 109.
       time:
           start:
               year: 2011
           end:
               year: 2120
           step_size: 365

Simulation outputs
------------------

Running this simulation will produce two output files:

+ Person-years and life expectancy data will be saved to ``out_py_le.csv``.

+ Morbidity and mortality rates will be saved to ``out_morb_mort.csv``.

Outputs will be written for each :math:`\mathrm{age} \times \mathrm{sex}`
cohort at each year of the simulation.
Cohorts will be identified by their *year of birth* and their *sex*.

Person-years and life expectancy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This output will include the following columns:

population
  The number of people in this cohort.

bau_PYadj
  The adjusted person-years (DALYs or HALYs) for this cohort in the
  business-as-usual (BAU) scenario.

PYadj
  The adjusted person-years (DALYs or HALYs) for this cohort in the
  intervention scenario.

diff_PYadj
  The difference in adjusted person-years between the intervention scenario
  and the BAU scenario.

bau_LEadj
  The adjusted life expectancy for this cohort in the business-as-usual (BAU)
  scenario.

LEadj
  The adjusted life expectancy for this cohort in the intervention scenario.

diff_LEadj
  The difference in adjusted life expectancy between the intervention scenario
  and the BAU scenario.

The first few lines of output produced by the simulation defined in
``example1.yaml`` are shown below:

.. csv-table::
   :header-rows: 1

   year_of_birth,sex,age,year,population,PYadj,LEadj,bau_PYadj,bau_LEadj,diff_LEadj,diff_PYadj
   1901,female,109,2010,207.0,99.12582282825697,0.47886870931525105,99.12582282825697,0.47886870931525105,0.0,0.0
   1901,male,109,2010,86.625,41.79833696929562,0.4825204844940331,41.79833696929562,0.4825204844940331,0.0,0.0
   1902,female,108,2010,207.0,99.12582282825697,0.7233742736210527,99.12582282825697,0.7233742736210527,0.0,0.0
   1902,female,109,2011,105.69212568445639,50.61265181130093,0.4788687093152511,50.61265181130093,0.4788687093152511,0.0,0.0
   1902,male,108,2010,86.625,41.79833696929562,0.7249870374203383,41.79833696929562,0.7249870374203383,0.0,0.0
   1902,male,109,2011,43.52906419976232,21.003665147241186,0.4825204844940331,21.003665147241186,0.4825204844940331,0.0,0.0

   ...

.. note::

   Because ``example1.yaml`` does not include an intervention, the outputs for
   the business-as-usual scenario and for the intervention scenario are
   **identical**.

Morbidity and mortality rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This output will include the following columns:

bau_yld_rate
  The disability rate in the business-as-usual (BAU) scenario.

int_yld_rate
  The disability rate in the intervention scenario.

bau_mortality_rate
  The mortality rate in the business-as-usual (BAU) scenario.

int_mortality_rate
  The mortality rate in the intervention scenario.

The first few lines of output produced by the simulation defined in
``example1.yaml`` are shown below:

.. csv-table::
   :header-rows: 1

   year_of_birth,age,sex,year,bau_yld_rate,bau_mortality_rate,int_yld_rate,int_mortality_rate
   1902,109,female,2011,0.365984528,0.6721884,0.365984528,0.6721884
   1902,109,male,2011,0.357709846,0.6881596,0.357709846,0.6881596
   1903,108,female,2011,0.365984528,0.6721884,0.365984528,0.6721884
   1903,109,female,2012,0.365984528,0.6721884,0.365984528,0.6721884
   1903,108,male,2011,0.357709846,0.6881596,0.357709846,0.6881596
   1903,109,male,2012,0.357709846,0.6881596,0.357709846,0.6881596
   ...

.. note::

   Because ``example1.yaml`` does not include an intervention, the outputs for
   the business-as-usual scenario and for the intervention scenario are
   **identical**.
