.. _simulations:

Defining simulations
====================

A simulation is defined in terms of the **components** that are included in
the model, and the **configuration settings** that define key model
parameters.

Overview
--------

The components and configuration settings are collected into a single text
file, as shown in the example file below.

.. literalinclude:: ../../../test_adjPY.yaml
   :caption: An example simulation file.
   :language: yaml
   :lines: 11-21,23-
   :linenos:

In order to define your own simulation, you will need to provide data to
characterise:

- The :ref:`target population <popn_comp>` (lines 4--6 in the above example).
- Each :ref:`chronic disease <chronic_comp>` that you want to include (line 8
  in the above example).
- Any :ref:`acute diseases or events <acute_comp>` that you want to include
  (not shown in the above example).
- The :ref:`intervention(s) <interv_comp>` that will act on one or more
  disease rates (not shown in the above example).

You can then use the provided MSLT_ components (shown in the figure below) to
incorporate these features into a simulation.

.. note:: You may also want to adjust the :ref:`number of cohorts <popn_conf>`
   or the :ref:`simulation timescale <time_conf>`.

To run a simulation (e.g., as defined in the file ``simulation_1.yaml``), open
a terminal (Linux and OS X) or command prompt (Windows) and run the following
command:

.. code-block:: shell

   simulate run simulation_1.yaml

.. graphviz::
   :caption: The structure of an MSLT_ simulation, where each box identifies a
       component in the MSLT_ framework.
       The Life Table component (dashed box) is automatically included in
       every simulation, and its contents are defined by the other simulation
       components.
       User-provided inputs are shown in blue, simulation outputs are shown in
       green.

   digraph mslt_structure {
     graph [bgcolor="transparent", fontsize=8, scale=0.5];
     input [label="All-cause Mortality Rate\nDisability Rate\nDemographics",
            color="#1f78b4", margin="0.11"];
     popn [label="Population",
           href="../tutorial/simulation.html#population", target="_top",
           shape=box, margin="0.22,0.011"];
     mslt [label="Life Table",
           style="dashed",
           shape=box, margin="0.22,0.011"];
     obs [label="Observer",
          href="../tutorial/simulation.html#observers", target="_top",
          shape=box, margin="0.22,0.011"];
     output [label="Adjusted life-years\nAdjusted life-expectancy",
             color="#33a02c", margin="0.11"];
     input -> popn [color="#1f78b4"];
     popn -> mslt -> obs;
     obs -> output [color="#33a02c"];
     delta [label="ΔACMR\nΔYLDR", color="transparent"];
     chronic [label="Chronic Disease",
              href="../tutorial/simulation.html#chronic-diseases",
              target="_top",
              shape=box, margin="0.22,0.011"];
     chronic_input [
       label="Incidence\nRemission\nPrevalence\nMortality\nDisability",
       color="#1f78b4", margin="0.11"];
     chronic_input -> chronic [color="#1f78b4"];
     chronic -> delta -> mslt;
     acute [label="Acute Disease",
            href="../tutorial/simulation.html#acute-diseases-and-events",
            target="_top",
            shape=box, margin="0.22,0.011"];
     acute_input [label="Mortality\nDisability",
                  color="#1f78b4", margin="0.11"];
     acute_input -> acute [color="#1f78b4"];
     acute -> delta;
     pif [label="PIF", color="#1f78b4", margin="0.11"];
     intervention [label="Intervention",
                   href="../tutorial/simulation.html#interventions",
                   target="_top",
                   shape=box, margin="0.22,0.011"];
     pif -> intervention [color="#1f78b4"];
     intervention -> chronic;
     intervention -> acute;
   }

Components
----------

.. _popn_comp:

Population
^^^^^^^^^^

This will generally comprise:

1. The baseline population demographics (age, sex, cohort size);
2. An all-cause mortality rate (ACMR); and
3. A years lost due to disability (YLD) rate.

These data are currently loaded using three separate components:

.. code-block:: yaml

   components:
       mslt_port:
           population:
               - BasePopulation()
               - Mortality()
               - Disability()

.. _chronic_comp:

Chronic diseases
^^^^^^^^^^^^^^^^

Chronic diseases are characterised in terms of:

- Incidence rate;
- Remission rate;
- Case fatality rate;
- Initial prevalence; and
- Disability rate.

.. note:: These rates must be specified for each cohort (i.e., age and sex).

Any number of chronic diseases can be included in a simulation.
Each disease must be given a unique name, which should correspond to the name
used to identify the disease in the input data file.

.. code-block:: yaml

   components:
       mslt_port:
           disease:
               - Disease('chd')
               - Disease('stroke')
               - Disease('lungC')
               - Disease('colorectC')

.. _acute_comp:

Acute diseases and events
^^^^^^^^^^^^^^^^^^^^^^^^^

Acute diseases (such as lower respiratory tract infections) and acute events
(such as road traffic accidents) are characterised in terms of:

- Excess mortality rate; and
- Disability rate.

.. note:: These rates must be specified for each cohort (i.e., age and sex),
   and can also be specified separately for each year of the simulation.

Any number of acute diseases/events can be included in a simulation.
Each disease must be given a unique name, which should correspond to the name
used to identify the disease in the input data file.

.. code-block:: yaml

   components:
       mslt_port:
           disease:
               - AcuteDisease('lrti')

.. _interv_comp:

Interventions
^^^^^^^^^^^^^

Interventions are characterised as one or more Population Impact Fractions
(PIFs), each of which modifies an existing rate (e.g., the incidence rate of a
specific chronic disease).

.. note:: Interventions are a work in progress, we are still deciding exactly
   what form they will take.

.. _obs_comp:

Observers
^^^^^^^^^

Observers are used to record salient details during a simulation.
The MSLT_ framework currently provides one observer,
:class:`~mslt_port.observer.AdjustedPYandLE`.

Any number of observers can be included in a simulation.
Each observer should write to a different output file.

.. code-block:: yaml

   components:
       mslt_port:
           observer:
               - AdjustedPYandLE('adjusted-py-le.csv')

At the end of a simulation, it will write these data to the specified output
file:

.. csv-table::
   :header-rows: 1

   year_of_birth,sex,age,year,population,PYadj,LEadj
   1901,female,109,2010,207.0,99.12582282825697,0.47886870931525105
   1901,male,109,2010,86.625,41.79833696929562,0.4825204844940331
   1902,female,108,2010,207.0,99.12582282825697,0.7233742736210527
   1902,female,109,2011,105.69212568445639,50.61265181130093,0.4788687093152511
   1902,male,108,2010,86.625,41.79833696929562,0.7249870374203383
   1902,male,109,2011,43.52906419976232,21.003665147241186,0.4825204844940331
   ...

Configuration
-------------

.. _popn_conf:

Population
^^^^^^^^^^

The number of cohorts must be defined.
In the example simulations presented here, this is the number of distinct age
groups (110, for single-year cohorts aged from 0 to 109 years) multiplied by 2
(separating the population into male and female cohorts).

.. code-block:: yaml

   population:
       population_size: 220

.. _time_conf:

Simulation timescale
^^^^^^^^^^^^^^^^^^^^

The simulation timescale is defined in terms of the start and end years, and
the time-step size.
In the example simulations presented here, simulations are run for 109 years,
starting from 2011, using a time-step of 1 year:

.. code-block:: yaml

   time:
       start:
           year: 2011
       end:
           year: 2120
       step_size: 365
