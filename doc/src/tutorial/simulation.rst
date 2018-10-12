.. _simulations:

Defining simulations
====================

A simulation is defined in terms of the **components** that are included in
the model, and the **configuration settings** that define key model
parameters.

Overview
--------

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
- Excess mortality rate;
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

.. autoclass:: mslt_port.observer.AdjustedPYandLE

Any number of observers can be included in a simulation.
Each observer should write to a different output file.

.. code-block:: yaml

   components:
       mslt_port:
           observer:
               - AdjustedPYandLE('adjusted-py-le.csv')

At the end of a simulation, it will write these data to the specified output
file:

.. literalinclude:: ../../../adjusted-py-le.csv
   :language: text
   :lines: -6
   :append: ...

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
