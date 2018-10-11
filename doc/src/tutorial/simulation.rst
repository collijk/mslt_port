.. _simulations:

Defining simulations
====================

A simulation is defined in terms of the **components** that are included in
the model, and the **configuration settings** that define key model
parameters.

Components
----------

.. _popn_comp:

Population
^^^^^^^^^^

This will generally comprise:

1. The baseline population demographics;
2. An all-cause mortality rate (ACMR); and
3. A years lost due to disability (YLD) rate.

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

.. _acute_comp:

Acute diseases and events
^^^^^^^^^^^^^^^^^^^^^^^^^

Acute diseases (such as lower respiratory tract infections) and acute events
(such as road traffic accidents) are characterised in terms of:

- Excess mortality rate; and
- Disability rate.

.. note:: These rates must be specified for each cohort (i.e., age and sex),
   and can also be specified separately for each year of the simulation.

.. _interv_comp:

Interventions
^^^^^^^^^^^^^

Interventions are characterised as one or more Population Impact Fractions
(PIFs), each of which modifies an existing rate (e.g., the incidence rate of a
specific chronic disease).

.. _obs_comp:

Observers
^^^^^^^^^

Observers are used to record salient details during a simulation.
The MSLT_ framework currently provides one observer,
:class:`~mslt_port.observer.AdjustedPYandLE`.

.. autoclass:: mslt_port.observer.AdjustedPYandLE

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
