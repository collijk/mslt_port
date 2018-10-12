.. _inputs:

Input data
==========

Population
----------

Population data are provided as a CSV file that identifies each cohort by sex
and age (at the simulation start time).
The mortality rate and disability rate must be defined for each cohort, while
the population size is defined for each 5-year cohort (0--4, 5--9, etc).
For example, in the :ref:`example data below <input_popn>` there are 114,970
males aged 0--4.

.. note:: The population sizes should be provided for the **central** cohort
   (i.e., age 2 for the 0--4 cohorts, age 7 for the 5--9 cohorts, etc).

.. literalinclude:: ../../../data/inputs_take2.csv
   :caption: An example of population input data.
   :name: input_popn
   :language: text
   :lines: -11
   :append: ...

.. note:: The mortality and disability rates are indexed by **current age**
   and sex for the duration of each simulation.
   The mortality and disability rates for, e.g., 3-year-old males will be
   applied in each year of the simulation to males that are 3 years old at
   that time-step.

Chronic diseases
----------------

Chronic disease data are provided as a CSV file that defines the following
quantities for each population cohort:

- The incidence rate :math:`i`;
- The remission rate :math:`r`;
- The mortality rate :math:`f`;
- The disability rate :math:`DR`; and
- The initial prevalence :math:`prev`.

Data for multiple diseases can be included in a single CSV file, because the
column names must include the name of the disease (see the example data,
below).

.. literalinclude:: ../../../data/inputs_chd_only.csv
   :caption: An example of chronic disease input data (not showing rows for
      males aged less than 18 years).
   :name: input_chronic
   :language: text
   :lines: 1,20-25
   :append: ...

.. note:: The mortality and disability rates are indexed by **current age**
   and sex for the duration of each simulation.

Acute diseases and events
-------------------------

Acute diseases and events data are provided as a CSV file that defines an
excess mortality rate and a disability rate for each cohort.
Shown below are example data for lower respiratory tract infections (LRTIs):

.. literalinclude:: ../../../data/inputs_lrti.csv
   :caption: An example of acute disease input data.
   :name: input_acute
   :language: text
   :lines: -11
   :append: ...

.. note:: The mortality and disability rates are indexed by **current age**
   and sex for the duration of each simulation.

Interventions
-------------

.. note:: At this stage, we expect intervention input data to take the form of
   ``(sex, age, PIF)`` tables.

Data processing pipeline
------------------------

We anticipate building a data pipeline that will automate some, if not all, of
the steps required to transform data obtained from common sources into
suitable formats to use as input for MSLT_ simulations.
