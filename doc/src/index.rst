Multi-state Life Table macro-simulations using Vivarium
=======================================================

This repository implements a Multi-state Life Table (MSLT) macro-simulation
framework built on top of `Vivarium`_.

Quick start
-----------

1. Ensure that `Vivarium`_ and the `MSLT`_ framework are both
   :ref:`installed <installation>`.

2. Define a simulation, which comprises the :ref:`population <popn_comp>`,
   :ref:`time period <time_conf>`, :ref:`chronic <chronic_comp>` and
   :ref:`acute <acute_comp>` diseases, and
   :ref:`interventions <interv_comp>`.
   An example simulation definition is shown below.

   .. literalinclude:: ../../test_adjPY.yaml
      :caption: The contents of ``test_adjPY.yaml``.
      :language: yaml
      :lines: 11-21,23-

3. Run the simulation:

   .. code-block:: sh

      simulate run test_adjPY.yaml

4. Examine the simulation outputs.
   In the example simulation shown above, the ``AdjustedPYandLE`` observer was
   used to record adjusted person-years and adjusted life-expectancy for each
   population cohort at each year of the simulation.
   These data are written to the CSV file ``adjusted-py-le.csv``.
   An extract of this data file is shown below.

   .. literalinclude:: ../../adjusted-py-le.csv
      :caption: The contents of ``adjusted-py-le.csv``.
      :language: text
      :lines: -11
      :append: ...

.. toctree::
   :hidden:

   Home <self>

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   tutorial/installation
   tutorial/simulation
   tutorial/inputs

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/equations

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
