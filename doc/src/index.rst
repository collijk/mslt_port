Multi-state Life Table macro-simulations using Vivarium
=======================================================

Multi-state life table (MSLT) modelling is a tool that can be used to predict
the impact of preventative interventions on chronic disease morbidity and
mortality, in terms of standard metrics such as health-adjusted life years
(HALYs).

This documentation describes how to use this `MSLT`_ framework, as introduced
in the paper "Multistate lifetable modelling of preventive interventions:
Concept and Python code".
This macro-simulation framework is built on top of the `Vivarium`_
micro-simulation framework.

How to use this documentation
-----------------------------

First, ensure that `Vivarium`_ and the `MSLT`_ framework are both
:ref:`installed <installation>`.
Then refer to the documentation that best suits your learning style:

- To learn how to reproduce each of the simulations presented in the paper,
  and how to adapt these simulations to your own scenarios, see the
  :ref:`tutorial <tutorial_overview>`, which presents a bottom-up
  introduction to the MSLT framework.

- For a top-down introduction to the MSLT framework, see the :ref:`getting
  started guide <getting_started>`.

- To learn about the underlying equations and the Python classes that
  implement these equations, see the :ref:`reference documentation
  <reference>`.

.. toctree::
   :hidden:

   Home <self>
   installation

.. toctree::
   :maxdepth: 2
   :caption: Tutorial

   tutorial/index
   tutorial/first_simulation

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   getting-started/index
   getting-started/simulation
   getting-started/inputs

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index
   reference/equations
   reference/api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
