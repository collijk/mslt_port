.. _equations:

.. math::

   \DeclareMathOperator{\APC}{APC}
   \DeclareMathOperator{\ACMR}{ACMR}
   \DeclareMathOperator{\PD}{PD}
   \DeclareMathOperator{\PY}{PY}
   \DeclareMathOperator{\LE}{LE}
   \DeclareMathOperator{\Pop}{Pop}
   \DeclareMathOperator{\Deaths}{Deaths}
   \DeclareMathOperator{\YLDrate}{YLDrate}
   \DeclareMathOperator{\PIF}{PIF}
   \DeclareMathOperator{\CPIF}{CPIF}
   \DeclareMathOperator{\IPIF}{InterventionPIF}
   \newcommand{\FPIF}{\ensuremath{1 - \PIF_{\mathrm{Final}}}}
   \DeclareMathOperator{\Prev}{Prev}


MSLT Equations
==============

The inputs are the all-cause mortality rate, :math:`ACMR(a, t_0)`, and the
annual percent change in ACMR, :math:`APC(a, t)`.
The outputs are disability-adjusted person-years, :math:`PY_{adj}(a,t)`, and
the disability-adjusted life expectancy, :math:`LE_{adj}(a,t)`.
Interventions are evaluated by comparing the business-as-usual outputs to
their intervention-specific values.

.. math::

   \begin{align}
     \ACMR(a, t+1) &= \ACMR(a, t) \times \left[ 1 + \frac{\APC(a, t)}{100}\right] \\
     \ACMR(a, t_0+n) &= \ACMR(a, t_0) \times
     \prod_{k=0}^{n-1} \left[ 1 + \frac{\APC(a, t_0 + k)}{100} \right] \\
     \PD(a, t) &= P(t < T_{death,a} < t + 1 | T_{death,a} > t) = 1 - e^{-\ACMR(a, t)} \\
     \PD_{cum}(a, t_0, n) &= \prod_{k = 0}^n \left[ 1 - e^{-\ACMR(a+k, t_0+k)} \right] \\
     \Deaths(a, t_0) &= \Pop(a, t_0) \times \PD(a, t_0) \\
     \Deaths_{cum}(a+n, t_0+n) &= \Pop(a, t_0) \times \PD_{cum}(a, t_0, n) \\
     \Pop(a+1, t_0+1) &= \Pop(a, t_0) - \Deaths(a, t_0) \\
     &= \Pop(a, t_0) \times (1 - \PD(a, t_0)) \\
     \Pop(a+n, t_0+n) &= \Pop(a, t_0) \prod_{k=0}^{n-1}\left[1 - \PD(a+k, t_0+k)\right]\\
     \PY(a, t_0) &= \Pop(a, t_0) \times \left(1 - \frac{\PD(a,t_0)}{2}\right) \\
     \PY(a+n, t_0+n) &= \Pop(a, t_0) \prod_{k=0}^{n-1} \left[ 1 - \PD(a+k, t_0+k) \right]
     \times \left(1 - \frac{\PD(a+k+1, t_0+k+1}{2}\right) \\
     \LE(a, t) &= \frac{\sum_{k=0}^{a_{\max}-a} \PY(a+k,t+k)}{\Pop(a, t)} \\
     \YLDrate &: a \mapsto \YLDrate \\
     \PY_{adj}(a, t) &= \PY(a, t) \times \left[1 - \YLDrate(a)\right] \\
     \LE_{adj}(a, t) &= \frac{\sum_{k=0}^{a_{\max}-a} \PY_{adj}(a+k, t+k)}{\Pop(a, t)}
   \end{align}

.. table:: Definition of symbols used in equations.

   =================  =============================================================
   Symbol             Definition
   =================  =============================================================
   :math:`\ACMR`      All-cause mortality rate
   :math:`\APC`       Annual percent change in :math:`\ACMR`
   :math:`\PD`        Probability of death in a cohort over a single year
   :math:`\Pop`       Number of individuals in a cohort
   :math:`\PY`        Person-years in a cohort over a single year
   :math:`\LE`        Life expectancy, relative to current age
   :math:`\YLDrate`   Year-life disability discount rate
   :math:`\PY_{adj}`  Person-years, adjusted for YLD
   :math:`\LE_{adj}`  Life expectancy, relative to current age and adjusted for YLD
   =================  =============================================================
