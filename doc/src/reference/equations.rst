.. _equations:

.. math::

   \DeclareMathOperator{\d}{d\!}
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


MSLT equations
==============

Life table
----------

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

Chronic disease
---------------

The equations for chronic disease prevalence, remission, and mortality come
from `Barendregt et al., 2003 <https://doi.org/10.1186/1478-7954-1-4>`_.
A key assumption in their derivation is the independence of mortality from all
causes:

   If it is assumed that mortality from all other causes is independent of the
   disease, i.e., that it is the same for healthy and diseased people, this
   implies that the transition hazards for incidence, remission and case
   fatality are not affected by the value of the 'all other causes' mortality.
   Therefore we can set the value of mortality from all other causes to 0
   (i.e., leave it out of the equations) and still derive the right values for
   the disease rates.

.. graphviz::
   :align: center
   :caption: The conceptual disease model.
       :math:`S`: number of healthy people (i.e., without the disease under
       consideration); :math:`C`: number of diseased people; :math:`D`: number
       of people dead from the disease; and :math:`M`: number of people dead
       from all other causes.
       There are four transition hazards: :math:`i`: incidence, :math:`r`:
       remission, :math:`f`: case fatality, and :math:`m`: all other
       mortality.

    digraph dis_mod {
      {rank=same; Sa Ma }
      {rank=same; Ca Da }
      Sa [label="S", shape=box];
      Ca [label="C", shape=box];
      Ma [label="M", shape=box];
      Da [label="D", shape=box];
      Sa -> Ca [label="i"];
      Ca -> Sa [label="r"];
      Sa -> Ma [label="m"];
      Ca -> Ma [label="m"];
      Ca -> Da [label="f"];
    }

With this simplifying assumption, the system of equations becomes:

.. math::

   \begin{align}
     \frac{\d{}S_a}{\d{}a} &= -i_a S_a + r_a C_a \\
     \frac{\d{}C_a}{\d{}a} &= -(f_a + r_a) C_a + i_a S_a \\
     \frac{\d{}D_a}{\d{}a} &= f_a C_a
   \end{align}

The eigenvalues :math:`\lambda` can then be calculated:

.. math::

   \begin{align}
     M &= \begin{pmatrix}
       -i & r & 0 \\
       i & -f - r & 0 \\
       0 & f & 0 \\
     \end{pmatrix} \\
     \det(M - \lambda I) &= - \lambda \cdot \left[
       \lambda^2 + \lambda \cdot (i + f + r) + if \right ] \\
     \lambda \ne 0 \implies \lambda &=
       \frac{- i - f - r \pm \sqrt{(i + f + r)^2 - 4if}}{2}
   \end{align}

.. note:: This is equivalent to the derivation in
   `Barendregt et al., 2003 <https://doi.org/10.1186/1478-7954-1-4>`_.

The resulting difference equations for :math:`S_a`, :math:`C_a`, and
:math:`D_a` are then:

.. math::

   \begin{align}
   S_a &= \frac{%
     2 (v_a − w_a) \left[ S_{a−1} (f_a + r_a) + C_{a−1} r_a \right]
     + S_{a−1} \left[ v_a (q_a − l_a) + w_a (q_a + l_a) \right]
     }{2q_a} \\
   C_a &= \frac{%
     (v_a - w_a) \left( 2\left[ (f_a + r_a)(S_{a-1} + C_{a-1}) - l_a S_{a-1}
       \right] - C_{a-1} l_a \right) - C_{a_1} q_a (v_a + w_a)
     }{2q_a} \\
   D_a &= \frac{%
     (v_a - w_a) \left[ 2 f_a C_{a-1} - l_a(S_{a-1} + C_{a-1}) \right]
     - q_a (S_{a-1} + C_{a-1})(v_a + w_a) + 2q_a (S_{a-1} + C_{a-1} + D_{a-1})
     }{2q_a}
   \end{align}

where we define the following convenience variables:

.. math::

   \begin{align}
   l_a &= i_a + r_a + f_a \\
   q_a &= \sqrt{(i + f + r)^2 - 4if} \\
   w_a &= \exp\left( \frac{-l_a + q_a}{2} \right) \\
   v_a &= \exp\left( \frac{-l_a - q_a}{2} \right)
   \end{align}
