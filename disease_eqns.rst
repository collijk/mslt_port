Interventions affect the ACMR and the life-year disability rate indirectly, because they directly affect disease prevalence.
The equations for disease incidence and mortality come from Barendregt et al., as circulated in my previous email.
 
This allows you to calculate the number of people alive after reaching age a, and the person-years (mean of numbers alive at a and at a+1).
 
The prevalence rate is then calculated (equivalent to dC/da).
 
The mortality risk is the number of people who died during the year (of being aged a) divided by the number of people alive at the start of that year.

The intervention scenario is identical to the business-as-usual scenario, except that incidence rates are different. Remission and case fatality rates remain the same.
 
To evaluate the effect of the intervention, calculate:
 
+ The difference in prevalence rate; i.e., the intervention value of dC/da minus the baseline value of dC/da.
 
+ The change in mortality rate (not mortality risk). Mortality rate is: - log(1 - mortality risk).
 
+ The change in PLYD rate: change in prevalence multiplied by the original PLYD rate for that disease.
 
+ The change in mortality rate and PLYD rate are fed back into the life table.

There are also acute events, such as injuries or lower respiratory tract infections, which are handled differently. It is possible, if not sensible, to implement them using the same approach and selecting a very large value for the remission rate.
