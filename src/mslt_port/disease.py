import numpy as np
import pandas as pd

from mslt_port.data import (get_incidence, get_remission, get_excess_mortality,
                            get_prevalence, get_disability, AGE_GROUP_END,
                            get_acute_excess_mortality, get_acute_disability)


class AcuteDisease:
    """
    An acute disease has a sufficiently short duration, relative to the
    time-step size, that it is not meaningful to talk about prevalence.
    Instead, it simply contributes an excess mortality rate, and/or a
    disability rate.

    Interventions may affect these rates:

    - `<disease>_intervention.excess_mortality`
    - `<disease>_intervention.yld_rate`

    where `<disease>` is the name given to an acute disease object.
    """

    def __init__(self, name):
        self.name = name

    def setup(self, builder):
        mty_rate = builder.lookup.build_table(get_acute_excess_mortality(self.name))
        yld_rate = builder.lookup.build_table(get_acute_disability(self.name))
        self.excess_mortality = builder.value.register_rate_producer(
            f'{self.name}.excess_mortality',
            source=mty_rate)
        self.int_excess_mortality = builder.value.register_rate_producer(
            f'{self.name}_intervention.excess_mortality',
            source=mty_rate)
        self.disability_rate = builder.value.register_rate_producer(
            f'{self.name}.yld_rate',
            source=yld_rate)
        self.int_disability_rate = builder.value.register_rate_producer(
            f'{self.name}_intervention.yld_rate',
            source=yld_rate)
        builder.value.register_value_modifier('mortality_rate', self.mortality_adjustment)
        builder.value.register_value_modifier('yld_rate', self.disability_adjustment)

    def mortality_adjustment(self, index, mortality_rate):
        delta = self.int_excess_mortality(index) - self.excess_mortality(index)
        return mortality_rate + delta

    def disability_adjustment(self, index, yld_rate):
        delta = self.int_disability_rate(index) - self.disability_rate(index)
        return yld_rate + delta


class Disease:

    def __init__(self, name):
        self.name = name

    def setup(self, builder):
        i = builder.lookup.build_table(get_incidence(self.name))
        r = builder.lookup.build_table(get_remission(self.name))
        f = builder.lookup.build_table(get_excess_mortality(self.name))
        yld_rate = builder.lookup.build_table(get_disability(self.name))

        self.initial_prevalence = builder.lookup.build_table(get_prevalence(self.name))

        self.incidence = builder.value.register_rate_producer(f'{self.name}.incidence', source=i)
        self.incidence_intervention = builder.value.register_rate_producer(f'{self.name}_intervention.incidence',
                                                                           source=i)
        self.remission = builder.value.register_rate_producer(f'{self.name}.remission', source=r)
        self.excess_mortality = builder.value.register_rate_producer(f'{self.name}.excess_mortality', source=f)
        self.disability_rate = builder.value.register_rate_producer(f'{self.name}.yld_rate', source=yld_rate)

        builder.value.register_value_modifier('mortality_rate', self.mortality_adjustment)
        builder.value.register_value_modifier('yld_rate', self.disability_adjustment)

        columns = [f'{self.name}_S', f'{self.name}_S_previous',
                   f'{self.name}_C', f'{self.name}_C_previous',
                   f'{self.name}_S_intervention', f'{self.name}_S_intervention_previous',
                   f'{self.name}_C_intervention', f'{self.name}_C_intervention_previous']
        builder.population.initializes_simulants(self.on_initialize_simulants,
                                                 creates_columns=columns,
                                                 requires_columns=['age', 'sex'])
        self.population_view = builder.population.get_view(columns + ['age', 'sex', 'population'])

        builder.event.register_listener('time_step__prepare', self.on_time_step_prepare)
        self.clock = builder.time.clock()
        builder.event.register_listener('collect_metrics', self.after_step)

    def on_initialize_simulants(self, pop_data):
        C = 1000 * self.initial_prevalence(pop_data.index)
        S = 1000 - C

        pop = pd.DataFrame({f'{self.name}_S': S,
                            f'{self.name}_C': C,
                            f'{self.name}_S_previous': S,
                            f'{self.name}_C_previous': C,
                            f'{self.name}_S_intervention': S,
                            f'{self.name}_C_intervention': C,
                            f'{self.name}_S_intervention_previous': S,
                            f'{self.name}_C_intervention_previous': C},
                           index=pop_data.index)

        self.population_view.update(pop)
        print('\nFollowing males aged 2 in 2010 ...\n')
        print('{:>3}       {:>15}  {:>15}  {:>15}  {:>9}'.format(
            'age', 'Healthy', 'Diseased', 'Dead', 'Mort Risk'))

    def on_time_step_prepare(self, event):
        idx = event.index
        pop = self.population_view.get(idx)
        S, C = pop[f'{self.name}_S'], pop[f'{self.name}_C']
        S_int, C_int = pop[f'{self.name}_S_intervention'], pop[f'{self.name}_C_intervention']

        new_S = self.update_S(
            S, C, self.incidence(idx), self.remission(idx), self.excess_mortality(idx),
            self.l(idx), self.q(idx), self.w(idx), self.v(idx)
        )
        new_C = self.update_C(
            S, C, self.incidence(idx), self.remission(idx), self.excess_mortality(idx),
            self.l(idx), self.q(idx), self.w(idx), self.v(idx)
        )

        new_S_intervention = self.update_S(
            S_int, C_int, self.incidence_intervention(idx), self.remission(idx), self.excess_mortality(idx),
            self.l_intervention(idx), self.q_intervention(idx), self.w_intervention(idx), self.v_intervention(idx)
        )
        new_C_intervention = self.update_C(
            S_int, C_int, self.incidence_intervention(idx), self.remission(idx), self.excess_mortality(idx),
            self.l_intervention(idx), self.q_intervention(idx), self.w_intervention(idx), self.v_intervention(idx)
        )
        pop_update = pd.DataFrame({f'{self.name}_S': new_S,
                                   f'{self.name}_C': new_C,
                                   f'{self.name}_S_previous': S,
                                   f'{self.name}_C_previous': C,
                                   f'{self.name}_S_intervention': new_S_intervention,
                                   f'{self.name}_C_intervention': new_C_intervention,
                                   f'{self.name}_S_intervention_previous': S_int,
                                   f'{self.name}_C_intervention_previous': C_int},
                                  index=pop.index)
        self.population_view.update(pop_update)
        # The index for males aged 2 years in 2010.
        cix = 5
        if pop.shape[0] <= cix:
            return
        age = self.clock().year - 2008 - 1
        mortality_risk = ((- new_S[cix] - new_C[cix] + S[cix] + C[cix])
                          / (S[cix] + C[cix]))
        mortality_risk_int = (- new_S_intervention[cix] - new_C_intervention[cix]
                              + S_int[cix] + C_int[cix]) / (S_int[cix] + C_int[cix])
        print('{:3d}  BAU  {: 15.9f}  {: 15.9f}  {: 15.9f}  {:0.4e}'.format(
            age,
            new_S[cix], new_C[cix], 1000 - new_S[cix] - new_C[cix],
            mortality_risk))
        print('     INT  {: 15.9f}  {: 15.9f}  {: 15.9f}  {:0.4e}'.format(
            new_S_intervention[cix], new_C_intervention[cix],
            1000 - new_S_intervention[cix] - new_C_intervention[cix],
            mortality_risk_int))
        # NOTE: See 'CHD' worksheet, column AG: delta mortality
        #       =(-LN(1-AD22))-(-LN(1-R22))
        #       AD22 is the intervention mortality risk
        #       R22 is the BAU mortality risk

    def after_step(self, event):
        idx = event.index
        pop = self.population_view.get(idx)
        age = self.clock().year - 2008
        # The index for males aged 2 years in 2010.
        cix = 5
        if pop.shape[0] <= cix:
            return
        prev_S = pop.loc[cix, f'{self.name}_S_previous']
        prev_C = pop.loc[cix, f'{self.name}_C_previous']
        prev_S_int = pop.loc[cix, f'{self.name}_S_intervention_previous']
        prev_C_int = pop.loc[cix, f'{self.name}_C_intervention_previous']
        new_S = pop.loc[cix, f'{self.name}_S']
        new_C = pop.loc[cix, f'{self.name}_C']
        new_S_int = pop.loc[cix, f'{self.name}_S_intervention']
        new_C_int = pop.loc[cix, f'{self.name}_C_intervention']
        mortality_risk = ((- new_S - new_C + prev_S + prev_C)
                          / (prev_S + prev_C))
        mortality_risk_int = ((- new_S_int - new_C_int + prev_S_int
                               + prev_C_int)
                              / (prev_S_int + prev_C_int))
        print('{:3d}  BAU  {: 15.9f}  {: 15.9f}  {: 15.9f}  {:0.4e}'.format(
            age,
            new_S, new_C, 1000 - new_S - new_C,
            mortality_risk))
        print('     INT  {: 15.9f}  {: 15.9f}  {: 15.9f}  {:0.4e}'.format(
            new_S_int, new_C_int, 1000 - new_S_int - new_C_int,
            mortality_risk_int))

    def mortality_adjustment(self, index, mortality_rate):
        pop = self.population_view.get(index)

        S, C = pop[f'{self.name}_S'], pop[f'{self.name}_C']
        S_prev, C_prev = pop[f'{self.name}_S_previous'], pop[f'{self.name}_C_previous']
        D, D_prev = 1000 - S - C, 1000 - S_prev - C_prev

        S_int, C_int = pop[f'{self.name}_S_intervention'], pop[f'{self.name}_C_intervention']
        S_int_prev, C_int_prev = pop[f'{self.name}_S_intervention_previous'], pop[f'{self.name}_C_intervention_previous']
        D_int, D_int_prev = 1000 - S_int - C_int, 1000 - S_int_prev - C_int_prev

        mortality_risk = (D - D_prev) / (S + C)
        mortality_risk_int = (D_int - D_int_prev) / (S_int + C_int)

        cause_mortality_rate = -np.log(1 - mortality_risk)
        cause_mortality_rate_int = -np.log(1 - mortality_risk_int)

        delta = cause_mortality_rate_int - cause_mortality_rate
        # The index for males aged 2 years in 2010.
        cix = 5
        if pop.shape[0] > cix:
            bau = mortality_rate[cix] + cause_mortality_rate[cix]
            print('    ACMR  {:15.9f} +{:15.9f} ={:15.9f}'.format(
                bau, delta[cix], bau + delta[cix]))

        return mortality_rate + (cause_mortality_rate_int - cause_mortality_rate)

    def disability_adjustment(self, index, yld_rate):
        pop = self.population_view.get(index)

        S, S_prev = pop[f'{self.name}_S'], pop[f'{self.name}_S_previous']
        C, C_prev = pop[f'{self.name}_C'], pop[f'{self.name}_C_previous']
        S_int, S_int_prev = pop[f'{self.name}_S_intervention'], pop[f'{self.name}_S_intervention_previous']
        C_int, C_int_prev = pop[f'{self.name}_C_intervention'], pop[f'{self.name}_C_intervention_previous']

        # The prevalence rate is the mean number of diseased people over the
        # year, divided by the mean number of alive people over the year.
        # The 0.5 multipliers in the numerator and denominator therefore cancel
        # each other out, and can be removed.
        prevalence_rate = (C + C_prev) / (S + C + S_prev + C_prev)
        prevalence_rate_int = (C_int + C_int_prev) / (S_int + C_int + S_int_prev + C_int_prev)

        return yld_rate + self.disability_rate(index) * (prevalence_rate_int - prevalence_rate)

    def l(self, index):
        i = self.incidence(index)
        r = self.remission(index)
        f = self.excess_mortality(index)

        return i + r + f

    def l_intervention(self, index):
        i = self.incidence_intervention(index)
        r = self.remission(index)
        f = self.excess_mortality(index)

        return i + r + f

    def q(self, index):
        i = self.incidence(index)
        r = self.remission(index)
        f = self.excess_mortality(index)

        return np.sqrt(i**2 + r**2 + f**2 + i*r + f*r - i*f)

    def q_intervention(self, index):
        i = self.incidence_intervention(index)
        r = self.remission(index)
        f = self.excess_mortality(index)

        return np.sqrt(i ** 2 + r ** 2 + f ** 2 + i * r + f * r - i * f)

    def w(self, index):
        l = self.l(index)
        q = self.q(index)

        return np.exp(-(l + q) / 2)

    def w_intervention(self, index):
        l = self.l_intervention(index)
        q = self.q_intervention(index)

        return np.exp(-(l + q) / 2)

    def v(self, index):
        l = self.l(index)
        q = self.q(index)

        return np.exp(-(l - q) / 2)

    def v_intervention(self, index):
        l = self.l_intervention(index)
        q = self.q_intervention(index)

        return np.exp(-(l - q) / 2)

    @staticmethod
    def update_S(S, C, i, r, f, l, q, w, v):
        new_S = (2*(v - w)*(S*(f + r) + C*r) + S*(v*(q - l) + w*(q + l))) / (2 * q)
        new_S[q == 0] = S[q == 0]
        # NOTE: try using simpler no-remission equations.
        new_S = S * np.exp(-i)
        return new_S

    @staticmethod
    def update_C(S, C, i, r, f, l, q, w, v):
        new_C = -((v - w)*(2*((f + r)*(S + C) - l*S) - l*C) - (v + w)*q*C) / (2 * q)
        new_C[q == 0] = C[q == 0]
        # NOTE: try using simpler no-remission equations.
        # new_C = S * (- np.expm1(-i)) + C * np.exp(-f)
        new_C = S * (1 - np.exp(-i)) + C * np.exp(-f)
        return new_C
