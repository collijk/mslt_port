import numpy as np
import pandas as pd

from mslt_port.data import (get_incidence, get_remission, get_excess_mortality,
                            get_prevalence, get_disability, AGE_GROUP_END)


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
        self.population_view = builder.population.get_view(columns)

        builder.event.register_listener('time_step__prepare', self.on_time_step_prepare)

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

        return mortality_rate + (cause_mortality_rate_int - cause_mortality_rate)

    def disability_adjustment(self, index, yld_rate):
        pop = self.population_view.get(index)

        C, C_prev = pop[f'{self.name}_C'], pop[f'{self.name}_C_previous']
        C_int, C_int_prev = pop[f'{self.name}_C_intervention'], pop[f'{self.name}_C_intervention_previous']

        prevalence_change = C - C_prev
        prevalence_change_int = C_int - C_int_prev

        return yld_rate + self.disability_rate(index) * (prevalence_change_int - prevalence_change)

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
        return new_S

    @staticmethod
    def update_C(S, C, i, r, f, l, q, w, v):
        new_C = -(v - w)*(2*((f + r)*(S + C) - l*S) - l*C) + (v + w)*q*C / (2 * q)
        new_C[q == 0] = C[q == 0]
        return new_C
