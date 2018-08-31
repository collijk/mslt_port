
class MortalityShift:

    def setup(self, builder):
        builder.value.register_value_modifier('mortality_rate', self.mortality_adjustment)

    def mortality_adjustment(self, index, rates):
        return rates * .5


class YLDShift:

    def setup(self, builder):
        builder.value.register_value_modifier('yld_rate', self.disability_adjustment)

    def disability_adjustment(self, index, rates):
        return rates * .5


class IncidenceShift:

    def __init__(self, name):
        self.name = name

    def setup(self, builder):
        builder.value.register_value_modifier(f'{self.name}_intervention.incidence', self.incidence_adjustment)

    def incidence_adjustment(self, index, rates):
        return rates * .5
