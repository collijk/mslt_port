import pandas as pd

from mslt_port.data import AGE_GROUP_END, YEAR_START, YEAR_END

class DiseaseObserver:

    def __init__(self, name, output_path):
        self.name = name
        self.output_path = output_path

    def setup(self, builder):
        columns = [f'{self.name}_S', f'{self.name}_S_previous',
                   f'{self.name}_C', f'{self.name}_C_previous',
                   f'{self.name}_S_intervention', f'{self.name}_S_intervention_previous',
                   f'{self.name}_C_intervention', f'{self.name}_C_intervention_previous']
        self.population_view = builder.population.get_view(columns)

        builder.event.register_listener('collect_metrics', self.on_collect_metrics)
        builder.event.register_listener('simulation_end', self.write_output)

        idx = pd.MultiIndex.from_product([range(AGE_GROUP_END + 1), range(YEAR_START, YEAR_END + 1), ['male', 'female']])
        self.data = pd.DataFrame(columns=['S', 'C', 'S_int', 'C_int'], index=idx)

    def on_collect_metrics(self, event):
        pass

    def write_output(self, event):
        pass

