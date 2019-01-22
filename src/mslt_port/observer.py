import pandas as pd

class TobaccoPrevalence:

    def __init__(self, output_file):
        self.output_file = output_file
        self.name = 'tobacco'

    def setup(self, builder):
        self.config = builder.configuration
        self.clock = builder.time.clock()
        self.bin_years = int(self.config[self.name]['delay'])

        view_columns = ['age', 'sex', 'bau_population', 'population'] + self.get_bin_names()
        self.population_view = builder.population.get_view(view_columns)

        self.tables = []
        self.table_cols = ['age', 'sex', 'year',
                           'bau_no', 'bau_yes', 'bau_previously', 'bau_population',
                           'int_no', 'int_yes', 'int_previously', 'int_population']

        builder.event.register_listener('collect_metrics',
                                        self.on_collect_metrics)
        builder.event.register_listener('simulation_end',
                                        self.write_output)

    def get_bin_names(self):
        """
        Return the bin names for both the BAU and the intervention scenario.

        These names take the following forms:

        - ``"name.no"``: The number of people who have never been exposed.
        - ``"name.yes"``: The number of people currently exposed.
        - ``"name.N"``: The number of people N years post-exposure.

          - The final bin is the number of people :math:`\ge N` years
            post-exposure.

        The intervention bin names take the form ``"name_intervention.X"``.
        """
        if self.bin_years == 0:
            delay_bins = [str(0)]
        else:
            delay_bins = [str(s) for s in range(self.bin_years + 2)]
        bins = ['no', 'yes'] + delay_bins
        bau_bins = ['{}.{}'.format(self.name, bin) for bin in bins]
        int_bins = ['{}_intervention.{}'.format(self.name, bin) for bin in bins]
        all_bins = bau_bins + int_bins
        return all_bins

    def on_collect_metrics(self, event):
        pop = self.population_view.get(event.index)
        if len(pop.index) == 0:
            # No tracked population remains.
            return

        bau_cols = [c for c in pop.columns.values
                    if c.startswith('{}.'.format(self.name))]
        int_cols = [c for c in pop.columns.values
                    if c.startswith('{}_intervention.'.format(self.name))]

        bau_denom = pop.reindex(columns=bau_cols).sum(axis=1)
        int_denom = pop.reindex(columns=int_cols).sum(axis=1)

        # print('{} {} {}'.format(self.clock().year, min(bau_denom), min(int_denom)))

        # pop['bau_no'] = pop['{}.no'.format(self.name)] / bau_denom
        # pop['bau_yes'] = pop['{}.yes'.format(self.name)] / bau_denom
        # pop['bau_previously'] = 1 - pop['bau_no'] - pop['bau_yes']
        # pop['int_no'] = pop['{}_intervention.no'.format(self.name)] / int_denom
        # pop['int_yes'] = pop['{}_intervention.yes'.format(self.name)] / int_denom
        # pop['int_previously'] = 1 - pop['int_no'] - pop['int_yes']

        pop['bau_no'] = pop['{}.no'.format(self.name)]
        pop['bau_yes'] = pop['{}.yes'.format(self.name)]
        pop['bau_previously'] = pop.reindex(columns=bau_cols).sum(axis=1) - pop['bau_no'] - pop['bau_yes']
        pop['int_no'] = pop['{}_intervention.no'.format(self.name)]
        pop['int_yes'] = pop['{}_intervention.yes'.format(self.name)]
        pop['int_previously'] = pop.reindex(columns=int_cols).sum(axis=1) - pop['int_no'] - pop['int_yes']

        # TODO: where the denominator is zero, set these columns to zero
        # bau_zero = bau_denom == 0.0
        # int_zero = int_denom == 0.0
        # pop.loc[bau_zero, ['bau_yes', 'bau_no', 'bau_previously']] = 0.0
        # pop.loc[int_zero, ['int_yes', 'int_no', 'int_previously']] = 0.0

        pop = pop.rename(columns={'population': 'int_population'})

        pop['year'] = self.clock().year
        # print(pop.reindex(columns=self.table_cols).reset_index(drop=True).head().to_string())
        self.tables.append(pop.reindex(columns=self.table_cols).reset_index(drop=True))


    def write_output(self, event):
        data = pd.concat(self.tables, ignore_index=True)
        data['year_of_birth'] = data['year'] - data['age']
        # Sort the table by cohort (i.e., generation and sex), and then by
        # calendar year, so that results are output in the same order as in
        # the spreadsheet models.
        data = data.sort_values(by=['year_of_birth', 'sex', 'age'], axis=0)
        data = data.reset_index(drop=True)
        # Re-order the table columns.
        cols = ['year_of_birth'] + self.table_cols
        data = data.reindex(columns=cols)
        data.to_csv(self.output_file, index=False)
