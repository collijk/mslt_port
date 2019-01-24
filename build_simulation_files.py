#!/usr/bin/env python

"""
This script generates the simulation definition files for each experiment.
"""

import jinja2
import os
import sys



def main(args=None):
    """
    Construct the data artifacts requires for the tobacco simulations.
    """

    template_file = 'yaml_template.in'
    with open(template_file, 'r') as f:
        template_contents = f.read()

    template = jinja2.Template(template_contents,
                               trim_blocks=True,
                               lstrip_blocks=True)

    out_format = 'exp_{}_bau-{}_delay-{}_{}.yaml'

    populations = ['non-maori', 'maori']
    bau_labels = ['const', 'decr']
    delay_labels = {0: '0yrs', 20: '20yrs'}
    intervention = {'erad': 'TobaccoEradication',
                    'tfg': 'TobaccoFreeGeneration',
                    'tax': None}

    for population in populations:
        for bau_label in bau_labels:
            for delay, delay_label in delay_labels.items():
                for interv_label, interv_class in intervention.items():
                    out_file = out_format.format(population, bau_label,
                                                 delay_label, interv_label)
                    basename = out_file[:-5]
                    template_args = {
                        'basename': basename,
                        'population': population,
                        'constant_prevalence': bau_label == 'const',
                        'delay': delay,
                        'intervention_class': interv_class,
                        'tobacco_tax': interv_label == 'tax',
                    }
                    if interv_class is None:
                        del template_args['intervention_class']
                    out_content = template.render(template_args)
                    with open(out_file, 'w') as f:
                        f.write(out_content)

    return 0


if __name__ == "__main__":
    sys.exit(main())
