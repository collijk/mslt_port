#!/usr/bin/env python
"""
Use this script to inspect profiling data for a Vivarium simulation.
"""

from pstats import Stats
import sys


def main(args=None):
    prof_file = 'tobacco_decr_20yrs.stats'
    stats = Stats(prof_file)
    stats.sort_stats('cumulative', 'calls')
    stats.print_stats()
    return 0


if __name__ == "__main__":
    sys.exit(main())
