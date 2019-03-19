#! /usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import logging
import numpy
import os

from collections import defaultdict

from downward.experiment import FastDownwardExperiment
from downward.reports import PlanningReport
from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

from lab.reports import Attribute, geometric_mean

from lab import tools

from common_setup import IssueExperiment

from domain_comparison import DomainComparisonReport
from oracle import OracleReport

exp = FastDownwardExperiment()

REVISIONS = [
    'fts-search-base-v2',
    'ab305ba7fa1f',
]

def remove_revision(run):
    algo = run['algorithm']
    for rev in REVISIONS:
        algo = algo.replace('{}-'.format(rev), '')
    run['algorithm'] = algo
    return run

exp.add_fetcher('data/2019-02-18-astar-baseline-eval',filter=[remove_revision],filter_algorithm=[
    'astar-hmax',
    'astar-masdfpbisim50k',
    'astar-masmiasmbisim50k',
])
exp.add_fetcher('data/2019-02-26-astar-hmax-eval',filter=[remove_revision],merge=True)
exp.add_fetcher('data/2019-02-26-astar-masdfp-eval',filter=[remove_revision],merge=True)
exp.add_fetcher('data/2019-02-26-astar-masmiasm-eval',filter=[remove_revision],merge=True)

ms_algorithm_time = Attribute('ms_algorithm_time', absolute=False, min_wins=True, functions=[geometric_mean])
ms_atomic_algorithm_time = Attribute('ms_atomic_algorithm_time', absolute=False, min_wins=True, functions=[geometric_mean])
ms_memory_delta = Attribute('ms_memory_delta', absolute=False, min_wins=True)
fts_transformation_time = Attribute('fts_transformation_time', absolute=False, min_wins=True, functions=[geometric_mean])
fts_search_task_construction_time = Attribute('fts_search_task_construction_time', absolute=False, min_wins=True, functions=[geometric_mean])
fts_plan_reconstruction_time = Attribute('fts_plan_reconstruction_time', absolute=False, min_wins=True, functions=[geometric_mean])
atomic_task_constructed = Attribute('atomic_task_constructed', absolute=True, min_wins=False)
extra_attributes = [
    ms_algorithm_time,
    ms_atomic_algorithm_time,
    ms_memory_delta,
    fts_transformation_time,
    fts_search_task_construction_time,
    fts_plan_reconstruction_time,
    atomic_task_constructed,
]

attributes = list(IssueExperiment.DEFAULT_TABLE_ATTRIBUTES)
attributes.extend(extra_attributes)

## Configs

all_configs = [
'astar-hmax',
'astar-hmax-transform-atomic',
'astar-hmax-transform-atomic-labelreduction',
'astar-hmax-transform-atomic-bisim-labelreduction',
'astar-hmax-transform-full-bisim-labelreduction-dfp100-t900',
'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900',
'astar-hmax-transform-full-bisim-labelreduction-dfp10000-t900',
'astar-hmax-transform-full-bisim-labelreduction-miasm100-t900',
'astar-hmax-transform-full-bisim-labelreduction-miasm1000-t900',
'astar-hmax-transform-full-bisim-labelreduction-miasm10000-t900',

'astar-masdfpbisim50k',
'astar-masdfpbisim50k-transform-atomic',
'astar-masdfpbisim50k-transform-atomic-labelreduction',
'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp100-t900',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp10000-t900',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-miasm100-t900',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-miasm1000-t900',
'astar-masdfpbisim50k-transform-full-bisim-labelreduction-miasm10000-t900',

'astar-masmiasmbisim50k',
'astar-masmiasmbisim50k-transform-atomic',
'astar-masmiasmbisim50k-transform-atomic-labelreduction',
'astar-masmiasmbisim50k-transform-atomic-bisim-labelreduction',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-dfp100-t900',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-dfp1000-t900',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-dfp10000-t900',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-miasm100-t900',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-miasm1000-t900',
'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-miasm10000-t900',
]

## HTML reports

exp.add_report(AbsoluteReport(attributes=['coverage'],filter_algorithm=all_configs))

## Latex reports

algo_to_print = {
    'astar-hmax': '\\SAS',
    'astar-hmax-transform-atomic': '\\atomic',
    'astar-hmax-transform-atomic-bisim-labelreduction': '\\atomicshrink',
    'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900': '\\fulldfp',
    'astar-hmax-transform-full-bisim-labelreduction-miasm1000-t900': '\\fullmiasm',
    'astar-masdfpbisim50k': '\\SAS',
    'astar-masdfpbisim50k-transform-atomic': '\\atomic',
    'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction': '\\atomicshrink',
    'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900': '\\fulldfp',
    'astar-masdfpbisim50k-transform-full-bisim-labelreduction-miasm1000-t900': '\\fullmiasm',
    'astar-masmiasmbisim50k': '\\SAS',
    'astar-masmiasmbisim50k-transform-atomic': '\\atomic',
    'astar-masmiasmbisim50k-transform-atomic-bisim-labelreduction': '\\atomicshrink',
    'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-dfp1000-t900': '\\fulldfp',
    'astar-masmiasmbisim50k-transform-full-bisim-labelreduction-miasm1000-t900': '\\fullmiasm',
}

exp.add_report(
    DomainComparisonReport(
        filter_algorithm=[
            'astar-hmax',
            'astar-hmax-transform-atomic',
            'astar-hmax-transform-atomic-bisim-labelreduction',
            'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900',
            'astar-hmax-transform-full-bisim-labelreduction-miasm1000-t900',
        ],
        algo_to_print=algo_to_print,
        format='tex',
        # attributes=attributes,
        attributes=['coverage'],
    ),
    outfile=os.path.join(exp.eval_dir, 'domain-comparison-coverage-astar-hmax.tex'),
)

exp.add_report(
    DomainComparisonReport(
        filter_algorithm=[
            'astar-masdfpbisim50k',
            'astar-masdfpbisim50k-transform-atomic',
            'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction',
            'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900',
            'astar-masdfpbisim50k-transform-full-bisim-labelreduction-miasm1000-t900',
        ],
        algo_to_print=algo_to_print,
        format='tex',
        # attributes=attributes,
        attributes=['coverage'],
    ),
    outfile=os.path.join(exp.eval_dir, 'domain-comparison-coverage-astar-masdfpbisim50k.tex'),
)


# plots

comparison_algo_pairs = [
    ('astar-hmax', 'astar-hmax-transform-atomic'),
    ('astar-hmax', 'astar-hmax-transform-atomic-bisim-labelreduction'),
    ('astar-hmax', 'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900'),
    ('astar-hmax-transform-atomic', 'astar-hmax-transform-atomic-bisim-labelreduction'),
    ('astar-hmax-transform-atomic', 'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900'),

    ('astar-masdfpbisim50k', 'astar-masdfpbisim50k-transform-atomic'),
    ('astar-masdfpbisim50k', 'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction'),
    ('astar-masdfpbisim50k', 'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900'),
    ('astar-masdfpbisim50k-transform-atomic', 'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction'),
    ('astar-masdfpbisim50k-transform-atomic', 'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900'),
]

tex_comparison_algo_pairs = [
    ('astar-hmax-transform-atomic', 'astar-hmax-transform-atomic-bisim-labelreduction'),
    ('astar-hmax-transform-atomic-bisim-labelreduction', 'astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900'),
    ('astar-masdfpbisim50k-transform-atomic', 'astar-masdfpbisim50k-transform-atomic-bisim-labelreduction'),
    ('astar-masdfpbisim50k-transform-atomic-bisim-labelreduction', 'astar-masdfpbisim50k-transform-full-bisim-labelreduction-dfp1000-t900'),
]

comparison_attributes = [
    'expansions_until_last_jump',
    'search_time',
    'total_time',
]

step_name = "make-absolute-scatter-plots"
scatter_dir = os.path.join(exp.eval_dir, "scatter-plots")
def make_scatter_plot(algo1, algo2, attribute):
    name = "-".join([attribute, algo1, 'vs', algo2])
    print "Make scatter plots for", name

    report = ScatterPlotReport(
        filter_algorithm=[algo1, algo2],
        attributes=[attribute],
        format='tex',
    )
    report(
        exp.eval_dir,
        os.path.join(scatter_dir, name))

def make_scatter_plots():
    for algo_pair in tex_comparison_algo_pairs:
        for attribute in comparison_attributes:
            make_scatter_plot(algo_pair[0], algo_pair[1], attribute)

exp.add_step(step_name, make_scatter_plots)

exp.run_steps()
