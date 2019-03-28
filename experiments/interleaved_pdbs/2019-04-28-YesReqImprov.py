#! /usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import os
import subprocess

from lab.environments import LocalEnvironment
from lab.reports import Attribute, geometric_mean

from downward.reports.compare import ComparativeReport
from domain_comparison import DomainComparisonReport
from downward import machines
from downward.experiment import FastDownwardExperiment
from downward.reports import PlanningReport
from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

from lab.reports import Attribute, geometric_mean

from lab import tools

from common_setup import IssueExperiment

from domain_comparison import DomainComparisonReport

from common_setup import IssueConfig, IssueExperiment, DEFAULT_OPTIMAL_SUITE, is_test_run, get_experiment_name

REVISION = '87632ae9671d'

BENCHMARKS_DIR = "/mnt/data_server/franco/benchmarks/"
REVISION_CACHE = "/mnt/data_server/franco/lab-data/cache/"

REPO = "/mnt/data_server/franco/nonagnosticpdbs_interleaved/"

def main(revisions=None):
    benchmarks_dir=os.path.expanduser('/mnt/data_server/franco/benchmarks')
    suite = DEFAULT_OPTIMAL_SUITE
    environment = machines.OLD_SERVERS

    if is_test_run():
        suite = ['gripper:prob01.pddl', 'depot:p01.pddl', 'mystery:prob07.pddl']
        environment = LocalEnvironment(processes=4)

    configs = {
        IssueConfig('astar-hmax-transform-atomic', ["--search", "astar(hmax)"]),
    }

    exp = IssueExperiment(
        revisions=revisions,
        configs=configs,
        environment=environment,
    )
    exp.add_suite(benchmarks_dir, suite)

    exp.add_parser(exp.EXITCODE_PARSER)
    exp.add_parser(exp.TRANSLATOR_PARSER)
    exp.add_parser(exp.SINGLE_SEARCH_PARSER)
    exp.add_parser(exp.PLANNER_PARSER)
    exp.add_parser('fts-parser.py')

    ms_algorithm_time = Attribute('ms_algorithm_time', absolute=False, min_wins=True, functions=[geometric_mean])
    ms_atomic_algorithm_time = Attribute('ms_atomic_algorithm_time', absolute=False, min_wins=True, functions=[geometric_mean])
    ms_memory_delta = Attribute('ms_memory_delta', absolute=False, min_wins=True)
    fts_transformation_time = Attribute('fts_transformation_time', absolute=False, min_wins=True, functions=[geometric_mean])
    transformed_task_variables = Attribute('transformed_task_variables', absolute=False, min_wins=True, functions=[sum])
    transformed_task_labels = Attribute('transformed_task_labels', absolute=False, min_wins=True, functions=[sum])
    transformed_task_facts = Attribute('transformed_task_facts', absolute=False, min_wins=True, functions=[sum])
    transformed_task_transitions = Attribute('transformed_task_transitions', absolute=False, min_wins=True, functions=[sum])
    fts_search_task_construction_time = Attribute('fts_search_task_construction_time', absolute=False, min_wins=True, functions=[geometric_mean])
    search_task_variables = Attribute('search_task_variables', absolute=False, min_wins=True, functions=[sum])
    search_task_labels = Attribute('search_task_labels', absolute=False, min_wins=True, functions=[sum])
    search_task_facts = Attribute('search_task_facts', absolute=False, min_wins=True, functions=[sum])
    search_task_transitions = Attribute('search_task_transitions', absolute=False, min_wins=True, functions=[sum])
    fts_plan_reconstruction_time = Attribute('fts_plan_reconstruction_time', absolute=False, min_wins=True, functions=[geometric_mean])
    atomic_task_constructed = Attribute('atomic_task_constructed', absolute=True, min_wins=False)
    #extra_attributes = [
    #    ms_algorithm_time,
    #    ms_atomic_algorithm_time,
    #    ms_memory_delta,
    #    fts_transformation_time,
    #    transformed_task_variables,
    #    transformed_task_labels,
    #    transformed_task_facts,
    #    transformed_task_transitions,
    #    fts_search_task_construction_time,
    #    search_task_variables,
    #    search_task_labels,
    #    search_task_facts,
    #    search_task_transitions,
    #    fts_plan_reconstruction_time,
    #    atomic_task_constructed,
    #]

    attributes = list(exp.DEFAULT_TABLE_ATTRIBUTES)
    #attributes.extend(extra_attributes)

    exp.add_step('build', exp.build)
    exp.add_step('start', exp.start_runs)
    exp.add_fetcher(name='fetch')

   # exp.add_absolute_report_step(attributes=attributes, filter_algorithm=[
   #     '{}-astar-hmax-transform-atomic'.format(REVISION),
   #     '{}-astar-hmax-transform-atomic-labelreduction'.format(REVISION),
   #     '{}-astar-hmax-transform-atomic-bisim-labelreduction'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-dfp100-t900'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-dfp1000-t900'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-dfp10000-t900'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-miasm100-t900'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-miasm1000-t900'.format(REVISION),
   #     '{}-astar-hmax-transform-full-bisim-labelreduction-miasm10000-t900'.format(REVISION),
   # ])

    #BASELINE_REV = 'fts-search-base-v2'
    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/GamerLocal-3GB',filter_algorithm=['LocalGamer-3GB'], merge=True)
    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-3GB-AvgH',filter_algorithm=['Interleaved-LocalGamer-3GB-AvgH'], merge=True)
    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-3GB-OpenList-AvgH',filter_algorithm=['Interleaved-LocalGamer-3GB-OpenList-AvgH'], merge=True)

    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-3GB-AvgH-Online', filter_algorithm=[
    #    'Interleaved-LocalGamer-3GB-AvgH-Online',
    #],merge=True)
    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/GamerLocal-RandWalk', filter_algorithm=[
    #    'LocalGamer-RandWalk',
    #],merge=True)
    #exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-RandWalk-AvgH', filter_algorithm=[
    #    'Interleaved-LocalGamer-RandWalk-AvgH',
    #],merge=True)


#CGAMER
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/CGamer',filter_algorithm=['CGamer'], merge=True)
#LG-RW-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/GamerLocal-RandWalk-DynLims-YesReqImprov', filter_algorithm=[
        'LG-RW-DL-YI',
    ],merge=True)
#LG-AD-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/GamerLocal-AvgDist-DynLims-YesReqImprov', filter_algorithm=[
        'LG-AD-DL-YI',
    ],merge=True)

#ILG-AD-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-GamerLocal-AvgDist-DynLims-YesReqImprov', filter_algorithm=[
        'ILG-AD-DL-YI',
    ],merge=True)
#ILG-RW-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-RandWalk-AvgH-DynLims-YesReqImprov', filter_algorithm=[
        'ILG-RW-DL-YI',
    ],merge=True)
#ILG-OL-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-LocalGamer-OpenListNoRandWalk-DynLims-YesReqImprov', filter_algorithm=[
        'ILG-OL-DL-YI',
    ],merge=True)
#ILG-OL+RW-DL-YI
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-GamerLocal-OpenListWithRandomWalk-DynLims-YesReqImprov-NoTieBreak', filter_algorithm=[
        'ILG-OL+RW-DL-YI',
    ],merge=True)

#ILG-OL+RW+TB-DL-YI    
    exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Interleaved-GamerLocal-OpenListWithRandomWalk-DynLims-YesReqImprov-YesTieBreak', filter_algorithm=[
        'ILG-OL+RW+TB-DL-YI',
    ],merge=True)

    outfile1 = os.path.join(exp.eval_dir, 'Gamer-RequireImprov.html')
    
#    exp.add_report(AbsoluteReport(attributes=attributes,
#        format='html',
#        ),
#        outfile=outfile1,
#    )
    exp.add_report(AbsoluteReport(attributes=attributes,
        filter_algorithm=[
         'LG-RW-DL-YI',
         'LG-AD-DL-YI',
         'ILG-AD-DL-YI',
         'ILG-RW-DL-YI',
         'ILG-OL-DL-YI',
         ],
        format='html',
        ),
        outfile=outfile1,
    )

    outfile = os.path.join(exp.eval_dir, get_experiment_name() + '-vs-regular-GamerLocal.html')
    exp.add_report(
        ComparativeReport(
            algorithm_pairs=[
                ('CGamer','Interleaved-LocalGamer-RandWalk-AvgH-DynLims-YesReqImprov'),
                ('CGamer','Interleaved-LocalGamer-RandWalk-AvgH-DynLims-YesReqImprov'),
                ('CGamer','Interleaved-LocalGamer-RandWalk-AvgH-DynLims-NoReqImprov'),
                #('LocalGamer-3GB', 'Interleaved-LocalGamer-3GB-AvgH'),
                #('LocalGamer-3GB', 'Interleaved-LocalGamer-3GB-OpenList-AvgH'),
                #('LocalGamer-3GB-AvgH', 'Interleaved-LocalGamer-3GB-OpenList-AvgH'),
                #('LocalGamer-3GB', 'LocalGamer-RandWalk'),
                #('LocalGamer-RandWalk','Interleaved-LocalGamer-RandWalk-AvgH'),
                #('Interleaved-LocalGamer-3GB-OpenList-AvgH','Interleaved-LocalGamer-RandWalk-AvgH'),
            ],
            format='html',
            attributes=attributes,
        ),
        outfile=outfile,
    )

    exp.add_report(
        DomainComparisonReport(
        filter_algorithm=[
         'LG-AD-DL-YI',
         'LG-RW-DL-YI',
         'ILG-AD-DL-YI',
         'ILG-RW-DL-YI',
         'ILG-OL-DL-YI',
         'ILG-OL+RW+TB-DL-YI',
         ],
        format='tex',
        attributes=['coverage'],
        ),
        outfile=os.path.join(exp.eval_dir, 'domain-comparison-coverage-Gamers.tex'),
    )

    tex_comparison_algo_pairs = [
        ('LG-AD-DL-YI', 'ILG-AD-DL-YI'),
        ('LG-AD-DL-YI', 'ILG-OL-DL-YI'),
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
            format='png',
        )
        report(
            exp.eval_dir,
            os.path.join(scatter_dir, name))

    def make_scatter_plots():
        for algo_pair in tex_comparison_algo_pairs:
            for attribute in comparison_attributes:
                make_scatter_plot(algo_pair[0], algo_pair[1], attribute)

    exp.add_step(step_name, make_scatter_plots)


    exp.add_step('publish-{}'.format(outfile), subprocess.call, ['publish', outfile])

    exp.run_steps()

main(revisions=[REVISION])
