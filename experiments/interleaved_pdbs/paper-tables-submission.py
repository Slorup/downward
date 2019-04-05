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

from common_setup import IssueExperiment,SANTI_OPTIMAL_SUITE

from domain_comparison import DomainComparisonReport
#from oracle import OracleReport
from plot_cumulative_solved import PlotCumulativeReport

exp = FastDownwardExperiment()

#D    
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/D', filter_algorithm=[
    'D',
],merge=True)
#R    
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/R', filter_algorithm=[
    'R',
],merge=True)
#RT    
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT', filter_algorithm=[
    'RT',
],merge=True)

#ID    
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/ID-BP', filter_algorithm=[
    'ID-BP',
],merge=True)
##IR    
#exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/IR-BP', filter_algorithm=[
#    'IR-BP',
#],merge=True)
##IO    
#exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/IO-BP', filter_algorithm=[
#    'IO-BP',
#],merge=True)
#IRT-BP  
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/IRT-BP/', filter_algorithm=[
    'IRT-BP',
],merge=True)
#IOT
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/IOT-BP-BeforePrintingRe_Evaluated/', filter_algorithm=[
    'IOT-BP',
],merge=True)


re_eval_nodes = Attribute('re_eval_nodes', absolute=False, min_wins=True, functions=[geometric_mean])
re_insert_nodes = Attribute('re_insert_nodes', absolute=False, min_wins=True, functions=[geometric_mean])

attributes = list(IssueExperiment.DEFAULT_TABLE_ATTRIBUTES)
#attributes.extend(extra_attributes)


all_configs = [
        'D',
#        'R',
        'RT',
        'ID-BP',
#        'IR-BP',
#        'IO-BP',
        'IRT-BP',
        'IOT-BP',
]



## HTML reports

exp.add_report(AbsoluteReport(attributes=['coverage'],filter_algorithm=all_configs))

## Latex reports

algo_to_print = {
    'D': 'G',
    'R': 'R',
    'RT': '900',
    'ID-BP': 'IG',
    'IR-BP': 'IR',
    'IO-BP': 'IO',
    'IRT-BP': 'IR',
    'IOT-BP': 'IO',
    'LmCut': 'LC',
    'Complementary2': 'C2',
    'Scorpion': 'Sc',
    'Cartesian-Online': 'CO',
    'RT-10sec' : '10s',
    'RT-100secs' : '100s',
    'RT-300secs': '300s',
    'RT-600secs' : '600s',
    'RT-1200secs': '1200s',
    'RT-1500secs': '1500s',
}

exp.add_report(
    DomainComparisonReport(
        filter_algorithm=all_configs,
        algo_to_print=algo_to_print,
        filter_domain=SANTI_OPTIMAL_SUITE,
        format='tex',
        attributes=['coverage'],
    ),
    outfile=os.path.join(exp.eval_dir, 'domain-comparison-our-configs.tex'),
)


# plots

comparison_algo_pairs = [
    ('RT', 'IRT'),
    ('IRT', 'IOT'),
]

tex_comparison_algo_pairs = [
    ('RT', 'IRT-BP'),
    ('IRT-BP', 'IOT-BP'),
]

comparison_attributes = [
    'expansions_until_last_jump',
#    'search_time',
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

def make_cumulative_plot(attribute, alg_list_name,alg_list):
    name = "-".join(["cumulative", attribute, alg_list_name])
    print "Make cumulative plots for", name
    
    report = PlotCumulativeReport(
        filter_algorithm=alg_list,
        algo_to_print=algo_to_print,
        attributes=[attribute],
        filter_domain=SANTI_OPTIMAL_SUITE,
        format='tex',
    )
    report(
        exp.eval_dir,
        os.path.join(scatter_dir, name))
    

def make_cumulative_plots():
    make_cumulative_plot("total_time", "gamer", ["D", "RT","IRT-BP","IOT-BP"])

exp.add_step("make-cumulative-plots", make_cumulative_plots)
#########################################################################################################################
#STATE OF THE ART
#LmCut
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/LmCut', filter_algorithm=[
    'LmCut',
],merge=True)
#C2
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/ipc-2018-seq-opt/Complementary2', filter_algorithm=[
    'Complementary2',
],merge=True)
#Sc
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/Scorpion', filter_algorithm=[
    'Scorpion',
],merge=True)
#CO
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/h2-preprocessor/Cartesian-Online', filter_algorithm=[
    'Cartesian-Online',
],merge=True)

soa_configs = [
        'IOT-BP',
        'LmCut',
        'Complementary2',
        'Scorpion',
        'Cartesian-Online',
        ]

exp.add_report(
    DomainComparisonReport(
        filter_algorithm=soa_configs,
        algo_to_print=algo_to_print,
        filter_domain=SANTI_OPTIMAL_SUITE,
        format='tex',
        attributes=['coverage'],
    ),
    outfile=os.path.join(exp.eval_dir, 'domain-comparison-SOA.tex'),
)


# plots


def make_cumulative_plots_SOA():
    make_cumulative_plot("total_time", "SOA", ["IOT-BP","LmCut","Complementary2","Scorpion","Cartesian-Online"])

exp.add_step("make-cumulative-plots_SOA", make_cumulative_plots_SOA)

######################################################################
#Fixed times

exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT-10secs', filter_algorithm=[
'RT-10sec',
],merge=True)
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT-100secss', filter_algorithm=[
    'RT-100secs',
],merge=True)
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT-300secss', filter_algorithm=[
    'RT-300secs',
],merge=True)
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT-600secss', filter_algorithm=[
    'RT-600secs',
],merge=True)
exp.add_fetcher('/mnt/data_server/franco/lab-data/interleaved-search/lemmy-experiments/RT-1200secss', filter_algorithm=[
    'RT-1200secs',
],merge=True)

timing_configs = [
        'RT-10sec',
        'RT-100secs',
        'RT-300secs',
        'RT-600secs',
        'RT',
        'RT-1200secs',
        'IOT-BP',
        ]

exp.add_report(
    DomainComparisonReport(
        filter_algorithm=timing_configs,
        algo_to_print=algo_to_print,
        filter_domain=SANTI_OPTIMAL_SUITE,
        format='tex',
        attributes=['coverage'],
    ),
    outfile=os.path.join(exp.eval_dir, 'domain-comparison-timings.tex'),
    )

def make_cumulative_plots_timings():
    make_cumulative_plot("total_time", "Timings", 
        ["IOT-BP", "RT-10sec", "RT-100secs", "RT-300secs", "RT-600secs", "RT", "RT-1200secs",]
        )
exp.add_step("make-cumulative-plots_timings", make_cumulative_plots_timings)

exp.run_steps()
