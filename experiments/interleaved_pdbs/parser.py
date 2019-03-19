#! /usr/bin/env python

from lab.parser import Parser

parser = Parser()
parser.add_pattern('ms_algorithm_time', 'Merge-and-shrink algorithm runtime: (.+)s', required=False, type=float)
parser.add_pattern('ms_atomic_algorithm_time', 'Merge-and-shrink atomic construction runtime: (.+)s', required=False, type=float)
parser.add_pattern('ms_memory_delta', 'Final peak memory increase of merge-and-shrink algorithm: (\d+) KB', required=False, type=int)
parser.add_pattern('fts_transformation_time', 'Transform time: (.+)s', required=False, type=float)
parser.add_pattern('transformed_task_variables', 'Main task: FTSTask with (\d+) variables, \d+ labels, \d+ facts, \d+ transitions', required=False, type=int)
parser.add_pattern('transformed_task_labels', 'Main task: FTSTask with \d+ variables, (\d+) labels, \d+ facts, \d+ transitions', required=False, type=float)
parser.add_pattern('transformed_task_facts', 'Main task: FTSTask with \d+ variables, \d+ labels, (\d+) facts, \d+ transitions', required=False, type=float)
parser.add_pattern('transformed_task_transitions', 'Main task: FTSTask with \d+ variables, \d+ labels, \d+ facts, (\d+) transitions', required=False, type=float)
parser.add_pattern('fts_search_task_construction_time', 'Done building search task wrapper for FTS task: (.+)s', required=False, type=float)
parser.add_pattern('search_task_variables', 'Heuristic task: FTSTask with (\d+) variables, \d+ labels, \d+ facts, \d+ transitions', required=False, type=int)
parser.add_pattern('search_task_labels', 'Heuristic task: FTSTask with \d+ variables, (\d+) labels, \d+ facts, \d+ transitions', required=False, type=float)
parser.add_pattern('search_task_facts', 'Heuristic task: FTSTask with \d+ variables, \d+ labels, (\d+) facts, \d+ transitions', required=False, type=float)
parser.add_pattern('search_task_transitions', 'Heuristic task: FTSTask with \d+ variables, \d+ labels, \d+ facts, (\d+) transitions', required=False, type=float)
parser.add_pattern('fts_plan_reconstruction_time', 'Plan reconstruction time: (.+)s', required=False, type=float)

def check_flags(content, props):
    atomic_task_constructed = False
    for line in content.split('\n'):
        if line == 'Main task constructed':
            atomic_task_constructed = True
        if line == 'Task solved without search':
            assert 'expansions_until_last_jump' not in props
            props['search_time'] = 0
            props['expansions'] = 0
            props['expansions_until_last_jump'] = 0
    props['atomic_task_constructed'] = atomic_task_constructed

parser.add_function(check_flags)

parser.parse()
