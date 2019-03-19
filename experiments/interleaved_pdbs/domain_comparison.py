from downward.reports import PlanningReport

from collections import defaultdict

def domain_mapping(domain):
    if domain == 'openstacks':
        return 'openstacks-06'
    elif 'openstacks' in domain:
        return 'openstacks-08-11-14'
    elif 'logistics' in domain:
        return 'logistics'
    else:
        result = domain.replace('-adl', '')
        result = domain.replace('-strips', '')
        result = result.replace('-08', '')
        result = result.replace('-opt08', '')
        result = result.replace('-opt11', '')
        result = result.replace('-opt14', '')
        result = result.replace('-opt18', '')
        result = result.replace('-sat08', '')
        result = result.replace('-sat11', '')
        result = result.replace('-sat14', '')
        result = result.replace('-sat18', '')
        return result

# TODO: this currently uses coverage and sum as the hard-coded attribute and
# aggregation function.
# NOTE: hacked to consider same domains of different years as one domain.
class DomainComparisonReport(PlanningReport):
    def __init__(self, algo_to_print = {}, **kwargs):
        PlanningReport.__init__(self, **kwargs)
        self.algo_to_print = algo_to_print

    def get_text(self):
        domains = set()
        domain_and_algo_to_coverage = defaultdict(int)
        algo_to_coverage = defaultdict(int)
        for (domain, problem), runs in self.problem_runs.items():
            mapped_domain = domain_mapping(domain)
            domains.add(mapped_domain)
            for run in runs:
                assert domain == run['domain']
                domain_and_algo_to_coverage[(mapped_domain, run['algorithm'])] += run['coverage']
                algo_to_coverage[run['algorithm']] += run['coverage']

        def turn_list_into_table_row(line):
            result = ''
            for index, value in enumerate(line):
                result += '{}'.format(value)
                if index == len(line) - 1:
                    result += ' \\\\'
                else:
                    result += ' & '
            return result

        def format_algo(algo):
            if algo in self.algo_to_print:
                return self.algo_to_print[algo]
            return algo

        lines = []

        header_line = ['']
        for algo in self.algorithms:
            header_line.append(format_algo(algo))
        header_line.append('tot')
        lines.append(turn_list_into_table_row(header_line))

        for algo1 in self.algorithms:
            line = []
            line.append(format_algo(algo1))
            for algo2 in self.algorithms:
                num_algo1_better = 0
                num_algo2_better = 0
                for domain in domains:
                    coverage1 = domain_and_algo_to_coverage[(domain, algo1)]
                    coverage2 = domain_and_algo_to_coverage[(domain, algo2)]
                    if coverage1 > coverage2:
                        num_algo1_better += 1
                    elif coverage2 > coverage1:
                        num_algo2_better += 1

                if algo1 == algo2:
                    entry = "--"
                elif num_algo1_better >= num_algo2_better:
                    entry = "\\textbf{{{}}}".format(num_algo1_better)
                else:
                    entry = str(num_algo1_better)
                line.append(entry)
            # aggregated value column
            line.append(algo_to_coverage[algo1])
            lines.append(turn_list_into_table_row(line))
        return '\n'.join(lines)
