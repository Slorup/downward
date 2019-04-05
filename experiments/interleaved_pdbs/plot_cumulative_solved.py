# -*- coding: utf-8 -*-

from collections import defaultdict

import os
from lab import tools

from matplotlib import ticker

from downward.reports.plot import MatplotlibPlot, Matplotlib, PgfPlots, PlotReport, MIN_AXIS

from downward.reports.plot import PlotReport, Matplotlib, MatplotlibPlot



class CumulativePgfPlots(PgfPlots):
    @classmethod
    def _format_coord(cls, coord):
        def format_value(v):
            return str(v) if isinstance(v, int) else '%f' % v
        return '(%s, %s)' % (format_value(coord[0]), format_value(coord[1]))

    @classmethod
    def _get_plot(cls, report):
        lines = []
        options = cls._get_axis_options(report)
        lines.append('\\begin{axis}[%s]' % cls._format_options(options))
        for category, coords in sorted(report.categories.items()):
            plot = {'no marks' : True}#{'only marks': True}
            lines.append(
                '\\addplot+[%s] coordinates {\n%s\n};' % (
                    cls._format_options(plot),
                    ' '.join(cls._format_coord(c) for c in coords)))
            if category:
                lines.append('\\addlegendentry{%s}' % category)
            elif report.has_multiple_categories:
                # None is treated as the default category if using multiple
                # categories. Add a corresponding entry to the legend.
                lines.append('\\addlegendentry{default}')
        lines.append('\\end{axis}')
        return lines

    @classmethod
    def _get_axis_options(cls, report):
        return {'xmode' : 'log', 'legend pos':'outer north east', 'xmin' : 1,  'ymin' : '0', 'ymax' : str(report.max_coverage), 'cycle list name' : 'color list'}



class PlotCumulativeReport(PlotReport):
    """
    Generate a cumulative coverage plot for a specific attribute.
    """
    def __init__(self, show_missing=True, get_category=None, algo_to_print={}, **kwargs):
        """
        See :class:`.PlotReport` for inherited arguments.

        Use the *filter_algorithm* keyword argument to select algorithms.

        If only one of the two algorithms has a value for a run, only
        add a coordinate if *show_missing* is True.

        """
        self.algo_to_print = algo_to_print
        # If the size has not been set explicitly, make it a square.
        matplotlib_options = kwargs.get('matplotlib_options', {})
        matplotlib_options.setdefault('figure.figsize', [8, 8])
        kwargs['matplotlib_options'] = matplotlib_options
        PlotReport.__init__(self, **kwargs)
        if not self.attribute:
            logging.critical('CumulativePlotReport needs exactly one attribute')
        # By default all values are in the same category.
        self.get_category = get_category or (lambda run1, run2: None)
        self.show_missing = show_missing
        self.xlim_left = self.xlim_left or MIN_AXIS
        self.ylim_bottom = self.ylim_bottom or MIN_AXIS
        if self.output_format == 'tex':
            self.writer = CumulativePgfPlots
        else:
            self.writer = ScatterMatplotlib

    def _set_scales(self, xscale, yscale):
        PlotReport._set_scales(self, xscale or self.attribute.scale or 'log', yscale)


    def _fill_categories(self, runs):
        # We discard the *runs* parameter.
        # Map category names to value tuples
        categories = defaultdict(list)
        xvalues = defaultdict(set)
        algorithms = set()
        for (domain, problem), runs in self.problem_runs.items():
            for run in runs:
                if run['coverage'] == 1:
                    xvalues[run["algorithm"]].add(int(run[self.attribute]))
                    algorithms.add(run["algorithm"])

        values = {}
        for a in algorithms:
            xvalues[a] = sorted(list(xvalues[a]))
            values [a] = {} 
            for x in xvalues[a]:
                values[a][x] = 0

        for (domain, problem), runs in self.problem_runs.items():
            for run in runs:
                if run['coverage'] == 1:
                    x = run[self.attribute]
                    a = run["algorithm"]

                    for xs in xvalues[a]:
                        if xs >= x:
                            values[a][xs] += 1

        self.max_coverage = 0
        for a in algorithms:
            if a not in self.algo_to_print:
                self.algo_to_print[a] = a
            self.max_coverage = max(self.max_coverage, values[a][xvalues[a][-1]])
            for x in xvalues[a]:
                categories[self.algo_to_print[a]].append((x, values[a][x]))

        return categories

    
    def _prepare_categories(self, categories):
        return PlotReport._prepare_categories(self, categories)

    def write(self):
        self.xlabel = self.xlabel or self.attribute
        self.ylabel = self.ylabel or "Instances Solved"

        suffix = '.' + self.output_format
        if not self.outfile.endswith(suffix):
            self.outfile += suffix
        tools.makedirs(os.path.dirname(self.outfile))
        self._write_plot(self.runs.values(), self.outfile)
