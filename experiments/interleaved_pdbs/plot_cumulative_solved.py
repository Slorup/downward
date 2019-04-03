# -*- coding: utf-8 -*-

from collections import defaultdict

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
            plot = {'only marks': True}
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
        # Add black line.
        start = min(report.min_x, report.min_y)
        if report.xlim_left is not None:
            start = min(start, report.xlim_left)
        if report.ylim_bottom is not None:
            start = min(start, report.ylim_bottom)
        end = max(report.max_x, report.max_y)
        if report.xlim_right:
            end = max(end, report.xlim_right)
        if report.ylim_top:
            end = max(end, report.ylim_top)
        if report.show_missing:
            end = max(end, report.missing_val)
        lines.append(
            '\\addplot[color=black] coordinates {(%f, %f) (%d, %d)};' %
            (start, start, end, end))
        lines.append('\\end{axis}')
        return lines

    @classmethod
    def _get_axis_options(cls, report):
        opts = PgfPlots._get_axis_options(report)
        # Add line for missing values.
        for axis in ['x', 'y']:
            opts['extra %s ticks' % axis] = report.missing_val
            opts['extra %s tick style' % axis] = 'grid=major'
        return opts



class PlotCumulativeReport(PlotReport):
    """
    Generate a cumulative coverage plot for a specific attribute.
    """
    def __init__(self, show_missing=True, get_category=None, **kwargs):
        """
        See :class:`.PlotReport` for inherited arguments.

        Use the *filter_algorithm* keyword argument to select algorithms.

        If only one of the two algorithms has a value for a run, only
        add a coordinate if *show_missing* is True.

        """
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
            self.writer = ScatterPgfPlots
        else:
            self.writer = ScatterMatplotlib

    def _set_scales(self, xscale, yscale):
        PlotReport._set_scales(self, xscale or self.attribute.scale or 'log', yscale)


    def _fill_categories(self, runs):
        # We discard the *runs* parameter.
        # Map category names to value tuples
        categories = defaultdict(list)
        xvalues = set()
        algorithms = set()
        for (domain, problem), runs in self.problem_runs.items():
            for run in runs:
                if run['coverage'] == 1:
                    xvalues.add(run[self.attribute])
                    algorithms.add(run[self.algorithm])

        
        values = {}
        for a in algorithms:
            values [a] = {} 
            for x in xvalues:
                values[a][x] = 0

        for (domain, problem), runs in self.problem_runs.items():
            for run in runs:
                if run['coverage'] == 1:
                    x = run[self.attribute]
                    a = run[self.algorithm]

                    for xs in xvalues:
                        if xs <= x:
                            values[a][xs] += 1

        for a in algorithms:
            for x in xvalues:
                categories[a].append((x, values[a][x]))
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
