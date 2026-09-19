#------------------------------------------------------------------------------
#     Copyright (C) 2026 Alberto Passalacqua
#------------------------------------------------------------------------------
# License
#     This file is part of OpenQBMM.
#
#     OpenQBMM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenQBMM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenQBMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     qbmm.plotting
#
# Description
#     Draw the figures of a validation case with matplotlib, in one style for
#     every case: the curves the case computes in the accent of the OpenQBMM
#     site, the references in its dark and grey neutrals after it, so that a
#     reader tells at a glance which curve is the result, and the figures sit
#     on the site in its own colours (the tokens of its industry.css).
#
#     The labels are written in the mathematics of matplotlib, $d_{43}$,
#     between dollar signs.
#
#------------------------------------------------------------------------------

import os

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt                               # noqa: E402

from . import validation                                      # noqa: E402


# --color-accent, and --color-text, --color-neutral-600, --color-accent-800
# and --color-neutral-400 of the site
computedColor = '#5980a6'

referenceColors = ['#1d1f20', '#7a7a7d', '#2c455d', '#b7b7ba']

lineStyles = {
    'line': '-',
    'dashed': '--',
    'dotted': ':',
    'dashDot': '-.',
    'points': 'none',
    'linesPoints': '-',
}

style = {
    'figure.figsize': (6.4, 4.4),
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.grid': True,
    'grid.color': '#d4d4d7',
    'grid.linewidth': 0.6,
    'axes.edgecolor': '#424244',
    'axes.linewidth': 0.8,
    'legend.frameon': False,
    'legend.fontsize': 9.5,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'lines.markeredgewidth': 1.3,
    'svg.fonttype': 'path',
    'pdf.fonttype': 42,
}


def _range(value, where):
    if value is None:
        return None

    if (
        not isinstance(value, list) or len(value) != 2
     or not all(isinstance(v, (int, float)) for v in value)
    ):
        raise validation.ValidationError(
            '%s: a range is written (min max)' % where
        )

    return value


def draw(figure, base, outputDir, formats):
    """Draw a figure into outputDir, once in every format. Returns the files
    written, and the curves that could not be drawn, with the reason"""
    missing = []

    with plt.rc_context(style):
        fig, ax = plt.subplots()

        nReference = 0

        for curve in figure.curves:
            try:
                points = curve.points(base, gaps=True)
            except FileNotFoundError as e:
                missing.append((curve.name, 'no such file: %s' % e))
                continue

            if not any(y is not None for _, y in points):
                missing.append((curve.name, 'no points'))
                continue

            points = points[::curve.every]

            xs = [x for x, _ in points]
            ys = [float('nan') if y is None else y for _, y in points]

            if curve.computed:
                color = computedColor
            else:
                color = referenceColors[nReference % len(referenceColors)]
                nReference += 1

            color = curve.color if curve.color is not None else color

            marker = curve.marker

            if marker is None and curve.style in ('points', 'linesPoints'):
                marker = 'o' if curve.computed else 's'

            ax.plot(
                xs, ys,
                linestyle=lineStyles[curve.style],
                marker=marker if marker is not None else 'None',
                markerfacecolor='none' if curve.computed else color,
                color=color,
                label=curve.label,
                zorder=3 if curve.computed else 2
            )

        where = figure.name

        xRange = _range(figure.get('xRange'), where + '/xRange')
        yRange = _range(figure.get('yRange'), where + '/yRange')

        if figure.get('xScale'):
            ax.set_xscale(figure.get('xScale'))

        if figure.get('yScale'):
            ax.set_yscale(figure.get('yScale'))

        if xRange:
            ax.set_xlim(*xRange)

        if yRange:
            ax.set_ylim(*yRange)

        ax.set_xlabel(figure.get('xLabel', ''))
        ax.set_ylabel(figure.get('yLabel', ''))

        if figure.get('title'):
            ax.set_title(figure.get('title'))

        legend = figure.get('legend', 'best')

        if legend != 'none' and ax.get_legend_handles_labels()[0]:
            ax.legend(loc=legend)

        fig.tight_layout()

        os.makedirs(outputDir, exist_ok=True)

        written = []

        for fmt in formats:
            path = os.path.join(outputDir, figure.name + '.' + fmt)
            # An SVG is drawn on no ground of its own, so that it sits on
            # the ground of the page that shows it
            fig.savefig(
                path, dpi=150 if fmt == 'png' else None,
                transparent=(fmt == 'svg')
            )
            written.append(path)

        plt.close(fig)

    return written, missing
