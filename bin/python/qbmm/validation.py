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
#     qbmm.validation
#
# Description
#     The figures a validation case draws and the comparisons it makes, read
#     from <case>/validation/validationDict.
#
#     A case declares each curve once, in the figure it is drawn in: the
#     file it is read from and the expressions of its abscissa and value. A
#     curve the case computes names, in its compare sub-dictionary, the
#     reference curves of the same figure it is held against and the
#     tolerance of each. The graph and the number that judges it are made
#     from the same declaration, so they cannot come to disagree.
#
#         figures
#         {
#             d43
#             {
#                 xLabel  "Time, $t$ [s]";
#                 yLabel  "$d_{43}$ [$\mu$m]";
#                 xRange  (1 200);
#
#                 curves
#                 {
#                     Vanni2000
#                     {
#                         file    "validationData/Vanni2000_d43.dat";
#                         label   "Rigorous solution (Vanni, 2000)";
#                     }
#
#                     OpenQBMM
#                     {
#                         file    ("../postProcessing/probes/0/moment.4.pb"
#                                  "../postProcessing/probes/0/moment.3.pb");
#                         y       "c2/c4";
#                         label   "QMOM (OpenQBMM)";
#                         style   points;
#
#                         compare
#                         {
#                             Vanni2000   0.02;
#                         }
#                     }
#                 }
#             }
#         }
#
#     Files hold columns of numbers, and are named relative to the directory
#     of validationDict. A list of files is read side by side, as paste does,
#     with the columns numbered across all of them. A file may be a pattern:
#     where it matches several files the last one is read, in the order of
#     the numbers in their paths, which is the latest time of a case that
#     sampled several.
#
#     The abscissa and the value of a curve are expressions over the columns,
#     written c1, c2, ... in the notation of Python, with the functions of
#     its math module (sqrt, exp, log, ...). They default to c1 and c2. A
#     value of None leaves the point out, and a gap in the line, as in
#
#         y   "sqrt((c7-c10)**2+(c8-c11)**2) if c2 > 1e-12 else None";
#
#     How a comparison is measured, and what its tolerance means, is
#     described in bin/qbmmCompareReference.
#
#------------------------------------------------------------------------------

import bisect
import glob
import math
import os
import re

from . import foamDict


dictName = 'validationDict'

figureKeys = {
    'title', 'xLabel', 'yLabel', 'xRange', 'yRange', 'xScale', 'yScale',
    'legend', 'curves'
}

curveKeys = {
    'file', 'x', 'y', 'label', 'style', 'color', 'marker', 'every',
    'computed', 'compare'
}

styles = {'line', 'dashed', 'dotted', 'dashDot', 'points', 'linesPoints'}


class ValidationError(Exception):
    """A validationDict that cannot be read or used"""


# * * * * * * * * * * * * * * * * Data files * * * * * * * * * * * * * * * * #

def _naturalKey(path):
    """A key that orders the numbers in paths by value: 0.5 < 2 < 10"""
    return [
        (0, float(part)) if _isNumber(part) else (1, part)
        for part in re.split(r'([0-9][0-9.eE+-]*)', path)
    ]


def _isNumber(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def resolve(pattern, base):
    """The file a name or a pattern stands for, or None"""
    paths = sorted(glob.glob(os.path.join(base, pattern)), key=_naturalKey)

    return paths[-1] if paths else None


def readRows(path):
    """The rows of numbers of a file, the lines of text left out"""
    rows = []

    with open(path, errors='replace') as f:
        for line in f:
            line = line.split('#')[0].strip()

            if not line:
                continue

            # A comma separates the columns of the csv files, a space or a
            # tab those of everything else
            parts = line.replace(',', ' ').split()

            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue                    # a header, or a line of text

    return rows


def readColumns(files, base):
    """The rows of a list of files, side by side as paste does"""
    tables = []

    for pattern in files:
        path = resolve(pattern, base)

        if path is None:
            raise FileNotFoundError(pattern)

        tables.append(readRows(path))

    if len(tables) == 1:
        return tables[0]

    return [sum(parts, []) for parts in zip(*tables)]


# * * * * * * * * * * * * * * * * Expressions * * * * * * * * * * * * * * * * #

_functions = {
    name: getattr(math, name)
    for name in (
        'sqrt', 'exp', 'log', 'log10', 'sin', 'cos', 'tan', 'atan', 'atan2',
        'pi', 'fabs', 'erf'
    )
}
_functions.update({'abs': abs, 'min': min, 'max': max})


class Expression:
    """An expression over the columns of a row"""

    def __init__(self, text, where):
        self.text = str(text)

        try:
            self.code = compile(self.text, where, 'eval')
        except SyntaxError as e:
            raise ValidationError(
                '%s: cannot read the expression %r: %s' % (where, self.text, e)
            )

    def __call__(self, row):
        names = {'c%d' % (i + 1): v for i, v in enumerate(row)}
        names.update(_functions)

        try:
            value = eval(self.code, {'__builtins__': {}}, names)  # noqa: S307
        except (ZeroDivisionError, IndexError, KeyError, NameError,
                ValueError, OverflowError):
            return None

        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return None

        value = float(value)

        return value if math.isfinite(value) else None


# * * * * * * * * * * * * * * * * * Curves  * * * * * * * * * * * * * * * * * #

class Curve:
    """A curve of a figure: where it is read from and how it is drawn"""

    def __init__(self, name, entries, where):
        unknown = set(entries) - curveKeys

        if unknown:
            raise ValidationError(
                '%s: unknown entries %s' % (where, ', '.join(sorted(unknown)))
            )

        if 'file' not in entries:
            raise ValidationError('%s: no file' % where)

        self.name = name
        self.where = where

        files = entries['file']
        self.files = [str(f) for f in files] if isinstance(files, list) \
            else [str(files)]

        self.x = Expression(entries.get('x', 'c1'), where + '/x')
        self.y = Expression(entries.get('y', 'c2'), where + '/y')
        self.label = str(entries.get('label', name))

        self.style = str(entries.get('style', 'line'))

        if self.style not in styles:
            raise ValidationError(
                '%s: style %s is not one of %s'
                % (where, self.style, ', '.join(sorted(styles)))
            )

        self.color = entries.get('color')
        self.marker = entries.get('marker')
        self.every = int(entries.get('every', 1))

        if self.every < 1:
            raise ValidationError('%s: every must be at least 1' % where)

        compare = entries.get('compare', {})

        if not isinstance(compare, dict):
            raise ValidationError('%s: compare is not a sub-dictionary' % where)

        #  A tolerance of report measures the curves without judging them
        self.compare = {}

        for reference, tolerance in compare.items():
            if tolerance == 'report':
                self.compare[reference] = None
            elif isinstance(tolerance, (int, float)) and tolerance > 0:
                self.compare[reference] = float(tolerance)
            else:
                raise ValidationError(
                    '%s: the tolerance of %s is neither a positive number '
                    'nor report' % (where, reference)
                )

        # A curve compared with a reference is one the case computes
        self.computed = (
            foamDict.switch(entries['computed']) if 'computed' in entries
            else bool(self.compare)
        )

    def points(self, base, gaps=False):
        """The (x, y) points of the curve, in the order of the file.

        A point whose value is not defined is left out, or, with gaps, kept
        with a value of None, so that a line drawn through the points
        breaks there."""
        points = []

        for row in readColumns(self.files, base):
            x = self.x(row)

            if x is None:
                continue

            y = self.y(row)

            if y is None and not gaps:
                continue

            points.append((x, y))

        return points


class Figure:
    """A figure of a case, and the curves it draws"""

    def __init__(self, name, entries, where):
        unknown = set(entries) - figureKeys

        if unknown:
            raise ValidationError(
                '%s: unknown entries %s' % (where, ', '.join(sorted(unknown)))
            )

        self.name = name
        self.entries = entries

        curves = entries.get('curves')

        if not isinstance(curves, dict) or not curves:
            raise ValidationError('%s: no curves' % where)

        self.curves = [
            Curve(n, e, where + '/curves/' + n) for n, e in curves.items()
            if _isDict(e, where + '/curves/' + n)
        ]

        names = {c.name for c in self.curves}

        for c in self.curves:
            for reference in c.compare:
                if reference not in names:
                    raise ValidationError(
                        '%s: %s compares with %s, which is not a curve of '
                        'the figure' % (where, c.name, reference)
                    )

                if reference == c.name:
                    raise ValidationError(
                        '%s: %s is compared with itself' % (where, c.name)
                    )

    def get(self, key, default=None):
        return self.entries.get(key, default)

    def curve(self, name):
        return next(c for c in self.curves if c.name == name)


def _isDict(entries, where):
    if not isinstance(entries, dict):
        raise ValidationError('%s is not a sub-dictionary' % where)

    return True


class Comparison:
    """A curve the case computes against a reference curve"""

    def __init__(self, figure, computed, reference, tolerance, single):
        self.figure = figure
        self.computed = computed
        self.reference = reference
        self.tolerance = tolerance

        # Named by the figure and the reference, and by the computed curve
        # too where the figure compares more than one
        self.name = (
            figure.name + '/' + reference.name if single
            else figure.name + '/' + computed.name + '/' + reference.name
        )


class ValidationDict:
    """The figures and comparisons of a case"""

    def __init__(self, case):
        self.case = case
        self.base = os.path.join(case, 'validation')
        self.path = os.path.join(self.base, dictName)

        try:
            entries = foamDict.read(self.path)
        except foamDict.FoamDictError as e:
            raise ValidationError(str(e))

        unknown = set(entries) - {'figures'}

        if unknown:
            raise ValidationError(
                '%s: unknown entries %s'
                % (self.path, ', '.join(sorted(unknown)))
            )

        figures = entries.get('figures', {})

        self.figures = [
            Figure(n, e, self.path + ': ' + n) for n, e in figures.items()
            if _isDict(e, self.path + ': ' + n)
        ]

        self.comparisons = []

        for figure in self.figures:
            compared = [c for c in figure.curves if c.compare]

            for c in compared:
                for reference, tolerance in c.compare.items():
                    self.comparisons.append(
                        Comparison(
                            figure, c, figure.curve(reference), tolerance,
                            len(compared) == 1
                        )
                    )

        names = [c.name for c in self.comparisons]

        if len(set(names)) != len(names):
            raise ValidationError('%s: two comparisons share a name' % self.path)


def exists(case):
    return os.path.isfile(os.path.join(case, 'validation', dictName))


# * * * * * * * * * * * * * * * *  Comparing  * * * * * * * * * * * * * * * * #

def interpolate(points, abscissae, x):
    """The value of a curve at x, or None outside the curve"""
    if not points or x < abscissae[0] or x > abscissae[-1]:
        return None

    i = bisect.bisect_left(abscissae, x)

    if i == 0:
        return points[0][1]

    x0, y0 = points[i - 1]
    x1, y1 = points[i]

    if x1 == x0:
        return y0

    return y0 + (y1 - y0)*(x - x0)/(x1 - x0)


def measure(reference, computed):
    """The rms and the largest difference between a computed curve and a
    reference, over the points of the reference the computed curve covers,
    as fractions of the scale of the reference, and the number of those
    points; None where the curves do not overlap"""
    reference = sorted(reference)
    computed = sorted(computed)

    values = [y for _, y in reference]
    scale = max(max(values) - min(values), max(abs(y) for y in values))

    if scale <= 0.0:
        scale = 1.0

    errors = []
    abscissae = [x for x, _ in computed]

    for x, y in reference:
        value = interpolate(computed, abscissae, x)

        if value is not None:
            errors.append(abs(value - y)/scale)

    if not errors:
        return None

    rms = math.sqrt(sum(e*e for e in errors)/len(errors))

    return rms, max(errors), len(errors)
