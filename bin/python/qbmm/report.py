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
#     qbmm.report
#
# Description
#     Collect the state of the cases that were run into one record: whether
#     each case ran to its end and how long it took, how far it sits from its
#     reference data, whether it still matches its baseline, and the unit
#     tests of the same build. The record is what the validation dashboard
#     reads, and what history keeps from one commit to the next.
#
#     Nothing is run here: a case is reported in the state its last run left
#     it in, and a case that was not run is reported as such.
#
#------------------------------------------------------------------------------

import datetime
import json
import os
import re
import shutil
import subprocess
import sys

from . import validation


# * * * * * * * * * * * * * * * * *  Cases  * * * * * * * * * * * * * * * * * #

def findCases(roots):
    """The cases under roots that declare a validation or record a baseline,
    in the order of their paths, each with the directory above its root,
    which its path is reported from"""
    cases = []

    for root in roots:
        top = os.path.dirname(os.path.abspath(root))

        for directory, dirs, files in os.walk(root):
            dirs[:] = sorted(
                d for d in dirs
                if not d.startswith('processor')
                and d not in ('constant', 'system', 'postProcessing')
                and not _isTime(d)
            )

            if (
                os.path.isfile(os.path.join(directory, 'validation',
                                            validation.dictName))
             or os.path.isfile(os.path.join(directory, 'regression',
                                            'baseline'))
            ):
                cases.append((directory, top))

    return sorted(set(cases))


def _isTime(name):
    try:
        float(name)
        return True
    except ValueError:
        return name.endswith('.orig')


def application(case):
    """The application the controlDict of a case names, or None"""
    path = os.path.join(case, 'system', 'controlDict')

    if not os.path.isfile(path):
        return None

    m = re.search(r'^\s*application\s+(\w+)\s*;', open(path).read(), re.M)

    return m.group(1) if m else None


def solverRun(case, app):
    """What the log of the application says of the last run: whether it
    ran, whether it reached its end, its execution time in seconds and the
    date of the log"""
    log = os.path.join(case, 'log.' + app) if app else None

    if not log or not os.path.isfile(log):
        return {'ran': False}

    completed = False
    seconds = None

    with open(log, errors='replace') as f:
        for line in f:
            m = re.match(r'ExecutionTime = ([0-9.eE+-]+) s', line)

            if m:
                seconds = float(m.group(1))

            if line.strip() == 'End':
                completed = True

    date = datetime.datetime.fromtimestamp(
        os.path.getmtime(log), datetime.timezone.utc
    )

    return {
        'ran': True,
        'completed': completed,
        'seconds': seconds,
        'date': date.isoformat(timespec='seconds'),
    }


def comparisons(spec):
    """The comparisons of a case, measured"""
    results = []

    for c in spec.comparisons:
        entry = {
            'name': c.name,
            'figure': c.figure.name,
            'computed': c.computed.label,
            'reference': c.reference.label,
            'tolerance': c.tolerance,
        }

        try:
            reference = c.reference.points(spec.base)
            computed = c.computed.points(spec.base)
        except FileNotFoundError as e:
            entry.update(verdict='missing', reason='no such file: %s' % e)
            results.append(entry)
            continue

        measured = (
            validation.measure(reference, computed)
            if reference and computed else None
        )

        if measured is None:
            entry.update(verdict='missing', reason='no points to compare')
            results.append(entry)
            continue

        rms, largest, n = measured

        if c.tolerance is None:
            verdict = 'report'
        elif rms <= c.tolerance:
            verdict = 'pass'
        else:
            verdict = 'fail'

        entry.update(
            rms=rms, largest=largest, points=n,
            referencePoints=len(reference), verdict=verdict
        )
        results.append(entry)

    return results


def regression(case, bin):
    """The check of a case against its baseline, with the largest relative
    difference it found"""
    if not os.path.isfile(os.path.join(case, 'regression', 'baseline')):
        return {'status': 'none'}

    out = subprocess.run(
        [sys.executable, os.path.join(bin, 'qbmmRegression'), 'check', case],
        capture_output=True, text=True
    )

    spans = [
        float(x) for x in re.findall(r'\(([0-9.eE+-]+) of span\)', out.stdout)
    ]

    if out.returncode == 0:
        status = 'pass'
    elif 'produced nothing to check' in out.stderr:
        status = 'notRun'
    elif out.returncode == 1:
        status = 'fail'
    else:
        status = 'error'

    result = {'status': status}

    if spans:
        result['largest'] = max(spans)

    problems = [
        line.strip() for line in out.stdout.splitlines()
        if line.startswith('  ')
    ]

    if problems:
        result['problems'] = problems[:20]
        result['nProblems'] = len(problems)

    return result


def status(record):
    """What a case comes to: passed, failed, or notRun"""
    run = record['run']
    judged = [c for c in record['comparisons'] if c['verdict'] != 'report']
    regression = record['regression']['status']

    produced = (
        run.get('ran')
     or any(c['verdict'] != 'missing' for c in record['comparisons'])
     or regression in ('pass', 'fail')
    )

    if not produced:
        return 'notRun'

    if run.get('ran') and not run.get('completed'):
        return 'failed'

    if run.get('allrun', 0) != 0:
        return 'failed'

    if regression in ('fail', 'error', 'notRun'):
        return 'failed'

    if any(c['verdict'] in ('fail', 'missing') for c in judged):
        return 'failed'

    return 'passed'


def collectCase(case, top, bin, allrun=None):
    """The record of a case, and its validationDict where it has one.
    allrun is the exit status of the Allrun that ran it, where known."""
    app = application(case)

    record = {
        'path': os.path.relpath(case, top),
        'application': app,
        'run': solverRun(case, app),
        'comparisons': [],
        'regression': regression(case, bin),
        'figures': [],
    }

    if allrun is not None:
        record['run']['allrun'] = allrun

    if validation.exists(case):
        try:
            spec = validation.ValidationDict(case)
        except validation.ValidationError as e:
            record['error'] = str(e)
            record['status'] = 'failed'
            return record, None

        record['comparisons'] = comparisons(spec)
    else:
        spec = None

    record['status'] = status(record)

    return record, spec


# * * * * * * * * * * * * * * * * * *  Tests  * * * * * * * * * * * * * * * * * #

def unitTests(log):
    """The unit tests a log of test/Allrun reports, in their order, each
    with the directory it ran in where the log gives it"""
    tests = []
    directory = None

    for line in open(log, errors='replace'):
        line = line.strip()

        m = re.match(r'^Running (Test-\S+) in (\S+?),?(?: which .*)?$', line)

        if m:
            directory = m.group(2)
            continue

        m = re.match(r'^(Test-\S+?): (passed|FAILED)\b(.*)$', line)

        if m:
            tests.append({
                'name': m.group(1),
                'directory': directory if directory != '.' else None,
                'passed': m.group(2) == 'passed',
                'note': m.group(3).strip(' ()'),
            })
            directory = None

    return tests


# * * * * * * * * * * * * * * * * *  Build  * * * * * * * * * * * * * * * * * #

def build(top, tier=None):
    """The commit the sources are at, and whether they differ from it"""
    def git(*args):
        out = subprocess.run(
            ['git', '-C', top] + list(args), capture_output=True, text=True
        )
        return out.stdout.strip() if out.returncode == 0 else None

    commit = git('rev-parse', 'HEAD')
    dirty = git('status', '--porcelain', '--untracked-files=no')

    return {
        'commit': commit,
        'shortCommit': commit[:8] if commit else None,
        'branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
        'subject': git('log', '-1', '--format=%s'),
        'modified': bool(dirty),
        'openfoam': os.environ.get('WM_PROJECT_VERSION'),
        'tier': tier,
        'date': datetime.datetime.now(datetime.timezone.utc)
                .isoformat(timespec='seconds'),
    }


def summary(cases, tests):
    counts = {'passed': 0, 'failed': 0, 'notRun': 0}

    for c in cases:
        counts[c['status']] += 1

    ran = [c for c in cases if c['status'] != 'notRun']
    compared = [x for c in cases for x in c['comparisons']]
    judged = [x for x in compared if x['tolerance'] is not None]

    return {
        'cases': len(cases),
        'passed': counts['passed'],
        'failed': counts['failed'],
        'notRun': counts['notRun'],
        'comparisons': len(compared),
        'comparisonsJudged': len(judged),
        'comparisonsPassed': sum(x['verdict'] == 'pass' for x in judged),
        'comparisonsFailed': sum(
            x['verdict'] in ('fail', 'missing')
            for c in ran for x in c['comparisons']
            if x['tolerance'] is not None
        ),
        'comparisonsNotRun': sum(
            x['tolerance'] is not None
            for c in cases if c['status'] == 'notRun'
            for x in c['comparisons']
        ),
        'comparisonsReported': sum(x['verdict'] == 'report' for x in compared),
        'baselines': sum(c['regression']['status'] != 'none' for c in cases),
        'baselinesPassed': sum(
            c['regression']['status'] == 'pass' for c in cases
        ),
        'baselinesNotRun': sum(
            c['regression']['status'] != 'none' for c in cases
            if c['status'] == 'notRun'
        ),
        'tests': len(tests),
        'testsPassed': sum(t['passed'] for t in tests),
    }


def historyEntry(record):
    """What history keeps of a record: the status of every case and the rms
    of every comparison, by commit"""
    return {
        'commit': record['build']['commit'],
        'shortCommit': record['build']['shortCommit'],
        'modified': record['build']['modified'],
        'date': record['build']['date'],
        'summary': record['summary'],
        'cases': {
            c['path']: {
                'status': c['status'],
                'rms': {
                    x['name']: x['rms'] for x in c['comparisons'] if 'rms' in x
                },
            }
            for c in record['cases']
        },
    }


def mergeHistory(history, entry, keep=100):
    """History with the entry of a commit replacing any earlier one of the
    same commit, newest last"""
    history = [
        h for h in history
        if not (h['commit'] == entry['commit']
                and h['modified'] == entry['modified'])
    ]
    history.append(entry)

    return history[-keep:]


# * * * * * * * * * * * * * * * * * * Site  * * * * * * * * * * * * * * * * * #

currentStyle = (
    'color:var(--color-accent-900);border-bottom:2px solid '
    'var(--color-accent);padding-bottom:2px'
)


def _element(html, tag):
    m = re.search(r'[ \t]*<%s\b.*?</%s>' % (tag, tag), html, re.S)
    return m.group(0) if m else None


def _moveDown(match):
    """A link of a page of the site, as seen from a page one directory
    below it"""
    attr, target = match.group(1), match.group(2)

    if re.match(r'^[a-z][a-z0-9+.-]*:|^/|^\.\./', target):
        return match.group(0)

    if target.startswith('#'):
        target = 'index.html' + target

    return '%s="../%s"' % (attr, target)


def siteChrome(site):
    """The header and the footer of the home page of the site, rewritten for
    the validation page: every link one directory up, the style that marks
    the page a visitor is on taken off the link of the home page's own
    section and put on the Validation link. None where the site has no
    home page with both."""
    path = os.path.join(site, 'index.html')

    if not os.path.isfile(path):
        return None

    html = open(path, errors='replace').read()

    header = _element(html, 'header')
    footer = _element(html, 'footer')

    if header is None or footer is None:
        return None

    header = re.sub(r'\b(href|src)="([^"]*)"', _moveDown, header)
    footer = re.sub(r'\b(href|src)="([^"]*)"', _moveDown, footer)

    # The Validation link is the page itself, and marked as the current one
    def current(m):
        attrs = re.sub(r'\s*style="[^"]*"', '', m.group(1))
        return '<a%s style="%s">' % (
            attrs.replace('href="../validation/"', 'href="./"'), currentStyle
        )

    header, n = re.subn(
        r'<a((?:(?!>)[^>])*href="\.\./validation/?"[^>]*)>', current, header
    )

    if n == 0:
        return None

    return header, footer


def page(template, data, chrome=None, style='industry.css'):
    """The dashboard page: its template with the data, and with the header
    and the footer of the site where they are given"""
    if chrome:
        header, footer = chrome
        template = re.sub(
            r'<!--@HEADER@-->.*?<!--@/HEADER@-->',
            lambda m: header.strip(), template, flags=re.S
        )
        template = re.sub(
            r'<!--@FOOTER@-->.*?<!--@/FOOTER@-->',
            lambda m: footer.strip(), template, flags=re.S
        )

    template = template.replace('@STYLE@', style)

    return template.replace(
        '/*@DATA@*/null', data.replace('</', '<\\/')
    )


# * * * * * * * * * * * * * * * * * Writing * * * * * * * * * * * * * * * * * #

class ReportError(Exception):
    """A report that cannot be written where it is asked for"""


def destination(site=None, output=None, style=None):
    """Where a report goes: its directory, the header and the footer of the
    site where it goes into one, and the stylesheet the page links"""
    if site:
        if style:
            raise ReportError('a report into a site uses its stylesheet')

        chrome = siteChrome(site)

        if chrome is None:
            raise ReportError(
                '%s has no index.html with a header linking validation/ and '
                'a footer' % site
            )

        return output or os.path.join(site, 'validation'), chrome, \
            '../industry.css'

    if not output:
        raise ReportError('no directory to write the report into')

    return output, None, 'industry.css'


def write(cases, output, bin, tests=None, chrome=None, styleHref=None,
          style=None, allrun=None, tier=None, log=print):
    """Collect the cases, given as (directory, top) pairs, and write the
    dashboard of them into output: the page, results.json, history.json and
    the figures. allrun maps the directory of a case to the exit status of
    the Allrun that ran it. Returns the record written."""
    from . import plotting

    allrun = allrun or {}
    records = []

    for case, top in cases:
        record, spec = collectCase(case, top, bin, allrun.get(case))

        if spec is not None and record['status'] != 'notRun':
            outDir = os.path.join(output, 'figures', record['path'])

            for figure in spec.figures:
                computed = {c.name for c in figure.curves if c.computed}

                try:
                    written, missing = plotting.draw(
                        figure, spec.base, outDir, ['svg']
                    )
                except validation.ValidationError as e:
                    log('%s: %s' % (record['path'], e))
                    continue

                # A figure without the curves the case computed shows
                # nothing of the case
                if computed & {name for name, _ in missing}:
                    for path in written:
                        os.remove(path)
                    continue

                record['figures'].append({
                    'name': figure.name,
                    'file': os.path.relpath(written[0], output),
                })

        records.append(record)

    unitTestList = unitTests(tests) if tests else []

    result = {
        'build': build(os.path.dirname(bin), tier),
        'summary': summary(records, unitTestList),
        'cases': records,
        'tests': unitTestList,
    }

    os.makedirs(output, exist_ok=True)

    with open(os.path.join(output, 'results.json'), 'w') as f:
        json.dump(result, f, indent=1)

    historyPath = os.path.join(output, 'history.json')
    history = []

    if os.path.isfile(historyPath):
        with open(historyPath) as f:
            history = json.load(f)

    history = mergeHistory(history, historyEntry(result))

    with open(historyPath, 'w') as f:
        json.dump(history, f, indent=1)

    # The page carries its data as well, so that it opens from the disk,
    # where a browser does not let it fetch the files next to it
    with open(os.path.join(bin, 'python', 'qbmm', 'dashboard.html')) as f:
        template = f.read()

    with open(os.path.join(output, 'index.html'), 'w') as f:
        f.write(page(
            template, json.dumps({'results': result, 'history': history}),
            chrome, styleHref or 'industry.css'
        ))

    if style:
        shutil.copy(style, os.path.join(output, 'industry.css'))

    return result
