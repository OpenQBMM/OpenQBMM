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
#     qbmm.foamDict
#
# Description
#     Read a dictionary written in the format of OpenFOAM, so that the files
#     the scripts read look like the files of the case they sit in.
#
#     What is read is the part of the format a case file of the scripts
#     uses: sub-dictionaries, entries ended by a semicolon, lists in
#     parentheses, strings in double quotes, and comments of both kinds. An
#     entry is a word, a number, a string, or a list of them, and the
#     FoamFile header is read and dropped. The directives (#include, ...)
#     and the macros ($name) of OpenFOAM are not expanded: a file that uses
#     them is refused rather than read wrongly.
#
#     A sub-dictionary is returned as a dict, which keeps the order of its
#     entries, so that the curves of a figure are drawn in the order they
#     are written.
#
#------------------------------------------------------------------------------

import re


class FoamDictError(Exception):
    """A file that is not a dictionary this module reads"""


_token = re.compile(
    r'''
      (?P<space>\s+)
    | (?P<lineComment>//[^\n]*)
    | (?P<blockComment>/\*.*?\*/)
    | (?P<string>"(?:[^"\\\n]|\\.)*")
    | (?P<open>"[^\n]*)
    | (?P<punct>[{}();])
    | (?P<word>[^\s{}();"]+)
    ''',
    re.S | re.X
)


def _tokens(text, path):
    """The tokens of a text, each with the line it is on"""
    tokens = []
    pos = 0
    line = 1

    while pos < len(text):
        m = _token.match(text, pos)

        if not m:
            raise FoamDictError(
                '%s:%d: cannot read from %r' % (path, line, text[pos:pos + 20])
            )

        kind = m.lastgroup
        value = m.group(kind)

        if kind == 'open':
            raise FoamDictError(
                '%s:%d: a string is not closed on its line' % (path, line)
            )

        if kind == 'string':
            # Only a quote is escaped, so that the backslashes of the
            # mathematics of a label are kept as they are written
            tokens.append(('string', value[1:-1].replace('\\"', '"'), line))
        elif kind in ('punct', 'word'):
            tokens.append((kind, value, line))

        line += value.count('\n')
        pos = m.end()

    return tokens


def _number(word):
    try:
        return int(word)
    except ValueError:
        pass

    try:
        return float(word)
    except ValueError:
        return None


class _Parser:
    def __init__(self, tokens, path):
        self.tokens = tokens
        self.i = 0
        self.path = path

    def error(self, message):
        line = (
            self.tokens[self.i][2] if self.i < len(self.tokens)
            else (self.tokens[-1][2] if self.tokens else 1)
        )

        raise FoamDictError('%s:%d: %s' % (self.path, line, message))

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self):
        token = self.peek()

        if token is None:
            self.error('the file ends inside an entry')

        self.i += 1

        return token

    def dictionary(self, closing):
        entries = {}

        while True:
            token = self.peek()

            if token is None:
                if closing:
                    self.error('a sub-dictionary is not closed')
                return entries

            kind, value, _ = token

            if kind == 'punct' and value == '}' and closing:
                self.i += 1
                return entries

            if kind not in ('word', 'string'):
                self.error('expected a keyword, found %r' % value)

            if value.startswith('#') or value.startswith('$'):
                self.error(
                    '%s is not expanded by this reader; write the entry out'
                    % value
                )

            self.i += 1
            keyword = value

            if keyword in entries:
                self.error('%s is given twice' % keyword)

            nextToken = self.peek()

            if nextToken and nextToken[:2] == ('punct', '{'):
                self.i += 1
                entries[keyword] = self.dictionary(closing=True)
                continue

            values = []

            while True:
                kind, value, _ = self.take()

                if kind == 'punct' and value == ';':
                    break

                if kind == 'punct' and value == '(':
                    values.append(self.list())
                elif kind == 'punct':
                    self.error('unexpected %r in %s' % (value, keyword))
                else:
                    values.append(self.value(kind, value))

            if not values:
                self.error('%s has no value' % keyword)

            entries[keyword] = values[0] if len(values) == 1 else values

    def list(self):
        items = []

        while True:
            kind, value, _ = self.take()

            if kind == 'punct' and value == ')':
                return items

            if kind == 'punct' and value == '(':
                items.append(self.list())
            elif kind == 'punct' and value == '{':
                items.append(self.dictionary(closing=True))
            elif kind == 'punct':
                self.error('unexpected %r in a list' % value)
            else:
                items.append(self.value(kind, value))

    def value(self, kind, value):
        if kind == 'string':
            return value

        number = _number(value)

        return value if number is None else number


def parse(text, path='<string>'):
    """The entries of a dictionary given as text"""
    entries = _Parser(_tokens(text, path), path).dictionary(closing=False)
    entries.pop('FoamFile', None)

    return entries


def read(path):
    """The entries of the dictionary in a file"""
    with open(path, errors='replace') as f:
        return parse(f.read(), path)


def switch(value):
    """The value of an entry written as a switch of OpenFOAM"""
    words = {
        'yes': True, 'on': True, 'true': True, 'y': True, 't': True,
        'no': False, 'off': False, 'false': False, 'n': False, 'f': False,
        'none': False
    }

    if isinstance(value, bool):
        return value

    if isinstance(value, str) and value.lower() in words:
        return words[value.lower()]

    raise FoamDictError('%r is not a switch' % (value,))
