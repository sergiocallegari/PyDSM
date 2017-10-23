# -*- coding: utf-8 -*-

# Copyright (c) 2016, Sergio Callegari
# All rights reserved.

# This file is part of PyDSM.

# PyDSM is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyDSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

# This file includes code from Sphinx 1.3.5
# covered by the following copyright and permission notice

# Copyright (c) 2007-2016 by the Sphinx team
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Monkey patch Sphinx, in order to use the ltytable package

* This is to provide better support support for multi page tables.
  * The ltytable package is not an official ctan package, rather an attempt
    at mixing tabulary and longtable, based on a proposal by David Carlise.
    see http://tex.stackexchange.com/questions/78075/multi-page-with-tabulary
  * Namely use the longtabu environment in place of longable.
  * The new ltabulary environment supports tabulary features such as the
    possibility of having automatically sized paragraph type columns.  This
    is needed by autosummary tables to avoid overfull boxes in the
    function descriptions.

* Also have autosummary use the new column types
"""

from __future__ import unicode_literals
import os
from os import path
from sphinx.writers import latex
from sphinx.builders.latex import LaTeXBuilder
from sphinx.locale import _
from sphinx.ext.autosummary import autosummary_table, Autosummary
from sphinx import addnodes
from docutils import nodes
from docutils.statemachine import ViewList
from sphinx.util.console import bold
from sphinx.util.osutil import copyfile

latex_package_path = (path.dirname(path.realpath(__file__)) +
                      '/../texinputs/')

old_finish = LaTeXBuilder.finish


def new_finish(self):
    old_finish(self)
    self.info(bold('copying more TeX support files...'))
    for filename in os.listdir(latex_package_path):
        if not filename.startswith('.'):
            copyfile(path.join(latex_package_path, filename),
                     path.join(self.outdir, filename))


def depart_table(self, node):
    if self.table.rowcount > 30:
        self.table.longtable = True
    self.popbody()
    if not self.table.longtable and self.table.caption is not None:
        self.body.append('\n\n\\begin{threeparttable}\n'
                         '\\capstart\\caption{')
        for caption in self.table.caption:
            self.body.append(caption)
        self.body.append('}')
        for id in self.next_table_ids:
            self.body.append(self.hypertarget(id, anchor=False))
        if node['ids']:
            self.body.append(self.hypertarget(node['ids'][0],
                                              anchor=False))
        self.next_table_ids.clear()
    if self.table.longtable:
        self.body.append('\n\\begin{ltabulary}')
        endmacro = '\\end{ltabulary}\n\n'
    elif self.table.has_verbatim:
        self.body.append('\n\\begin{tabular}')
        endmacro = '\\end{tabular}\n\n'
    elif self.table.has_problematic and not self.table.colspec:
        # if the user has given us tabularcolumns, accept them and use
        # tabulary nevertheless
        self.body.append('\n\\begin{tabular}')
        endmacro = '\\end{tabular}\n\n'
    else:
        self.body.append('\n\\begin{tabulary}{\\linewidth}')
        endmacro = '\\end{tabulary}\n\n'
    if self.table.colspec:
        self.body.append(self.table.colspec)
    else:
        if self.table.has_problematic:
            colwidth = 0.95 / self.table.colcount
            colspec = ('p{%.3f\\linewidth}|' % colwidth) * \
                self.table.colcount
            self.body.append('{|' + colspec + '}\n')
        elif self.table.longtable:
            self.body.append('{|' + ('l|' * self.table.colcount) + '}\n')
        else:
            self.body.append('{|' + ('L|' * self.table.colcount) + '}\n')
    if self.table.longtable and self.table.caption is not None:
        self.body.append(u'\\caption{')
        for caption in self.table.caption:
            self.body.append(caption)
        self.body.append('}')
        for id in self.next_table_ids:
            self.body.append(self.hypertarget(id, anchor=False))
        self.next_table_ids.clear()
        self.body.append(u'\\\\\n')
    if self.table.longtable:
        self.body.append('\\hline\n')
        self.body.extend(self.tableheaders)
        self.body.append('\\endfirsthead\n\n')
        self.body.append('\\multicolumn{%s}{c}%%\n' % self.table.colcount)
        self.body.append(r'{{\textsf{\tablename\ \thetable{} -- %s}}} \\'
                         % _('continued from previous page'))
        self.body.append('\n\\hline\n')
        self.body.extend(self.tableheaders)
        self.body.append('\\endhead\n\n')
        self.body.append((r'\hline \multicolumn{%s}{|r|}{{\textsf{%s}}}' +
                          r'\\ \hline')
                         % (self.table.colcount,
                            _('Continued on next page')))
        self.body.append('\n\\endfoot\n\n')
        self.body.append('\\endlastfoot\n\n')
    else:
        self.body.append('\\hline\n')
        self.body.extend(self.tableheaders)
    self.body.extend(self.tablebody)
    self.body.append(endmacro)
    if not self.table.longtable and self.table.caption is not None:
        self.body.append('\\end{threeparttable}\n\n')
    self.unrestrict_footnote(node)
    self.table = None
    self.tablebody = None


def get_table(self, items):
    """Generate a proper list of table nodes for autosummary:: directive.

    *items* is a list produced by :meth:`get_items`.
    """
    table_spec = addnodes.tabular_col_spec()
    table_spec['spec'] = 'lL'

    table = autosummary_table('')
    real_table = nodes.table('', classes=['longtable'])
    table.append(real_table)
    group = nodes.tgroup('', cols=2)
    real_table.append(group)
    group.append(nodes.colspec('', colwidth=10))
    group.append(nodes.colspec('', colwidth=90))
    body = nodes.tbody('')
    group.append(body)

    def append_row(*column_texts):
        row = nodes.row('')
        for text in column_texts:
            node = nodes.paragraph('')
            vl = ViewList()
            vl.append(text, '<autosummary>')
            self.state.nested_parse(vl, 0, node)
            try:
                if isinstance(node[0], nodes.paragraph):
                    node = node[0]
            except IndexError:
                pass
            row.append(nodes.entry('', node))
        body.append(row)

    for name, sig, summary, real_name in items:
        qualifier = 'obj'
        if 'nosignatures' not in self.options:
            col1 = ':%s:`%s <%s>`\ %s' % (qualifier, name, real_name, sig)
        else:
            col1 = ':%s:`%s <%s>`' % (qualifier, name, real_name)
        col2 = summary
        append_row(col1, col2)

    return [table_spec, table]


def setup(app):
    # First of all, assure that the latex package files are copyed to the
    # latex build directory
#    LaTeXBuilder.finish = new_finish
#    # Then, assure that the ltytable package is loaded in place of the
#    # longtable package
#    latex.LaTeXTranslator.default_elements['ltytable'] = (
#        '\\usepackage{ltytable}')
#    latex.HEADER = (latex.HEADER
#                    .replace("%(longtable)s",
#                             "%(ltytable)s"))
#    # And use ltabulary in place of longtable
#    latex.LaTeXTranslator.depart_table = depart_table
#    # also redefine the get_table func in autosummary
#    Autosummary.get_table = get_table
    return {'version': '0.1'}
