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

"""
Monkey patch Sphinx, to use textgreek for Greek letters

* Restrict the tex_replacement table, to stop escaping Greek letters.
  * These are much better supported with the textgreek package.
  * Converting Greek text to math causes issues with hyperref when Greek
    letters occur in section headings.
  * Converting Greek text to math makes it impossible to make it follow
    the formatting of the surrounding text *eg. bold, italics, etc).
  * Note that this change also requires the unicode option to be passed to
    hyperref
"""

from __future__ import unicode_literals
from sphinx.util import texescape
from sphinx.writers import latex

greek_letters = 'αβγδεζηθικλμνξοπρςστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ'


def setup(app):
    tex_replacements_dict = dict(texescape.tex_replacements)
    for letter in greek_letters:
        if letter in tex_replacements_dict:
            del tex_replacements_dict[letter]
    texescape.tex_replacements = tex_replacements_dict.items()
#    latex.DEFAULT_SETTINGS['textgreek'] = (
#        '\\usepackage[artemisia]{textgreek}')
#    latex.HEADER = (latex.HEADER
#                    .replace("%(usepackages)s",
#                             "%(textgreek)s\n" +
#                             "%(usepackages)s")
#                    .replace("\\usepackage{sphinx}",
#                             "\\usepackage{sphinx}\n" +
#                             "\\hypersetup{unicode}"))
    return {'version': '0.1'}
