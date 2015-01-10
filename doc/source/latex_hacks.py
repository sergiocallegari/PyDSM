# -*- coding: utf-8 -*-


def setup(app):
    from sphinx.util import texescape
    texescape.tex_replacements = [
        # map TeX special chars
        (u'$', ur'\$'),
        (u'%', ur'\%'),
        (u'&', ur'\&'),
        (u'#', ur'\#'),
        (u'_', ur'\_'),
        (u'{', ur'\{'),
        (u'}', ur'\}'),
        (u'[', ur'{[}'),
        (u']', ur'{]}'),
        (u'`', ur'{}`'),
        (u'\\', ur'\textbackslash{}'),
        (u'~', ur'\textasciitilde{}'),
        (u'<', ur'\textless{}'),
        (u'>', ur'\textgreater{}'),
        (u'^', ur'\textasciicircum{}'),
        # map special Unicode characters to TeX commands
        (u'¶', ur'\P{}'),
        (u'§', ur'\S{}'),
        (u'€', ur'\texteuro{}'),
        (u'∞', ur'\(\infty\)'),
        (u'±', ur'\(\pm\)'),
        (u'→', ur'\(\rightarrow\)'),
        (u'‣', ur'\(\rightarrow\)'),
        # used to separate -- in options
        (u'\ufeff', ur'{}'),
        # map some special Unicode characters to similar ASCII ones
        (u'─', ur'-'),
        (u'⎽', ur'\_'),
        (u'╲', ur'\textbackslash{}'),
        (u'|', ur'\textbar{}'),
        (u'│', ur'\textbar{}'),
        (u'ℯ', ur'e'),
        (u'ⅈ', ur'i'),
        (u'₁', ur'1'),
        (u'₂', ur'2')
        ]
    from sphinx.writers.latex import LaTeXTranslator
    del LaTeXTranslator.default_elements['longtable']
    LaTeXTranslator.default_elements['tabu'] = '\usepackage{tabu}'
