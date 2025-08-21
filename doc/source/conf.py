# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

sys.path.append(str(Path('.').resolve()))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pydsm._version import __version__

project = 'PyDSM'
copyright = '2012-%Y, Sergio Callegari'
author = 'Sergio Callegari'
release = __version__
version = ".".join(release.split(".", 2)[:2])


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

nitpicky = True

extensions = [
    # 'latex_textgreek',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    "sphinx.ext.intersphinx",
    'numpydoc',
    'sphinx.ext.autosummary',
    # 'latex_ltytable',
]

templates_path = ['_templates']
exclude_patterns = []

modindex_common_prefix = ['pydsm.']

intersphinx_mapping = {
    #"python": ("https://docs.python.org/3", None),
    "scipy":  ("https://docs.scipy.org/doc/scipy", None),
    # "cvxopt": ("https://cvxopt.org/userguide", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

html_logo = "../Figures/pydsm_logo_small.png"

html_show_sourcelink = False

# -- Autosummary ---------------------------------------

autosummary_generate = True


def mk_docstring_trimmer(pre, post=0, what=None, name=None):
    '''Create a docstring processing function for autodoc that trims lines.

    Returns a function conforming to the signature of docstring
    processors in autodoc, that trims the initial `pre` lines and the
    final `post` lines of the docstring for items of ttype `what`
    named `name`.
    '''
    def process(app, what_, name_, obj, options, lines):
        if what and what_ not in what:
            return
        if name and name_ != name:
            return
        del lines[:pre]
        if post:
            # remove one trailing blank line.
            if lines and not lines[-1]:
                lines.pop(-1)
            del lines[-post:]
        # make sure there is a blank line at the end
        if lines and lines[-1]:
            lines.append('')
    return process


# Set up autodoc so that the heading of the PyDSM module docstring is
# trimmed away.
def setup(app):
    from sphinx.ext.autodoc import cut_lines
    app.connect('autodoc-process-docstring',
                mk_docstring_trimmer(3, what=['module'], name='pydsm'))
