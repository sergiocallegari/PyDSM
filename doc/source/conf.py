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

extensions = ['latex_textgreek', 'sphinx.ext.autodoc', 'sphinx.ext.imgmath',
              'numpydoc', 'sphinx.ext.autosummary', 'latex_ltytable']

templates_path = ['_templates']
exclude_patterns = []

modindex_common_prefix = ['pydsm.']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
        'donate.html',
    ]
}

html_logo = "../Figures/pydsm_logo_small.png"

# -- Autosummary ---------------------------------------
autosummary_generate = True


# Ignore header of main package
def better_cut_lines(pre, post=0, what=None, name=None):
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


def setup(app):
    from sphinx.ext.autodoc import cut_lines
    app.connect('autodoc-process-docstring',
                better_cut_lines(3, what=['module'], name='pydsm'))
