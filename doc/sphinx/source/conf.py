# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import subprocess


# -- Project information -----------------------------------------------------

project = 'cRacklet'
copyright = '2021, Fabian Barras, Thibault Roch, Philippe H Geubelle, Jean-François Molinari'
author = 'Fabian Barras, Thibault Roch, Philippe H Geubelle, Jean-François Molinari'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
    'breathe',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

html_logo = '../../logo.svg'


# -- Extension configuration -------------------------------------------------

# If on RTD build, run doxygen
on_read_the_docs = os.environ.get('READTHEDOCS') == 'True'

if on_read_the_docs:
    subprocess.call('cd ../../; mkdir -p build/doxygen; '
                    + 'doxygen doxygen/Doxyfile', shell=True)
    
breathe_projects = {
    'cRacklet': '../../build/doxygen/xml'
}

breathe_default_project = 'cRacklet'

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
}
