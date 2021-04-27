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
import shutil
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import subprocess

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.intersphinx',
              'sphinx.ext.todo',
              'sphinx.ext.coverage',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.ifconfig',
              'breathe']

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = '_static/logo.svg'

if os.environ.get('READTHEDOCS', None) is not None:
    print("${READTHEDOCS} = " + os.environ.get('READTHEDOCS', None))
if read_the_docs_build:
    cRacklet_path = "."
else:
    cRacklet_path = "@CMAKE_CURRENT_BINARY_DIR@"
    os.makedirs("@CMAKE_CURRENT_BINARY_DIR@/_static", exist_ok=True)
    shutil.copyfile(
        os.path.join('@CMAKE_CURRENT_SOURCE_DIR@', html_logo),
        os.path.join(cRacklet_path, html_logo))
    
print("cRacklet_path = '{}'".format(cRacklet_path))
subprocess.call('ls; pwd', shell=True)
subprocess.call("cd {} && doxygen".format(cRacklet_path), shell=True)

breathe_projects = {"cRacklet": os.path.join(cRacklet_path, "doxygenxml")}
breathe_default_project = "cRacklet"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# -- Project information -----------------------------------------------------

project = 'cRacklet'
copyright = '2012-2013 Fabian Barras, 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)'
author = 'Fabian Barras, Thibault Roch, Philippe H Geubelle, Jean-François Molinari'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = 'v0.1'
# The full version, including alpha/beta/rc tags.
release = 'v0.1 α'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["CmakeLists.txt"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    html_theme = 'default'
else:
    html_theme = 'sphinx_rtd_theme'

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
html_sidebars = {
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

# -- Options for HTMLHelp output ------------------------------------------
# Output file base name for HTML help builder.
htmlhelp_basename = 'cRackletdoc'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}


# -- Extension configuration -------------------------------------------------
    
intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
}
