# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
import os.path
import sys

# -- Project information -----------------------------------------------------
project = 'ASKIPY'
copyright = '2024, Wolfgang Friederich'
author = 'Wolfgang Friederich'
release = '1.1'

# -- General configuration ---------------------------------------------------

sys.path.insert(0, os.path.abspath("C:/Users/wolle/work_git/ASKI_SEM3D_CAR"))

extensions = ['sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------

#html_theme = 'pydata_sphinx_theme'
html_theme = 'classic' #'sphinxdoc'
html_static_path = ['_static']
html_theme_options = {
    "stickysidebar": "true"
}
