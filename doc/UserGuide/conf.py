# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Unified Post Processor Users Guide'
copyright = '2020'
author = ' '

# The full version, including alpha/beta/rc tags
release = ' '


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
]
autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 4

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = '.rst'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The master toctree document.
master_doc = 'index'

# -- Linkcheck options -------------------------------------------------

# Avoid a 403 Forbidden error when accessing certain links (e.g., noaa.gov)
# Can be found using navigator.userAgent inside a browser console.
user_agent = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36"

# Ignore working links that cause a linkcheck 403 error.
linkcheck_ignore = []

linkcheck_allowed_redirects = {r"https://github.com/JCSDA/crtm/wiki/.*/.*": 
                                 r"https://raw.githubusercontent.com/wiki/JCSDA/crtm/.*/.*",
                              }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_path = ["_themes", ]
html_logo = 'https://github.com/ufs-community/ufs/wiki/images/ufs.png'

# html_theme_options = {}
html_theme_options = {
    "body_max_width": "none",
    "navigation_depth": 6,
    }

# html_sidebar_options = {}
html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_context = {}

# Add custom style sheets
def setup(app):
    app.add_css_file('custom.css')  # may also be an URL
    app.add_css_file('theme_overrides.css')  # may also be an URL

# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
   'hpc-stack': ('https://hpc-stack-epic.readthedocs.io/en/develop/', None),
   'spack-stack': ('https://spack-stack.readthedocs.io/en/develop/', None),
   'ufs-wm': ('https://ufs-weather-model.readthedocs.io/en/develop/', None),
   'srw': ('https://ufs-srweather-app.readthedocs.io/en/develop', None),
}
