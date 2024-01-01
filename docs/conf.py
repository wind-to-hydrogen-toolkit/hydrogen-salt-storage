# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
sys.path.append('..')

project = 'hydrogen-salt-storage-optimisation'
copyright = '2023, Nithiya Streethran'
author = 'Nithiya Streethran'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc']

# disable sorting of functions by alphabetical order
autodoc_member_order = 'bysource'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

html_theme_options = {
    'icon_links': [
        {
            # Label for this link
            'name': 'GitHub',
            # URL where the link will redirect
            'url': 'https://github.com/nmstreethran/hydrogen-salt-storage-optimisation',  # required
            # Icon class (if 'type': 'fontawesome'), or path to local image (if 'type': 'local')
            'icon': 'fa-brands fa-github',
            # The type of image to be used (see below for details)
            'type': 'fontawesome',
        }
    ],
    'navbar_align': 'right'
}
