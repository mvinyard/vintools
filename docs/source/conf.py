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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'vintools'
copyright = '2021, Michael E. Vinyard'
author = 'Michael E. Vinyard'

# The full version, including alpha/beta/rc tags
release = '0.0.54'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
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
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = "../../documents/scdiffeq_logo_simple.png"

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/<your-org>/<your-repo>",
            "icon": "fab fa-github-square",
        },
        {
            "name": "GitLab",
            "url": "https://gitlab.com/<your-org>/<your-repo>",
            "icon": "fab fa-gitlab",
        },
        {
            "name": "Twitter",
            "url": "https://twitter.com/<your-handle>",
            "icon": "fab fa-twitter-square",
        },
    ],
}

html_theme_options = {
"navbar_start": ["navbar-logo"],
"navbar_center": ["navbar-nav"],
"navbar_end": ["navbar-icon-links"]
}

html_sidebars = {
    "<page_pattern>": ["list", "of", "templates"]
}

html_theme_options = {
  "footer_items": ["copyright", "sphinx-version"],
}
