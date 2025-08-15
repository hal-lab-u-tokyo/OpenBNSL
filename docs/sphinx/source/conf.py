# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OpenBNSL'
copyright = '2024, Ryota Miyagi'
author = 'Ryota Miyagi'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
]
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

breathe_projects = {
    "openbnsllib": "../../doxygen/build/xml",
}
breathe_default_project = "openbnsllib"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": False,
    "special-members": "__init__",
}

html_sidebars = {
    "**": [
        "globaltoc.html", 
        "relations.html", 
        "searchbox.html", 
    ],
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
