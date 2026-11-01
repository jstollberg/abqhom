# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Sphinx configuration for the abqhom documentation."""

import os
import sys


sys.path.insert(0, os.path.abspath("../../src"))


# -- Project information -----------------------------------------------------

project = "abqhom"
copyright = "2026, Jonathan Stollberg"
author = "Jonathan Stollberg"
version = "0.1"
release = "0.1.0"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
language = "en"

# Abaqus is unavailable on Read the Docs. Gmsh and Paraqus are not required to
# render their interfaces, so mocking keeps the documentation build isolated
# from all solver and meshing installations.
autodoc_mock_imports = [
    "abaqus",
    "abaqusConstants",
    "gmsh",
    "mesh",
    "paraqus",
    "regionToolset",
]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autosummary_generate = True

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = "abqhom documentation"
html_short_title = "abqhom"


# -- Options for EPUB output -------------------------------------------------

epub_show_urls = "footnote"
