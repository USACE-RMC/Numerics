"""Sphinx configuration for the pynumerics documentation.

Build from ``python/pynumerics`` with the ``docs`` extra installed::

    python -m sphinx -W --keep-going -b html docs docs/_build/html

The API reference imports every public module. Importing pynumerics does not start .NET, so the
build needs neither the runtime nor ``Numerics.dll``; the code blocks of the guide pages run in
``tests/test_docs.py`` instead.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pynumerics._version import __version__  # noqa: E402

project = "pynumerics"
author = "Haden Smith"
copyright = "U.S. Army Corps of Engineers Risk Management Center"
version = release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build"]
templates_path = ["_templates"]

myst_heading_anchors = 3

autosummary_generate = True
autosummary_imported_members = False
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "none"

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False

html_theme = "alabaster"
html_title = f"pynumerics {release}"
html_theme_options = {"description": "RMC.Numerics for Python", "fixed_sidebar": True}
