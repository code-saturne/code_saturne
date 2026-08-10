# Sphinx configuration

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

project = "TurbulenceModel"
release = "1.0.0"
language = "fr"

# Allow import of project package ---
# Gu up: .../code_saturne/model/docs/source -> repo/
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

# --- Extensions Sphinx ---
extensions = [
    "sphinx.ext.autodoc",         # API auto doc from docstrings/signatures
    "sphinx.ext.viewcode",        # Link to the code
    "sphinx.ext.napoleon",        # Docstrings Google/NumPy -> clean reST
    "sphinx_autodoc_typehints",   # Displays hints
    # "sphinxcontrib.mermaid",    # (optionnal) activate for Mermaid diagrams
]

# Docstrings: prefer Google style
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Autodoc: keep the code order
autodoc_member_order = "bysource"

# If the environment does not have all code_saturne deps, avoid import failure:
autodoc_mock_imports = [
    "code_saturne.model.Common",
    "code_saturne.model.XMLvariables",
    "code_saturne.model.XMLmodel",
    "code_saturne.model.NumericalParamGlobalModel",
    "code_saturne.model.DefineUserScalarsModel",
    # add/remove here based on environnement
]

# (Optionnal) Type hints: short names
typehints_fully_qualified = False

# --- Theme ---
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
}

# Optional exclusions
exclude_patterns = []

