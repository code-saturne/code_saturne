# Configuration Sphinx


# Configuration Sphinx
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

project = "TurbulenceModel"
release = "1.0.0"
language = "fr"

# --- Rendre importable le package du projet ---
# On remonte : .../code_saturne/model/docs/source -> repo/
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

# --- Extensions Sphinx ---
extensions = [
    "sphinx.ext.autodoc",         # Doc API auto depuis docstrings/signatures
    "sphinx.ext.viewcode",        # Lien vers le code
    "sphinx.ext.napoleon",        # Docstrings Google/NumPy -> reST propre
    "sphinx_autodoc_typehints",   # Affiche les annotations de type
    # "sphinxcontrib.mermaid",    # (optionnel) active si tu veux des diagrammes Mermaid
]

# Docstrings: prefer Google style
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Autodoc : garder l'ordre du code
autodoc_member_order = "bysource"

# If the environment does not have all Code_Saturne deps, avoid import failure:
autodoc_mock_imports = [
    "code_saturne.model.Common",
    "code_saturne.model.XMLvariables",
    "code_saturne.model.XMLmodel",
    "code_saturne.model.NumericalParamGlobalModel",
    "code_saturne.model.DefineUserScalarsModel",
    # ajoute/supprime ici selon ton environnement
]

# (Optionnel) Type hints : noms courts
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

