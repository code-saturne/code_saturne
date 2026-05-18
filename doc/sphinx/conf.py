import os

project = 'code_saturne'

# Breathe is optionnal - only if Doxygen XML is available
builddir = os.environ.get("builddir", "")
doxygen_xml = os.path.join(builddir, "docs/doxygen/src/xml/") if builddir else ""
doxygen_available = builddir and os.path.exists(
    os.path.join(doxygen_xml, "index.xml"))

if doxygen_available:
    extensions = ['breathe', 'myst_parser']

    myst_enable_extensions = ['colon_fence']
    breathe_projects = {"code_saturne": doxygen_xml}
    breathe_default_project = "code_saturne"
else:
    extensions = []
    if builddir:
        print(f"WARNING: Doxygen XML not found in {doxygen_xml}")
        print("         Building documentation without C++ API reference.")
    else:
        print("WARNING: builddir not set. Building model documentation only.")

html_theme = 'sphinx_rtd_theme'
nitpicky = False
cpp_id_attributes = []
cpp_paren_attributes = []
exclude_patterns = [
    '_build',
    'api/cs__math_8h.rst',
    'api/cs__math_8cpp.rst',
    'api/cs__math__hip_8h.rst',
]
