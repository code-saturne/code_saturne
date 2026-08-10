import os

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
        print(f"Warning: Doxygen XML not found in {doxygen_xml}")
        print("          Building documentation without C++ API reference.")
    else:
        print("Warning: builddir not set. Building model documentation only.")

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
