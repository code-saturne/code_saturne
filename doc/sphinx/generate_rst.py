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

"""
Generate one .rst file per Doxygen XML entry in the api/ folder.
Must be run from builddir, with srcdir and builddir passed as arguments:
    python3 <srcdir>/doc/sphinx/generate_rst.py <srcdir> <builddir>
"""

import os
import sys
import xml.etree.ElementTree as ET

if len(sys.argv) < 3:
    raise SystemExit(
        "Usage: python3 generate_rst.py <srcdir> <builddir>\n"
        "Example: python3 /path/to/src/doc/sphinx/generate_rst.py "
        "/path/to/src /path/to/build"
    )

srcdir   = sys.argv[1]
builddir = sys.argv[2]

# Doxygen XML is in builddir
xml_dir = os.path.join(builddir, "doc", "doxygen", "src", "xml")

# Output RST files go in builddir, not in sources
rst_dir = os.path.join(builddir, "doc", "sphinx", "api")
os.makedirs(rst_dir, exist_ok=True)

files = []
for f in sorted(os.listdir(xml_dir)):
    if not f.endswith(".xml") or f in ("index.xml", "Doxyfile.xml"):
        continue
    xml_path = os.path.join(xml_dir, f)
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        compounddef = root.find(".//compounddef")
        if compounddef is None:
            continue
        kind = compounddef.get("kind")
        compoundname = compounddef.find("compoundname")
        title = compounddef.find("title")
        if compoundname is None:
            continue
        real_name = compoundname.text
        display_title = title.text if title is not None else real_name
        safe_name = f.replace(".xml", "")
        if kind == "file":
            directive = f".. doxygenfile:: {real_name}"
        elif kind == "page":
            directive = f".. doxygenpage:: {real_name}"
        elif kind == "struct":
            directive = f".. doxygenstruct:: {real_name}\n   :members:"
        elif kind == "class":
            directive = f".. doxygenclass:: {real_name}\n   :members:"
        else:
            continue
        rst_path = os.path.join(rst_dir, safe_name + ".rst")
        with open(rst_path, "w") as rst:
            rst.write(f"{display_title}\n{'=' * len(display_title)}\n\n")
            rst.write(f"{directive}\n")
            rst.write(f"   :project: code_saturne\n")
        files.append((safe_name, display_title))
    except Exception as e:
        print(f"Error on {f}: {e}")

# Generate index
with open(os.path.join(rst_dir, "index.rst"), "w") as idx:
    idx.write("API Reference\n=============\n\n")
    idx.write(".. toctree::\n   :maxdepth: 1\n\n")
    for safe_name, _ in sorted(files):
        idx.write(f"   {safe_name}\n")

print(f"{len(files)} RST files generated in {rst_dir}")
