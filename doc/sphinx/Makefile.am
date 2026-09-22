## Process this file with automake to produce Makefile.in
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

sphinxdir = $(pkgdatadir)/doc/sphinx

# Generate RST files from Doxygen XML and copy static files
sphinx-doc:
	echo "Copying static Sphinx files to builddir..."; \
	cp $(srcdir)/index.rst $(top_builddir)/doc/sphinx/; \
	cp $(srcdir)/section_*.md $(top_builddir)/doc/sphinx/; \
	if [ -f "$(top_builddir)/doc/doxygen/src/xml/index.xml" ]; then \
	  echo "Generating API RST from Doxygen XML..."; \
	  python3 $(srcdir)/generate_rst.py $(top_srcdir) $(top_builddir); \
	else \
	  echo "Doxygen XML not found, skipping API RST generation."; \
	fi

# Generate Sphinx HTML documentation (RST from Doxygen XML + sphinx-build)
html-local: sphinx-doc
	@if test -x "$(SPHINXBUILD)" > /dev/null 2>&1; then \
	  echo "Generating Sphinx documentation..."; \
	  builddir="$(abs_top_builddir)" $(SPHINXBUILD) -b html $(top_builddir)/doc/sphinx \
	    $(top_builddir)/doc/sphinx/html \
	    -c $(srcdir) \
	    -d $(top_builddir)/doc/sphinx/.doctrees \
	    -D breathe_projects.code_saturne=$(abs_top_builddir)/doc/doxygen/src/xml; \
	  echo "Sphinx documentation generated."; \
	else \
	  echo "sphinx-build not found, skipping Sphinx documentation."; \
	fi

# Install only if HTML already generated
install-data-local:
	@if test -d "$(top_builddir)/doc/sphinx/html"; then \
	  $(MKDIR_P) "$(DESTDIR)$(pkgdatadir)/doc/sphinx"; \
	  cp -R $(top_builddir)/doc/sphinx/html \
	    $(DESTDIR)$(pkgdatadir)/doc/sphinx/; \
	  echo "Sphinx documentation installed."; \
	else \
	  echo "Sphinx HTML not found, skipping install. Run 'make doc' first."; \
	fi

# Clean generated files
clean-local:
	rm -rf $(top_builddir)/doc/sphinx/html \
	       $(top_builddir)/doc/sphinx/.doctrees \
	       $(top_builddir)/doc/sphinx/api \
	       $(top_builddir)/doc/sphinx/_build \
	       $(top_builddir)/doc/sphinx/index.rst \
	       $(top_builddir)/doc/sphinx/section_*.rst
