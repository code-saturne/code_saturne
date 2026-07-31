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

# Sphinx documentation sources (to be distributed)
sphinxdir = $(pkgdatadir)/doc/sphinx

# dist_sphinx_DATA = \
# conf.py \
# index.rst \
# section_advanced_modeling.rst \
# section_api.rst \
# section_basic_modeling.rst \
# section_developer_guide.rst \
# section_model.rst \
# section_running_computation.rst \
# section_studymanager.rst \
# section_user_examples.rst \
# generate_rst.py \
# generate_toctree.py

# Install Sphinx doc if sphinx-build is available
install-data-local:
	@if test -x "$(SPHINXBUILD)" > /dev/null 2>&1; then \
	  echo "Generating Sphinx documentation..."; \
	  if [ -f "$(top_builddir)/doc/doxygen/src/xml/index.xml" ]; then \
	    echo "Generating API RST from Doxygen XML..."; \
	    python3 $(srcdir)/generate_rst.py $(top_srcdir) $(top_builddir); \
	    python3 $(srcdir)/generate_toctree.py $(top_srcdir) $(top_builddir); \
	  else \
	    echo "Doxygen XML not found, skipping API generation."; \
	  fi && \
	  $(SPHINXBUILD) -b html $(top_builddir)/doc/sphinx $(DESTDIR)$(pkgdatadir)/doc/sphinx/html \
	    -d $(top_builddir)/doc/sphinx/.doctrees \
	    -D breathe_projects.code_saturne=$(top_builddir)/doc/doxygen/src/xml; \
	  echo "Sphinx documentation installed."; \
	else \
	  echo "sphinx-build not found, skipping Sphinx documentation."; \
	fi

clean-local:
	rm -rf _build

