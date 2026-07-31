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

# Model rules XML files
modelrulesdir = $(pkgdatadir)/model
dist_modelrules_DATA = \
TurbulenceRules.xml \
ThermalRules.xml \
TimeSteppingRules.xml \
BoundaryConditionsRules.xml \
FieldsRules.xml \
TurbulentFluxRules.xml \
LagrangianRules.xml \
ElectricalRules.xml

# Documentation HTML generee depuis les XML
modeldocsdir = $(pkgdatadir)/model/docs

install-data-local:
	@if command -v sphinx-build > /dev/null 2>&1 && sphinx-build --version > /dev/null 2>&1; then \
	  echo "Generating model rules documentation..."; \
	  $(MKDIR_P) $(builddir/docs) && cd $(builddir)/docs && python3 generate_doc.py && make html; \
	  $(MKDIR_P) $(DESTDIR)$(modeldocsdir); \
	  cp -r $(srcdir)/docs/build/html/. $(DESTDIR)$(modeldocsdir)/; \
	  echo "Documentation installed in $(modeldocsdir)"; \
	else \
	  echo "sphinx-build not found, skipping documentation generation"; \
	fi

uninstall-local:
	rm -rf $(DESTDIR)$(modeldocsdir)

clean-local:
	rm -rf $(srcdir)/docs/build
	rm -rf $(srcdir)/docs/source/turbulence
	rm -rf $(srcdir)/docs/source/thermal
	rm -rf $(srcdir)/docs/source/timestepping
	rm -rf $(srcdir)/docs/source/boundary_conditions
	rm -rf $(srcdir)/docs/source/fields
