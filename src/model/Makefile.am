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

# PLE dependency
if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part
AM_CPPFLAGS = \
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(CPPFLAGS) $(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)
AM_LDFLAGS =

# Public header files (to be installed)

csmodel_includedir=$(pkgincludedir)/model
csmodel_include_HEADERS = \
cs_turbulence_rules_manager.h \
cs_thermal_rules_manager.h \
cs_time_stepping_rules_manager.h \
cs_fields_rules_manager.h \
cs_turbulent_flux_rules_manager.h \
cs_boundary_conditions_rules_manager.h \
cs_lagr_rules_manager.h \
cs_elec_rules_manager.h

# Library source files

noinst_LIBRARIES = libcsmodel.a

libcsmodel_a_SOURCES = \
cs_turbulence_rules_manager.cpp \
cs_thermal_rules_manager.cpp \
cs_time_stepping_rules_manager.cpp \
cs_fields_rules_manager.cpp \
cs_turbulent_flux_rules_manager.cpp \
cs_boundary_conditions_rules_manager.cpp \
cs_lagr_rules_manager.cpp \
cs_elec_rules_manager.cpp
