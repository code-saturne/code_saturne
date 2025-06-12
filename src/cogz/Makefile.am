## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2025 EDF S.A.
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
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

cscogz_includedir=$(pkgincludedir)/cogz
cscogz_include_HEADERS = \
cs_combustion_boundary_conditions.h \
cs_combustion_bsh.h \
cs_combustion_d3p.h \
cs_combustion_ebu.h \
cs_combustion_slfm.h \
cs_combustion_gas.h \
cs_combustion_ht_convert.h \
cs_combustion_lw.h \
cs_combustion_physical_properties.h \
cs_combustion_read_data.h \
cs_soot_model.h \
cs_cogz_headers.h

# Library source files

noinst_LIBRARIES = libcscogz.a
libcscogz_a_SOURCES = \
cs_combustion_boundary_conditions.cpp \
cs_combustion_bsh.cpp \
cs_combustion_d3p.cpp \
cs_combustion_ebu.cpp \
cs_combustion_slfm.cpp \
cs_combustion_gas.cpp \
cs_combustion_ht_convert.cpp \
cs_combustion_lw.cpp \
cs_combustion_physical_properties.cpp \
cs_combustion_read_data.cpp \
cs_soot_model.cpp

# Rules for CUDA (not known by Automake)

if HAVE_CUDA
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
