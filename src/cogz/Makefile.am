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

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_OPT) $(FCFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

cscogz_includedir=$(pkgincludedir)/cogz
cscogz_include_HEADERS = \
cs_combustion_boundary_conditions.h \
cs_combustion_bsh.h \
cs_combustion_ebu.h \
cs_combustion_gas.h \
cs_combustion_ht_convert.h \
cs_combustion_lw.h \
cs_combustion_physical_properties.h \
cs_combustion_read_data.h \
cs_soot_model.h \
cs_steady_laminar_flamelet_source_terms.h \
cs_cogz_headers.h

# Library source files

noinst_LIBRARIES = libcscogz.a
libcscogz_a_SOURCES = \
cs_combustion_boundary_conditions.cpp \
cs_combustion_bsh.cpp \
cs_combustion_ebu.cpp \
cs_combustion_gas.cpp \
cs_combustion_ht_convert.cpp \
cs_combustion_lw.cpp \
cs_combustion_physical_properties.cpp \
cs_combustion_read_data.cpp \
cs_soot_model.cpp \
cs_steady_laminar_flamelet_source_terms.cpp \
coini1.f90 \
d3phst.f90 \
d3pini.f90 \
d3pint.f90 \
lwcgfu.f90 \
lwcphy.f90 \
lwcurl.f90 \
pdflwc.f90 \
pdfpp3.f90 \
pdfpp4.f90 \
cs_steady_laminar_flamelet_init.f90 \
cs_steady_laminar_flamelet_physical_prop.f90 \
cs_steady_laminar_flamelet_verify.f90 \
cs_steady_laminar_flamelet_read_base.f90

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -f *__genmod.f90 *.mod
	-rm -rf .deps
