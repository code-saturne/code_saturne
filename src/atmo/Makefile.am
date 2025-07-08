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
$(CPPFLAGS) $(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)
AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_OPT) $(FCFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

csatmo_includedir=$(pkgincludedir)/atmo
csatmo_include_HEADERS = \
cs_air_props.h \
cs_at_opt_interp.h \
cs_at_data_assim.h \
cs_atmo.h \
cs_atmo_chemistry.h \
cs_atmo_aerosol_ssh.h \
cs_atmo_headers.h \
cs_atmo_imbrication.h \
cs_atmo_profile_std.h \
cs_atmo_source_terms.h \
cs_atmo_variables.h \
cs_atmo_soliva.h \
cs_atmo_solcat.h \
cs_atmo_solmoy.h \
cs_atmo_sol.h \
cs_atprke.h \
cs_intprf.h

# Library source files

noinst_LIBRARIES = libcsatmo.a
libcsatmo_a_SOURCES = \
atr1vf.f90 \
cs_air_props.cpp \
cs_at_opt_interp.cpp \
cs_at_data_assim.cpp \
cs_atmo.cpp \
cs_atmo_chemistry.cpp \
cs_atmo_aerosol_ssh.cpp \
cs_atmo_imbrication.cpp \
cs_atmo_profile_std.cpp \
cs_atmo_source_terms.cpp \
cs_atmo_variables.cpp \
cs_atmo_soliva.cpp \
cs_atmo_solcat.cpp \
cs_atmo_solmoy.cpp \
cs_atmo_sol.cpp \
cs_atprke.cpp \
cs_intprf.cpp \
kinrates.f90 \
rayigc.f90 \
rayir.f90 \
rayive.f90 \
rayso.f90

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
	-rm -rf .deps

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif
