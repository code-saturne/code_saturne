## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft  \
-I$(top_srcdir)/src/fvm  \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/pprt \
$(CPPFLAGS_PLE) $(MPI_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_gwf.h \
cs_gwf_hydraulic_model.h \
cs_gwf_param.h \
cs_gwf_priv.h \
cs_gwf_soil.h \
cs_gwf_sspf.h \
cs_gwf_toolbox.h \
cs_gwf_tpf.h \
cs_gwf_tracer.h \
cs_gwf_uspf.h \
cs_gwf_headers.h

# Library source files

noinst_LIBRARIES = libcsgwf.a
libcsgwf_a_SOURCES = \
cs_gwf.cpp \
cs_gwf_priv.cpp \
cs_gwf_soil.cpp \
cs_gwf_sspf.cpp \
cs_gwf_toolbox.cpp \
cs_gwf_tpf.cpp \
cs_gwf_tracer.cpp \
cs_gwf_uspf.cpp

libcsgwf_a_LIBADD =

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
