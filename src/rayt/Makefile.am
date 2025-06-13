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
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)
AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

csrayt_includedir=$(pkgincludedir)/rayt
csrayt_include_HEADERS = \
cs_rad_transfer_absorption.h \
cs_rad_transfer_dir.h \
cs_rad_transfer_fields.h \
cs_rad_transfer_options.h \
cs_rad_transfer_wall_flux.h \
cs_rad_transfer_bcs.h \
cs_rad_transfer_fsck.h \
cs_rad_transfer_rcfsk.h \
cs_rad_transfer_adf_models.h \
cs_rad_transfer_pun.h \
cs_rad_transfer_solve.h \
cs_rad_transfer_modak.h \
cs_rad_transfer_source_terms.h \
cs_rad_transfer_restart.h \
cs_rad_transfer.h \
cs_rad_headers.h

# Library source files

noinst_LIBRARIES = libcsrayt.a
libcsrayt_a_SOURCES = \
cs_rad_transfer_absorption.cpp \
cs_rad_transfer.cpp \
cs_rad_transfer_dir.cpp \
cs_rad_transfer_fields.cpp \
cs_rad_transfer_options.cpp \
cs_rad_transfer_bcs.cpp \
cs_rad_transfer_wall_flux.cpp \
cs_rad_transfer_fsck.cpp \
cs_rad_transfer_rcfsk.cpp \
cs_rad_transfer_adf_models.cpp \
cs_rad_transfer_pun.cpp \
cs_rad_transfer_solve.cpp \
cs_rad_transfer_modak.cpp \
cs_rad_transfer_source_terms.cpp \
cs_rad_transfer_restart.cpp

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
