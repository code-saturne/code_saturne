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

AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Add cs_user_modules file to distribution as it is not automatically accounted
# for by automake
EXTRA_DIST =

# Public header files (to be installed)

csmeg_includedir=$(pkgincludedir)/meg
csmeg_include_HEADERS = \
cs_meg_prototypes.h \
cs_meg_xdef_wrapper.h \
cs_meg_headers.h

# Library source files

megfiles = \
cs_meg_boundary_function.cxx \
cs_meg_fsi_struct.cxx \
cs_meg_immersed_boundaries_inout.cxx \
cs_meg_initialization.cxx \
cs_meg_post_calculator.cxx \
cs_meg_post_output.cxx \
cs_meg_post_profile.cxx \
cs_meg_source_terms.cxx \
cs_meg_volume_function.cxx \
cs_meg_xdef_wrapper.cxx

noinst_LIBRARIES = libcsmeg.a
libcsmeg_a_SOURCES = $(megfiles)

clean-local:
	-rm -rf .deps
