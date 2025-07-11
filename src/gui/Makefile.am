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

AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

csgui_includedir=$(pkgincludedir)/gui
csgui_include_HEADERS = \
cs_gui.h \
cs_gui_boundary_conditions.h \
cs_gui_conjugate_heat_transfer.h \
cs_gui_mesh.h \
cs_gui_mobile_mesh.h \
cs_gui_output.h \
cs_gui_particles.h \
cs_gui_radiative_transfer.h \
cs_gui_specific_physics.h \
cs_gui_util.h \
cs_tree_xml.h \
cs_gui_headers.h

# Library source files

noinst_LIBRARIES = libcsgui.a

###

libcsgui_a_SOURCES = \
cs_gui.cpp \
cs_gui_boundary_conditions.cpp \
cs_gui_conjugate_heat_transfer.cpp \
cs_gui_mesh.cpp \
cs_gui_mobile_mesh.cpp \
cs_gui_output.cpp \
cs_gui_particles.cpp \
cs_gui_radiative_transfer.cpp \
cs_gui_specific_physics.cpp \
cs_gui_util.cpp \
cs_tree_xml.cpp

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
