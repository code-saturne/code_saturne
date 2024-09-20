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
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/comb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_geom.h \
cs_join.h \
cs_join_intersect.h \
cs_join_merge.h \
cs_join_mesh.h \
cs_join_perio.h \
cs_join_post.h \
cs_join_set.h \
cs_join_split.h \
cs_join_update.h \
cs_join_util.h \
cs_mesh.h \
cs_mesh_adjacencies.h \
cs_mesh_bad_cells.h \
cs_mesh_boundary.h \
cs_mesh_boundary_layer.h \
cs_mesh_builder.h \
cs_mesh_cartesian.h \
cs_mesh_coherency.h \
cs_mesh_coarsen.h \
cs_mesh_connect.h \
cs_mesh_extrude.h \
cs_mesh_from_builder.h \
cs_mesh_group.h \
cs_mesh_halo.h \
cs_mesh_headers.h \
cs_mesh_location.h \
cs_mesh_intersect.h \
cs_mesh_quality.h \
cs_mesh_quantities.h \
cs_mesh_refine.h \
cs_mesh_remove.h \
cs_mesh_save.h \
cs_mesh_to_builder.h \
cs_mesh_warping.h \
cs_mesh_smoother.h \
cs_partition.h \
cs_stl.h \
cs_symmetry_faces_filter.h

# Library source files

noinst_LIBRARIES = libcsmesh.a \
		   libcspartition.a

libcsmesh_a_SOURCES = \
cs_geom.cpp \
cs_join.cpp \
cs_join_intersect.cpp \
cs_join_merge.cpp \
cs_join_mesh.cpp \
cs_join_perio.cpp \
cs_join_post.cpp \
cs_join_set.cpp \
cs_join_split.cpp \
cs_join_update.cpp \
cs_join_util.cpp \
cs_mesh.cpp \
cs_mesh_adjacencies.cpp \
cs_mesh_bad_cells.cpp \
cs_mesh_boundary.cpp \
cs_mesh_boundary_layer.cpp \
cs_mesh_builder.cpp \
cs_mesh_cartesian.cpp \
cs_mesh_coarsen.cpp \
cs_mesh_coherency.cpp \
cs_mesh_connect.cpp \
cs_mesh_extrude.cpp \
cs_mesh_from_builder.cpp \
cs_mesh_group.cpp \
cs_mesh_halo.cpp \
cs_mesh_intersect.cpp \
cs_mesh_location.cpp \
cs_mesh_quality.cpp \
cs_mesh_quantities.cpp \
cs_mesh_refine.cpp \
cs_mesh_remove.cpp \
cs_mesh_save.cpp \
cs_mesh_smoother.cpp \
cs_mesh_to_builder.cpp \
cs_mesh_warping.cpp \
cs_stl.cpp \
cs_symmetry_faces_filter.cpp

# Partitioner (may require extra headers)

libcspartition_a_CPPFLAGS = $(AM_CPPFLAGS) \
$(METIS_CPPFLAGS) $(SCOTCH_CPPFLAGS)
libcspartition_a_SOURCES = cs_partition.cpp

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/mesh -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
