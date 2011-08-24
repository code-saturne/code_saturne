#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2011 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
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
cs_mesh_coherency.h \
cs_mesh_connect.h \
cs_mesh.h \
cs_mesh_halo.h \
cs_mesh_quality.h \
cs_mesh_quantities.h \
cs_mesh_save.h \
cs_mesh_thinwall.h \
cs_mesh_warping.h \
cs_mesh_smoother.h

# Library source files

noinst_LTLIBRARIES = libcsmesh.la

libcsmesh_la_SOURCES = \
cs_join.c \
cs_join_intersect.c \
cs_join_merge.c \
cs_join_mesh.c \
cs_join_perio.c \
cs_join_post.c \
cs_join_set.c \
cs_join_split.c \
cs_join_update.c \
cs_join_util.c \
cs_mesh.c \
cs_mesh_coherency.c \
cs_mesh_connect.c \
cs_mesh_halo.c \
cs_mesh_quality.c \
cs_mesh_quantities.c \
cs_mesh_save.c \
cs_mesh_thinwall.c \
cs_mesh_warping.c \
cs_mesh_smoother.c

libcsmesh_la_LDFLAGS = -no-undefined
