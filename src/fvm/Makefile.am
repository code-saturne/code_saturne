#-------------------------------------------------------------------------------
#   This file is part of the "Finite Volume Mesh" library, intended to provide
#   finite volume mesh and associated fields I/O and manipulation services.
#
#   Copyright (C) 2004-2010  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = -I$(top_srcdir)/src/bft @MPI_CPPFLAGS@
AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@
AM_LDFLAGS = @CGNS_LDFLAGS@ @MED_LDFLAGS@ @HDF5_LDFLAGS@ @MPI_LDFLAGS@

# Conditionally compiled extensions

libfvm_cgns_la_CPPFLAGS = -I$(top_srcdir)/src/bft @CGNS_CPPFLAGS@ @MPI_CPPFLAGS@
libfvm_med_la_CPPFLAGS = -I$(top_srcdir)/src/bft @HDF5_CPPFLAGS@ @MED_CPPFLAGS@ @MPI_CPPFLAGS@

# Public header files (to be installed)

include_HEADERS = \
fvm_block_to_part.h \
fvm_defs.h \
fvm_file.h \
fvm_group.h \
fvm_interface.h \
fvm_io_num.h \
fvm_neighborhood.h \
fvm_nodal.h \
fvm_nodal_append.h \
fvm_nodal_extrude.h \
fvm_nodal_from_desc.h \
fvm_nodal_order.h \
fvm_nodal_project.h \
fvm_nodal_triangulate.h \
fvm_nodal_extract.h \
fvm_order.h \
fvm_parall.h \
fvm_part_to_block.h \
fvm_periodicity.h \
fvm_point_location.h \
fvm_selector.h \
fvm_trace.h \
fvm_triangulate.h \
fvm_writer.h

noinst_HEADERS = \
fvm_block_dist.h \
fvm_box.h \
fvm_box_priv.h \
fvm_box_tree.h \
fvm_config_defs.h \
fvm_convert_array.h \
fvm_gather.h \
fvm_morton.h \
fvm_nodal_priv.h \
fvm_selector_postfix.h \
fvm_tesselation.h \
fvm_to_cgns.h \
fvm_to_med.h \
fvm_to_ensight.h \
fvm_to_ensight_case.h \
fvm_to_ensight_v1.h \
fvm_to_text.h \
fvm_writer_helper.h \
fvm_writer_priv.h

# Library source files

lib_LTLIBRARIES = libfvm.la
noinst_LTLIBRARIES =
libfvm_la_LIBADD =

libfvm_la_SOURCES = \
fvm_block_dist.c \
fvm_block_to_part.c \
fvm_box.c \
fvm_box_tree.c \
fvm_convert_array.c \
fvm_defs.c \
fvm_file.c \
fvm_gather.c \
fvm_group.c \
fvm_interface.c \
fvm_io_num.c \
fvm_morton.c \
fvm_neighborhood.c \
fvm_nodal.c \
fvm_nodal_append.c \
fvm_nodal_extrude.c \
fvm_nodal_from_desc.c \
fvm_nodal_order.c \
fvm_nodal_project.c \
fvm_nodal_triangulate.c \
fvm_nodal_extract.c \
fvm_order.c \
fvm_parall.c \
fvm_part_to_block.c \
fvm_periodicity.c \
fvm_point_location.c \
fvm_selector.c \
fvm_selector_postfix.c \
fvm_tesselation.c \
fvm_to_ensight.c \
fvm_to_ensight_case.c \
fvm_to_ensight_v1.c \
fvm_to_text.c \
fvm_trace.c \
fvm_triangulate.c \
fvm_writer.c \
fvm_writer_helper.c

if HAVE_CGNS
noinst_LTLIBRARIES += libfvm_cgns.la
libfvm_la_LIBADD += libfvm_cgns.la
libfvm_cgns_la_SOURCES = fvm_to_cgns.c
endif

if HAVE_MED
noinst_LTLIBRARIES += libfvm_med.la
libfvm_la_LIBADD += libfvm_med.la
libfvm_med_la_SOURCES = fvm_to_med.c
endif

libfvm_la_LDFLAGS = -no-undefined -version-info 0:0:0
