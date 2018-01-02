## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT)
AM_LDFLAGS = $(CCM_LDFLAGS) $(CGNS_LDFLAGS) $(MED_LDFLAGS) \
$(HDF5_LDFLAGS) $(MPI_LDFLAGS)

# Conditionally compiled extensions

libfvm_cgns_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(CGNS_CPPFLAGS) $(MPI_CPPFLAGS)
libfvm_med_la_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(HDF5_CPPFLAGS) $(MED_CPPFLAGS) $(MPI_CPPFLAGS)

# Public header files (to be installed)

pkginclude_HEADERS = \
fvm_defs.h \
fvm_group.h \
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
fvm_periodicity.h \
fvm_point_location.h \
fvm_selector.h \
fvm_trace.h \
fvm_triangulate.h \
fvm_writer.h

noinst_HEADERS = \
fvm_box.h \
fvm_box_priv.h \
fvm_box_tree.h \
fvm_convert_array.h \
fvm_hilbert.h \
fvm_morton.h \
fvm_nodal_priv.h \
fvm_selector_postfix.h \
fvm_tesselation.h \
fvm_to_ccm.h \
fvm_to_cgns.h \
fvm_to_med.h \
fvm_to_catalyst.h \
fvm_to_ensight.h \
fvm_to_ensight_case.h \
fvm_to_medcoupling.h \
fvm_to_plot.h \
fvm_to_time_plot.h \
fvm_writer_helper.h \
fvm_writer_priv.h

# Library source files

lib_LTLIBRARIES =
pkglib_LTLIBRARIES = 
noinst_LTLIBRARIES = libfvm.la libfvm_filters.la
libfvm_la_LIBADD = $(LTLIBINTL)
libfvm_filters_la_LIBADD =

libfvm_la_SOURCES = \
fvm_box.c \
fvm_box_tree.c \
fvm_convert_array.c \
fvm_defs.c \
fvm_group.c \
fvm_hilbert.c \
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
fvm_periodicity.c \
fvm_point_location.c \
fvm_selector.c \
fvm_selector_postfix.c \
fvm_tesselation.c \
fvm_trace.c \
fvm_triangulate.c

libfvm_filters_la_SOURCES = \
fvm_to_ensight.c \
fvm_to_ensight_case.c \
fvm_to_plot.c \
fvm_to_time_plot.c \
fvm_writer.c \
fvm_writer_helper.c

if HAVE_CCM
noinst_LTLIBRARIES += libfvm_ccm.la
libfvm_ccm_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/mesh \
$(CFLAGS_DBG) $(CFLAGS_OPT) \
$(CCM_CPPFLAGS) $(MPI_CPPFLAGS)
libfvm_filters_la_LIBADD += libfvm_ccm.la
libfvm_ccm_la_SOURCES = fvm_to_ccm.c
endif

if HAVE_CGNS
noinst_LTLIBRARIES += libfvm_cgns.la
libfvm_filters_la_LIBADD += libfvm_cgns.la
libfvm_cgns_la_SOURCES = fvm_to_cgns.c
endif

if HAVE_MED
noinst_LTLIBRARIES += libfvm_med.la
libfvm_filters_la_LIBADD += libfvm_med.la
libfvm_med_la_SOURCES = fvm_to_med.c
endif

if HAVE_MEDCOUPLING

fvm_medcoupling_la_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) \
$(MEDCOUPLING_CPPFLAGS) $(MPI_CPPFLAGS)
pkglib_LTLIBRARIES += fvm_medcoupling.la

fvm_medcoupling_la_SOURCES = fvm_to_medcoupling.cxx
fvm_medcoupling_la_LDFLAGS = -module -avoid-version
fvm_medcoupling_la_LIBADD = \
$(MEDCOUPLING_LDFLAGS) $(MEDCOUPLING_LIBS) \
$(MEDCOUPLINGRUNPATH)

endif

# Catalyst (plugin or linked)

if HAVE_CATALYST

if HAVE_PLUGIN_CATALYST

fvm_catalyst_la_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(CATALYST_CPPFLAGS) $(MPI_CPPFLAGS)
fvm_catalyst_la_CXXFLAGS = \
$(CXXFLAGS_DBG) $(CXXFLAGS_OPT) \
$(CATALYST_CXXFLAGS)
fvm_catalyst_la_LIBADD = \
$(CATALYST_LDFLAGS) $(CATALYST_LIBS) \
$(CATALYSTRUNPATH)
fvm_catalyst_la_SOURCES = fvm_to_catalyst.cxx
pkglib_LTLIBRARIES += fvm_catalyst.la
fvm_catalyst_la_LDFLAGS = -module -avoid-version

else

libfvm_catalyst_la_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
$(CATALYST_CPPFLAGS) $(MPI_CPPFLAGS)
libfvm_catalyst_la_CXXFLAGS = \
$(CXXFLAGS_DBG) $(CXXFLAGS_OPT) \
$(CATALYST_CXXFLAGS)
libfvm_catalyst_la_LIBADD = \
$(CATALYST_LDFLAGS) $(CATALYST_LIBS) \
$(CATALYSTRUNPATH)
libfvm_catalyst_la_SOURCES = fvm_to_catalyst.cxx
noinst_LTLIBRARIES += libfvm_catalyst.la
libfvm_filters_la_LIBADD += libfvm_catalyst.la

endif

endif

libfvm_la_LDFLAGS = -no-undefined
libfvm_filters_la_LDFLAGS = -no-undefined
