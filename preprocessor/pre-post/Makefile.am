## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

AM_CPPFLAGS = -I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
$(HDF5_CPPFLAGS) $(HDF5_CPPFLAGS_MPI) \
$(MED_CPPFLAGS) $(CGNS_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = \
$(MED_LDFLAGS) $(CGNS_LDFLAGS) $(HDF5_LDFLAGS) $(HDF5_LDFLAGS_MPI)

# Headers

noinst_HEADERS = \
ecs_table_comm.h \
ecs_table_post_cgns.h \
ecs_table_post_ens.h \
ecs_table_post.h \
ecs_table_post_med.h \
ecs_maillage_ncs.h \
ecs_maillage_post.h \
ecs_maillage_pre.h \
ecs_med.h \
ecs_med_priv.h \
ecs_post_cgns.h \
ecs_post_cgns_priv.h \
ecs_post_ens.h \
ecs_post_ens_priv.h \
ecs_post.h \
ecs_post_med.h \
ecs_pre_ccm.h \
ecs_pre_cgns.h \
ecs_pre_ens.h \
ecs_pre_gambit.h \
ecs_pre_gmsh.h \
ecs_pre.h \
ecs_pre_ideas.h \
ecs_pre_med.h \
ecs_pre_nopo.h

# Library source files

lib_LTLIBRARIES =
pkglib_LTLIBRARIES = 
noinst_LTLIBRARIES = libprepostecs.la
libprepostecs_la_LIBADD =

libprepostecs_la_CPPFLAGS = \
-DPKGLIBDIR=\"'$(pkglibdir)'\" \
$(AM_CPPFLAGS)

libprepostecs_la_SOURCES = \
ecs_table_comm.c \
ecs_table_post.c \
ecs_table_post_cgns.c \
ecs_table_post_ens.c \
ecs_table_post_med.c \
ecs_maillage_ncs.c \
ecs_maillage_post.c \
ecs_maillage_pre.c \
ecs_med.c \
ecs_post.c \
ecs_post_cgns.c \
ecs_post_ens.c \
ecs_post_med.c \
ecs_pre.c \
ecs_pre_cgns.c \
ecs_pre_ens.c \
ecs_pre_gambit.c \
ecs_pre_gmsh.c \
ecs_pre_ideas.c \
ecs_pre_med.c \
ecs_pre_nopo.c

# CCMIO (plugin or linked)

if HAVE_CCM

if HAVE_PLUGINS

ecs_ccm_la_CPPFLAGS = \
-I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
$(CCM_CPPFLAGS) $(CCM_CXXFLAGS) $(MPI_CPPFLAGS)
ecs_ccm_la_LIBADD = $(CCM_LDFLAGS) $(CCM_LIBS) $(CCMRUNPATH)
ecs_ccm_la_SOURCES = ecs_pre_ccm.c
pkglib_LTLIBRARIES += ecs_ccm.la
ecs_ccm_la_LDFLAGS = -module -avoid-version

else

libecs_ccm_la_CPPFLAGS = \
-I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
$(CCM_CPPFLAGS) $(CCM_CXXFLAGS) $(MPI_CPPFLAGS)
libecs_ccm_la_LIBADD = $(CCM_LDFLAGS) $(CCM_LIBS) $(CCMRUNPATH)
libecs_ccm_la_SOURCES = ecs_pre_ccm.c
noinst_LTLIBRARIES += libecs_ccm.la
libprepostecs_la_LIBADD += libecs_ccm.la

endif

endif

libprepostecs_la_LDFLAGS = -no-undefined

