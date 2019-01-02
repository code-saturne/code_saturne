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

# Relocatable installation

if HAVE_RELOCATABLE
ORIGINRUNPATH = -R '\$$ORIGIN/../../lib'
endif

if HAVE_PLUGINS
  PLUGINRUNPATH = $(PLUGINRUNPATH_R)$(PLUGINRUNPATH_N)
endif

# CCM

if HAVE_CCM
if HAVE_PLUGINS
else
  LDADD_CCM = $(CCM_LDFLAGS) $(CCMRUNPATH) $(CCM_LIBS)
endif
endif

# Main part

AM_CPPFLAGS = -DLOCALEDIR=\"'$(localedir)'\" \
-I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/pre-post -I$(top_srcdir)/preprocessor/util

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

# Global LDFLAGS must be used so as to handle the possibility of a CXX link.
AM_LDFLAGS = $(ORIGINRUNPATH) $(PLUGINRUNPATH)

# Main program

pkglibexec_PROGRAMS = cs_preprocess

if HAVE_LINK_CXX
# Dummy C++ source to cause C++ linking.
nodist_EXTRA_cs_preprocess_SOURCES = dummy.cxx
cs_preprocess_LINK = $(CXXLINK)
else
cs_preprocess_LINK = $(LINK)
endif

cs_preprocess_CPPFLAGS = \
-DLOCALEDIR=\"'$(localedir)'\" \
$(AM_CPPFLAGS) \
$(HDF5_CPPFLAGS) $(MED_CPPFLAGS) $(CGNS_CPPFLAGS)

cs_preprocess_SOURCES = \
ecs_cmd.c \
ecs_main.c

noinst_HEADERS = \
ecs_cmd.h

cs_preprocess_LDADD = \
$(top_builddir)/preprocessor/pre-post/libprepostecs.la \
$(top_builddir)/preprocessor/base/libbaseecs.la \
$(top_builddir)/preprocessor/util/libutilecs.la $(LTLIBINTL) \
$(LDADD_CCM) \
$(MED_LDFLAGS) $(MEDRUNPATH) $(MED_LIBS) \
$(CGNS_LDFLAGS) $(CGNSRUNPATH) $(CGNS_LIBS) \
$(HDF5_LDFLAGS) $(HDF5RUNPATH) $(HDF5_LIBS) $(HDF5_MPI_LIBS) \
-lm

