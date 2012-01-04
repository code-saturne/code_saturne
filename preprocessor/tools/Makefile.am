## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

AM_CPPFLAGS = -I$(top_srcdir)/preprocessor/util
AM_CFLAGS =
AM_LDFLAGS =

# Main program

pkglibexec_PROGRAMS =

if HAVE_MED
pkglibexec_PROGRAMS += igghexa_to_med
endif

# Conditionally compiled executables

if HAVE_MED

igghexa_to_med_CPPFLAGS = $(AM_CPPFLAGS) \
$(HDF5_CPPFLAGS) $(MED_CPPFLAGS)
igghexa_to_med_SOURCES = igghexa_to_med.c
igghexa_to_med_LDADD = \
$(top_builddir)/preprocessor/pre-post/libprepostecs.la \
$(top_builddir)/preprocessor/base/libbaseecs.la \
$(top_builddir)/preprocessor/util/libutilecs.la $(LTLIBINTL) \
$(HDF5_LDFLAGS) $(HDF5_LIBS) $(MED_LDFLAGS) $(MED_LIBS) \
$(CGNS_LDFLAGS) $(CGNS_LIBS) \
-lm

endif

