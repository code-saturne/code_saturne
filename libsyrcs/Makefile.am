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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = $(CPPFLAGS_PLE) $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS =

# Headers

noinst_HEADERS = \
syr_comm.h \
syr_coupling.h \
syr_defs.h

# Library source files

lib_LTLIBRARIES = libsyrcs.la
libsyrcs_la_LIBADD =

libsyrcs_la_SOURCES = \
syr_comm.c \
syr_coupling.c \
syr_cs.c \
syr_defs.c

libsyrcs_la_LDFLAGS = -no-undefined -version-info 0:0:0
