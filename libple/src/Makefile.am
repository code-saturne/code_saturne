#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2020  EDF S.A.
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

AM_CPPFLAGS = $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(MPI_LDFLAGS)

# Public header files (to be installed)

pleincludedir = $(includedir)
pleinclude_HEADERS = \
ple_coupling.h \
ple_defs.h \
ple_locator.h

noinst_HEADERS = \
ple_config_defs.h

# Library source files

lib_LTLIBRARIES = libple.la
noinst_LTLIBRARIES =
libple_la_LIBADD = $(MPI_LIBS) -lm
libple_la_SOURCES = ple_defs.c ple_locator.c

if HAVE_MPI
libple_la_SOURCES += ple_coupling.c
endif

libple_la_LDFLAGS = $(MPI_LDFLAGS) -no-undefined -version-info 2:0:0
