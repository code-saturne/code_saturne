#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2004-2018  EDF S.A.
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

AM_CPPFLAGS = -I$(top_srcdir)/src $(MPI_CPPFLAGS)
AM_CFLAGS   =
AM_LDFLAGS  = $(MPI_LDFLAGS)
LDADD = $(top_builddir)/src/libple.la $(MPI_LIBS)

check_PROGRAMS =

if HAVE_MPI
check_PROGRAMS += ple_coupling_test
ple_coupling_test_SOURCES = ple_coupling_test.c
ple_coupling_test_CPPFLAGS = -I$(top_srcdir)/src $(MPI_CPPFLAGS)
ple_coupling_test_LDFLAGS  = $(MPI_LDFLAGS)
ple_coupling_test_LDADD = $(top_builddir)/src/libple.la $(MPI_LIBS) $(INTLLIBS) -lm
endif

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)

