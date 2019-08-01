#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2019  EDF S.A.
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

if HAVE_MPI4PY

AM_CPPFLAGS = $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(MPI_LDFLAGS)

MPI4PY_CFLAGS = -I${shell ${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'}

pkgpython_PYTHON = __init__.py Coupling.py Init.py
pkgpyexec_LTLIBRARIES = libpyplecoupling.la
libpyplecoupling_la_SOURCES = ple_coupling_wrapper.c
libpyplecoupling_la_CFLAGS = $(MPI4PY_CFLAGS) $(PYTHON_CPPFLAGS) -I$(srcdir)/../src $(MPI_CPPFLAGS)
libpyplecoupling_la_CPPFLAGS = $(MPI4PY_CFLAGS) $(PYTHON_CPPFLAGS) -I$(srcdir)/../src $(MPI_CPPFLAGS)
libpyplecoupling_la_LDFLAGS = $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS)
libpyplecoupling_la_LIBADD = $(top_builddir)/src/libple.la


endif
