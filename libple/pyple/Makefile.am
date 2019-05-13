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

if HAVE_SWIG
if HAVE_MPI4PY

AM_CPPFLAGS = $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(MPI_LDFLAGS)

# SWIG Library for python
BUILT_SOURCES = swig_wrap.c
SWIG_SOURCES = swig.i
#PYSWIG_SOURCES = pyple.py swig_wrap.c

pkgpython_PYTHON = __init__.py mpi_utilities.py
nodist_pkgpython_PYTHON = pyple.py
pkgpyexec_LTLIBRARIES = _pyple.la
_pyple_la_SOURCES = $(SWIG_SOURCES)
nodist__pyple_la_SOURCES = $(BUILT_SOURCES)
_pyple_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) $(SWIG_PYTHON_FLAGS) -I$(srcdir)/../src $(MPI_CPPFLAGS)
_pyple_la_LDFLAGS = -module
_pyple_la_LIBADD = ../src/libple.la

$(BUILT_SOURCES): $(SWIG_SOURCES)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I. $(AX_SWIG_PYTHON_CPPFLAGS) $(SWIG_PYTHON_FLAGS) -I$(srcdir)/../src -o $@ $<

endif
endif

# Distribution

EXTRA_DIST = \
swig.i
