## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = 

# Library source files

noinst_LTLIBRARIES = libcsatmo.la
libcsatmo_la_SOURCES = \
atini1.f90 \
atiniv.f90 \
atlecm.f90 \
atmcls.f90 \
atmsol.f90 \
atmstd.f90 \
atphyv.f90 \
atprke.f90 \
atprop.f90 \
attycl.f90 \
atr1vf.f90 \
attssc.f90 \
atvarp.f90 \
etheq.f90  \
intprf.f90 \
nuclea.f90 \
satfun.f90 \
rayigc.f90 \
rayir.f90 \
rayive.f90 \
rayso.f90 \
raysze.f90 \
solcat.f90 \
soliva.f90 \
solmoy.f90 \
solvar.f90
libcsatmo_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
