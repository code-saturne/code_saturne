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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_at_plugin.h \
cs_at_opt_interp.h \
cs_at_data_assim.h

# Library source files

noinst_LTLIBRARIES = libcsatmo.la
libcsatmo_la_SOURCES = \
atini0.f90 \
atini1.f90 \
atiniv.f90 \
atleca.f90 \
atlecc.f90 \
atlecm.f90 \
atmcls.f90 \
atmsol.f90 \
atmstd.f90 \
atphyv.f90 \
atprke.f90 \
atprop.f90 \
atr1vf.f90 \
attssc.f90 \
attycl.f90 \
atvarp.f90 \
chem_luscheme_siream.f90 \
chem_luscheme1.f90 \
chem_luscheme2.f90 \
chem_luscheme3.f90 \
chem_solvelu.f90 \
chem_source_terms.f90 \
compute_gaseous_chemistry.f90 \
compute_siream.f90 \
cs_at_plugin.c \
cs_at_opt_interp.c \
cs_at_data_assim.c \
cs_at_source_terms.f90 \
etheq.f90  \
intprf.f90 \
kinrates.f90 \
nuclea.f90 \
rayigc.f90 \
rayir.f90 \
rayive.f90 \
rayso.f90 \
raysze.f90 \
roschem.f90 \
satfun.f90 \
solcat.f90 \
soliva.f90 \
solmoy.f90 \
solvar.f90

libcsatmo_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
