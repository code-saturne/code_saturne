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
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCMODINCLUDE)$(top_builddir)/src/turb \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_benchmark.h \
cs_blas.h \
cs_gradient.h \
cs_gradient_perio.h \
cs_grid.h \
cs_matrix.h \
cs_matrix_priv.h \
cs_matrix_util.h \
cs_multigrid.h \
cs_sles.h

# Library source files

noinst_LTLIBRARIES = libcsalge.la

libcsalge_la_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_HOT)
libcsalge_la_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_HOT)
libcsalge_la_SOURCES = \
cs_benchmark.c \
cs_blas.c \
cs_gradient.c \
cs_gradient_perio.c \
cs_grid.c \
cs_matrix.c \
cs_matrix_util.c \
cs_multigrid.c \
cs_sles.c \
bilsca.f90 \
bilsct.f90 \
bilsc2.f90 \
bilscv.f90 \
bilsc4.f90 \
diften.f90 \
diftnv.f90 \
divmas.f90 \
divrij.f90 \
grdcel.f90 \
grdpre.f90 \
grdvec.f90 \
grdvni.f90 \
grdpot.f90 \
inimas.f90 \
inimav.f90 \
invers.f90 \
itrgrp.f90 \
itrgrv.f90 \
itrmas.f90 \
itrmav.f90 \
matrdt.f90 \
matrix.f90 \
matrvv.f90 \
matrxv.f90 \
prodsc.f90 \
projts.f90 \
projtv.f90 \
vectds.f90 \
viscfa.f90 \
visort.f90 \
vitens.f90 \
vistnv.f90 \
viortv.f90
libcsalge_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
