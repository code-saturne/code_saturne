## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_benchmark.h \
cs_blas.h \
cs_convection_diffusion.h \
cs_face_viscosity.h \
cs_divergence.h \
cs_gradient.h \
cs_gradient_perio.h \
cs_grid.h \
cs_matrix.h \
cs_matrix_building.h \
cs_matrix_default.h \
cs_matrix_priv.h \
cs_matrix_tuning.h \
cs_matrix_util.h \
cs_multigrid.h \
cs_sles.h \
cs_sles_default.h \
cs_sles_it.h

# Library source files

noinst_LTLIBRARIES = libcsalge.la

libcsalge_la_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_HOT)
libcsalge_la_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_HOT)
libcsalge_la_SOURCES = \
cs_benchmark.c \
cs_blas.c \
cs_bw_time_diff.c \
cs_convection_diffusion.c \
cs_face_viscosity.c \
cs_divergence.c \
cs_gradient.c \
cs_gradient_perio.c \
cs_grid.c \
cs_matrix.c \
cs_matrix_building.c \
cs_matrix_default.c \
cs_matrix_tuning.c \
cs_matrix_util.c \
cs_multigrid.c \
cs_sles.c \
cs_sles_default.c \
cs_sles_it.c \
bilsca.f90 \
bilscv.f90 \
grdcel.f90 \
grdpre.f90 \
grdvec.f90 \
grdpot.f90 \
prodsc.f90 \
visort.f90
libcsalge_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
