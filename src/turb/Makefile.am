## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
-I$(top_srcdir)/src/atmo \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/pprt \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/alge \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
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
cs_turbulence_rotation.h \
cs_les_balance.h \
cs_les_filter.h \
cs_les_inflow.h \
cs_les_mu_t.h \
cs_turbulence_bc.h \
cs_turbulence_inflow.h \
cs_turbulence_ke.h \
cs_turbulence_kw.h \
cs_turbulence_ml.h \
cs_turbulence_model.h \
cs_turbulence_rij.h \
cs_turbulence_sa.h \
cs_turbulence_v2f.h \
cs_turbulence_headers.h

# Library source files

noinst_LTLIBRARIES = libcsturb.la

libcsturb_la_SOURCES = \
cs_les_balance.c \
cs_les_filter.c \
cs_les_inflow.c \
cs_les_mu_t.c \
cs_turbulence_bc.c \
cs_turbulence_inflow.c \
cs_turbulence_ke.c \
cs_turbulence_kw.c \
cs_turbulence_ml.c \
cs_turbulence_model.c \
cs_turbulence_rij.c \
cs_turbulence_rotation.c \
cs_turbulence_sa.c \
cs_turbulence_v2f.c \
clprij.f90 \
divrit.f90 \
reseps.f90 \
resrij.f90 \
resrij2.f90 \
resrit.f90 \
resssg2.f90 \
rijech.f90 \
rijthe.f90 \
rijthe2.f90 \
turrij.f90
libcsturb_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
