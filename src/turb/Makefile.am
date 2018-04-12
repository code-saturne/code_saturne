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
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/cdo \
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
cs_les_filter.h \
cs_les_inflow.h \
cs_turbulence_bc.h \
cs_turbulence_model.h

# Library source files

noinst_LTLIBRARIES = libcsturb.la

libcsturb_la_SOURCES = \
cs_les_filter.c \
cs_les_inflow.c \
cs_turbulence_bc.c \
cs_turbulence_model.c \
clipke.f90 \
clipsa.f90 \
clpalp.f90 \
clprij.f90 \
clpv2f.f90 \
divrit.f90 \
resalp.f90 \
reseps.f90 \
resrij.f90 \
resrij2.f90 \
resrit.f90 \
resssg.f90 \
resssg2.f90 \
resv2f.f90 \
rijech.f90 \
rijthe.f90 \
rijthe2.f90 \
rotcor.f90 \
tsepls.f90 \
turbke.f90 \
turbkw.f90 \
turbsa.f90 \
turrij.f90 \
visdyn.f90 \
vislmg.f90 \
vissma.f90 \
vissst.f90 \
visv2f.f90 \
viswal.f90 \
vor2cl.f90 \
vordep.f90 \
vorimp.f90 \
vorin0.f90 \
vorini.f90 \
vorlgv.f90 \
vorpre.f90 \
vortex.f90 \
vorver.f90 \
vorvit.f90
libcsturb_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
