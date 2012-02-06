## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/cplv \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Add user_modules file to distribution as it is not automatically accounted
# for by automake

EXTRA_DIST = user_modules.f90

# Library source files

userfiles = \
cs_user_boundary_conditions.f90 \
cs_user_coupling.c \
cs_user_initialization.f90 \
cs_user_mesh.c \
cs_user_postprocess.c \
cs_user_solver.c \
cs_user_source_terms.f90 \
usalcl.f90 \
usaste.f90 \
uscfpv.f90 \
uscfth.f90 \
usctdz.f90 \
uselen.f90 \
uselph.f90 \
uselrc.f90 \
ushist.f90 \
usini1.f90 \
uskpdc.f90 \
uslabo.f90 \
usladp.f90 \
uslaed.f90 \
uslaen.f90 \
uslafe.f90 \
uslag1.f90 \
uslag2.f90 \
uslain.f90 \
uslapr.f90 \
uslaru.f90 \
uslast.f90 \
uslatc.f90 \
uslatp.f90 \
usphyv.f90 \
usproj.f90 \
uspt1d.f90 \
usray1.f90 \
usray2.f90 \
usray3.f90 \
usray4.f90 \
usray5.f90 \
ussmag.f90 \
usstru.f90 \
usthht.f90 \
ustmgr.f90 \
ustsma.f90 \
usvima.f90 \
usvist.f90 \
usvort.f90 \
usvosy.f90 \
usvpst.f90

noinst_LTLIBRARIES = libcsuser.la
libcsuser_la_SOURCES = $(userfiles)
libcsuser_la_LDFLAGS = -no-undefined

userdir = $(pkgdatadir)/user
user_DATA = $(userfiles) user_modules.f90

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
