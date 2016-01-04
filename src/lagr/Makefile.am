## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
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
cs_lagr_clogging.h \
cs_lagr_roughness.h \
cs_lagr_dlvo.h \
cs_lagr_extract.h \
cs_lagr_perio.h \
cs_lagr_restart.h \
cs_lagr_tracking.h \
cs_lagr_utils.h

# Library source files

noinst_LTLIBRARIES = libcslagr.la
libcslagr_la_SOURCES = \
cs_lagr_clogging.c \
cs_lagr_roughness.c \
cs_lagr_dlvo.c \
cs_lagr_extract.c \
cs_lagr_perio.c \
cs_lagr_restart.c \
cs_lagr_tracking.c \
cs_lagr_utils.c \
diverv.f90 \
lagadh.f90 \
lagaff.f90 \
lagbar.f90 \
lagcar.f90 \
lagcli.f90 \
lagcou.f90 \
lagdcl.f90 \
lagdep.f90 \
lagdif.f90 \
lageje.f90 \
lagent.f90 \
lageqp.f90 \
lages1.f90 \
lages2.f90 \
lagesd.f90 \
lagesp.f90 \
laggeo.f90 \
laggra.f90 \
laghis.f90 \
lagich.f90 \
lagidp.f90 \
lagimp.f90 \
lagini.f90 \
lagipn.f90 \
lagitf.f90 \
lagitg.f90 \
lagitp.f90 \
laglec.f90 \
laglis.f90 \
lagnew.f90 \
lagopt.f90 \
lagout.f90 \
lagphy.f90 \
lagpoi.f90 \
lagprj.f90 \
lagres.f90 \
lagrus.f90 \
lagsec.f90 \
lagsta.f90 \
lagstf.f90 \
lagswe.f90 \
lagtmp.f90 \
lagune.f90
libcslagr_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod
