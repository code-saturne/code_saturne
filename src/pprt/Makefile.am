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

pkginclude_HEADERS =

# Library source files

noinst_LTLIBRARIES = libcspprt.la
libcspprt_la_SOURCES = \
ppcabs.f90 \
ppini1.f90 \
ppinii.f90 \
ppiniv.f90 \
ppinv2.f90 \
pplecd.f90 \
pppdfr.f90 \
ppphyv.f90 \
ppprcl.f90 \
ppprop.f90 \
ppray4.f90 \
pptbht.f90 \
pptssc.f90 \
pptycl.f90 \
ppvarp.f90
libcspprt_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod
