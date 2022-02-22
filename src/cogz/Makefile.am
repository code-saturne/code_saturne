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

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

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
cs_combustion_gas.h \
cs_cogz_headers.h

# Library source files

noinst_LTLIBRARIES = libcscogz.la
libcscogz_la_SOURCES = \
cs_combustion_gas.c \
coini1.f90 \
colecd.f90 \
coprop.f90 \
cothht.f90 \
covarp.f90 \
d3phst.f90 \
d3pini.f90 \
d3pint.f90 \
d3pphy.f90 \
d3ptcl.f90 \
d3pver.f90 \
ebuini.f90 \
ebuphy.f90 \
ebutcl.f90 \
ebutss.f90 \
ebuver.f90 \
gauss.f90 \
lwcgfu.f90 \
lwcini.f90 \
lwcphy.f90 \
lwctcl.f90 \
lwctss.f90 \
lwcurl.f90 \
lwcver.f90 \
pdflwc.f90 \
pdfpp3.f90 \
pdfpp4.f90 \
sootsc.f90 \
cs_steady_laminar_flamelet_init.f90 \
cs_steady_laminar_flamelet_bcond.f90 \
cs_steady_laminar_flamelet_physical_prop.f90 \
cs_steady_laminar_flamelet_source_terms.f90 \
cs_steady_laminar_flamelet_verify.f90 \
cs_steady_laminar_flamelet_read_base.f90

libcscogz_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod
