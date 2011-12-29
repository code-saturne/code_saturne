

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS =

# Library source files

noinst_LTLIBRARIES = libcselec.la
libcselec_la_SOURCES = \
elflux.f90 \
elini1.f90 \
eliniv.f90 \
ellecd.f90 \
elphyv.f90 \
elprop.f90 \
elthht.f90 \
eltssc.f90 \
elvarp.f90 \
elveri.f90
libcselec_la_LDFLAGS = -no-undefined

elecusers = \
$(top_srcdir)/users/elec/uselen.f90 \
$(top_srcdir)/users/elec/uselph.f90 \
$(top_srcdir)/users/elec/uselrc.f90

libcselec_la_SOURCES += $(elecusers)

elecdir = $(pkgdatadir)/users/elec
elec_DATA = $(elecusers)

clean-local:
	-rm -f *__genmod.f90 *.mod
