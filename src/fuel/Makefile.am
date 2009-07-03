#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main part

AM_FCFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/pprt \
-I$(top_srcdir)/include/cogz \
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/fuel \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/fuel/fuincl.h

# Library source files

noinst_LTLIBRARIES = libcsfuel.la
libcsfuel_la_SOURCES = \
fucym1.f90 \
fucyno.f90 \
fuflux.f90 \
fuini1.f90 \
fuiniv.f90 \
fulecd.f90 \
fupdfr.f90 \
fuphy1.f90 \
fuphy2.f90 \
fuphyv.f90 \
fuprop.f90 \
fuptcl.f90 \
futeh1.f90 \
futeh2.f90 \
futhp1.f90 \
futhp2.f90 \
futssc.f90 \
futsvc.f90 \
futsvi.f90 \
fuvarp.f90 \
fuveri.f90 \
memfu1.f90
libcsfuel_la_LDFLAGS = -no-undefined

libcsfuel_la_SOURCES += \
$(top_srcdir)/users/fuel/usfucl.f90 \
$(top_srcdir)/users/fuel/usfui1.f90 \
$(top_srcdir)/users/fuel/usfuiv.f90

fueldir = $(pkgdatadir)/users/fuel
fuel_DATA = \
$(top_srcdir)/users/fuel/usfucl.f90 \
$(top_srcdir)/users/fuel/usfui1.f90 \
$(top_srcdir)/users/fuel/usfuiv.f90
