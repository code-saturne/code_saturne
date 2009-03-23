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
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS = @LDFLAGS_DBG@ @LDFLAGS_OPT@

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/cogz/coincl.h

# Library source files

noinst_LTLIBRARIES = libcscogz.la
libcscogz_la_SOURCES = \
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
lwcgfu.f90 \
lwcini.f90 \
lwcphy.f90 \
lwctcl.f90 \
lwctss.f90 \
lwcurl.f90 \
lwcver.f90 \
pdflwc.f90 \
pdfpp3.f90 \
pdfpp4.f90
libcscogz_la_LDFLAGS = -no-undefined

libcscogz_la_SOURCES += \
$(top_srcdir)/users/cogz/usd3p1.f90 \
$(top_srcdir)/users/cogz/usd3pc.f90 \
$(top_srcdir)/users/cogz/usd3pi.f90 \
$(top_srcdir)/users/cogz/usebu1.f90 \
$(top_srcdir)/users/cogz/usebuc.f90 \
$(top_srcdir)/users/cogz/usebui.f90 \
$(top_srcdir)/users/cogz/uslwc1.f90 \
$(top_srcdir)/users/cogz/uslwcc.f90 \
$(top_srcdir)/users/cogz/uslwci.f90

cogzdir = $(pkgdatadir)/users/cogz
cogz_DATA = \
$(top_srcdir)/users/cogz/usd3p1.f90 \
$(top_srcdir)/users/cogz/usd3pc.f90 \
$(top_srcdir)/users/cogz/usd3pi.f90 \
$(top_srcdir)/users/cogz/usebu1.f90 \
$(top_srcdir)/users/cogz/usebuc.f90 \
$(top_srcdir)/users/cogz/usebui.f90 \
$(top_srcdir)/users/cogz/uslwc1.f90 \
$(top_srcdir)/users/cogz/uslwcc.f90 \
$(top_srcdir)/users/cogz/uslwci.f90
