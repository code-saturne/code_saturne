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
-I$(top_srcdir)/include/atmo \
-I$(top_srcdir)/include/cfbl \
-I$(top_srcdir)/include/cogz \
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/ctwr \
-I$(top_srcdir)/include/elec \
-I$(top_srcdir)/include/fuel \
-I$(top_srcdir)/include/lagr \
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/pprt/ppcpfu.h \
$(top_srcdir)/include/pprt/ppincl.h \
$(top_srcdir)/include/pprt/ppppar.h \
$(top_srcdir)/include/pprt/ppthch.h

# Library source files

noinst_LTLIBRARIES = libcspprt.la
libcspprt_la_SOURCES = \
memppt.f90 \
ppcabs.f90 \
ppclim.f90 \
ppcsca.f90 \
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

libcspprt_la_SOURCES += \
$(top_srcdir)/users/pprt/usppmo.f90

pprtdir = $(pkgdatadir)/users/pprt
pprt_DATA = \
$(top_srcdir)/users/pprt/usppmo.f90
