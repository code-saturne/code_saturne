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
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/fuel \
-I$(top_srcdir)/include/lagr \
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS = @LDFLAGS_DBG@ @LDFLAGS_OPT@

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/rayt/radiat.h

# Library source files

noinst_LTLIBRARIES = libcsrayt.la
libcsrayt_la_SOURCES = \
memra1.f90 \
memra2.f90 \
memra3.f90 \
raybrd.f90 \
raycli.f90 \
raydak.f90 \
raydir.f90 \
raydom.f90 \
rayens.f90 \
raylec.f90 \
rayopt.f90 \
rayout.f90 \
raypar.f90 \
raypun.f90 \
raysca.f90 \
raysol.f90 \
rmodak.f90
libcsrayt_la_LDFLAGS = -no-undefined

libcsrayt_la_SOURCES += \
$(top_srcdir)/users/rayt/usray1.f90 \
$(top_srcdir)/users/rayt/usray2.f90 \
$(top_srcdir)/users/rayt/usray3.f90 \
$(top_srcdir)/users/rayt/usray4.f90 \
$(top_srcdir)/users/rayt/usray5.f90

raytdir = $(pkgdatadir)/users/rayt
rayt_DATA = \
$(top_srcdir)/users/rayt/usray1.f90 \
$(top_srcdir)/users/rayt/usray2.f90 \
$(top_srcdir)/users/rayt/usray3.f90 \
$(top_srcdir)/users/rayt/usray4.f90 \
$(top_srcdir)/users/rayt/usray5.f90
