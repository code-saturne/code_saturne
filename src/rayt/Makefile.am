#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2010 EDF S.A., France
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
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS =

# Library source files

noinst_LTLIBRARIES = libcsrayt.la
libcsrayt_la_SOURCES = \
memra1.f90 \
memra3.f90 \
raycli.f90 \
raydak.f90 \
raydir.f90 \
raydom.f90 \
raylec.f90 \
rayopt.f90 \
rayout.f90 \
raypar.f90 \
raypun.f90 \
raysca.f90 \
raysol.f90 \
rmodak.f90
libcsrayt_la_LDFLAGS = -no-undefined

raytusers = \
$(top_srcdir)/users/rayt/usray1.f90 \
$(top_srcdir)/users/rayt/usray2.f90 \
$(top_srcdir)/users/rayt/usray3.f90 \
$(top_srcdir)/users/rayt/usray4.f90 \
$(top_srcdir)/users/rayt/usray5.f90

libcsrayt_la_SOURCES += $(raytusers)

raytdir = $(pkgdatadir)/users/rayt
rayt_DATA = $(raytusers)

clean-local:
	-rm -f *__genmod.f90 *.mod
