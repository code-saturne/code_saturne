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
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS =

# Library source files

noinst_LTLIBRARIES = libcscplv.la
libcscplv_la_SOURCES = \
cpcym2.f90 \
cpflux.f90 \
cpini1.f90 \
cpiniv.f90 \
cplecd.f90 \
cplin1.f90 \
cplini.f90 \
cplph1.f90 \
cplphy.f90 \
cplpro.f90 \
cpltcl.f90 \
cpltss.f90 \
cpltsv.f90 \
cplvar.f90 \
cplver.f90 \
cplym1.f90 \
cppdf4.f90 \
cppdfr.f90 \
cpphy1.f90 \
cpphy2.f90 \
cpphyv.f90 \
cpprop.f90 \
cpptcl.f90 \
cprays.f90 \
cpteh1.f90 \
cpteh2.f90 \
cptehm.f90 \
cpthp1.f90 \
cpthp2.f90 \
cptssc.f90 \
cptsvc.f90 \
cptsvi.f90 \
cpvarp.f90 \
cpveri.f90 \
memcp1.f90
libcscplv_la_LDFLAGS = -no-undefined

cplvusers = \
$(top_srcdir)/users/cplv/uscpcl.f90 \
$(top_srcdir)/users/cplv/uscpi1.f90 \
$(top_srcdir)/users/cplv/uscpiv.f90 \
$(top_srcdir)/users/cplv/uscpl1.f90 \
$(top_srcdir)/users/cplv/uscplc.f90

libcscplv_la_SOURCES += $(cplvusers)

cplvdir = $(pkgdatadir)/users/cplv
cplv_DATA = $(cplvusers)

clean-local:
	-rm -f *__genmod.f90 *.mod
