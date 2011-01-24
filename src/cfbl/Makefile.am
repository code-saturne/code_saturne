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

noinst_LTLIBRARIES = libcscfbl.la
libcscfbl_la_SOURCES = \
cfbsc2.f90 \
cfbsc3.f90 \
cfcdts.f90 \
cfdivs.f90 \
cfdttv.f90 \
cfener.f90 \
cfini1.f90 \
cfiniv.f90 \
cfmsfl.f90 \
cfmsgs.f90 \
cfmsvl.f90 \
cfmsvs.f90 \
cfphyv.f90 \
cfprop.f90 \
cfqdmv.f90 \
cfrusb.f90 \
cfvarp.f90 \
cfxtcl.f90 \
memcfe.f90 \
memcfm.f90 \
memcft.f90 \
memcfv.f90
libcscfbl_la_LDFLAGS = -no-undefined

cfblusers = \
$(top_srcdir)/users/cfbl/uscfcl.f90 \
$(top_srcdir)/users/cfbl/uscfpv.f90 \
$(top_srcdir)/users/cfbl/uscfth.f90 \
$(top_srcdir)/users/cfbl/uscfx1.f90 \
$(top_srcdir)/users/cfbl/uscfx2.f90 \
$(top_srcdir)/users/cfbl/uscfxi.f90

libcscfbl_la_SOURCES += $(cfblusers)

cfbldir = $(pkgdatadir)/users/cfbl
cfbl_DATA = $(cfblusers)

clean-local:
	-rm -f *__genmod.f90 *.mod
