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

AM_CPPFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/lagr \
@FVM_CPPFLAGS@ @BFT_CPPFLAGS@ @MPI_CPPFLAGS@

AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@ @CFLAGS_EXT@

AM_FCFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/pprt \
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/lagr \
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/lagr/cs_lagr.h \
$(top_srcdir)/include/lagr/cs_lagr_perio.h \
$(top_srcdir)/include/lagr/lagdim.h \
$(top_srcdir)/include/lagr/lagpar.h \
$(top_srcdir)/include/lagr/lagran.h

# Library source files

noinst_LTLIBRARIES = libcslagr.la
libcslagr_la_SOURCES = \
cs_lagr.c \
cs_lagr_perio.c \
diverv.f90 \
enslag.f90 \
enswaf.f90 \
lagaff.f90 \
lagcar.f90 \
lagcel.f90 \
lagcou.f90 \
lagdeb.f90 \
lageli.f90 \
lagent.f90 \
lageqp.f90 \
lagerr.f90 \
lages1.f90 \
lages2.f90 \
lagesp.f90 \
lagfch.f90 \
laggra.f90 \
laghis.f90 \
lagich.f90 \
lagidp.f90 \
lagimp.f90 \
lagini.f90 \
lagipn.f90 \
lagitf.f90 \
lagitg.f90 \
lagitp.f90 \
laglec.f90 \
laglis.f90 \
lagnew.f90 \
lagnpr.f90 \
lagnwc.f90 \
lagopt.f90 \
lagout.f90 \
lagphy.f90 \
lagpoi.f90 \
lagrus.f90 \
lagsta.f90 \
lagstf.f90 \
lagtri.f90 \
lagune.f90 \
memla1.f90 \
memla2.f90 \
ouestu.f90
libcslagr_la_LDFLAGS = -no-undefined

libcslagr_la_SOURCES += \
$(top_srcdir)/users/lagr/uslabo.f90 \
$(top_srcdir)/users/lagr/usladp.f90 \
$(top_srcdir)/users/lagr/uslaed.f90 \
$(top_srcdir)/users/lagr/uslaen.f90 \
$(top_srcdir)/users/lagr/uslafe.f90 \
$(top_srcdir)/users/lagr/uslag1.f90 \
$(top_srcdir)/users/lagr/uslag2.f90 \
$(top_srcdir)/users/lagr/uslain.f90 \
$(top_srcdir)/users/lagr/uslapr.f90 \
$(top_srcdir)/users/lagr/uslaru.f90 \
$(top_srcdir)/users/lagr/uslast.f90 \
$(top_srcdir)/users/lagr/uslatc.f90 \
$(top_srcdir)/users/lagr/uslatp.f90

lagrdir = $(pkgdatadir)/users/lagr
lagr_DATA = \
$(top_srcdir)/users/lagr/uslabo.f90 \
$(top_srcdir)/users/lagr/usladp.f90 \
$(top_srcdir)/users/lagr/uslaed.f90 \
$(top_srcdir)/users/lagr/uslaen.f90 \
$(top_srcdir)/users/lagr/uslafe.f90 \
$(top_srcdir)/users/lagr/uslag1.f90 \
$(top_srcdir)/users/lagr/uslag2.f90 \
$(top_srcdir)/users/lagr/uslain.f90 \
$(top_srcdir)/users/lagr/uslapr.f90 \
$(top_srcdir)/users/lagr/uslaru.f90 \
$(top_srcdir)/users/lagr/uslast.f90 \
$(top_srcdir)/users/lagr/uslatc.f90 \
$(top_srcdir)/users/lagr/uslatp.f90
