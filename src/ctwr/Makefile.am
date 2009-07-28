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
-I$(top_srcdir)/include/ctwr \
@BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@

AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@

AM_FCFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/pprt \
-I$(top_srcdir)/include/ctwr \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/ctwr/cs_ctwr_air_props.h \
$(top_srcdir)/include/ctwr/cs_ctwr.h \
$(top_srcdir)/include/ctwr/cs_ctwr_halo.h \
$(top_srcdir)/include/ctwr/cs_ctwr_mesh.h \
$(top_srcdir)/include/ctwr/ctincl.h

# Library source files

noinst_LTLIBRARIES = libcsctwr.la
libcsctwr_la_SOURCES = \
cs_ctwr.c \
cs_ctwr_air_props.c \
cs_ctwr_f2c.f90 \
cs_ctwr_halo.c \
cs_ctwr_mesh.c \
ctini1.f90 \
ctiniv.f90 \
ctphyv.f90 \
cttssc.f90 \
ctvarp.f90
libcsctwr_la_LDFLAGS = -no-undefined

libcsctwr_la_SOURCES += \
$(top_srcdir)/users/ctwr/usctcl.f90 \
$(top_srcdir)/users/ctwr/usctdz.f90 \
$(top_srcdir)/users/ctwr/uscti1.f90 \
$(top_srcdir)/users/ctwr/usctiv.f90

ctwrdir = $(pkgdatadir)/users/ctwr
ctwr_DATA = \
$(top_srcdir)/users/ctwr/usctcl.f90 \
$(top_srcdir)/users/ctwr/usctdz.f90 \
$(top_srcdir)/users/ctwr/uscti1.f90 \
$(top_srcdir)/users/ctwr/usctiv.f90
