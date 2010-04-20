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

AM_CPPFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/ctwr \
-I$(top_srcdir)/include/lagr \
@FVM_CPPFLAGS@ @MEI_CPPFLAGS@ @BFT_CPPFLAGS@ @MPI_CPPFLAGS@ @LIBXML2_CPPFLAGS@ @BLAS_CPPFLAGS@

AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@

AM_LDFLAGS =

# Library source files

lib_LTLIBRARIES = libsaturne.la
libsaturne_la_SOURCES = cs_solver.c
libsaturne_la_LIBADD = \
$(top_builddir)/src/base/libcshot.la \
$(top_builddir)/src/base/libcsbase.la \
$(top_builddir)/src/base/libcscore.la \
$(top_builddir)/src/pprt/libcspprt.la \
$(top_builddir)/src/atmo/libcsatmo.la \
$(top_builddir)/src/cfbl/libcscfbl.la \
$(top_builddir)/src/cogz/libcscogz.la \
$(top_builddir)/src/cplv/libcscplv.la \
$(top_builddir)/src/ctwr/libcsctwr.la \
$(top_builddir)/src/elec/libcselec.la \
$(top_builddir)/src/fuel/libcsfuel.la \
$(top_builddir)/src/lagr/libcslagr.la \
$(top_builddir)/src/mati/libcsmati.la \
$(top_builddir)/src/rayt/libcsrayt.la
libsaturne_la_LDFLAGS = -no-undefined -version 0:0:0

# Main program

bin_PROGRAMS = cs_solver

# Code_Saturne IO utility program

bin_PROGRAMS += cs_io_dump

# MPI partitioning program

if HAVE_PARTITIONER
bin_PROGRAMS += cs_partition
endif

cs_solver_SOURCES =
cs_solver_LDADD = libsaturne.la $(LTLIBINTL) \
@FVM_LDFLAGS@ @FVM_LIBS@ \
@MEI_LDFLAGS@ @MEI_LIBS@ \
@BFT_LDFLAGS@ @BFT_LIBS@ \
@MPI_LDFLAGS@ @MPI_LIBS@ \
@LIBXML2_LDFLAGS@ @LIBXML2_LIBS@ \
@BLAS_LDFLAGS@ @BLAS_LIBS@ \
@FCLIBS@

# Code_Saturne IO utility (minimal dependencies)

cs_io_dump_SOURCES = cs_io_dump.c

# Conditionally compiled executables

if HAVE_PARTITIONER

cs_partition_CPPFLAGS = $(AM_CPPFLAGS) \
@BFT_CPPFLAGS@ @METIS_CPPFLAGS@ @SCOTCH_CPPFLAGS@
cs_partition_SOURCES = cs_partition.c
cs_partition_LDADD = $(top_builddir)/src/base/libcscore.la \
@FVM_LDFLAGS@ @FVM_LIBS@ \
@BFT_LDFLAGS@ @BFT_LIBS@ \
@METIS_LDFLAGS@ @METIS_LIBS@ \
@SCOTCH_LDFLAGS@ @SCOTCH_LIBS@ \
@MPI_LDFLAGS@ @MPI_LIBS@

endif
