#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Preprocessor.
#
#   Copyright (C) 2005-2010  EDF
#
#   The Code_Saturne Preprocessor is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Preprocessor is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = -I$(top_srcdir)/prepro/base -I$(top_srcdir)/prepro/pre-post -I$(top_srcdir)/prepro/util
AM_CFLAGS =
AM_LDFLAGS =

# Main program

bin_PROGRAMS = cs_preprocess

cs_preprocess_CPPFLAGS = $(AM_CPPFLAGS) \
@HDF5_CPPFLAGS@ @MED_CPPFLAGS@ @CGNS_CPPFLAGS@

cs_preprocess_SOURCES = \
ecs_cmd.c \
ecs_main.c

cs_preprocess_LDADD = \
$(top_builddir)/prepro/pre-post/libprepostecs.la \
$(top_builddir)/prepro/base/libbaseecs.la \
$(top_builddir)/prepro/util/libutilecs.la $(LTLIBINTL) \
@HDF5_LDFLAGS@ @HDF5_LIBS@ @ADF_LDFLAGS@ @ADF_LIBS@ \
@CCM_LDFLAGS@ @CCM_LIBS@ @MED_LDFLAGS@ @MED_LIBS@ \
@CGNS_LDFLAGS@ @CGNS_LIBS@ \
-lm

