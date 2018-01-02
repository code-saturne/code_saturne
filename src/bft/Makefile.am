## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/base \
$(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
bft_backtrace.h \
bft_error.h \
bft_mem.h \
bft_mem_usage.h \
bft_printf.h

# Library source files

noinst_LTLIBRARIES = libbft.la
libbft_la_SOURCES  = \
bft_backtrace.c \
bft_error.c \
bft_mem.c \
bft_mem_usage.c \
bft_printf.c

libbft_la_LIBADD = $(LTLIBINTL)
libbft_la_LDFLAGS = -no-undefined
