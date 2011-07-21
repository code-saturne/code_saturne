#-------------------------------------------------------------------------------
#   This file is part of the "Base Functions and Types" library, intended to
#   simplify and enhance portability, memory and I/O use for scientific codes.
#
#   Copyright (C) 2004-2011  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
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
bft_file.h \
bft_mem.h \
bft_mem_usage.h \
bft_printf.h \
bft_sys_info.h \
bft_timer.h

# Library source files

lib_LTLIBRARIES = libbft.la
libbft_la_SOURCES   = \
bft_backtrace.c \
bft_error.c \
bft_file.c \
bft_mem.c \
bft_mem_usage.c \
bft_printf.c \
bft_sys_info.c \
bft_timer.c

libbft_la_LDFLAGS = -no-undefined -version-info 0:0:0
