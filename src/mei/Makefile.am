# Process this file with automake to produce Makefile.in

#------------------------------------------------------------------------------
#   This file is part of the "Mathematical Expression Interpreter" library.
#
#   Copyright (C) 2008-2011  EDF
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
#------------------------------------------------------------------------------

AM_CPPFLAGS = -I$(top_srcdir)/src/bft
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS =

AM_LFLAGS = -l
AM_YFLAGS = -d

BUILT_SOURCES = mei_parser.h

pkginclude_HEADERS = \
mei_evaluate.h \
mei_hash_table.h \
mei_node.h

noinst_HEADERS = \
mei_math_util.h \
mei_parser_glob.h

lib_LTLIBRARIES = libmei.la
libmei_la_SOURCES = \
mei_evaluate.c \
mei_hash_table.c \
mei_math_util.c \
mei_node.c \
mei_parser_glob.c \
mei_parser.y \
mei_scanner.l

libmei_la_LDFLAGS = -no-undefined -version-info 0:0:0
