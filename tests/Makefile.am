# Process this file with automake to produce Makefile.in

#------------------------------------------------------------------------------
#   This file is part of the "Mathematical Expression Interpreter" library.
#
#   Copyright (C) 2008-2009  EDF
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
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = -I$(top_srcdir)/src/mei @BFT_CPPFLAGS@
AM_CFLAGS   =
AM_LDFLAGS  = @BFT_LDFLAGS@
LDADD = $(top_builddir)/src/mei/libmei.la \
        @BFT_LIBS@ -lm

check_PROGRAMS = mei_test

mei_test_SOURCES = mei_test_main.c mei_test_graph.c

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)

