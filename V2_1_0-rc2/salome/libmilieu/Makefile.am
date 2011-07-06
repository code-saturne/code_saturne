#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne CFD code.
#
#   Copyright (C) 2011  EDF S.A.
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

AM_CPPFLAGS = $(SALOME_KERNEL_CPPFLAGS)
AM_CFLAGS =
AM_LDFLAGS = $(SALOME_KERNEL_LDFLAGS)

# Public header files (to be installed)

include_HEADERS = \
communication.h \
donnees.h \
runmilieu.h \
utilitaires.h

# Library source files

lib_LTLIBRARIES = libmilieu.la
noinst_LTLIBRARIES =
libmilieu_la_LIBADD =
libmilieu_la_SOURCES = \
communication.c \
donnees.c \
runmilieu.c \
utilitaires.c

libmilieu_la_LDFLAGS = -no-undefined -version-info 0:0:0
