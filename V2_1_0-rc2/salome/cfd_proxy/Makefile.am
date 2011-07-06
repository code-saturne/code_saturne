#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne CFD code.
#
#   Copyright (C) 2006-2011  EDF
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

AM_CPPFLAGS = $(MPI_CPPFLAGS) $(SALOME_KERNEL_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(MPI_LDFLAGS) $(SALOME_KERNEL_LDFLAGS)

# Headers

include_HEADERS = \
cfd_proxy_api.h

noinst_HEADERS = \
cfd_proxy_defs.h \
cfd_proxy_child.h \
cfd_proxy_comm.h \
cfd_proxy_forward.h \
cfd_proxy_loader.h

# Library source files

lib_LTLIBRARIES = libcfdproxy.la
libcfdproxy_la_LIBADD =

libcfdproxy_la_SOURCES = \
cfd_proxy_defs.c \
cfd_proxy_child.c \
cfd_proxy_comm.c  \
cfd_proxy_forward.c  \
cfd_proxy_loader.c  \
cfd_proxy_api.c

libcfdproxy_la_LDFLAGS = -no-undefined -version-info 0:0:0
