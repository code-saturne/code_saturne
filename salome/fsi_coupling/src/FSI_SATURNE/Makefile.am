#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2011 EDF S.A., France
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

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = FSI

AM_CFLAGS = $(SALOME_INCLUDES) -fexceptions

salomelib_LTLIBRARIES = libFSI_SATURNEEngine.la
salomeinclude_HEADERS = FSI_SATURNE.hxx

libFSI_SATURNEEngine_la_SOURCES = FSI_SATURNE.cxx 
nodist_libFSI_SATURNEEngine_la_SOURCES =
libFSI_SATURNEEngine_la_CXXFLAGS = -I$(top_builddir)/salome/fsi_coupling/idl  $(SALOME_INCLUDES) -I$(top_srcdir)/salome/cfd_proxy
libFSI_SATURNEEngine_la_LIBADD = $(top_builddir)/salome/cfd_proxy/libcfdproxy.la $(top_builddir)/salome/fsi_coupling/idl/libSalomeIDLFSI.la $(SALOME_LIBS)
libFSI_SATURNEEngine_la_LDFLAGS = 

