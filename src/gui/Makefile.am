#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2010 EDF S.A., France
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
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/base \
$(MPI_CPPFLAGS) $(LIBXML2_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
cs_gui.h \
cs_gui_boundary_conditions.h \
cs_gui_conjugate_heat_transfer.h \
cs_gui_mesh.h \
cs_gui_mobile_mesh.h \
cs_gui_particles.h \
cs_gui_radiative_transfer.h \
cs_gui_specific_physics.h \
cs_gui_util.h \
cs_gui_variables.h

# Library source files

noinst_LTLIBRARIES = libcsgui.la

###

libcsgui_la_SOURCES = \
cs_gui.c \
cs_gui_boundary_conditions.c \
cs_gui_conjugate_heat_transfer.c \
cs_gui_mesh.c \
cs_gui_mobile_mesh.c \
cs_gui_particles.c \
cs_gui_radiative_transfer.c \
cs_gui_specific_physics.c \
cs_gui_util.c
