## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

AM_CPPFLAGS =
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS =

# Headers

noinst_HEADERS = \
ecs_backtrace.h \
ecs_comm.h \
ecs_def.h \
ecs_elt_typ_liste.h \
ecs_fic.h \
ecs_file.h \
ecs_mem.h \
ecs_mem_usage.h \
ecs_tab.h \
ecs_timer.h

# Library source files

noinst_LTLIBRARIES = libutilecs.la

libutilecs_la_SOURCES = \
ecs_backtrace.c \
ecs_comm.c \
ecs_def.c \
ecs_fic.c \
ecs_file.c \
ecs_mem.c \
ecs_mem_usage.c \
ecs_tab.c \
ecs_timer.c

libutilecs_la_LDFLAGS = -no-undefined

