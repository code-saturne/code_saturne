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

AM_CPPFLAGS = -I$(top_srcdir)/preprocessor/pre-post -I$(top_srcdir)/preprocessor/util
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS =

# Headers

noinst_HEADERS = \
ecs_table_att.h \
ecs_table_def.h \
ecs_table.h \
ecs_table_priv.h \
ecs_descr_chaine.h \
ecs_descr.h \
ecs_descr_priv.h \
ecs_famille_chaine.h \
ecs_famille.h \
ecs_famille_priv.h \
ecs_maillage.h \
ecs_maillage_priv.h

# Library source files

noinst_LTLIBRARIES = libbaseecs.la

libbaseecs_la_SOURCES = \
ecs_table.c \
ecs_table_att.c \
ecs_table_def.c \
ecs_descr.c \
ecs_descr_chaine.c \
ecs_famille.c \
ecs_famille_chaine.c \
ecs_maillage.c

libbaseecs_la_LDFLAGS = -no-undefined

