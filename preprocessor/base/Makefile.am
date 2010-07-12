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
AM_CFLAGS =
AM_LDFLAGS =

# Library source files

noinst_LTLIBRARIES = libbaseecs.la

libbaseecs_la_SOURCES = \
ecs_champ.c \
ecs_champ_att.c \
ecs_champ_def.c \
ecs_descr.c \
ecs_descr_chaine.c \
ecs_famille.c \
ecs_famille_chaine.c \
ecs_maillage.c

libbaseecs_la_LDFLAGS = -no-undefined

