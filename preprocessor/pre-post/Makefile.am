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

AM_CPPFLAGS = -I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
@HDF5_CPPFLAGS@ @CCM_CPPFLAGS@ @MED_CPPFLAGS@ @CGNS_CPPFLAGS@
AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@
AM_LDFLAGS = \
@HDF5_LDFLAGS@ @ADF_LDFLAGS@ @CCM_LDFLAGS@ @MED_LDFLAGS@ @CGNS_LDFLAGS@

# Library source files

noinst_LTLIBRARIES = libprepostecs.la

libprepostecs_la_SOURCES = \
ecs_cgns.c \
ecs_champ_comm.c \
ecs_champ_post.c \
ecs_champ_post_cgns.c \
ecs_champ_post_ens.c \
ecs_champ_post_med.c \
ecs_maillage_ncs.c \
ecs_maillage_post.c \
ecs_maillage_pre.c \
ecs_med.c \
ecs_post.c \
ecs_post_cgns.c \
ecs_post_ens.c \
ecs_post_med.c \
ecs_pre.c \
ecs_pre_ccm.c \
ecs_pre_cgns.c \
ecs_pre_comet.c \
ecs_pre_ens.c \
ecs_pre_gambit.c \
ecs_pre_gmsh.c \
ecs_pre_ideas.c \
ecs_pre_med.c \
ecs_pre_nopo.c

libprepostecs_la_LDFLAGS = -no-undefined

