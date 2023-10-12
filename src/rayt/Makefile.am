## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/atmo \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/comb \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/pprt \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/alge \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCMODINCLUDE)$(top_builddir)/src/turb \
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT) $(FCFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_rad_transfer_absorption.h \
cs_rad_transfer_dir.h \
cs_rad_transfer_fields.h \
cs_rad_transfer_options.h \
cs_rad_transfer_wall_flux.h \
cs_rad_transfer_bcs.h \
cs_rad_transfer_fsck.h \
cs_rad_transfer_adf_models.h \
cs_rad_transfer_pun.h \
cs_rad_transfer_solve.h \
cs_rad_transfer_modak.h \
cs_rad_transfer_source_terms.h \
cs_rad_transfer_restart.h \
cs_rad_transfer.h \
cs_rad_headers.h

# Library source files

noinst_LIBRARIES = libcsrayt.a
libcsrayt_a_SOURCES = \
cs_rad_transfer_absorption.c \
cs_rad_transfer.c \
cs_rad_transfer_dir.c \
cs_rad_transfer_fields.c \
cs_rad_transfer_options.c \
cs_rad_transfer_bcs.c \
cs_rad_transfer_wall_flux.c \
cs_rad_transfer_fsck.c \
cs_rad_transfer_adf_models.c \
cs_rad_transfer_pun.c \
cs_rad_transfer_solve.c \
cs_rad_transfer_modak.c \
cs_rad_transfer_source_terms.c \
cs_rad_transfer_restart.c

clean-local:
	-rm -f *__genmod.f90 *.mod
