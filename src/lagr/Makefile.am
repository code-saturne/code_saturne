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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/user \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_lagr.h \
cs_lagr_clogging.h \
cs_lagr_agglo.h \
cs_lagr_fragmentation.h \
cs_lagr_head_losses.h \
cs_lagr_roughness.h \
cs_lagr_dlvo.h \
cs_lagr_event.h \
cs_lagr_extract.h \
cs_lagr_injection.h \
cs_lagr_geom.h \
cs_lagr_post.h \
cs_lagr_restart.h \
cs_lagr_query.h \
cs_lagr_tracking.h \
cs_lagr_stat.h \
cs_lagr_new.h \
cs_lagr_precipitation_model.h \
cs_lagr_deposition_model.h \
cs_lagr_car.h \
cs_lagr_coupling.h \
cs_lagr_adh.h \
cs_lagr_resuspension.h \
cs_lagr_poisson.h \
cs_lagr_options.h \
cs_lagr_orientation.h \
cs_lagr_aux_mean_fluid_quantities.h \
cs_lagr_particle.h \
cs_lagr_print.h \
cs_lagr_sde.h \
cs_lagr_sde_model.h \
cs_lagr_lec.h \
cs_lagr_log.h \
cs_lagr_prototypes.h \
cs_lagr_headers.h

# Library source files

noinst_LTLIBRARIES = libcslagr.la
libcslagr_la_SOURCES = \
cs_lagr_clogging.c \
cs_lagr_agglo.c \
cs_lagr_fragmentation.c \
cs_lagr_roughness.c \
cs_lagr.c \
cs_lagr_dlvo.c \
cs_lagr_event.c \
cs_lagr_extract.c \
cs_lagr_injection.c \
cs_lagr_post.c \
cs_lagr_restart.c \
cs_lagr_query.c \
cs_lagr_tracking.c \
cs_lagr_stat.c \
cs_lagr_adh.c \
cs_lagr_print.c \
cs_lagr_car.c \
cs_lagr_coupling.c \
cs_lagr_deposition_model.c \
cs_lagr_geom.c \
cs_lagr_poisson.c \
cs_lagr_aux_mean_fluid_quantities.c \
cs_lagr_head_losses.c \
cs_lagr_new.c \
cs_lagr_options.c \
cs_lagr_orientation.c \
cs_lagr_particle.c \
cs_lagr_resuspension.c \
cs_lagr_precipitation_model.c \
cs_lagr_sde.c \
cs_lagr_sde_model.c \
cs_lagr_lec.c \
cs_lagr_log.c

libcslagr_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod
