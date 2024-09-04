## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
-I$(top_srcdir)/src/atmo \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/comb \
-I$(top_srcdir)/src/user \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_lagr.h \
cs_lagr_adh.h \
cs_lagr_agglo.h \
cs_lagr_aux_mean_fluid_quantities.h \
cs_lagr_car.h \
cs_lagr_clogging.h \
cs_lagr_coupling.h \
cs_lagr_deposition_model.h \
cs_lagr_dlvo.h \
cs_lagr_event.h \
cs_lagr_extract.h \
cs_lagr_fragmentation.h \
cs_lagr_geom.h \
cs_lagr_head_losses.h \
cs_lagr_injection.h \
cs_lagr_lec.h \
cs_lagr_log.h \
cs_lagr_new.h \
cs_lagr_options.h \
cs_lagr_orientation.h \
cs_lagr_particle.h \
cs_lagr_poisson.h \
cs_lagr_porosity.h \
cs_lagr_post.h \
cs_lagr_precipitation_model.h \
cs_lagr_print.h \
cs_lagr_query.h \
cs_lagr_restart.h \
cs_lagr_resuspension.h \
cs_lagr_roughness.h \
cs_lagr_sde.h \
cs_lagr_sde_model.h \
cs_lagr_stat.h \
cs_lagr_tracking.h \
cs_lagr_prototypes.h \
cs_lagr_headers.h

# Library source files

noinst_LIBRARIES = libcslagr.a
libcslagr_a_SOURCES = \
cs_lagr.cpp \
cs_lagr_adh.cpp \
cs_lagr_agglo.cpp \
cs_lagr_aux_mean_fluid_quantities.cpp \
cs_lagr_car.cpp \
cs_lagr_clogging.cpp \
cs_lagr_coupling.cpp \
cs_lagr_deposition_model.cpp \
cs_lagr_dlvo.cpp \
cs_lagr_event.cpp \
cs_lagr_extract.cpp \
cs_lagr_fragmentation.cpp \
cs_lagr_geom.cpp \
cs_lagr_head_losses.cpp \
cs_lagr_injection.cpp \
cs_lagr_lec.cpp \
cs_lagr_log.cpp \
cs_lagr_new.cpp \
cs_lagr_options.cpp \
cs_lagr_orientation.cpp \
cs_lagr_particle.cpp \
cs_lagr_poisson.cpp \
cs_lagr_porosity.cpp \
cs_lagr_post.cpp \
cs_lagr_precipitation_model.cpp \
cs_lagr_print.cpp \
cs_lagr_query.cpp \
cs_lagr_restart.cpp \
cs_lagr_resuspension.cpp \
cs_lagr_roughness.cpp \
cs_lagr_sde.cpp \
cs_lagr_sde_model.cpp \
cs_lagr_stat.cpp \
cs_lagr_tracking.cpp

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

# Additional rules

clean-local:
	-rm -rf .deps
