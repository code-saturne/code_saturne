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
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/atmo \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cfbl \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/darc \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCFLAGS_DBG) $(FCFLAGS_OPT) $(FCFLAGS_SHARED)

AM_LDFLAGS =

# Add cs_user_modules file to distribution as it is not automatically accounted
# for by automake

EXTRA_DIST = cs_user_modules.f90 \
cs_user_boundary_conditions.f90 cs_user_extra_operations.f90 \
cs_user_initialization.f90  \
cs_user_parameters.f90 \
cs_user_physical_properties.f90 \
cs_user_source_terms.f90 \
cs_user_porosity.f90

# Library source files

userfiles = \
cs_user_1d_wall_thermal.c \
cs_user_atmospheric_model.f90 \
cs_user_boundary_conditions.c \
cs_user_wall_condensation.c \
cs_user_coupling.c \
cs_user_electric_scaling.c \
cs_user_extra_operations.c \
cs_user_fluid_structure_interaction.f90 \
cs_user_head_losses.c \
cs_user_hgn.c \
cs_user_initialization.c \
cs_user_lagr_model.c \
cs_user_lagr_particle.c \
cs_user_lagr_boundary_conditions.c \
cs_user_lagr_volume_conditions.c \
cs_user_les_inflow.c \
cs_user_mesh.c \
cs_user_metal_structures_source_terms.f90 \
cs_user_paramedmem_coupling.c \
cs_user_parameters.c \
cs_user_performance_tuning.c \
cs_user_physical_properties.c \
cs_user_postprocess.c \
cs_user_porosity.c \
cs_user_radiative_transfer.c \
cs_user_radiative_transfer_bcs.c \
cs_user_solver.c \
cs_user_source_terms.c \
cs_user_turbomachinery.c \
cs_user_zones.c \
cs_user_boundary_conditions_ale.f90 \
usatch.f90

# Special treatment for files existing both in f90 and C:
cs_user_f_boundary_conditions = $(top_srcdir)/src/user/cs_user_boundary_conditions.f90

cs_f_user_boundary_conditions.$(OBJEXT): $(cs_user_f_boundary_conditions)
	$(FCCOMPILE) -c -o cs_f_user_boundary_conditions.$(OBJEXT)  $(cs_user_f_boundary_conditions)

cs_user_f_parameters = $(top_srcdir)/src/user/cs_user_parameters.f90

cs_f_user_parameters.$(OBJEXT): $(cs_user_f_parameters)
	$(FCCOMPILE) -c -o cs_f_user_parameters.$(OBJEXT)  $(cs_user_f_parameters)

cs_user_f_extra_operations = $(top_srcdir)/src/user/cs_user_extra_operations.f90

cs_f_user_extra_operations.$(OBJEXT): $(cs_user_f_extra_operations)
	$(FCCOMPILE) -c -o cs_f_user_extra_operations.$(OBJEXT)  $(cs_user_f_extra_operations)

cs_user_f_initialization = $(top_srcdir)/src/user/cs_user_initialization.f90

cs_f_user_initialization.$(OBJEXT): $(cs_user_f_initialization)
	$(FCCOMPILE) -c -o cs_f_user_initialization.$(OBJEXT)  $(cs_user_f_initialization)

cs_user_f_physical_properties = $(top_srcdir)/src/user/cs_user_physical_properties.f90

cs_f_user_physical_properties.$(OBJEXT): $(cs_user_f_physical_properties)
	$(FCCOMPILE) -c -o cs_f_user_physical_properties.$(OBJEXT)  $(cs_user_f_physical_properties)

cs_user_f_source_terms = $(top_srcdir)/src/user/cs_user_source_terms.f90

cs_f_user_source_terms.$(OBJEXT): $(cs_user_f_source_terms)
	$(FCCOMPILE) -c -o cs_f_user_source_terms.$(OBJEXT)  $(cs_user_f_source_terms)

cs_user_f_porosity = $(top_srcdir)/src/user/cs_user_porosity.f90

cs_f_user_porosity.$(OBJEXT): $(cs_user_f_porosity)
	$(FCCOMPILE) -c -o cs_f_user_porosity.$(OBJEXT)  $(cs_user_f_porosity)

noinst_LIBRARIES = libcsuser.a
libcsuser_a_SOURCES = $(userfiles)
libcsuser_a_LIBADD = cs_f_user_parameters.$(OBJEXT) cs_f_user_extra_operations.$(OBJEXT) \
cs_f_user_initialization.$(OBJEXT) cs_f_user_physical_properties.$(OBJEXT) \
cs_f_user_source_terms.$(OBJEXT) cs_f_user_boundary_conditions.$(OBJEXT) \
cs_f_user_porosity.$(OBJEXT)

userdir = $(pkgdatadir)/user_sources/REFERENCE
user_DATA = $(userfiles) cs_user_modules.f90 \
cs_user_boundary_conditions.f90 \
cs_user_parameters.f90 \
cs_user_extra_operations.f90 \
cs_user_initialization.f90 \
cs_user_physical_properties.f90 \
cs_user_source_terms.f90 \
cs_user_porosity.f90

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

# Finish user files installation
# - tag with the package version number
# - remove the test preventing from using a source file
# Done with two sed passes (does not work otherwise)

install-data-hook:
	c_user_files=`ls ${top_srcdir}/src/user/*.c 2>/dev/null`; \
	for f in $${c_user_files}; do \
	  file=`basename "$${f}"` && \
	  $(PYTHON) "${top_srcdir}/build-aux/cs_user_src_filter.py" "$(cs_version_short)" < "${top_srcdir}/src/user/$${file}" > "$(DESTDIR)$(userdir)/$${file}"; \
	done ; \
	cd $(DESTDIR)$(userdir) && \
	f90_user_files=`ls *.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!                      code_saturne version $(cs_version_short)\n"; \
	version2="\!                      ------------------------"; \
	comment1="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START"; \
	comment2="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END"; \
	for file in $${f90_user_files}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" | \
	  $(SED) -e "/$${comment1}/,/$${comment2}/d" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done
