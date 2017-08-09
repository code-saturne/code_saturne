## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(PETSC_CPPFLAGS) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

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
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

EXTRA_DIST =

# Library source files

examplefiles = \
cs_user_1d_wall_thermal.c \
cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c \
cs_user_cdo-condif.c \
cs_user_cdo-groundwater.c \
cs_user_cdo_extra_op-verif-diffusion.c \
cs_user_cdo_numerics-examples.c \
cs_user_coupling-saturne.c \
cs_user_coupling-syrthes.c \
cs_user_coupling.c \
cs_user_electric_scaling.c \
cs_user_extra_operations-balance_by_zone.c \
cs_user_extra_operations-scalar_balance.c \
cs_user_extra_operations-turbomachinery.c \
cs_user_head_losses.c \
cs_user_lagr_boundary_conditions.c \
cs_user_lagr_model.c \
cs_user_lagr_particle.c \
cs_user_mesh-boundary.c \
cs_user_mesh-groups_from_free_faces.c \
cs_user_mesh-input-save.c \
cs_user_mesh-joining.c \
cs_user_mesh-modify.c \
cs_user_mesh-periodicity.c \
cs_user_mesh-quality.c \
cs_user_parameters-base.c \
cs_user_parameters-ctwr.c \
cs_user_parameters-fans.c \
cs_user_parameters-linear_solvers.c \
cs_user_parameters-time_moments.c \
cs_user_performance_tuning-matrix.c \
cs_user_performance_tuning-numbering.c \
cs_user_performance_tuning-parallel-io.c \
cs_user_performance_tuning-partition.c \
cs_user_postprocess.c \
cs_user_postprocess-profiles.c \
cs_user_postprocess-sfc.c \
cs_user_radiative_transfer.c \
cs_user_radiative_transfer_bcs.c \
cs_user_solver-heat-equation.c \
cs_user_turbomachinery.c \
cs_user_zones.c \
cs_user_atmospheric_model.f90 \
cs_user_boundary_conditions-advanced.f90 \
cs_user_boundary_conditions-atmospheric.f90 \
cs_user_boundary_conditions-auto_inlet_profile.f90 \
cs_user_boundary_conditions-base.f90 \
cs_user_boundary_conditions-compressible.f90 \
cs_user_boundary_conditions-cooling_towers.f90 \
cs_user_boundary_conditions-electric_arcs.f90 \
cs_user_boundary_conditions-fuel.f90 \
cs_user_boundary_conditions-gas_3ptchem.f90 \
cs_user_boundary_conditions-gas_ebu.f90 \
cs_user_boundary_conditions-gas_libby_williams.f90 \
cs_user_boundary_conditions-mapped_inlet.f90 \
cs_user_boundary_conditions-pulverized_coal.f90 \
cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 \
cs_user_boundary_conditions-richards.f90 \
cs_user_boundary_conditions_ale-base.f90 \
cs_user_boundary_mass_source_terms-condensation.f90 \
cs_user_boundary_mass_source_terms-nzones_condensation.f90 \
cs_user_extra_operations-energy_balance.f90 \
cs_user_extra_operations-force_temperature.f90 \
cs_user_extra_operations-global_efforts.f90 \
cs_user_extra_operations-nusselt_calculation.f90 \
cs_user_extra_operations-parallel_operations.f90 \
cs_user_extra_operations-richards_flux.f90 \
cs_user_fluid_structure_interaction.f90 \
cs_user_initialization-atmospheric.f90 \
cs_user_initialization-base.f90 \
cs_user_initialization-compressible.f90 \
cs_user_initialization-cooling_towers.f90 \
cs_user_initialization-electric_arcs.f90 \
cs_user_initialization-fuel.f90 \
cs_user_initialization-gas_3ptchem.f90 \
cs_user_initialization-gas_ebu.f90 \
cs_user_initialization-gas_libby_williams.f90 \
cs_user_initialization-pulverized_coal.f90 \
cs_user_initialization-richards.f90 \
cs_user_initialization-time_step.f90 \
cs_user_initialization-unified_combustion_coal.f90 \
cs_user_les_inflow-base.f90 \
cs_user_mass_source_terms.f90 \
cs_user_metal_structures_source_terms-condensation.f90 \
cs_user_modules-user-arrays.f90 \
cs_user_parameters-cavitation.f90 \
cs_user_parameters-output.f90 \
cs_user_parameters-richards.f90 \
cs_user_parameters-rotation.f90 \
cs_user_parameters.f90 \
cs_user_physical_properties-coal_drift.f90 \
cs_user_physical_properties-compressible_flow.f90 \
cs_user_physical_properties-richards_sat.f90 \
cs_user_physical_properties-richards_unsat.f90 \
cs_user_physical_properties-scalar_drift.f90 \
cs_user_porosity.f90 \
cs_user_postprocess_var.f90 \
cs_user_source_terms-richards_decay.f90 \
cs_user_source_terms-scalar_in_a_channel.f90 \
cs_user_source_terms.f90 \
usatch_4spe5reac.f90 \
usthht.f90 \
usvort.f90 \
usvosy.f90

noinst_LTLIBRARIES = libcsexample.la
libcsexample_la_SOURCES = $(examplefiles)
libcsexample_la_LDFLAGS = -no-undefined

exampledir = $(pkgdatadir)/user_examples
example_DATA = $(examplefiles)

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

# Finish example files installation
# - tag with the package version number
# Done with two sed passes (does not work otherwise)

install-data-hook:
	cd $(DESTDIR)$(pkgdatadir)/user_examples && \
	f90_user_examples=`ls *.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!                      Code_Saturne version $(VERSION)\n"; \
	version2="\!                      --------------------------"; \
	for file in $${f90_user_examples}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done ;\
	c_user_examples=`ls *.c 2>/dev/null`; \
	stubvers="\/\* VERS \*\/"; \
	version="\/\* Code_Saturne version $(VERSION) \*\/"; \
	for file in $${c_user_examples}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version}/" < "$${file}" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done

