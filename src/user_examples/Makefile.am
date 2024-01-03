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
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cfbl \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/comb \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/meg \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/gwf \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(PETSC_CPPFLAGS) \
$(HYPRE_CPPFLAGS) \
$(MUMPS_CPPFLAGS) \
$(MEDCOUPLING_CPPFLAGS) \
$(MED_CPPFLAGS) \
$(HDF5_CPPFLAGS) \
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

EXTRA_DIST = \
cs_cad_intersect.h \
cs_cad_intersect.cxx \
cs_user_porosity-from_cad.c

# Library source files

examplefiles = \
cs_user_1d_wall_thermal.c \
cs_user_wall_condensation.c \
cs_user_boundary_conditions.c \
cs_user_boundary_conditions-advanced.c \
cs_user_boundary_conditions_ale-base.c \
cs_user_boundary_conditions-atmospheric.c \
cs_user_boundary_conditions-auto_inlet_profile.c \
cs_user_boundary_conditions-base.c \
cs_user_boundary_conditions-compressible.c \
cs_user_boundary_conditions-cooling_towers.c \
cs_user_boundary_conditions-electric_arcs.c \
cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c \
cs_user_boundary_conditions-mapped_inlet.c \
cs_user_boundary_conditions-medcoupling_remapper_2d.c \
cs_user_boundary_conditions-notebook.c \
cs_user_boundary_conditions-pulverized_coal.c \
cs_user_boundary_conditions-time_table.c \
cs_user_coupling-saturne.c \
cs_user_coupling-syrthes.c \
cs_user_electric_scaling.c \
cs_user_extra_operations-balance_by_zone.c \
cs_user_extra_operations-boundary_forces.c \
cs_user_extra_operations-force_temperature.c \
cs_user_extra_operations-mean_profiles.c \
cs_user_extra_operations-medcoupling_slice.c \
cs_user_extra_operations-nusselt_and_friction_calculation.c \
cs_user_extra_operations-nusselt_calculation.c \
cs_user_extra_operations-parallel_operations.c \
cs_user_extra_operations-scalar_balance.c \
cs_user_extra_operations-stopping_criterion.c \
cs_user_extra_operations-turbomachinery.c \
cs_user_extra_operations-verif_cdo_diffusion.c \
cs_user_extra_operations-vorticity_field.c \
cs_user_fluid_structure_interaction.c \
cs_user_fluid_structure_interaction-code_aster.c \
cs_user_head_losses.c \
cs_user_hgn-relax_time_scale.c \
cs_user_ibm.c \
cs_user_ibm_modify.c \
cs_user_initialization-atmospheric.c \
cs_user_initialization-base.c \
cs_user_initialization-compressible.c \
cs_user_initialization-cooling_towers.c \
cs_user_initialization-electric_arcs.c \
cs_user_initialization-gas_3ptchem.c \
cs_user_initialization-medcoupling_remapper_3d.c \
cs_user_initialization-pulverized_coal.c \
cs_user_initialization-time_step.c \
cs_user_lagr_boundary_conditions.c \
cs_user_lagr_volume_conditions.c \
cs_user_lagr_model.c \
cs_user_lagr_particle.c \
cs_user_lagr_particle-coal.c \
cs_user_les_inflow-base.c\
cs_user_mesh-boundary.c \
cs_user_mesh-cartesian.c \
cs_user_mesh-groups_from_free_faces.c \
cs_user_mesh-input-save.c \
cs_user_mesh-joining.c \
cs_user_mesh-modify.c \
cs_user_mesh-periodicity.c \
cs_user_mesh-quality.c \
cs_user_paramedmem_coupling-base.c \
cs_user_parameters-base.c \
cs_user_parameters-cdo-condif.c \
cs_user_parameters-cdo-gwf.c \
cs_user_parameters-cdo-gwf_user_soil.c \
cs_user_parameters-cdo-laplacian.c \
cs_user_parameters-cdo-solidification.c \
cs_user_parameters-ctwr.c \
cs_user_parameters-coupling.c \
cs_user_parameters-fans.c \
cs_user_parameters-cdo-linear_solvers.c \
cs_user_parameters-linear_solvers.c \
cs_user_parameters-rotation.c \
cs_user_parameters-volume_mass_injection.c \
cs_user_parameters-time_moments.c \
cs_user_performance_tuning-matrix.c \
cs_user_performance_tuning-numbering.c \
cs_user_performance_tuning-parallel-io.c \
cs_user_performance_tuning-partition.c \
cs_user_physical_properties.c \
cs_user_physical_properties-compressible_flow.c \
cs_user_physical_properties-coal_drift.c \
cs_user_physical_properties-scalar_drift.c \
cs_user_porosity.c \
cs_user_postprocess.c \
cs_user_postprocess-boundary_pressure_profile.c \
cs_user_postprocess-profiles.c \
cs_user_postprocess-probes-interpolate.c \
cs_user_postprocess-sfc.c \
cs_user_radiative_transfer.c \
cs_user_radiative_transfer_bcs.c \
cs_user_solver-heat-equation.c \
cs_user_source_terms-along_line.c \
cs_user_source_terms-base.c \
cs_user_source_terms-momentum.c \
cs_user_source_terms-scalar_in_a_channel.c \
cs_user_source_terms-turbulence.c \
cs_user_time_table.c \
cs_user_turbomachinery.c \
cs_user_zones.c \
cs_user_atmospheric_model.f90 \
cs_user_boundary_conditions-atmo.f90 \
cs_user_boundary_conditions-gas_3ptchem.c \
cs_user_boundary_conditions-gas_ebu.f90 \
cs_user_boundary_conditions-gas_libby_williams.f90 \
cs_user_initialization-gas_ebu.f90 \
cs_user_initialization-gas_libby_williams.f90 \
cs_user_parameters.f90 \
cs_user_profile.cxx \
cs_user_profile.h \
cs_user_profile_plot.py

noinst_LIBRARIES = libcsexample.a
libcsexample_a_SOURCES = $(examplefiles)

exampledir = $(pkgdatadir)/user_sources/EXAMPLES
example_DATA = $(examplefiles)

clean-local:
	-rm -f *__genmod.f90 *.mod

# Finish example files installation
# - tag with the package version number
# Done with two sed passes (does not work otherwise)

install-data-hook:
	c_user_files=`ls ${top_srcdir}/src/user_examples/*.c* 2>/dev/null`; \
	for f in $${c_user_files}; do \
	  file=`basename "$${f}"` && \
	  $(PYTHON) "${top_srcdir}/build-aux/cs_user_src_filter.py" "$(cs_version_short)" < "${top_srcdir}/src/user_examples/$${file}" > "$(DESTDIR)$(exampledir)/$${file}"; \
	done ; \
	h_user_files=`ls ${top_srcdir}/src/user_examples/*.h 2>/dev/null`; \
	for f in $${h_user_files}; do \
	  file=`basename "$${f}"` && \
	  $(PYTHON) "${top_srcdir}/build-aux/cs_user_src_filter.py" "$(cs_version_short)" < "${top_srcdir}/src/user_examples/$${file}" > "$(DESTDIR)$(exampledir)/$${file}"; \
	done ; \
	cd $(DESTDIR)$(exampledir) && \
	f90_user_examples=`ls *.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!                      code_saturne version $(cs_version_short)\n"; \
	version2="\!                      ------------------------"; \
	for file in $${f90_user_examples}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done
