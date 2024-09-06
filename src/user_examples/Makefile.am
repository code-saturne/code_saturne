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

AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

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
cs_user_porosity-from_cad.cpp

# Library source files

examplefiles = \
cs_user_1d_wall_thermal.cpp \
cs_user_wall_condensation.cpp \
cs_user_boundary_conditions.cpp \
cs_user_boundary_conditions-advanced.cpp \
cs_user_boundary_conditions_ale-base.cpp \
cs_user_boundary_conditions-atmospheric.cpp \
cs_user_boundary_conditions-auto_inlet_profile.cpp \
cs_user_boundary_conditions-base.cpp \
cs_user_boundary_conditions-compressible.cpp \
cs_user_boundary_conditions-cooling_towers.cpp \
cs_user_boundary_conditions-electric_arcs.cpp \
cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.cpp \
cs_user_boundary_conditions-mapped_inlet.cpp \
cs_user_boundary_conditions-medcoupling_remapper_2d.cpp \
cs_user_boundary_conditions-notebook.cpp \
cs_user_boundary_conditions-pulverized_coal.cpp \
cs_user_boundary_conditions-time_table.cpp \
cs_user_coupling-saturne.cpp \
cs_user_coupling-syrthes.cpp \
cs_user_electric_scaling.cpp \
cs_user_extra_operations-balance_by_zone.cpp \
cs_user_extra_operations-boundary_forces.cpp \
cs_user_extra_operations-force_temperature.cpp \
cs_user_extra_operations-mean_profiles.cpp \
cs_user_extra_operations-medcoupling_slice.cpp \
cs_user_extra_operations-nusselt_and_friction_calculation.cpp \
cs_user_extra_operations-nusselt_calculation.cpp \
cs_user_extra_operations-parallel_operations.cpp \
cs_user_extra_operations-scalar_balance.cpp \
cs_user_extra_operations-stopping_criterion.cpp \
cs_user_extra_operations-turbomachinery.cpp \
cs_user_extra_operations-verif_cdo_diffusion.cpp \
cs_user_extra_operations-vorticity_field.cpp \
cs_user_fluid_structure_interaction.cpp \
cs_user_fluid_structure_interaction-code_aster.cpp \
cs_user_head_losses.cpp \
cs_user_hgn-relax_time_scale.cpp \
cs_user_ibm.cpp \
cs_user_ibm_modify.cpp \
cs_user_initialization-atmospheric.cpp \
cs_user_initialization-base.cpp \
cs_user_initialization-compressible.cpp \
cs_user_initialization-cooling_towers.cpp \
cs_user_initialization-electric_arcs.cpp \
cs_user_initialization-gas_3ptchem.cpp \
cs_user_initialization-medcoupling_remapper_3d.cpp \
cs_user_initialization-pulverized_coal.cpp \
cs_user_initialization-time_step.cpp \
cs_user_lagr_boundary_conditions.cpp \
cs_user_lagr_volume_conditions.cpp \
cs_user_lagr_model.cpp \
cs_user_lagr_particle.cpp \
cs_user_lagr_particle-coal.cpp \
cs_user_les_inflow-base.cpp \
cs_user_mesh-boundary.cpp \
cs_user_mesh-cartesian.cpp \
cs_user_mesh-groups_from_free_faces.cpp \
cs_user_mesh-input-save.cpp \
cs_user_mesh-joining.cpp \
cs_user_mesh-modify.cpp \
cs_user_mesh-periodicity.cpp \
cs_user_mesh-quality.cpp \
cs_user_paramedmem_coupling-base.cpp \
cs_user_parameters-base.cpp \
cs_user_parameters-cdo-condif.cpp \
cs_user_parameters-cdo-gwf.cpp \
cs_user_parameters-cdo-gwf_user_soil.cpp \
cs_user_parameters-cdo-laplacian.cpp \
cs_user_parameters-cdo-navsto.cpp \
cs_user_parameters-cdo-solidification.cpp \
cs_user_parameters-ctwr.cpp \
cs_user_parameters-coupling.cpp \
cs_user_parameters-fans.cpp \
cs_user_parameters-cdo-linear_solvers.cpp \
cs_user_parameters-linear_solvers.cpp \
cs_user_parameters-rotation.cpp \
cs_user_parameters-volume_mass_injection.cpp \
cs_user_parameters-time_moments.cpp \
cs_user_performance_tuning-matrix.cpp \
cs_user_performance_tuning-numbering.cpp \
cs_user_performance_tuning-parallel-io.cpp \
cs_user_performance_tuning-partition.cpp \
cs_user_physical_properties.cpp \
cs_user_physical_properties-compressible_flow.cpp \
cs_user_physical_properties-coal_drift.cpp \
cs_user_physical_properties-scalar_drift.cpp \
cs_user_porosity.cpp \
cs_user_postprocess.cpp \
cs_user_postprocess-boundary_pressure_profile.cpp \
cs_user_postprocess-profiles.cpp \
cs_user_postprocess-probes-interpolate.cpp \
cs_user_postprocess-sfc.cpp \
cs_user_radiative_transfer.cpp \
cs_user_radiative_transfer_bcs.cpp \
cs_user_solver-heat-equation.cpp \
cs_user_source_terms-along_line.cpp \
cs_user_source_terms-base.cpp \
cs_user_source_terms-momentum.cpp \
cs_user_source_terms-scalar_in_a_channel.cpp \
cs_user_source_terms-turbulence.cpp \
cs_user_time_table.cpp \
cs_user_turbomachinery.cpp \
cs_user_zones.cpp \
cs_user_atmospheric_model.f90 \
cs_user_boundary_conditions-gas_3ptchem.cpp \
cs_user_profile.cxx \
cs_user_profile.h \
cs_user_profile_plot.py

noinst_LIBRARIES = libcsexample.a
libcsexample_a_SOURCES = $(examplefiles)

exampledir = $(pkgdatadir)/user_sources/EXAMPLES
example_DATA = $(examplefiles)

clean-local:
	-rm -f *__genmod.f90 *.mod
	-rm -rf .deps

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/user_examples -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

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
