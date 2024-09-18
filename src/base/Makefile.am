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
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/cfbl \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/comb \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/gwf \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/meg \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
-I$(top_srcdir)/src/atmo \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/alge \
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

# Files not built through classical Automake rules need to be explicitely
# added to distribution.

EXTRA_DIST = \
cs_base.c \
cs_coolprop.cxx \
cs_version.h.in

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_1d_wall_thermal.h \
cs_1d_wall_thermal_check.h \
cs_ale.h \
cs_all_to_all.h \
cs_array.h \
cs_array_reduce.h \
cs_assert.h \
cs_ast_coupling.h \
cs_base.h \
cs_base_accel.h \
cs_base_fortran.h \
cs_block_dist.h \
cs_block_to_part.h \
cs_boundary.h \
cs_boundary_conditions.h \
cs_boundary_conditions_check.h \
cs_boundary_conditions_coupling.h\
cs_boundary_conditions_set_coeffs.h \
cs_boundary_conditions_set_coeffs_symmetry.h \
cs_boundary_conditions_set_coeffs_turb.h \
cs_boundary_conditions_type.h \
cs_boundary_zone.h \
cs_calcium.h \
cs_compute_thermo_pressure_density.h \
cs_time_step_compute.h \
cs_control.h \
cs_coupling.h \
cs_crystal_router.h \
cs_defs.h \
cs_dilatable_scalar_diff_st.h \
cs_dispatch.h \
cs_drift_convective_flux.h \
cs_equation_iterative_solve.h \
cs_execution_context.h \
cs_ext_library_info.h \
cs_ext_neighborhood.h \
cs_fan.h \
cs_field.h \
cs_field_default.h \
cs_field_pointer.h \
cs_field_operator.h \
cs_file.h \
cs_file_csv_parser.h \
cs_flag_check.h \
cs_fp_exception.h \
cs_function.h \
cs_function_default.h \
cs_gas_mix.h \
cs_halo.h \
cs_halo_perio.h \
cs_head_losses.h \
cs_ht_convert.h \
cs_ibm.h \
cs_initialize_fields.h \
cs_interface.h \
cs_interpolate.h \
cs_internal_coupling.h \
cs_io.h \
cs_log.h \
cs_log_iteration.h \
cs_log_setup.h \
cs_map.h \
cs_mass_source_terms.h \
cs_math.h \
cs_measures_util.h \
cs_mobile_structures.h \
cs_rank_neighbors.h \
cs_notebook.h \
cs_numbering.h \
cs_opts.h \
cs_order.h \
cs_param_types.h \
cs_parameters.h \
cs_parameters_check.h \
cs_parall.h \
cs_part_to_block.h \
cs_physical_constants.h \
cs_physical_properties.h \
cs_physical_properties_default.h \
cs_porosity_from_scan.h \
cs_porous_model.h \
cs_post.h \
cs_post_default.h \
cs_post_util.h \
cs_preprocess.h \
cs_preprocessor_data.h \
cs_pressure_correction.h \
cs_probe.h \
cs_prototypes.h \
cs_random.h \
cs_range_set.h \
cs_renumber.h \
cs_resource.h \
cs_restart.h \
cs_restart_default.h \
cs_restart_main_and_aux.h \
cs_restart_map.h \
cs_rotation.h \
cs_runaway_check.h \
cs_sat_coupling.h \
cs_scalar_clipping.h \
cs_search.h \
cs_selector.h \
cs_setup.h \
cs_solid_zone.h \
cs_solve_all.h \
cs_solve_equation.h \
cs_solve_navier_stokes.h \
cs_solve_transported_variables.h \
cs_sort.h \
cs_sort_partition.h \
cs_syr_coupling.h \
cs_sys_coupling.h \
cs_system_info.h \
cs_thermal_model.h \
cs_theta_scheme.h \
cs_time_moment.h \
cs_time_control.h \
cs_time_plot.h \
cs_time_step.h \
cs_time_stepping.h \
cs_time_table.h \
cs_timer.h \
cs_timer_stats.h \
cs_tree.h \
cs_turbomachinery.h \
cs_utilities.h \
cs_velocity_pressure.h \
cs_vof.h \
cs_volume_zone.h \
cs_volume_mass_injection.h \
cs_wall_condensation.h \
cs_wall_condensation_1d_thermal.h \
cs_wall_distance.h \
cs_wall_functions.h \
cs_xdef_eval_at_zone.h \
cs_zone.h \
cs_base_headers.h \
cs_headers.h

noinst_HEADERS = \
cs_boundary_conditions_priv.h

if HAVE_CUDA
pkginclude_HEADERS += \
cs_base_cuda.h \
cs_math_cuda.cuh \
cs_array_cuda.h \
cs_cuda_contrib.h

noinst_HEADERS += \
cs_halo_cuda.h
endif

# Library source files

noinst_LIBRARIES = \
libcscore.a \
libcsbase.a libcsrenumber.a

# Core library.
# We use a specific rule for cs_base.c, so as to pass
# it localedir and pkgdatadir info.

cs_base_CPPFLAGS = \
-DLOCALEDIR=\"'$(localedir)'\" \
-DPKGDATADIR=\"'$(pkgdatadir)'\" \
-DPKGLIBDIR=\"'$(pkglibdir)'\" \
$(AM_CPPFLAGS)

cs_base = $(top_srcdir)/src/base/cs_base.c
cs_base.$(OBJEXT): $(cs_base)
	$(CC) -c $(DEFS) $(DEFAULT_INCLUDES) $(cs_base_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) $(cs_base)

libcscore_a_LIBADD = cs_base.$(OBJEXT)

libcscore_a_SOURCES = \
cs_all_to_all.cpp \
cs_block_dist.cpp \
cs_block_to_part.cpp \
cs_crystal_router.cpp \
cs_defs.c \
cs_file.cpp \
cs_file_csv_parser.cpp \
cs_fp_exception.cpp \
cs_ht_convert.cpp \
cs_interface.cpp \
cs_io.cpp \
cs_log.cpp \
cs_map.cpp \
cs_math.cpp \
cs_order.cpp \
cs_part_to_block.cpp \
cs_rank_neighbors.cpp \
cs_system_info.cxx \
cs_timer.c \
cs_tree.cpp

if HAVE_ACCEL
libcscore_a_SOURCES += cs_base_accel.cxx
endif

if HAVE_CUDA
libcscore_a_SOURCES += \
cs_base_cuda.cu \
cs_halo_cuda.cu
endif

libcsbase_a_SOURCES = \
condli.f90 \
cs_1d_wall_thermal.cpp \
cs_1d_wall_thermal_check.cpp \
cs_ale.cpp \
cs_array.cpp \
cs_array_reduce.cpp \
cs_ast_coupling.cpp \
cs_base_fortran.c \
cs_boundary.cpp \
cs_boundary_conditions_check.cpp \
cs_boundary_conditions_coupling.cpp\
cs_boundary_conditions_set_coeffs_symmetry.cpp \
cs_boundary_conditions_set_coeffs_turb.cpp \
cs_boundary_conditions_set_coeffs.cpp \
cs_boundary_conditions_type.cpp \
cs_boundary_conditions.cpp \
cs_boundary_zone.cpp \
cs_calcium.cpp \
cs_compute_thermo_pressure_density.cpp \
cs_control.cpp \
cs_coupling.cpp \
cs_dilatable_scalar_diff_st.cpp \
cs_drift_convective_flux.cpp \
cs_equation_iterative_solve.cpp \
cs_ext_library_info.c \
cs_ext_neighborhood.cpp \
cs_fan.cpp \
cs_field.c \
cs_field_default.c \
cs_field_operator.cpp \
cs_field_pointer.c \
cs_flag_check.cpp \
cs_function.cpp \
cs_function_default.cpp \
cs_gas_mix.cpp \
cs_halo_perio.cpp \
cs_halo.cpp \
cs_head_losses.cpp \
cs_ibm.cpp \
cs_initialize_fields.cpp \
cs_internal_coupling.cpp \
cs_interpolate.c \
cs_log_iteration.cpp \
cs_log_setup.cpp \
cs_mass_source_terms.cpp \
cs_measures_util.cpp \
cs_mobile_structures.cpp \
cs_notebook.cpp \
cs_numbering.cpp \
cs_opts.cpp \
cs_parall.cpp \
cs_param_types.cpp \
cs_parameters.c \
cs_parameters_check.cpp \
cs_physical_constants.c \
cs_physical_properties.c \
cs_physical_properties_default.c \
cs_porosity_from_scan.cpp \
cs_porous_model.cpp \
cs_post.c \
cs_post_default.cpp \
cs_post_util.cpp \
cs_preprocess.c \
cs_preprocessor_data.cpp \
cs_pressure_correction.cpp \
cs_probe.cpp \
cs_random.cpp \
cs_range_set.cpp \
cs_resource.cpp \
cs_restart.cpp \
cs_restart_default.cpp \
cs_restart_main_and_aux.cpp \
cs_restart_map.cpp \
cs_rotation.c \
cs_runaway_check.cpp \
cs_sat_coupling.cpp \
cs_scalar_clipping.c\
cs_search.cpp \
cs_selector.cpp \
cs_selector_f2c.f90 \
cs_setup.c \
cs_solid_zone.cpp \
cs_solve_all.cpp \
cs_solve_equation.cpp \
cs_solve_navier_stokes.cpp \
cs_solve_transported_variables.cpp \
cs_sort.cpp \
cs_sort_partition.cpp \
cs_syr_coupling.cpp \
cs_sys_coupling.cpp \
cs_thermal_model.cpp \
cs_theta_scheme.cpp \
cs_time_control.cpp \
cs_time_moment.c \
cs_time_plot.cpp \
cs_time_step_compute.cpp \
cs_time_step.c \
cs_time_stepping.cpp \
cs_time_table.cpp \
cs_timer_stats.cpp \
cs_turbomachinery.c \
cs_utilities.cpp \
cs_velocity_pressure.cpp \
cs_vof.cpp \
cs_volume_mass_injection.cpp \
cs_volume_zone.cpp \
cs_wall_condensation.c \
cs_wall_condensation_1d_thermal.cpp \
cs_wall_distance.cpp \
cs_wall_functions.c \
cs_xdef_eval_at_zone.cpp \
csclli.f90 \
csinit.f90 \
csopli.f90 \
csprnt.f90 \
findpt.f90 \
fldprp.f90 \
fldvar.f90 \
impini.f90 \
iniini.f90 \
inivar.f90 \
majgeo.f90 \
varpos.f90

libcsbase_a_LIBADD =
local_libraries_add =

# Renumbering (may require extra headers)

libcsrenumber_a_CPPFLAGS = $(AM_CPPFLAGS) \
$(METIS_CPPFLAGS) $(SCOTCH_CPPFLAGS)
libcsrenumber_a_SOURCES = cs_renumber.cpp

# Conditionally compiled extensions
lib_LIBRARIES =
pkglib_LIBRARIES =

# MEDCoupling wrappers, plugin or linked-in.

pkginclude_HEADERS += \
cs_medcoupling_mesh.hxx \
cs_medcoupling_utils.h \
cs_medcoupling_remapper.h \
cs_medcoupling_intersector.h \
cs_medcoupling_postprocess.h \
cs_paramedmem_coupling.h \
cs_paramedmem_remapper.h

libcs_medcoupling_a_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(CPPFLAGS_PLE) \
$(MEDCOUPLING_CPPFLAGS) \
$(MED_CPPFLAGS) \
$(HDF5_CPPFLAGS) \
$(MPI_CPPFLAGS)

libcs_medcoupling_a_SOURCES = \
cs_medcoupling_utils.cxx \
cs_medcoupling_mesh.cxx \
cs_medcoupling_remapper.cxx \
cs_medcoupling_intersector.cxx \
cs_medcoupling_postprocess.cxx \
cs_paramedmem_coupling.cxx \
cs_paramedmem_remapper.cxx

noinst_LIBRARIES += libcs_medcoupling.a

if HAVE_PLUGIN_MEDCOUPLING

local_libraries_add += cs_medcoupling.${SHARED_LIB_EXT}
LDADD_MEDCOUPLING = \
$(MEDCOUPLING_LDFLAGS) $(MEDCOUPLING_LIBS) \
$(MEDCOUPLINGRUNPATH)

cs_medcoupling.${SHARED_LIB_EXT}: libcs_medcoupling.a
	@echo "Building cs_medcoupling plugin..."
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py --echo --linker ${CS_LD} -o $@ --whole-archive-start libcs_medcoupling.a --whole-archive-end ${LDFLAGS_SHARED} ${LDFLAGS} ${LDADD_MEDCOUPLING} $(AM_LIBS)

local_libraries_add += cs_medcoupling.${SHARED_LIB_EXT}

else

libcsbase_a_LIBADD += \
libcs_medcoupling_a-cs_medcoupling_intersector.$(OBJEXT) \
libcs_medcoupling_a-cs_medcoupling_mesh.$(OBJEXT) \
libcs_medcoupling_a-cs_medcoupling_postprocess.$(OBJEXT) \
libcs_medcoupling_a-cs_medcoupling_remapper.$(OBJEXT) \
libcs_medcoupling_a-cs_medcoupling_utils.$(OBJEXT) \
libcs_medcoupling_a-cs_paramedmem_coupling.$(OBJEXT) \
libcs_medcoupling_a-cs_paramedmem_remapper.$(OBJEXT)

endif

# EOS (Equations of State) support, always as plugin.

if HAVE_EOS

libcs_eos_a_CPPFLAGS = \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(EOS_CPPFLAGS) $(MPI_CPPFLAGS)

pkginclude_HEADERS += cs_eos.hxx
noinst_LIBRARIES += libcs_eos.a
libcs_eos_a_SOURCES = cs_eos.cxx

if ENABLE_SHARED
local_libraries_add += cs_eos.${SHARED_LIB_EXT}
LDADD_EOS = \
$(EOS_LDFLAGS) $(EOS_LIBS) $(EOSRUNPATH)

cs_eos.${SHARED_LIB_EXT}: libcs_eos.a
	@echo "Building cs_eos plugin..."
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py --echo --linker ${CS_LD} -o $@ libcs_eos_a-cs_eos.$(OBJEXT) ${LDFLAGS_SHARED} ${LDFLAGS} ${LDADD_EOS} $(AM_LIBS)

local_libraries_add += cs_eos.${SHARED_LIB_EXT}

else

libcsbase_a_LIBADD += \
libcs_eos_a-cs_eos.$(OBJEXT)

endif

endif

# Coolprop library, plugin or linked-in.

if HAVE_COOLPROP

cs_coolprop_CPPFLAGS = \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(COOLPROP_CPPFLAGS) $(MPI_CPPFLAGS)
pkginclude_HEADERS += cs_coolprop.hxx

cs_coolprop = $(top_srcdir)/src/base/cs_coolprop.cxx
cs_coolprop.$(OBJEXT): $(cs_coolprop)
	$(CXX) -c $(DEFS) $(DEFAULT_INCLUDES) $(cs_coolprop_CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) ${cs_coolprop}

if HAVE_PLUGINS

local_libraries_add += cs_coolprop.${SHARED_LIB_EXT}
LDADD_COOLPROP = \
$(COOLPROP_LDFLAGS) $(COOLPROP_LIBS) \
$(COOLPROPRUNPATH)

cs_coolprop.${SHARED_LIB_EXT}: cs_coolprop.$(OBJEXT)
	@echo "Building cs_coolprop plugin..."
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py --echo --linker ${CS_LD} -o $@ cs_coolprop.$(OBJEXT) ${LDFLAGS_SHARED} ${CXXFLAGS_SHARED} ${LDFLAGS} ${LDADD_COOLPROP} $(AM_LIBS)

local_libraries_add += cs_coolprop.${SHARED_LIB_EXT}

else

libcsbase_a_LIBADD += cs_coolprop.$(OBJEXT)

endif

endif

# Rules for CUDA (not known by Automake)

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

.cu.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

# Additional rules

all-local: ${local_libraries_add}

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    echo "rm -f $$p"; \
	    rm -f $$p; \
	  else :; fi; \
	done
	-rm -rf .deps

# Remark:
# When installing shared libraries, remove library then copy it, rather than
# simply overwriting it, so that if a running instance is using that library,
# it will still access the old copy (not really removed until closed based
# on POSIX file semantics), instead of a mofified copy that can lead to
# SIGBUS errors due to an incomplete (in progress) copy or possibly other
# address changes.

install-exec-local: ${local_libraries_add}
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list" -a -n "$(pkglibdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(pkglibdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(pkglibdir)" || exit 1; \
	  for p in $$list; do \
	    echo " rm -f '$(DESTDIR)$(pkglibdir)/$$p'"; \
	    rm -f "$(DESTDIR)$(pkglibdir)/$$p"; \
	    echo " cp -fp $$p '$(DESTDIR)$(pkglibdir)'"; \
	    cp -fpP $$p "$(DESTDIR)$(pkglibdir)"; \
	  done; \
	}

uninstall-local:
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	for p in $$list; do \
	  echo " rm -f '$(DESTDIR)$(pkglibdir)/$$p'"; \
	  rm -f "$(DESTDIR)$(pkglibdir)/$$p"; \
	done; \
	rmdir "$(DESTDIR)$(pkglibdir)" || exit 0
