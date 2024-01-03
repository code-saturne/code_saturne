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
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

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
cs_boundary_conditions_set_coeffs.h \
cs_boundary_conditions_set_coeffs_symmetry.h \
cs_boundary_conditions_set_coeffs_turb.h \
cs_boundary_conditions_type.h \
cs_boundary_zone.h \
cs_calcium.h \
cs_time_step_compute.h \
cs_control.h \
cs_coupling.h \
cs_crystal_router.h \
cs_defs.h \
cs_drift_convective_flux.h \
cs_equation_iterative_solve.h \
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
cs_restart_map.h \
cs_rotation.h \
cs_runaway_check.h \
cs_sat_coupling.h \
cs_scalar_clipping.h \
cs_search.h \
cs_selector.h \
cs_solid_zone.h \
cs_solve_equation.h \
cs_solve_navier_stokes.h \
cs_sort.h \
cs_sort_partition.h \
cs_syr_coupling.h \
cs_sys_coupling.h \
cs_system_info.h \
cs_thermal_model.h \
cs_time_moment.h \
cs_time_control.h \
cs_time_plot.h \
cs_time_step.h \
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
cs_halo_cuda.h \
cs_cuda_contrib.h
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
cs_all_to_all.c \
cs_block_dist.c \
cs_block_to_part.c \
cs_crystal_router.c \
cs_defs.c \
cs_file.c \
cs_file_csv_parser.c \
cs_fp_exception.c \
cs_ht_convert.c \
cs_interface.c \
cs_io.c \
cs_log.c \
cs_math.c \
cs_rank_neighbors.c \
cs_map.c \
cs_order.c \
cs_part_to_block.c \
cs_system_info.c \
cs_timer.c \
cs_tree.c

if HAVE_ACCEL
libcscore_a_SOURCES += cs_base_accel.cxx
endif

if HAVE_CUDA
libcscore_a_SOURCES += \
cs_base_cuda.cu \
cs_halo_cuda.cu
endif

libcsbase_a_SOURCES = \
addfld.f90 \
caltri.f90 \
cs_mass_source_terms.c \
condli.f90 \
cou1di.f90 \
cou1do.f90 \
cs_wall_condensation.c \
cs_wall_condensation_1d_thermal.c \
cs_1d_wall_thermal.c \
cs_1d_wall_thermal_check.c \
cs_ale.c \
cs_array.c \
cs_array_reduce.c \
cs_ast_coupling.c \
cs_base_fortran.c \
cs_boundary.c \
cs_boundary_conditions.c \
cs_boundary_conditions_check.c \
cs_boundary_conditions_set_coeffs.c \
cs_boundary_conditions_set_coeffs_symmetry.c \
cs_boundary_conditions_set_coeffs_turb.c \
cs_boundary_conditions_type.c \
cs_boundary_zone.c \
cs_calcium.c \
csc2cl.f90 \
csc2ts.f90 \
csccel.f90 \
cscfbr.f90 \
cscini.f90 \
csclli.f90 \
cscloc.f90 \
cscpce.f90 \
cscpfb.f90 \
csinit.f90 \
csopli.f90 \
csprnt.f90 \
cs_time_step_compute.c \
cs_control.c \
cs_coupling.c \
cs_drift_convective_flux.c \
cs_ext_library_info.c \
cs_ext_neighborhood.c \
cs_equation_iterative_solve.c \
cs_fan.c \
cs_field.c \
cs_field_default.c \
cs_field_pointer.c \
cs_field_operator.c \
cs_flag_check.c \
cs_function.c \
cs_function_default.c \
cs_gas_mix.c \
cs_gas_mix_initialization.f90 \
cs_gas_mix_physical_properties.f90 \
cs_halo.c \
cs_halo_perio.c \
cs_head_losses.c \
cs_ibm.c \
cs_initialize_fields.c \
cs_internal_coupling.c \
cs_interpolate.c \
cs_log_iteration.c \
cs_log_setup.c \
cs_mobile_structures.c \
cs_measures_util.c \
cs_mesh_tagmr.f90 \
cs_notebook.c \
cs_numbering.c \
cs_opts.c \
cs_parall.c \
cs_param_types.c \
cs_parameters.c \
cs_parameters_check.c \
cs_physical_constants.c \
cs_physical_properties.c \
cs_porosity_from_scan.c \
cs_porous_model.c \
cs_post.c \
cs_post_default.c \
cs_post_util.c \
cs_preprocess.c \
cs_preprocessor_data.c \
cs_pressure_correction.c \
cs_probe.c \
cs_random.c \
cs_range_set.c \
cs_resource.c \
cs_restart.c \
cs_restart_default.c \
cs_restart_map.c \
cs_rotation.c \
cs_runaway_check.c \
cs_sat_coupling.c \
cs_scalar_clipping.c\
cs_search.c \
cs_selector.c \
cs_selector_f2c.f90 \
cs_solid_zone.c \
cs_solve_equation.c \
cs_solve_navier_stokes.c \
cs_sort.c \
cs_sort_partition.c \
cs_syr_coupling.c \
cs_sys_coupling.c \
cs_thermal_model.c \
cs_time_table.c \
cs_tagmri.f90 \
cs_tagmro.f90 \
cs_time_control.c \
cs_time_moment.c \
cs_time_plot.c \
cs_time_step.c \
cs_timer_stats.c \
cs_turbomachinery.c \
cs_utilities.c \
cs_velocity_pressure.c \
cs_vof.cxx \
cs_volume_mass_injection.c \
cs_volume_zone.c \
cs_wall_distance.c \
cs_wall_functions.c \
cs_xdef_eval_at_zone.c \
diffst.f90 \
distpr2.f90 \
ecrava.f90 \
ecrlis.f90 \
findpt.f90 \
fldini.f90 \
fldvar.f90 \
fldprp.f90 \
idrbla.f90 \
impini.f90 \
iniini.f90 \
initi1.f90 \
iniusi.f90 \
iniva0.f90 \
inivar.f90 \
iprbla.f90 \
lecamo.f90 \
lecamp.f90 \
lecamx.f90 \
majgeo.f90 \
modini.f90 \
mmtycl.f90 \
phyvar.f90 \
pthrbm.f90 \
scalai.f90 \
schtmp.f90 \
tdesi1.f90 \
tridim.f90 \
varpos.f90 \
verini.f90

libcsbase_a_LIBADD =
local_libraries_add =

# Renumbering (may require extra headers)

libcsrenumber_a_CPPFLAGS = $(AM_CPPFLAGS) \
$(METIS_CPPFLAGS) $(SCOTCH_CPPFLAGS)
libcsrenumber_a_SOURCES = cs_renumber.c

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

local_libraries_add += cs_eos.${SHARED_LIB_EXT}
LDADD_EOS = \
$(EOS_LDFLAGS) $(EOS_LIBS) $(EOSRUNPATH)

cs_eos.${SHARED_LIB_EXT}: libcs_eos.a
	@echo "Building cs_eos plugin..."
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py --echo --linker ${CS_LD} -o $@ libcs_eos_a-cs_eos.$(OBJEXT) ${LDFLAGS_SHARED} ${LDFLAGS} ${LDADD_EOS} $(AM_LIBS)

local_libraries_add += cs_eos.${SHARED_LIB_EXT}

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

.cu.$(OBJEXT):
	$(NVCC) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

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
