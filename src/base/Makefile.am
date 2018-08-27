## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/rayt \
-I$(top_srcdir)/src/turb \
-I$(top_srcdir)/src/atmo \
$(CPPFLAGS_PLE) \
$(FREESTEAM_CPPFLAGS) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/alge \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/darc \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCMODINCLUDE)$(top_builddir)/src/turb \
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_1d_wall_thermal.h \
cs_1d_wall_thermal_check.h \
cs_ale.h \
cs_all_to_all.h \
cs_array_reduce.h \
cs_assert.h \
cs_ast_coupling.h \
cs_base.h \
cs_base_fortran.h \
cs_block_dist.h \
cs_block_to_part.h \
cs_boundary_conditions.h \
cs_boundary_zone.h \
cs_calcium.h \
cs_control.h \
cs_coupling.h \
cs_crystal_router.h \
cs_defs.h \
cs_equation_iterative_solve.h \
cs_ext_neighborhood.h \
cs_fan.h \
cs_file.h \
cs_field.h \
cs_field_default.h \
cs_field_pointer.h \
cs_field_operator.h \
cs_flag_check.h \
cs_fp_exception.h \
cs_gas_mix.h \
cs_halo.h \
cs_halo_perio.h \
cs_head_losses.h \
cs_interface.h \
cs_interpolate.h \
cs_internal_coupling.h \
cs_io.h \
cs_log.h \
cs_log_iteration.h \
cs_log_setup.h \
cs_map.h \
cs_math.h \
cs_measures_util.h \
cs_rank_neighbors.h \
cs_numbering.h \
cs_opts.h \
cs_order.h \
cs_parameters.h \
cs_parameters_check.h \
cs_parall.h \
cs_part_to_block.h \
cs_physical_constants.h \
cs_physical_properties.h \
cs_post.h \
cs_post_default.h \
cs_post_util.h \
cs_preprocess.h \
cs_preprocessor_data.h \
cs_probe.h \
cs_prototypes.h \
cs_random.h \
cs_range_set.h \
cs_renumber.h \
cs_resource.h \
cs_restart.h \
cs_restart_default.h \
cs_rotation.h \
cs_sat_coupling.h \
cs_search.h \
cs_selector.h \
cs_sort.h \
cs_sort_partition.h \
cs_stokes_model.h \
cs_syr4_coupling.h \
cs_syr_coupling.h \
cs_system_info.h \
cs_thermal_model.h \
cs_time_moment.h \
cs_time_plot.h \
cs_time_step.h \
cs_timer.h \
cs_timer_stats.h \
cs_tree.h \
cs_turbomachinery.h \
cs_volume_zone.h \
cs_wall_functions.h \
cs_zone.h

# Library source files

noinst_LTLIBRARIES = libcscorep.la libcscore.la \
		     libcsbase.la libcsrenumber.la

# Core library.
# We use a helper library containing only cs_base.c, so as to pass
# it localedir and pkgdatadir info.

libcscorep_la_CPPFLAGS = \
-DLOCALEDIR=\"'$(localedir)'\" \
-DPKGDATADIR=\"'$(pkgdatadir)'\" \
-DPKGLIBDIR=\"'$(pkglibdir)'\" \
$(AM_CPPFLAGS)

libcscorep_la_SOURCES = cs_base.c
libcscorep_la_LDFLAGS = -no-undefined

libcscore_la_SOURCES = \
cs_all_to_all.c \
cs_block_dist.c \
cs_block_to_part.c \
cs_crystal_router.c \
cs_defs.c \
cs_file.c \
cs_fp_exception.c \
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
libcscore_la_LDFLAGS = -no-undefined
libcscore_la_LIBADD = libcscorep.la $(top_builddir)/src/fvm/libfvm.la

###

libcsbase_la_SOURCES = \
addfld.f90 \
alelav.f90 \
alemav.f90 \
altycl.f90 \
b_t_to_h.f90 \
b_h_to_t.f90 \
c_h_to_t.f90 \
calhyd.f90 \
caltri.f90 \
catsma.f90 \
catsmv.f90 \
catsmt.f90 \
cdomod.f90 \
clca66.f90 \
clpsca.f90 \
clptrg.f90 \
clptur.f90 \
clsyvt.f90 \
copain_model.f90 \
condensation_source_terms.f90 \
condli.f90 \
cou1di.f90 \
cou1do.f90 \
coupbi.f90 \
coupbo.f90 \
covofi.f90 \
covofv.f90 \
cptssy.f90 \
cpvosy.f90 \
cs_1d_wall_thermal.c \
cs_1d_wall_thermal_check.c \
cs_ale.c \
cs_array_reduce.c \
cs_ast_coupling.c \
cs_base_fortran.c \
cs_boundary_conditions.c \
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
cs_control.c \
cs_coupling.c \
cs_ext_neighborhood.c \
cs_equation_iterative_solve.c \
cs_fan.c \
cs_field.c \
cs_field_default.c \
cs_field_pointer.c \
cs_field_operator.c \
cs_flag_check.c \
cs_gas_mix.c \
cs_halo.c \
cs_halo_perio.c \
cs_head_losses.c \
cs_interpolate.c \
csinit.f90 \
cs_log_iteration.c \
cs_log_setup.c \
cs_numbering.c \
cs_measures_util.c \
cs_mesh_tagmr.f90 \
cs_metal_structures_tag.f90 \
cs_gas_mix_initialization.f90 \
cs_gas_mix_physical_properties.f90 \
csopli.f90 \
cs_opts.c \
cs_parall.c \
cs_parameters.c \
cs_parameters_check.c \
cs_physical_constants.c \
cs_physical_properties.c \
cs_post.c \
cs_post_default.c \
cs_post_util.c \
cs_preprocess.c \
cs_preprocessor_data.c \
cs_probe.c \
cs_random.c \
cs_range_set.c \
cs_resource.c \
cs_restart.c \
cs_restart_default.c \
cs_rotation.c \
cs_sat_coupling.c \
cs_search.c \
cs_selector.c \
cs_selector_f2c.f90 \
cs_sort.c \
cs_sort_partition.c \
cs_stokes_model.c \
cs_syr4_coupling.c \
cs_syr_coupling.c \
cs_thermal_model.c \
cs_tagmri.f90 \
cs_tagmro.f90 \
cs_time_moment.c \
cs_time_plot.c \
cs_time_step.c \
cs_timer_stats.c \
cs_turbomachinery.c \
cs_volume_zone.c \
cs_wall_functions.c \
cs_internal_coupling.c \
csprnt.f90 \
diffst.f90 \
distpr.f90 \
distpr2.f90 \
distyp.f90 \
driflu.f90 \
dttvar.f90 \
dvvpst.f90 \
ecrava.f90 \
ecrlis.f90 \
findpt.f90 \
fldini.f90 \
fldtri.f90 \
fldvar.f90 \
fldprp.f90 \
haltyp.f90 \
idrbla.f90 \
impini.f90 \
iniini.f90 \
initi1.f90 \
initi2.f90 \
iniusi.f90 \
iniva0.f90 \
inivar.f90 \
iprbla.f90 \
lecamo.f90 \
lecamp.f90 \
lecamx.f90 \
majgeo.f90 \
metal_structures_copain_model.f90 \
modini.f90 \
mmtycl.f90 \
navstv.f90 \
newmrk.f90 \
phyvar.f90 \
post_util.f90 \
precli.f90 \
predfl.f90 \
predvv.f90 \
prehyd.f90 \
prmoy0.f90 \
pthrbm.f90 \
resopv.f90 \
resvoi.f90 \
scalai.f90 \
schtmp.f90 \
stdtcl.f90 \
strdep.f90 \
strhis.f90 \
strini.f90 \
strpre.f90 \
tdesi1.f90 \
tridim.f90 \
tspdcv.f90 \
typecl.f90 \
varpos.f90 \
vericl.f90 \
verini.f90 \
verlon.f90 \
visecv.f90
libcsbase_la_LDFLAGS = -no-undefined
libcsbase_la_LIBADD =

# Renumbering (may require extra headers)

libcsrenumber_la_CPPFLAGS = $(AM_CPPFLAGS) \
$(METIS_CPPFLAGS) $(SCOTCH_CPPFLAGS)
libcsrenumber_la_SOURCES = cs_renumber.c
libcsrenumber_la_LDFLAGS = -no-undefined

# Conditionally compiled extensions
lib_LTLIBRARIES =
pkglib_LTLIBRARIES =

if HAVE_PARAMEDMEM

pkginclude_HEADERS += \
cs_paramedmem_coupling.hxx \
cs_medcoupling_remapper.hxx \
cs_medcoupling_utils.hxx

if HAVE_PLUGIN_MEDCOUPLING

cs_paramedmem_la_CPPFLAGS = \
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
pkglib_LTLIBRARIES += cs_paramedmem.la
cs_paramedmem_la_SOURCES = \
cs_paramedmem_coupling.cxx \
cs_medcoupling_remapper.cxx \
cs_medcoupling_utils.cxx
cs_paramedmem_la_LDFLAGS = -module -avoid-version
cs_paramedmem_la_LIBADD = \
$(MEDCOUPLING_LDFLAGS) $(MEDCOUPLING_LIBS) \
$(MEDCOUPLINGRUNPATH)

else

libcs_paramedmem_la_CPPFLAGS = \
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
libcs_paramedmem_la_SOURCES = \
cs_paramedmem_coupling.cxx \
cs_medcoupling_remapper.cxx \
cs_medcoupling_utils.cxx
libcs_paramedmem_la_LIBADD = \
$(MEDCOUPLING_LDFLAGS) $(MEDCOUPLING_LIBS) \
$(MEDCOUPLINGRUNPATH)
noinst_LTLIBRARIES += libcs_paramedmem.la
libcsbase_la_LIBADD += libcs_paramedmem.la

endif

endif

if HAVE_EOS

cs_eos_la_CPPFLAGS = \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(EOS_CPPFLAGS) $(MPI_CPPFLAGS)
pkginclude_HEADERS += cs_eos.hxx
pkglib_LTLIBRARIES += cs_eos.la
cs_eos_la_SOURCES = cs_eos.cxx
cs_eos_la_LDFLAGS = -module -avoid-version

endif

if HAVE_COOLPROP

if HAVE_PLUGINS

cs_coolprop_la_CPPFLAGS = \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(COOLPROP_CPPFLAGS) $(MPI_CPPFLAGS)
cs_coolprop_la_LIBADD = \
$(COOLPROP_LDFLAGS) $(COOLPROP_LIBS) \
$(COOLPROPRUNPATH)
pkginclude_HEADERS += cs_coolprop.hxx
pkglib_LTLIBRARIES += cs_coolprop.la
cs_coolprop_la_SOURCES = cs_coolprop.cxx
cs_coolprop_la_LDFLAGS = -module -avoid-version

else

libcs_coolprop_la_CPPFLAGS = \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
$(COOLPROP_CPPFLAGS) $(MPI_CPPFLAGS)
pkginclude_HEADERS += cs_coolprop.hxx
libcs_coolprop_la_LIBADD = \
$(COOLPROP_LDFLAGS) $(COOLPROP_LIBS) \
$(COOLPROPRUNPATH)
libcs_coolprop_la_SOURCES = cs_coolprop.cxx
noinst_LTLIBRARIES += libcs_coolprop.la
libcsbase_la_LIBADD += libcs_coolprop.la

endif

endif

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod
