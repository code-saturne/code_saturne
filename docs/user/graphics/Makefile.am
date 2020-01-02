## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

# Xfig files to process

FIGFILES = \
cs_components.fig \
join_periodic.fig \
join_tolerance.fig \
vortex.fig

# Other files

OTHERFILES = \
gui_ale_external.png \
gui_ale_internal.png \
gui_ale_mei.png \
gui_atmo_read.png \
gui_atmospheric_user_s_guide_v91.png \
gui_atmospheric_user_s_guide_v92.png \
gui_atmospheric_user_s_guide_v93.png \
gui_atmospheric_user_s_guide_v95.png \
gui_atmospheric_user_s_guide_v96.png \
gui_atmospheric_user_s_guide_v97.png \
gui_atmospheric_user_s_guide_v98.png \
gui_bc_parameters.png \
gui_bc_regions.png \
gui_calculation_features.png \
gui_case_dir.png \
gui_coal_bc.png \
gui_coal_classes.png \
gui_coal_composition.png \
gui_coal_oxydant.png \
gui_coal_reaction.png \
gui_coal_solution_control.png \
gui_coal_time_average.png \
gui_body_forces.png \
gui_density_law.png \
gui_fluid_props.png \
gui_formula_button.png \
gui_global_res_parameters.png \
gui_head_loss_coeffs.png \
gui_head_loss_regions.png \
gui_initialisation.png \
gui_lagr_bc.png \
gui_lagr_global_advanced.png \
gui_lagr_global_settings.png \
gui_lagr_output.png \
gui_lagr_statistics.png \
gui_mesh_join.png \
gui_mesh_periodicity.png \
gui_numerical_parameters.png \
gui_output_log.png \
gui_output_meshes.png \
gui_output_profiles.png \
gui_output_writers.png \
gui_phys_prop_reference_values.png \
gui_prepare_execution.png \
gui_rad_transf_do_advanced.png \
gui_rad_transf_do_params.png \
gui_rad_transf_p1_advanced.png \
gui_rad_transf_p1_params.png \
gui_rad_transf_post_output.png \
gui_rad_transf_wall_model.png \
gui_rad_transf_wall_params.png \
gui_syrthes_coupling.png \
gui_thermal_scalar.png \
gui_time_averages.png \
gui_time_step.png \
gui_turbulence_models.png \
gui_user_scal_def_init.png \
halo.pdf \
rota_perio_parall.jpg

# Distributed files

EXTRA_DIST = $(FIGFILES) $(OTHERFILES)

# Generated documents (pdf files)

PDF_FIGFILES = $(FIGFILES:.fig=.pdf)

# New suffixes and targets

SUFFIXES = .fig .pdf

# On Windows, fig2dev may be broken and generate a file named "-" instead
# of the target. Thus, we try to renamed such files if they exist.

.fig.pdf:
	"$(FIG2DEV)" -L pdf $< > $@
	@if test -f - ; then \mv -f - $@ ; fi

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local: $(PDF_FIGFILES)

clean-local:
	-rm -f $(PDF_FIGFILES)

distclean-local: clean-local
