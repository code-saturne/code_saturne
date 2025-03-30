#ifndef __CS_BASE_HEADERS_H__
#define __CS_BASE_HEADERS_H__

/*============================================================================
 * Global code_saturne headers file for easier include
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_1d_wall_thermal.h"
#include "base/cs_1d_wall_thermal_check.h"
#include "base/cs_ale.h"
#include "base/cs_all_to_all.h"
#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_array_2dspan.h"
#include "base/cs_assert.h"
#include "base/cs_ast_coupling.h"
#include "base/cs_base.h"
#include "base/cs_base_accel.h"
#include "base/cs_block_dist.h"
#include "base/cs_block_to_part.h"
#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#ifdef __cplusplus
#include "base/cs_boundary_conditions_set_coeffs.h"
#endif
#include "base/cs_boundary_zone.h"
#include "base/cs_coupling.h"
#include "base/cs_defs.h"
#include "base/cs_drift_convective_flux.h"
#include "base/cs_equation_iterative_solve.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_fan.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_file.h"
#include "base/cs_file_csv_parser.h"
#include "base/cs_flag_check.h"
#include "base/cs_fp_exception.h"
#include "base/cs_function.h"
#include "base/cs_function_default.h"
#include "base/cs_gas_mix.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_head_losses.h"
#include "base/cs_ht_convert.h"
#include "base/cs_ibm.h"
#include "base/cs_interface.h"
#include "base/cs_interpolate.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_log.h"
#include "base/cs_log_iteration.h"
#include "base/cs_map.h"
#include "base/cs_mass_source_terms.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_measures_util.h"
#include "base/cs_mobile_structures.h"
#include "base/cs_notebook.h"
#include "base/cs_numbering.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_param_types.h"
#include "base/cs_parameters.h"
#include "base/cs_part_to_block.h"
#include "base/cs_physical_constants.h"
#include "base/cs_physical_properties.h"
#include "base/cs_porosity_from_scan.h"
#include "base/cs_porous_model.h"
#include "base/cs_post.h"
#include "base/cs_post_util.h"
#include "base/cs_preprocess.h"
#include "base/cs_preprocessor_data.h"
#include "base/cs_probe.h"
#include "base/cs_prototypes.h"
#include "base/cs_random.h"
#include "base/cs_renumber.h"
#include "base/cs_restart.h"
#include "base/cs_restart_map.h"
#include "base/cs_rotation.h"
#include "base/cs_runaway_check.h"
#include "base/cs_sat_coupling.h"
#include "base/cs_scalar_clipping.h"
#include "base/cs_search.h"
#include "base/cs_selector.h"
#include "base/cs_syr_coupling.h"
#include "base/cs_sys_coupling.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_moment.h"
#include "base/cs_time_plot.h"
#include "base/cs_time_step.h"
#include "base/cs_time_table.h"
#include "base/cs_timer.h"
#include "base/cs_timer_stats.h"
#include "base/cs_tree.h"
#include "base/cs_turbomachinery.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_volume_mass_injection.h"
#include "base/cs_volume_zone.h"
#include "base/cs_vof.h"
#include "base/cs_wall_condensation.h"
#include "base/cs_wall_condensation_1d_thermal.h"
#include "base/cs_wall_functions.h"
#include "base/cs_xdef_eval_at_zone.h"
#include "base/cs_zone.h"

#include "base/cs_medcoupling_remapper.h"
#include "base/cs_medcoupling_intersector.h"
#include "base/cs_medcoupling_postprocess.h"
#include "base/cs_paramedmem_remapper.h"
#include "base/cs_paramedmem_coupling.h"

/*----------------------------------------------------------------------------*/

#endif /* __CS_BASE_HEADERS_H__ */
