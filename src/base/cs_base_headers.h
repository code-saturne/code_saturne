#ifndef __CS_BASE_HEADERS_H__
#define __CS_BASE_HEADERS_H__

/*============================================================================
 * Global Code_Saturne headers file for easier include
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_1d_wall_thermal.h"
#include "cs_1d_wall_thermal_check.h"
#include "cs_ale.h"
#include "cs_all_to_all.h"
#include "cs_array_reduce.h"
#include "cs_assert.h"
#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_block_to_part.h"
#include "cs_boundary.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_coupling.h"
#include "cs_defs.h"
#include "cs_equation_iterative_solve.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_flag_check.h"
#include "cs_fp_exception.h"
#include "cs_gas_mix.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_head_losses.h"
#include "cs_interface.h"
#include "cs_interpolate.h"
#include "cs_internal_coupling.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_math.h"
#include "cs_notebook.h"
#include "cs_numbering.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_properties.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_preprocess.h"
#include "cs_preprocessor_data.h"
#include "cs_probe.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_restart.h"
#include "cs_restart_map.h"
#include "cs_rotation.h"
#include "cs_sat_coupling.h"
#include "cs_selector.h"
#include "cs_stokes_model.h"
#include "cs_time_moment.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"
#include "cs_tree.h"
#include "cs_turbomachinery.h"
#include "cs_volume_zone.h"
#include "cs_wall_functions.h"
#include "cs_zone.h"

#if defined(HAVE_MEDCOUPLING_LOADER)
#include "cs_medcoupling_remapper.hxx"
#endif

#if defined(HAVE_PARAMEDMEM)
#include "cs_paramedmem_coupling.hxx"
#endif

/*----------------------------------------------------------------------------*/

#endif /* __CS_BASE_HEADERS_H__ */
