/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* Code_Saturne version 6.0-alpha */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_internal_coupling.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_param_cdo.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(void)
{
  /* Add boundary values for all scalars */
  {
    int n_fields = cs_field_n_fields();
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE)
        cs_parameters_add_boundary_values(f);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
