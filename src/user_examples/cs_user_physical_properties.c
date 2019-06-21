/*============================================================================
 * This function is called each time step to define physical properties
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

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
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t *domain)
{

  /* Check fields exists */
  if (CS_F_(lambda) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error lambda not variable\n"));
  if (CS_F_(rho) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error rho not variable\n"));
  if (CS_F_(cp) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error cp not variable\n"));

  cs_real_t *cpro_lambda = CS_F_(lambda)->val;
  cs_real_t *cpro_rho = CS_F_(rho)->val;
  cs_real_t *cpro_cp = CS_F_(cp)->val;

  /* Impose thermal conductivity, density and specific heat for solid zones */
  {
    /* Volume zone "CM" must be defined in the GUI or in cs_user_zones.c */
    const cs_zone_t *z = cs_volume_zone_by_name("CM");

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      cs_lnum_t cell_id = z->elt_ids[i];
      cpro_lambda[cell_id] = 3.3;
      cpro_rho[cell_id] = 7800.;
      cpro_cp[cell_id] = 444.;
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
