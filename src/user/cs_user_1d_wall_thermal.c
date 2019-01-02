/*============================================================================
 * Data Entry of the 1D wall thermal module.
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
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_1d_wall_thermal.h"
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
 * \file cs_user_1d_wall_thermal.c
 *
 * \brief Data Entry of the 1D wall thermal module.
 *
 * \param[in]   iappel   Call number:
 *                       - 1: first call during the initialization (called once)
 *                       Setting the number of cells nfpt1d.
 *                       - 2: second call during the initialization (called once)
 *                       Filling ifpt1d, nppt1d, eppt1d and rgpt1d arrays.
 *                       - 3: called at each time step
 *                       Filling iclt1d, xlmbt1, rcpt1d and dtpt1d arrays.

 * \param[in]   isuit1   Rereading of the restart file:
 *                       - 0: No rereading
 *                            (meshing and wall temperature reinitialization)
 *                       - 1: Rereading of the restart file for the 1-Dimension
 *                            thermal module
 *                       - isuite: Rereading only if the computational fluid
 *                            dynamic is a continuation of the computation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

void
cs_user_1d_wall_thermal(int iappel,
                        int isuit1)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
