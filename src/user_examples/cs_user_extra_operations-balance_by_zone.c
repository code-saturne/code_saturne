/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

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
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  //!< [example_1]

  cs_balance_by_zone("all[]",
                     "temperature");

  //!< [example_1]

  //!< [example_2]

  cs_balance_by_zone("box[-0.5, 1.3, 0.0, 1.0, 1.9, 1.0]",
                     "scalar1");

  //!< [example_2]

  //!< [example_3]
  cs_real_t normal[3] = {0., 0., 1.,};

  cs_surface_balance("selection_criterion", "scalar1", normal);
  //!< [example_3]

  /* More advanced usage for pressure drop */

  {
    /*< [example_4] */
    const char criteria[] = "zone_group";

    cs_lnum_t   n_selected_cells = 0;
    cs_lnum_t  *selected_cells = NULL;

    cs_real_t balance[CS_BALANCE_P_N_TERMS];

    BFT_MALLOC(selected_cells, domain->mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_cells,
                              selected_cells);

    cs_balance_by_zone_compute("scalar1",
                               n_selected_cells,
                               selected_cells,
                               balance);

    BFT_FREE(selected_cells);

    cs_balance_term_t  mass_in_idx = CS_BALANCE_MASS_IN;
    cs_balance_term_t  mass_out_idx = CS_BALANCE_MASS_OUT;

    bft_printf("inlet mass flow  (scalar 1): %g\n"
               "outlet mass flow (scalar 1): %g\n",
               balance[mass_in_idx],
               balance[mass_out_idx]);
    /*< [example_4] */
  }
  //!< [example_5]

  cs_pressure_drop_by_zone("zone_group");

  //!< [example_5]

  /* More advanced usage for pressure drop */

  //!< [example_6]
  {
    const char criteria[] = "zone_group";

    cs_lnum_t   n_selected_cells = 0;
    cs_lnum_t  *selected_cells = NULL;

    cs_real_t balance[CS_BALANCE_P_N_TERMS];

    BFT_MALLOC(selected_cells, domain->mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_cells,
                              selected_cells);

    cs_pressure_drop_by_zone_compute(n_selected_cells,
                                     selected_cells,
                                     balance);

    BFT_FREE(selected_cells);

    cs_balance_p_term_t  rhou_in_idx = CS_BALANCE_P_RHOU_IN;
    cs_balance_p_term_t  rhou_out_idx = CS_BALANCE_P_RHOU_OUT;

    bft_printf("inlet mass flow  (rho.u): %g\n"
               "outlet mass flow (rho.u): %g\n",
               balance[rhou_in_idx],
               balance[rhou_out_idx]);
  }
  //!< [example_6]
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
