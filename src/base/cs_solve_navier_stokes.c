/*============================================================================
 * Solve the Navier-Stokes equations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_assert.h"
#include "cs_atmo.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solve_navier_stokes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

void
cs_f_navier_stokes_total_pressure(void);

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_navier_stokes_total_pressure(void)
{
  cs_solve_navier_stokes_update_total_pressure(cs_glob_mesh,
                                               cs_glob_mesh_quantities,
                                               cs_glob_fluid_properties);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update total pressure (defined as a post-processed property).
 *
 * For the compressible module, the solved pressure is already
 * the total pressure.
 *
 * Note: for Eddy Viscosity Models, the TKE may be included in the
 * solved pressure.
 *
 * \param[in]     m   pointer to mesh structure
 * \param[in]     mq  pointer to mesh quantities structure
 * \param[in]     fp  pointer to fluid properties structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_navier_stokes_update_total_pressure(const cs_mesh_t              *m,
                                             const cs_mesh_quantities_t   *mq,
                                             const cs_fluid_properties_t  *fp)
{
  /* TODO: use a function pointer here to adapt to different cases */

  cs_field_t *f = cs_field_by_name_try("total_pressure");

  if ((CS_F_(p) == NULL) || (f == NULL))
    return;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  const cs_real_t *xyzp0 = fp->xyzp0;
  const cs_real_t p0 = fp->p0, pred0 = fp->pred0, ro0 = fp->ro0;

  cs_real_t *cpro_prtot = f->val;
  const cs_real_t *cvar_pr = CS_F_(p)->val;

  const cs_real_3_t *cpro_momst = NULL;

  if (cs_glob_atmo_option->open_bcs_treatment != 0)
    cpro_momst
      = (const cs_real_3_t *)cs_field_by_name("momentum_source_terms")->val;

  /* Update cell values */

# pragma omp parallel if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

    if (cpro_momst == NULL) {
      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {
        cpro_prtot[c_id] =   cvar_pr[c_id]
                           + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                  cell_cen[c_id],
                                                                  gxyz)
                           + p0 - pred0;
      }
    }
    else {
      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
        cpro_prtot[c_id] =   cvar_pr[c_id]
                           + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                  cell_cen[c_id],
                                                                  gxyz)
                           + p0 - pred0
                           - cs_math_3_distance_dot_product(xyzp0,
                                                            cell_cen[c_id],
                                                            cpro_momst[c_id]);
    }

    /* For Eddy Viscosity Models, "2/3 rho k"
       is included in the solved pressure */

    if (  (   cs_glob_turb_model->itytur == 2
           || cs_glob_turb_model->itytur == 5
           || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA)
        && cs_glob_turb_rans_model->igrhok != 1) {

      const cs_real_t *cvar_k = CS_F_(k)->val;
      const cs_real_t *cpro_rho = CS_F_(rho)->val;

      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
        cpro_prtot[c_id] -= 2.0/3 * cpro_rho[c_id]*cvar_k[c_id];
    }

  } /* cell values update */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
