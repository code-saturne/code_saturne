
/*============================================================================
 * Methods for head losses due to particle deposit
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

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Functions dealing with head losses due to particle deposit
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_gradient.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_parall.h"
#include "cs_prototypes.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_stat.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_head_losses.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the porosity in wall-normal cells based on the mean deposit height
 * (which is only known at the boundary faces)
 *----------------------------------------------------------------------------*/

static void
_porcel(cs_real_t         mdiam[],
        cs_real_t         porosi[],
        const cs_lnum_t   bc_type[])
{
  /* Initialization */

  const cs_mesh_t  *m = cs_glob_mesh;

  cs_lnum_t nfabor = m->n_b_faces;
  cs_lnum_t ncel   = m->n_cells;
  cs_lnum_t ncelet = m->n_cells_with_ghosts;

  cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_surf;
  cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;

  cs_lagr_boundary_interactions_t
    *lag_bdy_i = cs_glob_lagr_boundary_interactions;

  /*Compute distance to deposition wall, and check deposition exists on
    at least one face.*/

  int indic = 0;

  int  *itytmp;
  BFT_MALLOC(itytmp, nfabor, int);

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    if (   (bc_type[ifac] == CS_SMOOTHWALL || bc_type[ifac] == CS_ROUGHWALL)
        && bound_stat[ifac + lag_bdy_i->ihdepm * nfabor] > 0.0) {
      indic = 1;
      itytmp[ifac] = CS_INLET;
    }
    else
      itytmp[ifac] = CS_INDEF;
  }

  cs_parall_counter_max(&indic, 1);

  /* mean particle diameter and porosity value due to
     the deposit in each cell*/
  for (cs_lnum_t iel = 0; iel < ncelet; iel++) {
    porosi[iel] = 1.0;
    mdiam[iel]  = 0.0;
  }

  /* Nothing more to do if no deposition at this stage */
  if (indic == 0) {
    BFT_FREE(itytmp);
    return;
  }

  /* Allocate temporary arrays for the distance resolution */
  cs_real_t    *distpw, *coefap, *coefbp;
  cs_real_3_t  *q;
  cs_real_t    *masflu, *depvol;

  BFT_MALLOC(distpw, ncelet, cs_real_t);
  BFT_MALLOC(coefap, nfabor, cs_real_t);
  BFT_MALLOC(coefbp, nfabor, cs_real_t);
  BFT_MALLOC(q, ncelet, cs_real_3_t);
  BFT_MALLOC(masflu, ncelet, cs_real_t);
  BFT_MALLOC(depvol, ncelet, cs_real_t);

  if (indic > 0)
    CS_PROCF(distpr, DISTPR)(itytmp, distpw);

  /* Compute  n = Grad(DISTPW)/|Grad(DISTPW)|
   * ======================================== */

  /* Distance to wall is 0 at the wall by definition, zero flux elsewhere   */

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    if (itytmp[ifac] == CS_INLET) {
      /* Dirichlet Boundary Condition for gradients */
      coefap[ifac] = 0.0;
      coefbp[ifac] = 0.0;
    }
    else {
      /* Neumann Boundary Condition for gradients */
      coefap[ifac] = 0.0;
      coefbp[ifac] = 1.0;
    }

  }

  /* Compute the gradient of the distance to the wall */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_scalar("Work array",
                     gradient_type,
                     halo_type,
                     1,          /* inc */
                     1,          /* recompute_cocg */
                     100,        /* nswrgy */
                     0,          /* idimtr */
                     0,          /* hyd_p_flag */
                     1,          /* w_stride */
                     2,          /* iwarny */
                     -1,         /* imligy */
                     1.e-5,      /* epsrgy */
                     0.,         /* extray */
                     1.5,        /* climgy */
                     NULL,
                     coefap,
                     coefbp,
                     distpw,
                     NULL,
                     NULL, /* internal coupling */
                     q);

  /* Normalisation (caution, gradient can be zero sometimes) */

  for (cs_lnum_t iel = 0; iel < ncel; iel++) {

    cs_real_t xnorme = CS_MAX(cs_math_3_norm(q[iel]),
                              cs_math_epzero);

    for (cs_lnum_t isou = 0; isou < 3; isou++)
      q[iel][isou] /= xnorme;

  }

  /* Paralellism and periodicity */

  cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD, (cs_real_t *)q, 3);

  if (m->n_init_perio > 0)
    cs_halo_perio_sync_var_vect(m->halo, CS_HALO_STANDARD,
                                (cs_real_t *)q, 3);

  /* Compute  porosity
   * ================= */

  cs_real_t epsi = 1e-06;

  /* time step  */
  cs_real_t dt = 0.0001;

  /* volume of the deposited particles
     (calculated from mean deposit height)*/
  for (cs_lnum_t iel = 0; iel < ncelet; iel++)
    depvol[iel]    = 0.0;

  indic = 0;
  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    cs_lnum_t iel = ifabor[ifac];
    depvol[iel] += bound_stat[lag_bdy_i->ihdepm * nfabor + ifac] * surfbn[ifac];
    mdiam[iel]  += bound_stat[lag_bdy_i->ihdiam * nfabor + ifac];

  }

  const cs_real_t fmult = (1.0 - cs_glob_lagr_clogging_model->mporos);

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    cs_lnum_t iel = ifabor[ifac];
    porosi[iel] =   (volume[iel] - fmult * depvol[iel]) / volume[iel];

    if (porosi[iel] < cs_glob_lagr_clogging_model->mporos)
      indic = 1;

  }

  /* Paralellism and periodicity */

  cs_parall_counter_max(&indic, 1);

  cs_lnum_t nn = 0;
  while (indic > 0) {

    for (cs_lnum_t iel = 0; iel < ncelet; iel++)
      masflu[iel]  = 0.0;

    for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_i_faces; ifac++) {

      cs_real_t prod1 = 0.0;
      cs_real_t prod2 = 0.0;
      cs_lnum_t iel1  = cs_glob_mesh->i_face_cells[ifac][0];
      cs_lnum_t iel2  = cs_glob_mesh->i_face_cells[ifac][1];

      for (cs_lnum_t ii = 0; ii < 3; ii++) {

        prod1 += q[iel1][ii] * cs_glob_mesh_quantities->i_face_normal[ifac * 3 + ii];
        prod2 -= q[iel2][ii] * cs_glob_mesh_quantities->i_face_normal[ifac * 3 + ii];

      }

      if (porosi[iel1] >= cs_glob_lagr_clogging_model->mporos || prod1 <= epsi)
        masflu[iel1] = masflu[iel1];

      else {

        masflu[iel1] -=   (porosi[iel1] - cs_glob_lagr_clogging_model->mporos)
                        / dt * volume[iel1];
        masflu[iel2] +=   (porosi[iel1] - cs_glob_lagr_clogging_model->mporos)
                        / dt * volume[iel1];
        mdiam[iel2]   = mdiam[iel1];

      }

      if (porosi[iel2] >= cs_glob_lagr_clogging_model->mporos || prod2 <= epsi)
        masflu[iel2] = masflu[iel2];

      else {

        masflu[iel2] -=   (porosi[iel2] - cs_glob_lagr_clogging_model->mporos)
                        / dt * volume[iel2];
        masflu[iel1] +=   (porosi[iel2] - cs_glob_lagr_clogging_model->mporos)
                        / dt * volume[iel2];
        mdiam[iel1]   = mdiam[iel2];

      }

    }
    indic     = 0;

    for (cs_lnum_t iel = 0; iel < ncel; iel++)
      porosi[iel]  = porosi[iel] + (dt / volume[iel]) * masflu[iel];


    /* Paralellism and periodicity    */
    if (cs_glob_rank_id >= 0 || cs_glob_mesh->n_init_perio > 0) {

      cs_mesh_sync_var_scal(porosi);
      cs_mesh_sync_var_scal(mdiam);

    }

    for (cs_lnum_t iel = 0; iel < ncel; iel++) {

      if (porosi[iel] < cs_glob_lagr_clogging_model->mporos)
        indic = 1;

    }

    cs_parall_counter_max(&indic, 1);

    nn++;

    if (nn >= 100)
      bft_error(__FILE__, __LINE__, 0,
                "==========================================\n"
                " Error nn > 100\n"
                "\n"
                " Stop inside the porcel subroutine\n"
                "==========================================\n");

  }

  /* Free memory */

  BFT_FREE(masflu);
  BFT_FREE(depvol);
  BFT_FREE(q);
  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(distpw);
  BFT_FREE(itytmp);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Define Head losses to take into account deposit in the flow
 *
 * \param[in]   n_hl_cells  number of cells on which to apply head losses
 * \param[in]   cell_ids    ids of cells on which to apply head losses
 * \param[in]   bc_type     boundary face type
 * \param[out]  cku         head loss coefficients at matchin cells
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_head_losses(cs_lnum_t        n_hl_cells,
                    const cs_lnum_t  cell_ids[],
                    const cs_lnum_t  bc_type[],
                    cs_real_t        cku[][6])
{
  cs_lnum_t ncel   = cs_glob_mesh->n_cells;
  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;

  cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  /* Check that head loss zone definitions are consistent */
  if (n_hl_cells != ncel)
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("in Lagrangian module"),
       _("The number of cells in the head loss zones must cover\n"
         "the whole mesh (though the local head loss may be zero).\n"));

  /* ==================================================================
   * cku: compute head loss coefficients in the calculation coordinates,
   *            organized in order k11, k22, k33, k12, k13, k23
   * Note:
   *     - make sure diagonal coefficients are positive. The calculation
   *      may crash if this is not the case, and no further check will
   *      be done
   * ==================================================================*/

  /* ===================================================================
   * Porosity calculation for the influence of the deposit on the flow
   * by head losses
   * ====================================================================*/

  cs_real_t *mdiam;
  BFT_MALLOC(mdiam, ncelet, cs_real_t);

  /* cs_lnum_t poro_id; */
  /* field_get_id_try ("clogging_porosity", &poro_id); */

  /* TODO this field is never built */
  cs_field_t *f_poro = cs_field_by_name_try("clogging_porosity");

  cs_real_t *lporo;
  if (f_poro == NULL)
    BFT_MALLOC(lporo, ncelet, cs_real_t);
  else
    lporo = f_poro->val;

  _porcel(mdiam, lporo, bc_type);

  /* Calculation of the head loss term with the Ergun law
   * mdiam :  mean diameter of deposited particles
   * lcell :  characteristic length in the flow direction */

  for (cs_lnum_t ielpdc = 0; ielpdc < n_hl_cells; ielpdc++) {

    cs_lnum_t iel = cell_ids[ielpdc];

    if (mdiam[iel] > 0.0) {

      cs_real_t lcell = pow(volume[iel], 1.0 / 3.0);
      cs_real_t romf  = extra->cromf->val[iel];
      cs_real_t visccf = extra->viscl->val[iel] / romf;
      cs_real_t v      = sqrt (  pow(extra->vel->vals[1][iel * 3 + 0], 2)
                               + pow(extra->vel->vals[1][iel * 3 + 1], 2)
                               + pow(extra->vel->vals[1][iel * 3 + 2], 2));

      cs_real_t ck =  v * 1.75 * (1 - lporo[iel]) / pow (lporo[iel], 3.0)
                    * lcell / mdiam[iel]
                    +  (lcell * 150.0 * visccf) / (romf * pow(mdiam[iel], 2))
                     * pow((1 - lporo[iel]), 2) / lporo[iel] * 3;
      cku[iel][0] = ck;
      cku[iel][1] = ck;
      cku[iel][2] = ck;
      cku[iel][3] = 0.0;
      cku[iel][4] = 0.0;
      cku[iel][5] = 0.0;

    }

  }

  if (f_poro == NULL)
    BFT_FREE(lporo);

  BFT_FREE(mdiam);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
