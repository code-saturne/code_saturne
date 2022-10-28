/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  /*! [loc_var_dec] */
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)domain->mesh_quantities->b_face_normal;
  /*! [loc_var_dec] */

  /* Assign boundary conditions to boundary faces here
   *
   * For each subset:
   * - find the matchig zone to filter boundary faces of a given subset
   * - loop on faces from a zone
   *   - set the boundary condition for each face
   *---------------------------------------------------------------------------*/

  /*! [init] */
  int nbelec = cs_glob_transformer->nbelec;
  int nbtrf  = cs_glob_transformer->nbtrf;

  cs_data_joule_effect_t *transfo = cs_get_glob_transformer();

  cs_real_t *sir, *sii, *sirt, *siit;
  cs_real_6_t *sirb, *siib, *ur, *ui;
  int *nborne;
  char name[8];

  BFT_MALLOC(sir, nbelec, cs_real_t);
  BFT_MALLOC(sii, nbelec, cs_real_t);

  BFT_MALLOC(sirt, nbtrf, cs_real_t);
  BFT_MALLOC(siit, nbtrf, cs_real_t);

  BFT_MALLOC(siib, nbtrf, cs_real_6_t);
  BFT_MALLOC(sirb, nbtrf, cs_real_6_t);
  BFT_MALLOC(ui, nbtrf, cs_real_6_t);
  BFT_MALLOC(ur, nbtrf, cs_real_6_t);

  BFT_MALLOC(nborne, nbtrf, int);
  /*! [init] */

  /* 1 - Computation of intensity (A/m2) for each electrode */

  /*! [pre_init] */
  for (int i = 0; i < nbelec; i++) {
    sir[i] = 0.;
    sii[i] = 0.;
  }

  for (int i = 0; i < nbtrf; i++) {
    sirt[i] = 0.;
    siit[i] = 0.;
    nborne[i] = 0;
  }

  if (cs_glob_time_step->nt_cur < (cs_glob_time_step->nt_prev + 2)) {
    for (int i = 0; i < nbtrf; i++) {
      transfo->uroff[i] = 0.;
      transfo->uioff[i] = 0.;
    }
  }
  /*! [pre_init] */

  /*! [step_1] */
  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];

  for (int i = 0; i < nbelec; i++) {

    sprintf(name, "%07d", transfo->ielecc[i]);

    const cs_zone_t *z = cs_boundary_zone_by_name(name);

    cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
    cs_real_3_t *cpro_curim = NULL;
    if (ieljou == 4)
      cpro_curim = (cs_real_3_t *)(CS_F_(curim)->val);

    for (cs_lnum_t ilelt = 0; ilelt < z->n_elts; ilelt++) {

      cs_lnum_t face_id = z->elt_ids[ilelt];
      cs_lnum_t cell_id = b_face_cells[face_id];

      for (cs_lnum_t id = 0; id < 3; id++)
        sir[i] += cpro_curre[cell_id][id] * b_face_normal[id][face_id];

      if (ieljou == 4)
        for (cs_lnum_t id = 0; id < 3; id++)
          sii[i] += cpro_curim[cell_id][id] * b_face_normal[id][face_id];
    }

  }
  /*! [step_1] */

  /* 2 - Definition of Voltage on each termin of transformers */
  /* 2.1 Computation of Intensity on each termin of transformers */

  /*! [step_2_1] */
  for (int i = 0; i < nbelec; i++) {
    sirb[transfo->ielect[i]][transfo->ielecb[i]] = 0.;
    if (ieljou == 4)
      siib[transfo->ielect[i]][transfo->ielecb[i]] = 0.;
  }

  for (int i = 0; i < nbelec; i++) {
    if (transfo->ielect[i] != 0) {
      sirb[transfo->ielect[i]][transfo->ielecb[i]] += sir[i];
      if (ieljou == 4)
        siib[transfo->ielect[i]][transfo->ielecb[i]] += sii[i];
    }
  }
  /*! [step_2_1] */

  /* 2.2 RVoltage on each termin */

  /*! [step_2_2] */
  for (int ntf = 0; ntf < nbtrf; ntf++) {
    /* Primary and Secondary in Triangle */
    if (transfo->ibrpr[ntf] == 0 &&
        transfo->ibrsec[ntf] == 0) {
      nborne[ntf] = 3;
      cs_real_t rnbs2 = 3. * transfo->rnbs[ntf]
                           * transfo->rnbs[ntf];
      ur[ntf][0] = 1.154675 * transfo->tenspr[ntf]
                            / transfo->rnbs[ntf]
                 + (transfo->zr[ntf] * sirb[ntf][0]
                 -  transfo->zi[ntf] * siib[ntf][0]) / rnbs2;
      ur[ntf][1] = -0.5773 * transfo->tenspr[ntf]
                           / transfo->rnbs[ntf]
                 + (transfo->zr[ntf] * sirb[ntf][1]
                 -  transfo->zi[ntf] * siib[ntf][1]) / rnbs2;
      ur[ntf][2] =-0.5773 * transfo->tenspr[ntf]
                           / transfo->rnbs[ntf]
                 + (transfo->zr[ntf] * sirb[ntf][2]
                 -  transfo->zi[ntf] * siib[ntf][2]) / rnbs2;

      ui[ntf][0] = (transfo->zi[ntf] * sirb[ntf][0]
                 -  transfo->zr[ntf] * siib[ntf][0]) / rnbs2;
      ui[ntf][1] = (transfo->zi[ntf] * sirb[ntf][1]
                 -  transfo->zr[ntf] * siib[ntf][1]) / rnbs2;
      ui[ntf][2] = (transfo->zi[ntf] * sirb[ntf][2]
                 -  transfo->zr[ntf] * siib[ntf][2]) / rnbs2;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
              _("electric module : \n"
                "transformer matrix not defined\n"));
  }
  /*! [step_2_2] */

  /* 2.3 Total intensity for a transformer
   *     (zero valued WHEN Offset established) */

  /*! [step_2_3] */
  for (int ntf = 0; ntf < nbtrf; ntf++) {
    sirt[ntf] = 0.;
    if (ieljou == 4)
      siit[ntf] = 0.;
  }

  for (int i = 0; i < nbelec; i++) {
    if (transfo->ielect[i] != 0) {
      sirt[i] += sir[i];
      if (ieljou == 4)
        siit[i] += sii[i];
    }
  }
  /*! [step_2_3] */

  /* 2.4 Take in account of Offset */

  /*! [step_2_4] */
  cs_real_t capaeq = 3.;

  for (int ntf = 0; ntf < nbtrf; ntf++) {
    transfo->uroff[ntf] += sirt[ntf] / capaeq;
    if (ieljou == 4)
      transfo->uioff[ntf] += siit[ntf] / capaeq;
  }

  /* A reference transformer is assumed to have an Offset zero valued */
  if (transfo->ntfref > 0) {
    transfo->uroff[transfo->ntfref] = 0.;
    transfo->uioff[transfo->ntfref] = 0.;
  }

  for (int ntf = 0; ntf < nbtrf; ntf++) {
    for (int nb = 0; nb < nborne[ntf]; nb++) {
      ur[ntf][nb] += transfo->uroff[ntf];
      if (ieljou == 4)
        ui[ntf][nb] += transfo->uioff[ntf];
    }
  }

  /* Print of UROFF (real part of offset potential) */
  bft_printf(" ** INFORMATION ON TRANSFORMERS\n"
             "    ---------------------------------------\n"
             "\n"
             "      ---------------------------------\n"
             "      Number of Transformers   UROFF\n"
             "      ---------------------------------\n");
  for (int ntf = 0; ntf < nbtrf; ntf++)
    bft_printf("          %6i            %12.5E\n", ntf, transfo->uroff[ntf]);
  bft_printf("    ---------------------------------------\n");
  /*! [step_2_4] */

  /* 2.5 Take in account of Boundary Conditions */

  /*! [step_2_5] */

  int       *potr_icodcl  = CS_F_(potr)->bc_coeffs->icodcl;
  cs_real_t *potr_rcodcl1 = CS_F_(potr)->bc_coeffs->rcodcl1;
  cs_real_t *potr_rcodcl3 = CS_F_(potr)->bc_coeffs->rcodcl3;

  int       *poti_icodcl  = NULL;
  cs_real_t *poti_rcodcl1 = NULL;
  cs_real_t *poti_rcodcl3 = NULL;

  if (ieljou == 4) {
    poti_icodcl  = CS_F_(potr)->bc_coeffs->icodcl;
    poti_rcodcl1 = CS_F_(potr)->bc_coeffs->rcodcl1;
    poti_rcodcl3 = CS_F_(potr)->bc_coeffs->rcodcl3;
  }

  for (int i = 0; i < nbelec; i++) {

    sprintf(name, "%07d", transfo->ielecc[i]);

    const cs_zone_t *z = cs_boundary_zone_by_name(name);

    for (cs_lnum_t ilelt = 0; ilelt < z->n_elts; ilelt++) {
      cs_lnum_t face_id = z->elt_ids[ilelt];

      bc_type[face_id] = CS_SMOOTHWALL;

      if (transfo->ielect[i] != 0) {
        potr_icodcl[face_id] = 1;
        potr_rcodcl1[face_id] = ur[transfo->ielect[i]][transfo->ielecb[i]];

        if (ieljou == 4) {
          poti_icodcl[face_id] = 1;
          poti_rcodcl1[face_id] = ur[transfo->ielect[i]][transfo->ielecb[i]];
        }
      }
      else {
        potr_icodcl[face_id] = 3;
        potr_rcodcl3[face_id] = 0.;

        if (ieljou == 4) {
          poti_icodcl[face_id] = 3;
          poti_rcodcl3[face_id] = 0.;
        }
      }
    }

  }

  /* Test, if not any reference transformer
   *       a piece of wall may be at ground. */
  if (transfo->ntfref == 0) {
    int found = 0;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (bc_type[face_id] == CS_SMOOTHWALL) {
        if (potr_icodcl[face_id] == 1) {
          if (ieljou == 3) {
            if (fabs(potr_rcodcl1[face_id]) < 1.e-20)
              found = 1;
          }
          else if (ieljou == 4) {
            cs_real_t val = fabs(potr_rcodcl1[face_id]);
            if (fabs(poti_rcodcl1[face_id]) < 1.e-20 && val < 1.e-20)
              found = 1;
          }
        }
      }
    }
    cs_parall_max(1, CS_INT_TYPE, &found);
    if (!found)
      bft_error(__FILE__, __LINE__, 0,
              _("ERROR in JOULE : \n"
                "Lack of reference: choose a transformer for which\n"
                "offset is assumed zero or a face at ground on the\n"
                "boundary."));
  }
  /*! [step_2_5] */

  /*! [step_3] */
  BFT_FREE(sir);
  BFT_FREE(sii);
  BFT_FREE(sirb);
  BFT_FREE(siib);
  BFT_FREE(ur);
  BFT_FREE(ui);
  BFT_FREE(sirt);
  BFT_FREE(siit);
  BFT_FREE(nborne);
  /*! [step_3] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
