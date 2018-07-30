/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_elec_model.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_time_step.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in]     nvar          total number of variable BC's
 * \param[in]     bc_type       boundary face types
 * \param[in]     icodcl        boundary face code
 *                                - 1  -> Dirichlet
 *                                - 2  -> convective outlet
 *                                - 3  -> flux density
 *                                - 4  -> sliding wall and u.n=0 (velocity)
 *                                - 5  -> friction and u.n=0 (velocity)
 *                                - 6  -> roughness and u.n=0 (velocity)
 *                                - 9  -> free inlet/outlet (velocity)
 *                                inflowing possibly blocked
 * \param[in]     rcodcl        boundary condition values
 *                                rcodcl(3) = flux density value
 *                                (negative for gain) in W/m2
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         bc_type[],
                            int         icodcl[],
                            cs_real_t   rcodcl[])
{
  /*! [loc_var_dec] */
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *) cs_glob_mesh_quantities->b_face_normal;

  const int keyvar = cs_field_key_id("variable_id");
  /*! [loc_var_dec] */

  /* ===============================================================================
   * Assign boundary conditions to boundary faces here
   *
   * For each subset:
   * - use selection criteria to filter boundary faces of a given subset
   * - loop on faces from a subset
   *   - set the boundary condition for each face
   * =============================================================================== */

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

    cs_lnum_t nelts;
    cs_lnum_t *lstelt = NULL;

    sprintf(name, "%07d", transfo->ielecc[i]);
    cs_selector_get_b_face_list(name, &nelts, lstelt);

    cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
    cs_real_3_t *cpro_curim = NULL;
      if (ieljou == 4)
        cpro_curim = (cs_real_3_t *)(CS_F_(curim)->val);

    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {

      cs_lnum_t face_id = lstelt[ilelt];
      cs_lnum_t cell_id = b_face_cells[face_id];

      for (cs_lnum_t id = 0; id < 3; id++)
        sir[i] += cpro_curre[cell_id][id] * b_face_normal[id][face_id];

      if (ieljou == 4)
        for (cs_lnum_t id = 0; id < 3; id++)
          sii[i] += cpro_curim[cell_id][id] * b_face_normal[id][face_id];
    }

    BFT_FREE(lstelt);
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
      double rnbs2 = 3. * transfo->rnbs[ntf]
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
  double capaeq = 3.;

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
  bft_printf(" ** INFORMATIONS ON TRANSFOMERS\n");
  bft_printf("    ---------------------------------------\n");
  bft_printf("      ---------------------------------\n");
  bft_printf("      Number of Transfo        UROFF\n");
  bft_printf("      ---------------------------------\n");
  for (int ntf = 0; ntf < nbtrf; ntf++)
    bft_printf("          %6i            %12.5E\n", ntf, transfo->uroff[ntf]);
  bft_printf("    ---------------------------------------\n");
  /*! [step_2_4] */

  /* 2.5 Take in account of Boundary Conditions */

  /*! [step_2_5] */
  for (int i = 0; i < nbelec; i++) {

    cs_lnum_t nelts;
    cs_lnum_t *lstelt = NULL;

    sprintf(name, "%07d", transfo->ielecc[i]);
    cs_selector_get_b_face_list(name, &nelts, lstelt);

    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
      cs_lnum_t face_id = lstelt[ilelt];

      bc_type[face_id] = CS_SMOOTHWALL;

      if (transfo->ielect[i] != 0) {
        const cs_field_t *f = CS_F_(potr);
        int ivar = cs_field_get_key_int(f, keyvar) - 1;
        icodcl[ivar * n_b_faces + face_id] = 1;
        rcodcl[ivar * n_b_faces + face_id] = ur[transfo->ielect[i]][transfo->ielecb[i]];

        if (ieljou == 4) {
          f = CS_F_(poti);
          ivar = cs_field_get_key_int(f, keyvar) - 1;
          icodcl[ivar * n_b_faces + face_id] = 1;
          rcodcl[ivar * n_b_faces + face_id] = ur[transfo->ielect[i]][transfo->ielecb[i]];
        }
      }
      else {
        const cs_field_t *f = CS_F_(potr);
        int ivar = cs_field_get_key_int(f, keyvar) - 1;
        icodcl[ivar * n_b_faces + face_id] = 3;
        rcodcl[2 * n_b_faces * nvar + ivar * n_b_faces + face_id] = 0.;

        if (ieljou == 4) {
          f = CS_F_(poti);
          ivar = cs_field_get_key_int(f, keyvar) - 1;
          icodcl[ivar * n_b_faces + face_id] = 3;
          rcodcl[2 * n_b_faces * nvar + ivar * n_b_faces + face_id] = 0.;
        }
      }
    }

    BFT_FREE(lstelt);
  }

  /* Test, if not any reference transformer
   *       a piece of wall may be at ground. */
  if (transfo->ntfref == 0) {
    bool found = false;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (bc_type[face_id] == CS_SMOOTHWALL) {
        const cs_field_t *f = CS_F_(potr);
        int ivar = cs_field_get_key_int(f, keyvar) - 1;
        if (icodcl[ivar * n_b_faces + face_id] == 1) {
          if (ieljou == 3) {
            if (fabs(rcodcl[ivar * n_b_faces + face_id]) < 1.e-20)
              found = true;
          }
          else if (ieljou == 4) {
            double val = fabs(rcodcl[ivar * n_b_faces + face_id]);
            f = CS_F_(poti);
            ivar = cs_field_get_key_int(f, keyvar) - 1;
            if (fabs(rcodcl[ivar * n_b_faces + face_id]) < 1.e-20 &&
                val < 1.e-20)
              found = true;
          }
        }
      }
    }
    if (!found)
      bft_error(__FILE__, __LINE__, 0,
              _("ERROR in JOULE : \n"
                "Lack of reference : choose a transformer for wich\n"
                "offset is assumed zero or a face at ground on the\n"
                "boundary\n"));
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
