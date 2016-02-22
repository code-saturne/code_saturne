/*============================================================================
 * Boundary condition for transformers.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------*/

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

#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_elec_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


void
CS_PROCF (usclim_trans, UICLIM_TRANS)(cs_int_t   *icodcl,
                                      cs_real_t  *rcodcl,
                                      cs_int_t   *itypfb,
                                      cs_int_t   *izfppp,
                                      cs_int_t   *iparoi,
                                      cs_int_t   *nvarcl)
{
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;
  const cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;
  const cs_real_3_t *surfbo = (const cs_real_3_t *) cs_glob_mesh_quantities->b_face_normal;
  cs_lnum_t *lstelt = NULL;
  cs_lnum_t  nelts;
  int isca;
  cs_field_t *f;

  const int keyvar = cs_field_key_id("variable_id");

  /* ===============================================================================
   * Assign boundary conditions to boundary faces here
   *
   * For each subset:
   * - use selection criteria to filter boundary faces of a given subset
   * - loop on faces from a subset
   *   - set the boundary condition for each face
   * =============================================================================== */

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

  /* 1 - Computation of intensity (A/m2) for each electrode */
  for (int i = 0; i < nbelec; i++) {
    sir[i] = 0.;
    sii[i] = 0.;
  }

  for (int i = 0; i < nbtrf; i++) {
    sirt[i] = 0.;
    siit[i] = 0.;
    nborne[i] = 0;
  }

  if (cs_glob_time_step->nt_cur < (cs_glob_time_step->nt_prev + 2))
    for (int i = 0; i < nbtrf; i++) {
      transfo->uroff[i] = 0.;
      transfo->uioff[i] = 0.;
    }

  for (int i = 0; i < nbelec; i++) {
    sprintf(name, "%07d", transfo->ielecc[i]);
    cs_selector_get_b_face_num_list(name, &nelts, lstelt);

    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
      int ifac = lstelt[ilelt];
      int iel = ifabor[iel];

      for (int id = 0; id < 3; id++)
        sir[i] += CS_FI_(curre, id)->val[iel] * surfbo[id][ifac];

      if (cs_glob_elec_option->ieljou == 4)
        for (int id = 0; id < 3; id++)
          sii[i] += CS_FI_(curim, id)->val[iel] * surfbo[id][ifac];
    }

    BFT_FREE(lstelt);
  }

  /* 2 - Definition of Voltage on each termin of transformers */
  /* 2.1 Computation of Intensity on each termin of transformers */
  for (int i = 0; i < nbelec; i++) {
    sirb[transfo->ielect[i]][transfo->ielecb[i]] = 0.;
    if (cs_glob_elec_option->ieljou == 4)
      siib[transfo->ielect[i]][transfo->ielecb[i]] = 0.;
  }

  for (int i = 0; i < nbelec; i++) {
    if (transfo->ielect[i] != 0) {
      sirb[transfo->ielect[i]][transfo->ielecb[i]] += sir[i];
      if (cs_glob_elec_option->ieljou == 4)
        siib[transfo->ielect[i]][transfo->ielecb[i]] += sii[i];
    }
  }

  /* 2.2 RVoltage on each termin */
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

  /* 2.3 Total intensity for a transformer
   *     (zero valued WHEN Offset established) */
  for (int ntf = 0; ntf < nbtrf; ntf++) {
    sirt[ntf] = 0.;
    if (cs_glob_elec_option->ieljou == 4)
      siit[ntf] = 0.;
  }

  for (int i = 0; i < nbelec; i++) {
    if (transfo->ielect[i] != 0) {
      sirt[i] += sir[i];
      if (cs_glob_elec_option->ieljou == 4)
        siit[i] += sii[i];
    }
  }

  /* 2.4 Take in account of Offset */
  double capaeq = 3.;

  for (int ntf = 0; ntf < nbtrf; ntf++) {
    transfo->uroff[ntf] += sirt[ntf] / capaeq;
    if (cs_glob_elec_option->ieljou == 4)
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
      if (cs_glob_elec_option->ieljou == 4)
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

  /* 2.5 Take in account of Boundary Conditions */
  for (int i = 0; i < nbelec; i++) {
    sprintf(name, "%07d", transfo->ielecc[i]);
    cs_selector_get_b_face_num_list(name, &nelts, lstelt);

    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
      int ifac = lstelt[ilelt];
      int iel = ifabor[ifac];

      itypfb[ifac] = *iparoi;
      izfppp[ifac] = i;

      if (transfo->ielect[i] != 0) {
        f = CS_F_(potr);
        isca = cs_field_get_key_int(f, keyvar) - 1;
        icodcl[isca * nfabor + ifac] = 1;
        rcodcl[isca * nfabor + ifac] = ur[transfo->ielect[i]][transfo->ielecb[i]];

        if (cs_glob_elec_option->ieljou == 4) {
          f = CS_F_(poti);
          isca = cs_field_get_key_int(f, keyvar) - 1;
          icodcl[isca * nfabor + ifac] = 1;
          rcodcl[isca * nfabor + ifac] = ur[transfo->ielect[i]][transfo->ielecb[i]];
        }
      }
      else {
        f = CS_F_(potr);
        isca = cs_field_get_key_int(f, keyvar) - 1;
        icodcl[isca * nfabor + ifac] = 3;
        rcodcl[2 * nfabor * (*nvarcl) + isca * nfabor + ifac] = 0.;

        if (cs_glob_elec_option->ieljou == 4) {
          f = CS_F_(poti);
          isca = cs_field_get_key_int(f, keyvar) - 1;
          icodcl[isca * nfabor + ifac] = 3;
          rcodcl[2 * nfabor * (*nvarcl) + isca * nfabor + ifac] = 0.;
        }
      }
    }

    BFT_FREE(lstelt);
  }

  /* Test, if not any reference transformer
   *       a piece of wall may be at ground. */
  if (transfo->ntfref == 0) {
    bool itrouv = false;
    for (int ifac = 0; ifac < nfabor; ifac++) {
      if (itypfb[ifac] == *iparoi) {
        f = CS_F_(potr);
        isca = cs_field_get_key_int(f, keyvar) - 1;
        if (icodcl[isca * nfabor + ifac] == 1) {
          if (cs_glob_elec_option->ieljou == 3) {
            if (fabs(rcodcl[isca * nfabor + ifac]) < 1.e-20)
              itrouv = true;
          }
          else if (cs_glob_elec_option->ieljou == 4) {
            double val = fabs(rcodcl[isca * nfabor + ifac]);
            f = CS_F_(poti);
            isca = cs_field_get_key_int(f, keyvar) - 1;
            if (fabs(rcodcl[isca * nfabor + ifac]) < 1.e-20 &&
                val < 1.e-20)
              itrouv = true;
          }
        }
      }
    }
    if (!itrouv)
      bft_error(__FILE__, __LINE__, 0,
              _("ERROR in JOULE : \n"
                "Lack of reference : choose a transformer for wich\n"
                "offset is assumed zero or a face at ground on the\n"
                "boundary\n"));
  }

  BFT_FREE(sir);
  BFT_FREE(sii);
  BFT_FREE(sirb);
  BFT_FREE(siib);
  BFT_FREE(ur);
  BFT_FREE(ui);
  BFT_FREE(sirt);
  BFT_FREE(siit);
  BFT_FREE(nborne);
}


END_C_DECLS

