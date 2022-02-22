/*============================================================================
 * Compute rotation/curvature correction.
 *============================================================================*/

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

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gradient.h"
#include "cs_math.h"
#include "cs_physical_constants.h"
#include "cs_rotation.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_rotation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_turbulence_rotation.c
 *
 * Compute rotation/curvature correction for eddy-viscosity models.
 *
 * Two types of rotation/curvature correction can be computed, depending on
 * the specific eddy-viscosity model:
 *
 * - itycor = 1: - Cazalbou correction (variable Ce2 coefficient in the
 *                 destruction term of dissipation equation)
 *               - default correction for \f$ k - \epsilon \f$ type models,
 *                 including elliptic relaxation/blending models
 *                 (iturb = 20, 21, 50 or 51)
 *
 * - itycor = 2: - Spalart-Shur correction (production terms are multiplied
 *                 by a rotation function)
 *               - default correction for \f$ k - \omega \f$ SST or
 *                  Spalart-Allmaras (iturb = 60 or 70)
 */

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rotation/curvature correction for eddy-viscosity models.
 *
 * This function is called for the linear eddy viscosity RANS models,
 * when irccor = 1 is verified.
 *
 * \param[in]   dt      time step (per cell)
 * \param[out]  rotfct  rotation function of Spalart-Shur correction
 *                      at cell center
 * \param[out]  ce2rc   modified ce2 coeficient of Cazalbou correction
 *                      at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rotation_correction(const cs_real_t   dt[],
                                  cs_real_t         rotfct[],
                                  cs_real_t         ce2rc[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_real_t d1s2 =1./2.;

  /* Initialization
     ============== */

  /* Empty field at the first iteration */

  cs_real_6_t *cpro_straio
    = (cs_real_6_t *)cs_field_by_name("strain_rate_tensor")->val;

  /* Map field arrays */

  const cs_real_3_t *vela = (const cs_real_3_t *)CS_F_(vel)->val_pre;
  const cs_real_t *cvara_k = NULL;
  const cs_real_t *cvara_ep = NULL;
  const cs_real_t *cvara_omg = NULL;

  if (cs_glob_turb_rans_model->itycor == 1){
    cvara_k = (const cs_real_t *)CS_F_(k)->val_pre;
    cvara_ep = (const cs_real_t *)CS_F_(eps)->val_pre;
  }
  else if (cs_glob_turb_rans_model->itycor == 2){
    if (cs_glob_turb_model->iturb == 60){
      cvara_omg = (const cs_real_t *)CS_F_(omg)->val_pre;
    }
  }

  cs_real_t matrot[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  if (cs_glob_physical_constants->icorio == 1) {
    /* In case of a rotating frame, all cells belong to the same "rotor" */
    int r_num = 1;
    const cs_rotation_t *r = cs_glob_rotation + r_num;
    cs_rotation_add_coriolis_t(r, 1., matrot);
  }

  /* Preliminary calculations
     ======================== */

  /* Compute the strain rate and absolute vorticity tensor
     ----------------------------------------------------- */

  cs_real_6_t *strain = NULL;
  cs_real_3_t *vortab = NULL;
  cs_real_33_t *gradv = NULL;

  BFT_MALLOC(strain, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(vortab, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(CS_F_(vel),
                           true,   /* use_previous_t */
                           1 ,     /* inc */
                           gradv);

  /* Compute the strain rate tensor (symmetric)
   *          S_ij = 0.5(dU_i/dx_j+dU_j/dx_i)
   *
   * and the absolute vorticity tensor (anti-symmetric)
   *          W_ij = 0.5(dU_i/dx_j-dU_j/dx_i) + e_imj*Omega_m
   *
   * Only the non zero components in the upper triangle are stored */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    /* S11 */
    strain[c_id][0] = gradv[c_id][0][0];
    /* S22 */
    strain[c_id][1] = gradv[c_id][1][1];
    /* S33 */
    strain[c_id][2] = gradv[c_id][2][2];
    /* S12 */
    strain[c_id][3] = d1s2*(gradv[c_id][0][1] + gradv[c_id][1][0]);
    /* S13 */
    strain[c_id][4] = d1s2*(gradv[c_id][0][2] + gradv[c_id][2][0]);
    /* S23 */
    strain[c_id][5] = d1s2*(gradv[c_id][1][2] + gradv[c_id][2][1]);
    /* W12 */
    vortab[c_id][0] = d1s2*(gradv[c_id][0][1] - gradv[c_id][1][0]) + matrot[1][0];
    /* W13 */
    vortab[c_id][1] = d1s2*(gradv[c_id][0][2] - gradv[c_id][2][0]) + matrot[2][0];
    /* W23 */
    vortab[c_id][2] = d1s2*(gradv[c_id][1][2] - gradv[c_id][2][1]) + matrot[2][1];
  }

  /* Partially free memory (strain and vortab arrays are deallocated later) */

  BFT_FREE(gradv);

  /* Computation of:
   * --------------
   *
   *   brtild = 2.W_ik.S_jk(DS_ij/Dt + (e_imn.S_jn + e_jmn.S_in)*Omega_m)
   *     eta1 = S_ij.S_ij
   *     eta2 = W_ij.W_ij
   */

  cs_real_63_t *grdsij = NULL;
  cs_real_t *brtild = NULL;
  cs_real_t *eta1 = NULL;
  cs_real_t *eta2 = NULL;

  BFT_MALLOC(grdsij, n_cells_ext, cs_real_63_t);
  BFT_MALLOC(brtild, n_cells, cs_real_t);
  BFT_MALLOC(eta1, n_cells, cs_real_t);
  BFT_MALLOC(eta2, n_cells, cs_real_t);

  /* Index connectivity */

  /* istrai(i,j): position of the (i,j) component of the tensor */
  /*              in the strain and straio arrays */
  /* ivorab(i,j): position of the (i,j) component of the tensor */
  /*              in the vortab array */
  /* sigvor(i,j): sign of the (i,j) component of the absolute vorticity tensor */
  /*              = 1  if i > j */
  /*              = -1 if i < j */
  /*              = 0  if i = j */

  int istrai[3][3];
  int ivorab[3][3];
  cs_real_t sigvor[3][3];

  istrai[0][0] = 0;
  istrai[1][1] = 1;
  istrai[2][2] = 2;
  istrai[1][0] = 3;
  istrai[2][0] = 4;
  istrai[2][1] = 5;
  istrai[0][1] = istrai[1][0];
  istrai[0][2] = istrai[2][0];
  istrai[1][2] = istrai[2][1];

  ivorab[0][0] = 0;
  ivorab[1][1] = 0;
  ivorab[2][2] = 0;
  ivorab[1][0] = 0;
  ivorab[2][0] = 1;
  ivorab[2][1] = 2;
  ivorab[0][1] = ivorab[1][0];
  ivorab[0][2] = ivorab[2][0];
  ivorab[1][2] = ivorab[2][1];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i < j)
        sigvor[i][j] = -1.;
      else if (i == j)
        sigvor[i][j] = 0.;
      else
        sigvor[i][j] = 1.;
    }
  }

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    brtild[i] = 0.;
    eta1[i] = 0.;
    eta2[i] = 0.;
  }

  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 5
      || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {
    cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS){
    cs_field_get_key_struct(CS_F_(nusa), key_cal_opt_id, &var_cal_opt);
  }

  int nswrgp = var_cal_opt.nswrgr;
  int imligp = var_cal_opt.imligr;
  int iwarnp = var_cal_opt.verbosity;
  cs_real_t epsrgp = var_cal_opt.epsrgr;
  cs_real_t climgp = var_cal_opt.climgr;

  int inc = 1;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_tensor("Sij_gradient",
                     gradient_type,
                     halo_type,
                     inc,
                     nswrgp,
                     iwarnp,
                     imligp,
                     epsrgp,
                     climgp,
                     NULL,   /* Use default Neumann BC */
                     NULL,
                     strain,
                     grdsij);

  cs_real_t dsijdt, trrota, wiksjk;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {

        /* material derivative of S_ij */
        if (cs_glob_time_step_options->idtvar < 0)
          dsijdt = 0.;
        else
          dsijdt =    (strain[c_id][istrai[i][j]]
                    - cpro_straio[c_id][istrai[i][j]])/dt[c_id];

        dsijdt +=  cs_math_3_dot_product(vela[c_id],
                                         grdsij[c_id][istrai[i][j]]);

        /* (e_imn.S_jn+e_jmn.S_in)*Omega_m term */

        trrota = 0.;
        for (cs_lnum_t k = 0; k < 3; k++) {
          trrota +=   matrot[k][i] * strain[c_id][istrai[j][k]]
                    + matrot[k][j] * strain[c_id][istrai[i][k]];
        }

        /* W_ik.S_jk term */

        wiksjk = 0.;
        for (cs_lnum_t k = 0; k < 3; k++) {
          wiksjk +=   sigvor[i][k]*vortab[c_id][ivorab[i][k]]
                    * strain[c_id][istrai[j][k]];
        }

        /* brtild, eta1, eta2 (see the definitions above) */

        brtild[c_id] += 2. * wiksjk * (dsijdt + trrota);
        eta1[c_id] += cs_math_pow2(strain[c_id][istrai[i][j]]);
        eta2[c_id] += cs_math_pow2(sigvor[i][j] * vortab[c_id][ivorab[i][j]]);
      }
    }
  }

  /* Effective computation of the rotation correction
     ================================================ */

  cs_real_t stilde;
  cs_real_t wtilde;

  if (cs_glob_turb_rans_model->itycor == 1) {

    /* 2.1 Cazalbou correction
       ----------------------- */

    cs_real_t xk;
    cs_real_t xe;
    cs_real_t rotild;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Computation of stilde = sqrt(2.S_ij.S_ij)
         *          and wtilde = sqrt(W_ij.W_ij/2) */

      stilde = CS_MAX(sqrt(eta1[c_id] * 2.), 1.e-15);
      wtilde = CS_MAX(sqrt(eta2[c_id] * 0.5), 1.e-15);

      xk = CS_MAX(cvara_k[c_id], 1.e-15);
      xe = CS_MAX(cvara_ep[c_id], 1.e-15);
      rotild = xe/wtilde/xk;
      brtild[c_id] = -brtild[c_id]*xk/xe/pow(stilde,3.);

      /* Variable C_eps_2 coefficient of Cazalbou */

      ce2rc[c_id] =   cs_turb_ccaze2 + (cs_turb_ccaze2 - 1.)
                    / (1. + cs_turb_ccaza*pow(rotild,1.5))
                    + cs_turb_ccaze2*cs_turb_ccazsc*stilde*xk/xe
                    * (  tanh(cs_turb_ccazb*brtild[c_id] + cs_turb_ccazc)
                       - cs_turb_ccazd);

      ce2rc[c_id] = CS_MAX(ce2rc[c_id], 0.);
    }

  }
  else if (cs_glob_turb_rans_model->itycor == 2) {

    /* Spalart-Shur correction
     *------------------------
     * (including modifications of Smirnov & Menter, ASME, 2009) */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Computation of stilde = 2.S_ij.S_ij and wtilde = 2.W_ij.W_ij */

      stilde = CS_MAX(sqrt(eta1[c_id] * 2.), 1.e-15);
      wtilde = CS_MAX(sqrt(eta2[c_id] * 2.), 1.e-15);

      cs_real_t echtm2 = stilde;

      if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA)
        echtm2 = CS_MAX(echtm2, cs_turb_cmu*cs_math_pow2(cvara_omg[c_id]));

      brtild[c_id] = brtild[c_id] / sqrt(wtilde*cs_math_pow3(echtm2));

      cs_real_t rstar = sqrt(stilde)/sqrt(wtilde);

      /* Rotation function of Spalart & Shur */

      rotfct[c_id] =   (1. + cs_turb_cssr1)*2.*rstar/(1. + rstar)
                     * (1. - cs_turb_cssr3*atan(cs_turb_cssr2*brtild[c_id]))
                     - cs_turb_cssr1;

      rotfct[c_id] = CS_MIN(CS_MAX((double)rotfct[c_id],0.),1.25);

    }
  }

  /* Finalization
     ============ */

  /* Save the strain rate tensor for the next time step */

  if (cs_glob_time_step_options->idtvar >= 0) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 6; i++)
        cpro_straio[c_id][i] = strain[c_id][i];
    }
  }

  /* Free memory */

  BFT_FREE(strain);
  BFT_FREE(vortab);
  BFT_FREE(grdsij);
  BFT_FREE(brtild);
  BFT_FREE(eta1);
  BFT_FREE(eta2);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
