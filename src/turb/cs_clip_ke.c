/*============================================================================
 * Clipping of the turbulent kinetic energy and the turbulent
 * dissipation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_log.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_parall.h"

#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_parameters.h"
#include "cs_log_iteration.h"

#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_clip_ke.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_clip_ke.c

  Clipping of the turbulent kinetic energy and the turbulent
  dissipation.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of the turbulent kinetic energy and turbulent dissipation.
 *
 * \param[in]     n_cells  number of cells
 * \param[in]     iclip    indicator = 0 if viscl0 is used
 *                         otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_clip_ke(cs_lnum_t  n_cells,
           int        iclip)
{
  cs_turb_ref_values_t *turb_ref_values = cs_get_glob_turb_ref_values();
  cs_real_t almax = turb_ref_values->almax;

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t viscl0 = phys_pro->viscl0; /* reference pressure */
  cs_real_t ro0    = phys_pro->ro0; /* reference density */

  cs_real_t *crom    = (cs_real_t *)CS_F_(rho)->val;
  cs_real_t *cvar_k  = (cs_real_t *)CS_F_(k)->val;
  cs_real_t *cvar_ep = (cs_real_t *)CS_F_(eps)->val;
  cs_real_t *viscl   =  (cs_real_t *)CS_F_(mu)->val;

  cs_var_cal_opt_t vcopt;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &vcopt);

  int iwarnk = vcopt.verbosity;

  /* Une petite valeur pour eviter des valeurs exactement nulles. */

  const double epz2 = cs_math_pow2(cs_math_epzero);

  /* Postprocess clippings? */

  int key_clipping_id = cs_field_key_id("clipping_id");

  int clip_k_id = cs_field_get_key_int(CS_F_(k), key_clipping_id);
  cs_real_t *cpro_k_clipped = NULL;
  if (clip_k_id >= 0) {
    cpro_k_clipped = cs_field_by_id(clip_k_id)->val;
  }

  int clip_e_id = cs_field_get_key_int(CS_F_(eps), key_clipping_id);
  cs_real_t *cpro_e_clipped = NULL;
  if (clip_e_id >= 0) {
    cpro_e_clipped = cs_field_by_id(clip_e_id)->val;
  }

  /* Stockage Min et Max pour log
   * =============================*/

  const cs_real_t l_threshold = 1.e12;
  cs_real_t *cvar_var = NULL;
  cs_real_t var;
  cs_real_t vmax[2], vmin[2];

  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0)
      cvar_var = cvar_k;
    else if (ii == 1)
      cvar_var = cvar_ep;

    vmin[ii] =  l_threshold;
    vmax[ii] = -l_threshold;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      var = cvar_var[c_id];
      vmin[ii] = CS_MIN(vmin[ii], var);
      vmax[ii] = CS_MAX(vmax[ii], var);
    }
  }

  if (cpro_k_clipped != NULL) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_k_clipped[c_id] = 0.;
  }
  if (cpro_e_clipped != NULL) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_e_clipped[c_id] = 0.;
  }

  /* Detection des valeurs hors norme "physiques"
   * uniquement pour avertissement
   * ou dans le cas ICLKEP = 1
   * =============================================*/

  cs_gnum_t iclpke = 0;
  cs_lnum_t iclpmn[2] = {0, 0};
  cs_real_t xk, xe, xkmin, xepmin, xkm, xepm;

  if (iwarnk >= 2 || cs_glob_turb_rans_model->iclkep == 1) {

    if (iclip == 1) {

      xkm = 1296.*sqrt(cs_turb_cmu)/cs_math_pow2(almax);
      xepm = 46656.*cs_turb_cmu/cs_math_pow4(almax);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        xk = cvar_k[c_id];
        xe = cvar_ep[c_id];
        xkmin = xkm * cs_math_pow2(viscl[c_id] / crom[c_id]);
        xepmin = xepm * cs_math_pow3(viscl[c_id] / crom[c_id]);
        if (xk <= xkmin || xe <= xepmin) {
          if (cs_glob_turb_rans_model->iclkep == 1) {
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkmin - xk;
            cvar_k[c_id]  = xkmin;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepmin - xe;
            cvar_ep[c_id] = xepmin;
          }
          iclpke += 1;
        }
      }

    }
    else if (iclip == 0) {

      xkmin = 1296. * sqrt(cs_turb_cmu) / cs_math_pow2(almax)
                    * cs_math_pow2(viscl0/ro0);
      xepmin = 46656. * cs_turb_cmu/cs_math_pow4(almax)
                      * cs_math_pow3(viscl0/ro0);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++){
        xk = cvar_k[c_id];
        xe = cvar_ep[c_id];
        if (xk <= xkmin || xe <= xepmin) {
          if (cs_glob_turb_rans_model->iclkep == 1) {
            cvar_k[c_id]  = xkmin;
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkmin - xk;
            cvar_ep[c_id] = xepmin;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepmin - xe;
          }
          iclpke += 1;
        }
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Call of %s with option = %d"),
                __func__, iclip);

    /* save clip counts for log */

    if (cs_glob_turb_rans_model->iclkep == 1) {
      iclpmn[0] = iclpke;
      iclpmn[1] = iclpke;
    }

    /* logging */

    if (iwarnk >= 2) {
      cs_parall_sum(1, CS_GNUM_TYPE, &iclpke);

      cs_log_printf(CS_LOG_DEFAULT,
                    "\n "
                    "%llu k-epsilon values beyond the scales based on almax\n",
                    (unsigned long long)iclpke);
    }

  }

  /* "standard" clipping ICLKEP = 0
   * ==============================*/

  cs_lnum_t iclpk2 = 0, iclpe2 = 0;

  if (cs_glob_turb_rans_model->iclkep == 0) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++){
      xk = cvar_k[c_id];
      xe = cvar_ep[c_id];
      if (fabs(xk) <= epz2) {
        iclpk2 = iclpk2 + 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = epz2 - cvar_k[c_id];
        cvar_k[c_id] = CS_MAX(cvar_k[c_id],epz2);
      }
      else if(xk <= 0.) {
        iclpk2 = iclpk2 + 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = -xk;
        cvar_k[c_id] = -xk;
      }
      if (fabs(xe) <= epz2) {
        iclpe2 = iclpe2 + 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = epz2 - cvar_ep[c_id];
        cvar_ep[c_id] = CS_MAX(cvar_ep[c_id], epz2);
      }
      else if(xe <= 0.) {
        iclpe2 = iclpe2 + 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = - xe;
        cvar_ep[c_id] = - xe;
      }
    }

    /* save clip counts for log */

    iclpmn[0] = iclpk2;
    iclpmn[1] = iclpe2;
  }

  cs_lnum_t iclpmx[1] = {0};
  int id;

  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0)
      id = CS_F_(k)->id;
    else if (ii == 1)
      id = CS_F_(eps)->id;

    cs_log_iteration_clipping_field(id,
                                    iclpmn[ii],
                                    0,
                                    vmin + ii,
                                    vmax + ii,
                                    iclpmn + ii,
                                    iclpmx);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
