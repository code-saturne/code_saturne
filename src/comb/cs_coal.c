/*============================================================================
 * Coal combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal.c

  \brief Coal combustion model.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Coal combustion model parameters structure */

cs_coal_model_t  *cs_glob_coal_model = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_co_models_init(void);

void
cs_f_cp_models_init(void);

void
cs_f_cp_model_map_coal(void);

void
cs_f_ppcpfu_models_init(void);

void
cs_f_thch_models_init(void);

void
cs_f_coal_radst(int         id,
                cs_real_t  *smbrs,
                cs_real_t  *rovsdt);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coal model.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

static void
_coal_model_finalize(void)
{
  BFT_FREE(cs_glob_coal_model);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Take in account the radiative source terms in the particle equation
 * of a given class for pulverized coal flame.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id      <-- field id
 *   smbrs   <-- right side of the system
 *   rovsdt  <-- system diagonal
 *----------------------------------------------------------------------------*/

void
cs_f_coal_radst(int         id,
                cs_real_t  *smbrs,
                cs_real_t  *rovsdt)
{
  const cs_field_t *f = cs_field_by_id(id);

  cs_coal_rad_transfer_st(f, smbrs, rovsdt);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate coal combustion model.
 *
 * \return  pointer to coal combustion model structure.
 *
 * \param[in]  type  coal combustion model type
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_set_model(cs_coal_model_type_t  type)
{
  cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = type;

  if (type == CS_COMBUSTION_COAL_NONE) {
    BFT_FREE(cs_glob_coal_model);
    return NULL;
  }
  else if (cs_glob_coal_model != NULL) {
    cs_glob_coal_model->type = type;
    return cs_glob_coal_model;
  }

  /* Create and initialize model structure */

  cs_coal_model_t *cm = NULL;

  BFT_MALLOC(cm, 1, cs_coal_model_t);

  cs_glob_coal_model = cm;

  /* Members also present in gas combustion model */

  cm->n_gas_el_comp = 0;
  cm->n_gas_species = 0;
  cm->n_atomic_species = 0;
  cm->n_reactions = 0;
  cm->n_tab_points = 0;
  cm->pcigas = 0;
  cm->xco2 = 0;
  cm->xh2o = 0;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS; i++) {
    cm->wmole[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_OXYDANTS; i++) {
    cm->oxyo2[i] = 0;
    cm->oxyn2[i] = 0;
    cm->oxyh2o[i] = 0;
    cm->oxyco2[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_TABULATION_POINTS; i++) {
    cm->th[i] = 0;
  }

  /* Members specific to coal combustion model */

  cm->type = type;

  cm->n_coals = 0;
  cm->nclacp = 0;
  cm->nsolim = 0;
  cm->idrift = 0;

  cm->ieqnox = 1;
  cm->ieqco2 = 1;

  cm->ico = -1;
  cm->io2 = -1;
  cm->in2 = -1;
  cm->ico2 = -1;
  cm->ih2o = -1;

  cm->ckabs0 = 0;

  for (int i = 0; i < CS_COMBUSTION_MAX_COALS; i++) {
    cm->n_classes_per_coal[i] = 0;
    cm->ich[i] = 0;
    cm->ick[i] = 0;
    cm->iash[i] = 0;
    cm->iwat[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_TABULATION_POINTS; i++) {
    for (int j = 0; j < CS_COMBUSTION_COAL_MAX_SOLIDS; j++) {
      cm->ehsoli[i][j] = 0;
    }
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_SOLIDS; i++) {
    cm->wmols[i] = 0;
    cm->eh0sol[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++)
    cm->ichcor[i] = 0;

  for (int i = 0; i < CS_COMBUSTION_MAX_COALS; i++) {
    cm->cp2ch[i] = 0;
    cm->xashch[i] = 0;
    cm->xwatch[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++) {
    cm->diam20[i] = 0;
    cm->dia2mn[i] = 0;
    cm->rho20[i] = 0;
    cm->rho2mn[i] = 0;
    cm->xmp0[i] = 0;
    cm->xmasch[i] = 0;
  }

  /* Set finalization callback */

  cs_base_at_finalize(_coal_model_finalize);

  /* Set mappings with Fortran */

  cs_f_co_models_init();
  cs_f_cp_models_init();
  cs_f_cp_model_map_coal();
  cs_f_ppcpfu_models_init();
  cs_f_thch_models_init();

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take in account the radiative source terms in the particle equation
 *        of a given class for pulverized coal flame.
 *
 * \param[in]      f       pointer to scalar field
 * \param[in, out] smbrs   right and side (explicit ST contribution)
 * \param[in, out] rovsdt  system diagonal (implicit ST contribution)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_rad_transfer_st(const cs_field_t  *f,
                        cs_real_t         *smbrs,
                        cs_real_t         *rovsdt)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  /* Initialization
   * -------------- */

  const int keyccl = cs_field_key_id("scalar_class");
  const int numcla = cs_field_get_key_int(f, keyccl);
  const int ipcl   = 1 + numcla;

  /* Radiative source terms
   * ---------------------- */

  char f_name[64];

  snprintf(f_name, 63, "rad_st_implicit_%02d", ipcl);
  cs_real_t *cpro_tsri = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "rad_st_%02d", ipcl);
  cs_real_t *cpro_tsre = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "x_p_%02d", numcla);
  const cs_real_t *cval_x_p = cs_field_by_name(f_name)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cpro_tsri[c_id] = cs_math_fmax(-cpro_tsri[c_id], 0.);

    if (cval_x_p[c_id] > cs_math_epzero) {
      /* Explicit part */
      smbrs[c_id] += cpro_tsre[c_id]*cell_f_vol[c_id]*cval_x_p[c_id];

      /* Implicit part */
      rovsdt[c_id] += cpro_tsri[c_id]*cell_f_vol[c_id];
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
