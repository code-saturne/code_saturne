/*============================================================================
 * Coal combustion model.
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
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

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

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_coal_radst(int         id,
                cs_real_t  *smbrs,
                cs_real_t  *rovsdt);

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
 * \brief Create coal model.
 *
 * \return  pointer to coal model structure.
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_create(void)
{
  cs_coal_model_t *cm;

  BFT_MALLOC(cm, 1, cs_coal_model_t);

  cm->n_coals = 0;
  cm->nclacp = 0;
  cm->nsolim = 0;

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

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coal model.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_model_destroy(cs_coal_model_t **cm)
{
  if (cm != NULL) {
    cs_coal_model_t *_cm = *cm;
    BFT_FREE(_cm);
    *cm = _cm;
  }
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
