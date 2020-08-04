/*============================================================================
 * Porous model management
 *============================================================================*/

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

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_preprocess.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_porous_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_porous_model.c
 *
 * \brief Porous model management
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Choice of the porous model */
int cs_glob_porous_model = 0;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_porous_model_get_pointers(int   **iporos);

int
cs_f_porous_model_cell_is_active(cs_lnum_t  cell_id);

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to global variables.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iporos  --> pointer to cs_glob_porous_model
 *----------------------------------------------------------------------------*/

void
cs_f_porous_model_get_pointers(int   **iporos)
{
  *iporos = &(cs_glob_porous_model);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return 0 if cell is disabled, 1 otherwise
 *
 * \param[in]  cell_id
 *
 * \return  0  if cell is disabled
 */
/*----------------------------------------------------------------------------*/

int
cs_f_porous_model_cell_is_active(cs_lnum_t  cell_id)
{
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  return (1 - (mq->has_disable_flag
          * mq->c_disable_flag[mq->has_disable_flag * cell_id]));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set porous model option.
 *
 * parameters:
 *   porous_model <-- porous model option (> 0 for porosity)
 *----------------------------------------------------------------------------*/

void
cs_porous_model_set_model(int  porous_model)
{
  cs_glob_porous_model = porous_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize disable_flag
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_disable_flag(void)
{
  cs_mesh_t *m =cs_glob_mesh;
  cs_mesh_quantities_t *mq =cs_glob_mesh_quantities;

  if (cs_glob_porous_model > 0) {
    cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
    if (mq->c_disable_flag == NULL) {
      BFT_MALLOC(mq->c_disable_flag, n_cells_ext, int);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
        mq->c_disable_flag[cell_id] = 0;
    }
    else {
      cs_lnum_t n_cells = m->n_cells;
      BFT_REALLOC(mq->c_disable_flag, n_cells_ext, int);
      for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++)
        mq->c_disable_flag[cell_id] = 0;
    }
    if (m->halo != NULL)
      cs_halo_sync_untyped(m->halo, CS_HALO_STANDARD, sizeof(int),
                           mq->c_disable_flag);
  }
  else {
    if (mq->c_disable_flag == NULL)
      BFT_MALLOC(mq->c_disable_flag, 1, int);
    mq->c_disable_flag[0] = 0;
  }

  /* Update Fortran pointer quantities */
  cs_preprocess_mesh_update_fortran();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set (unset) has_disable_flag
 *
 * \param[in]  flag   1: on, 0: off
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_set_has_disable_flag(int  flag)
{
  cs_mesh_quantities_t *mq =cs_glob_mesh_quantities;

  mq->has_disable_flag = flag;

  /* if off, fluid surfaces point toward cell surfaces */
  /* Porous models */
  if (cs_glob_porous_model > 0) {
    if (flag == 0) {
      /* Set pointers of face quantities to the standard ones */
      if (cs_glob_porous_model == 3) {
        mq->i_f_face_normal = mq->i_face_normal;
        mq->b_f_face_normal = mq->b_face_normal;
        mq->i_f_face_surf   = mq->i_face_surf;
        mq->b_f_face_surf   = mq->b_face_surf;
        mq->i_f_face_factor = NULL;
        mq->b_f_face_factor = NULL;
      }
      mq->cell_f_vol = mq->cell_vol;
    }
    else {
      /* Use fluid surfaces and volumes */
      if (cs_glob_porous_model == 3) {
        mq->i_f_face_normal = cs_field_by_name("i_f_face_normal")->val;
        mq->b_f_face_normal = cs_field_by_name("b_f_face_normal")->val;
        mq->i_f_face_surf   = cs_field_by_name("i_f_face_surf")->val;
        mq->b_f_face_surf   = cs_field_by_name("b_f_face_surf")->val;
        mq->i_f_face_factor
          = (cs_real_2_t *)cs_field_by_name("i_f_face_factor")->val;
        mq->b_f_face_factor = cs_field_by_name("b_f_face_factor")->val;
      }
      mq->cell_f_vol        = cs_field_by_name("cell_f_vol")->val;
    }
  }

  /* Update Fortran pointer quantities */
  cs_preprocess_mesh_update_fortran();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Init fluid quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_fluid_quantities(void)
{
  if (cs_glob_porous_model == 3) {
    cs_mesh_init_fluid_sections(cs_glob_mesh, cs_glob_mesh_quantities);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
