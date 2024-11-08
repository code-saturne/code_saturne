/*============================================================================
 * Porous model management
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
#include "cs_porosity_from_scan.h"
#include "cs_preprocess.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_porous_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_porous_model.cpp
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

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

  static cs_lnum_t n_cells_prev = 0;

  if (cs_glob_porous_model > 0) {
    cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
    if (mq->c_disable_flag == nullptr || m->n_cells != n_cells_prev) {
      /* Currently should handle turbomachinery (ghost cells only can change),
         not mesh refinement */
      assert(m->n_cells == n_cells_prev || n_cells_prev == 0);

      CS_REALLOC_HD(mq->c_disable_flag, n_cells_ext, int, cs_alloc_mode);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
        mq->c_disable_flag[cell_id] = 0;

      n_cells_prev = m->n_cells;
    }
    else {
      if (mq->has_disable_flag != 0)
        CS_REALLOC_HD(mq->c_disable_flag, n_cells_ext, int, cs_alloc_mode);
      if (m->halo != nullptr)
        cs_halo_sync_untyped(m->halo, CS_HALO_STANDARD, sizeof(int),
                             mq->c_disable_flag);
    }
  }
  else {
    if (mq->c_disable_flag == nullptr)
      CS_MALLOC_HD(mq->c_disable_flag, 1, int, cs_alloc_mode);
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
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

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
        mq->i_f_face_cog = mq->i_face_cog;
        mq->b_f_face_cog = mq->b_face_cog;
        mq->i_f_face_factor = nullptr;
        mq->b_f_face_factor = nullptr;
        mq->cell_f_cen = mq->cell_cen;
        mq->cell_s_cen = nullptr;
        mq->c_w_face_normal = nullptr;
        mq->c_w_face_surf = nullptr;
        mq->c_w_face_cog = nullptr;
        mq->c_w_dist_inv = nullptr;
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
        mq->i_f_face_cog    = cs_field_by_name("i_f_face_cog")->val;
        mq->b_f_face_cog    = cs_field_by_name("b_f_face_cog")->val;
        mq->cell_f_cen      = cs_field_by_name("cell_f_cen")->val;
        mq->cell_s_cen      = cs_field_by_name("cell_s_cen")->val;
        mq->c_w_face_normal = cs_field_by_name("c_w_face_normal")->val;
        mq->c_w_face_surf   = cs_field_by_name("c_w_face_surf")->val;
        mq->c_w_face_cog    = cs_field_by_name("c_w_face_cog")->val;
        mq->c_w_dist_inv    = cs_field_by_name("c_w_dist_inv")->val;
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
/*!
 * \brief  Compute solid quantities and update porosity field
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_mesh_quantities_update(void)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_real_3_t *cen_points = nullptr;
  cs_field_t *f = cs_field_by_name_try("cell_scan_points_cog");
  if (f != nullptr)
    cen_points = (cs_real_3_t *)f->val;

  cs_mesh_quantities_solid_compute(m, cen_points,  mq);

  /* synchronize for use in fluid face factor calculation */
  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cs_field_by_name("porosity")->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic computation of the face porosity and factors.
 *
 * This is useful for the integral porous model.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_auto_face_porosity(void)
{
  if (cs_glob_porous_model < 3)
    return;

  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Get the cell porosity field value */
  cs_real_t *cpro_porosi = cs_field_by_name("porosity")->val;

  if (m->halo != nullptr)
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_porosi);

  /* Set interior face values */

  {
    const cs_lnum_t n_i_faces = m->n_i_faces;
    const cs_lnum_2_t *i_face_cells
      = (const cs_lnum_2_t *)m->i_face_cells;

    const cs_real_3_t *restrict i_face_normal
      = (const cs_real_3_t *)mq->i_face_normal;
    cs_real_3_t *restrict i_f_face_normal = (cs_real_3_t *)mq->i_f_face_normal;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      cs_lnum_t c_id0 = i_face_cells[face_id][0];
      cs_lnum_t c_id1 = i_face_cells[face_id][1];

      cs_real_t face_porosity = CS_MIN(cpro_porosi[c_id0], cpro_porosi[c_id1]);

      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

      mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

      if (mq->i_f_face_factor != nullptr) {
        if (face_porosity > cs_math_epzero) {
          mq->i_f_face_factor[face_id][0] = cpro_porosi[c_id0] / face_porosity;
          mq->i_f_face_factor[face_id][1] = cpro_porosi[c_id1] / face_porosity;
        }
        else {
          mq->i_f_face_factor[face_id][0] = 1.;
          mq->i_f_face_factor[face_id][1] = 1.;
        }
      }

    }
  }

  /* Set boundary face values */

  {
    const cs_lnum_t n_b_faces = m->n_b_faces;
    const cs_lnum_t *b_face_cells
      = (const cs_lnum_t *)m->b_face_cells;

    const cs_real_3_t *restrict b_face_normal
      = (const cs_real_3_t *)mq->b_face_normal;
    cs_real_3_t *restrict b_f_face_normal
      = (cs_real_3_t *)mq->b_f_face_normal;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_lnum_t c_id = b_face_cells[face_id];

      cs_real_t face_porosity = cpro_porosi[c_id];

      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

      mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

      if (mq->b_f_face_factor != nullptr) {
        if (face_porosity > cs_math_epzero) {
          mq->b_f_face_factor[face_id] = cpro_porosi[c_id] / face_porosity;
        }
        else {
          mq->b_f_face_factor[face_id] = 1.;
        }
      }

    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 *! \brief Penalize porosity and fluid surfaces.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_clip(void)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_halo_t *halo = m->halo;
  const cs_real_t *cell_vol = mq->cell_vol;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_real_t *porosi = cs_field_by_name("porosity")->val;
  cs_real_t *cell_f_vol = mq->cell_f_vol;
  int *c_disable_flag = mq->c_disable_flag;

  cs_halo_sync_var(halo, CS_HALO_STANDARD, porosi);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    /* Penalization of solid cells */
    if (porosi[c_id] < cs_math_epzero) {
      porosi[c_id] = 0.0;
      c_disable_flag[c_id] = 1;
    }

    cell_f_vol[c_id] = cell_vol[c_id] * porosi[c_id];
  }

  /* For integral formulation, in case of 0 fluid volume, clip fluid faces */
  if (cs_glob_porous_model == 3) {
    cs_real_3_t *i_f_face_normal = (cs_real_3_t *)mq->i_f_face_normal;
    cs_real_3_t *b_f_face_normal = (cs_real_3_t *)mq->b_f_face_normal;
    cs_real_t *i_f_face_surf = (cs_real_t *)mq->i_f_face_surf;
    cs_real_t *b_f_face_surf = mq->b_f_face_surf;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      /* TODO compute i_f_face_factor with porosi
         AND fluid surface and surface:
         epsilon_i*surface/f_surface */
      if (c_disable_flag[i_face_cells[f_id][0]] == 1) {
        i_f_face_normal[f_id][0] = 0.0;
        i_f_face_normal[f_id][1] = 0.0;
        i_f_face_normal[f_id][2] = 0.0;
        i_f_face_surf[f_id] = 0.;
      }
      else if (c_disable_flag[i_face_cells[f_id][1]] == 1) {
        i_f_face_normal[f_id][0] = 0.0;
        i_f_face_normal[f_id][1] = 0.0;
        i_f_face_normal[f_id][2] = 0.0;
        i_f_face_surf[f_id] = 0.0;
      }
    }

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      /* TODO
         compute i_f_face_factor with porosi AND fluid surface and surface:
         epsilon_i*surface/f_surface */
      if (c_disable_flag[b_face_cells[f_id]] == 1) {
        b_f_face_normal[f_id][0] = 0.0;
        b_f_face_normal[f_id][1] = 0.0;
        b_f_face_normal[f_id][2] = 0.0;
        b_f_face_surf[f_id] = 0.0;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Preprocess the fluid surfaces.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_fluid_surfaces_preprocessing(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_real_3_t *restrict i_face_normal = (cs_real_3_t *)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal = (cs_real_3_t *)mq->b_face_normal;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells = (const cs_lnum_t *)m->b_face_cells;

  cs_real_t *restrict cell_f_vol = mq->cell_f_vol;
  cs_real_3_t *restrict i_f_face_normal = (cs_real_3_t *)mq->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal = (cs_real_3_t *)mq->b_f_face_normal;
  cs_real_t *restrict i_f_face_surf = (cs_real_t *)mq->i_f_face_surf;
  cs_real_t *restrict b_f_face_surf = (cs_real_t *)mq->b_f_face_surf;

  /* Pointer to porosity field */
  cs_field_t *f = cs_field_by_name("porosity");

  /* fluid volume for restart (cs_mesh_quantities_solid_compute is skipped)
     (else cell_f_vol is computed in cs_mesh_quantities_solid_compute) */

  if (cs_glob_porosity_from_scan_opt->use_staircase) {
    for (cs_lnum_t c_id = 0; c_id < m->n_cells_with_ghosts; c_id++) {
      cell_f_vol[c_id] = f->val[c_id] * mq->cell_vol[c_id];
    }
  }

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_lnum_t c_id0 = i_face_cells[face_id][0];
    cs_lnum_t c_id1 = i_face_cells[face_id][1];

    if (mq->c_disable_flag[c_id0] == 1 || mq->c_disable_flag[c_id1] == 1) {
      i_f_face_normal[face_id][0] = 0.;
      i_f_face_normal[face_id][1] = 0.;
      i_f_face_normal[face_id][2] = 0.;
      i_f_face_surf[face_id] = 0.;
    }
    else {
      cs_real_t face_porosity = CS_MIN(f->val[c_id0], f->val[c_id1]);

      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

      mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);
    }

    if (mq->i_f_face_factor != nullptr) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mq->i_f_face_factor[face_id][0] = cpro_porosi[c_id0] / face_porosity;
      //  mq->i_f_face_factor[face_id][1] = cpro_porosi[c_id1] / face_porosity;
      //}
      //else {
        mq->i_f_face_factor[face_id][0] = 1.;
        mq->i_f_face_factor[face_id][1] = 1.;
      //}
    }
  }

  /* Set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    cs_lnum_t c_id = b_face_cells[face_id];

    if (mq->c_disable_flag[c_id] == 1) {
      b_f_face_normal[face_id][0] = 0.;
      b_f_face_normal[face_id][1] = 0.;
      b_f_face_normal[face_id][2] = 0.;
      b_f_face_surf[face_id] = 0.;
    }
    else {
      cs_real_t face_porosity = f->val[c_id];

      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

      mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);
    }

    if (mq->b_f_face_factor != nullptr) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mq->b_f_face_factor[face_id] = cpro_porosi[c_id] / face_porosity;
      //}
      //else {
      mq->b_f_face_factor[face_id] = 1.;
      //}
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS