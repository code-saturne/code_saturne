/*============================================================================
 * Porous model management
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "alge/cs_gradient.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "cdo/cs_domain.h"
#include "alge/cs_matrix_default.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parameters.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_porosity_from_scan.h"
#include "base/cs_post.h"
#include "base/cs_preprocess.h"
#include "base/cs_renumber.h"

#include "fvm/fvm_nodal_from_desc.h"
#include "fvm/fvm_nodal_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_porous_model.h"

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

static cs_porous_model_extra_faces _porous_model_extra_faces = {
  .ib_mesh = nullptr,
  .mesh_id = 0,
  .activate_post = false
};

/*============================================================================
 * Global variables
 *============================================================================*/

cs_porous_model_extra_faces *cs_glob_porous_model_extra_faces
= &_porous_model_extra_faces;

/*! Choice of the porous model */
int cs_glob_porous_model = 0;

/*! Specific mesh quantities associated with porous model */
cs_mesh_quantities_t  *cs_glob_mesh_quantities_f = nullptr;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free fluid mesh quantities
 */
/*----------------------------------------------------------------------------*/

static void
_porous_mesh_quantities_f_free(void)
{
  /* It is imperative to not use cs_mesh_quantities_destroy(mq_f),
   * since some members are fields, whose life cycle is managed
   * separately.
   *
   * Here, free the members which are not fields. */

  if (cs_glob_mesh_quantities_f != nullptr) {
    cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;
    cs_mesh_quantities_t *mq_f = cs_glob_mesh_quantities_f;
    CS_FREE(mq_f->c_disable_flag);
    CS_FREE(mq_f->b_sym_flag);
    CS_FREE(mq_f->djjpf);
    CS_FREE(mq_f->diipf);
    CS_FREE(mq_f->dofij);
    CS_FREE(mq_f->diipb);
    CS_FREE(mq_f->dijpf);
    CS_FREE(mq_f->weight);
    CS_FREE(mq_f->i_dist);
    CS_FREE(mq_f->b_dist);
    if (mq_f->i_face_u_normal != mq_g->i_face_u_normal)
      CS_FREE(mq_f->i_face_u_normal);
    if (mq_f->b_face_u_normal != mq_g->b_face_u_normal)
      CS_FREE(mq_f->b_face_u_normal);
    CS_FREE(mq_f->corr_grad_lin);
    CS_FREE(mq_f->corr_grad_lin_det);
    CS_FREE(mq_f);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Realloc boundary arrays for all defined fields.
 *
 * Location sized must thus be known.
 *
 * Fields that do not own their data should all have been mapped at this
 * stage, and are checked.
 *
 * \param[in]   n_ib_cells   immersed boundary cell number
 */
/*----------------------------------------------------------------------------*/

static void
_field_ibm_reallocate(cs_lnum_t  n_ib_cells)
{
  const int n_fields = cs_field_n_fields();

  cs_dispatch_context ctx;

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);
    assert(f != nullptr);

    if (f->is_owner) {

      if (f->location_id == CS_MESH_LOCATION_CELLS) {

        if (f->bc_coeffs != nullptr) {
          const int location_b_id = CS_MESH_LOCATION_BOUNDARY_FACES;
          const cs_lnum_t *n_b_elts = cs_mesh_location_get_n_elts(location_b_id);

          const cs_lnum_t n_i_end = f->dim*n_b_elts[0];
          const cs_lnum_t n_i_new = f->dim*n_ib_cells;
          const cs_lnum_t n_i_start = n_i_end - n_i_new;

          /* Realloc val_f for fields already computed (porosity) */

          bool has_limiter = false;
          if (f->bc_coeffs->val_f_lim != f->bc_coeffs->val_f) {
            CS_REALLOC_HD(f->bc_coeffs->val_f_lim, n_i_end,
                          cs_real_t, cs_alloc_mode);
            CS_REALLOC_HD(f->bc_coeffs->flux_lim, n_i_end,
                          cs_real_t, cs_alloc_mode);

            cs_real_t *val_f_lim = f->bc_coeffs->val_f_lim;
            cs_real_t *flux_lim = f->bc_coeffs->flux_lim;
            ctx.parallel_for(n_i_new, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
              val_f_lim[n_i_start + e_id] = 0.;
              flux_lim[n_i_start + e_id] = 0.;
            });
            has_limiter = true;
          }

          if (f->bc_coeffs->val_f != nullptr) {
            CS_REALLOC_HD(f->bc_coeffs->val_f, n_i_end, cs_real_t, cs_alloc_mode);
            CS_REALLOC_HD(f->bc_coeffs->flux, n_i_end, cs_real_t, cs_alloc_mode);

            cs_real_t *val_f = f->bc_coeffs->val_f;
            cs_real_t *flux = f->bc_coeffs->flux;
            ctx.parallel_for(n_i_new, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
              val_f[n_i_start + e_id] = 0.;
              flux[n_i_start + e_id] = 0.;
            });

            ctx.wait();

            if (!(has_limiter)) {
              f->bc_coeffs->val_f_lim = f->bc_coeffs->val_f;
              f->bc_coeffs->flux_lim = f->bc_coeffs->flux;
            }

          }
        }
      } /* End cell location */

      else if (f->location_id == CS_MESH_LOCATION_BOUNDARY_FACES) {

        const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);

        f->update_size();

        /* Initialization */
        const cs_lnum_t n_i_end = f->dim*n_elts[2];
        const cs_lnum_t n_i_new = f->dim*n_ib_cells;
        const cs_lnum_t n_i_start = n_i_end - n_i_new;

        for (int time_id = 0; time_id < f->n_time_vals; time_id++) {
          cs_real_t *vals = f->vals[time_id];

          ctx.parallel_for(n_i_new, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
            vals[n_i_start + e_id] = 0.;
          });

          ctx.wait();
        }
      }
    }
    else {
      if (f->val == nullptr)
        bft_error(__FILE__, __LINE__, 0,
                  _("Field \"%s\"\n"
                    " requires mapped values which have not been set."),
                  f->name);
    }
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fluid mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_map_mesh_quantites_f_and_compute(void)
{
  /* Make fluid surfaces of mesh quantity point to the created fields */
  cs_glob_mesh_quantities_f = cs_mesh_quantities_create();
  cs_mesh_quantities_t *mq_f = cs_glob_mesh_quantities_f;
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;
  mq_f->has_disable_flag = 1;

  /* Note that cs_glob_mesh_quantities_f is not owner of its members */
  if (cs_glob_porous_model == 3) {
    mq_f->i_face_normal = cs_field_by_name("i_f_face_normal")->val;
    mq_f->b_face_normal = cs_field_by_name("b_f_face_normal")->val;
    mq_f->i_face_surf   = cs_field_by_name("i_f_face_surf")->val;
    mq_f->b_face_surf   = cs_field_by_name("b_f_face_surf")->val;
    mq_f->i_f_face_factor
      = (cs_real_2_t *)cs_field_by_name("i_f_face_factor")->val;
    mq_f->b_f_face_factor = cs_field_by_name("b_f_face_factor")->val;
    mq_f->i_face_cog    = (cs_real_3_t *)cs_field_by_name("i_f_face_cog")->val;
    mq_f->b_face_cog    = (cs_real_3_t *)cs_field_by_name("b_f_face_cog")->val;
    mq_f->cell_cen      = (cs_real_3_t *)cs_field_by_name("cell_f_cen")->val;
    mq_f->cell_s_cen      = (cs_real_3_t *)cs_field_by_name("cell_s_cen")->val;
    mq_f->c_w_face_normal = cs_field_by_name("c_w_face_normal")->val;
    mq_f->c_w_face_surf   = cs_field_by_name("c_w_face_surf")->val;
    mq_f->c_w_face_cog    = cs_field_by_name("c_w_face_cog")->val;
    mq_f->c_w_dist_inv    = cs_field_by_name("c_w_dist_inv")->val;
  }
  else {
    mq_f->i_face_normal   = mq_g->i_face_normal;
    mq_f->b_face_normal   = mq_g->b_face_normal;
    mq_f->i_face_surf     = mq_g->i_face_surf;
    mq_f->b_face_surf     = mq_g->b_face_surf;
    mq_f->i_f_face_factor = mq_g->i_f_face_factor;
    mq_f->b_f_face_factor = mq_g->b_f_face_factor;
    mq_f->i_face_cog      = mq_g->i_face_cog;
    mq_f->b_face_cog      = mq_g->b_face_cog;
    mq_f->cell_cen        = mq_g->cell_cen;
    mq_f->cell_s_cen      = mq_g->cell_s_cen;
    mq_f->c_w_face_normal = mq_g->c_w_face_normal;
    mq_f->c_w_face_surf   = mq_g->c_w_face_surf;
    mq_f->c_w_face_cog    = mq_g->c_w_face_cog;
    mq_f->c_w_dist_inv    = mq_g->c_w_dist_inv;
  }
  mq_f->cell_vol        = cs_field_by_name("cell_f_vol")->val;

  cs_glob_mesh_quantities = mq_f;
  cs_glob_domain->mesh_quantities = mq_f;
  cs_mesh_quantities_compute(cs_glob_mesh, mq_f);
  cs_porous_model_init_fluid_quantities();

  /* Prepare the deallocation */
  cs_base_at_finalize(_porous_mesh_quantities_f_free);
}

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
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

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
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;
  cs_mesh_quantities_t *mq_f = cs_glob_mesh_quantities_f;

  /* if off, fluid surfaces point toward cell surfaces */
  /* Porous models */
  if (cs_glob_porous_model > 0) {
    if (flag == 0) {
      /* Use geometric quantities */
      cs_glob_mesh_quantities = mq_g;
    }
    else {
      /* Use fluid quantities */
      cs_glob_mesh_quantities = mq_f;
    }

    /* Update Fortran pointer quantities */
    cs_preprocess_mesh_update_fortran();
  }

  cs_glob_mesh_quantities->has_disable_flag = flag;
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
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  /* Get the cell porosity field value */
  cs_real_t *cpro_porosi = cs_field_by_name("porosity")->val;

  if (m->halo != nullptr)
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_porosi);

  /* Set interior face values */

  {
    const cs_lnum_t n_i_faces = m->n_i_faces;
    const cs_lnum_2_t *i_face_cells = m->i_face_cells;

    const cs_real_3_t *restrict i_face_normal
      = (const cs_real_3_t *)mq_g->i_face_normal;
    cs_real_3_t *restrict i_f_face_normal = (cs_real_3_t *)mq->i_face_normal;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      cs_lnum_t c_id0 = i_face_cells[face_id][0];
      cs_lnum_t c_id1 = i_face_cells[face_id][1];

      cs_real_t face_porosity = cs::min(cpro_porosi[c_id0], cpro_porosi[c_id1]);

      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

      mq->i_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

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
      = (const cs_real_3_t *)mq_g->b_face_normal;
    cs_real_3_t *restrict b_f_face_normal
      = (cs_real_3_t *)mq->b_face_normal;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_lnum_t c_id = b_face_cells[face_id];

      cs_real_t face_porosity = cpro_porosi[c_id];

      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

      mq->b_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

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
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_halo_t *halo = m->halo;
  const cs_real_t *cell_vol = mq_g->cell_vol;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_real_t *porosi = cs_field_by_name("porosity")->val;
  cs_real_t *cell_f_vol = mq->cell_vol;
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
    cs_real_3_t *i_f_face_normal = (cs_real_3_t *)mq->i_face_normal;
    cs_real_3_t *b_f_face_normal = (cs_real_3_t *)mq->b_face_normal;
    cs_real_t *i_f_face_surf = (cs_real_t *)mq->i_face_surf;
    cs_real_t *b_f_face_surf = mq->b_face_surf;

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
  const cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_3_t *restrict i_face_normal = (cs_real_3_t *)mq_g->i_face_normal;
  const cs_real_3_t *restrict b_face_normal = (cs_real_3_t *)mq_g->b_face_normal;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_real_t *restrict cell_f_vol = mq->cell_vol;
  cs_real_3_t *restrict i_f_face_normal = (cs_real_3_t *)mq->i_face_normal;
  cs_real_3_t *restrict b_f_face_normal = (cs_real_3_t *)mq->b_face_normal;
  cs_real_t *restrict i_f_face_surf = mq->i_face_surf;
  cs_real_t *restrict b_f_face_surf = mq->b_face_surf;

  /* Pointer to porosity field */
  cs_field_t *f = cs_field_by_name("porosity");

  /* fluid volume for restart (cs_mesh_quantities_solid_compute is skipped)
     (else cell_f_vol is computed in cs_mesh_quantities_solid_compute) */

  if (cs_glob_porosity_from_scan_opt->use_staircase) {
    for (cs_lnum_t c_id = 0; c_id < m->n_cells_with_ghosts; c_id++) {
      cell_f_vol[c_id] = f->val[c_id] * mq_g->cell_vol[c_id];
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
      cs_real_t face_porosity = cs::min(f->val[c_id0], f->val[c_id1]);

      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

      mq->i_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);
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

      mq->b_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);
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
/*!
 * \brief Convert cell array to boundary array
 *
 * \param[in]   n_ib_cells     immersed cell number
 * \param[in]   ibcell_cells   immersed cell to cell connectivity
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_convert_cell_to_boundary(const cs_lnum_t   n_ib_cells,
                                         const cs_lnum_t   ibcell_cells[])
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq_f = cs_glob_mesh_quantities;

  const cs_real_3_t *restrict c_w_face_normal
    = (const cs_real_3_t *)mq_f->c_w_face_normal;
  const cs_real_t *restrict c_w_face_surf
    = (const cs_real_t *)mq_f->c_w_face_surf;
  const cs_real_3_t *restrict c_w_face_cog
    = (const cs_real_3_t *)mq_f->c_w_face_cog;
  const cs_real_t *c_w_dist_inv = (const cs_real_t *)mq_f->c_w_dist_inv;
  const cs_real_3_t *cell_cen = mq_f->cell_cen;

  const cs_lnum_t n_b_faces_old = m->n_b_faces;
  const cs_lnum_t n_g_b_faces_all = m->n_g_b_faces_all;
  const cs_lnum_t n_b_faces_tot = n_b_faces_old + n_ib_cells;

  m->n_b_faces = n_b_faces_tot;
  cs_mesh_location_build(m,
                         CS_MESH_LOCATION_BOUNDARY_FACES);

  /* Update bc_type */
  cs_boundary_conditions_ibm_create(n_ib_cells);

  /* Realloc mesh quantities array wich are not fields */
  CS_REALLOC_HD(m->b_face_family, n_b_faces_tot, cs_lnum_t, cs_alloc_mode);
  CS_REALLOC_HD(m->b_face_cells, n_b_faces_tot, cs_lnum_t, cs_alloc_mode);
  /*CS_REALLOC_HD(m->b_face_vtx_idx, n_b_faces_tot+1, cs_lnum_t, cs_alloc_mode);
  CS_REALLOC_HD(m->b_face_vtx_lst, m->b_face_vtx_idx[n_b_faces_old]+n_glob_vtx,
                cs_lnum_t, cs_alloc_mode);*/

  CS_REALLOC_HD(mq_f->b_face_u_normal, n_b_faces_tot, cs_nreal_3_t, cs_alloc_mode);
  CS_REALLOC_HD(mq_f->b_dist, n_b_faces_tot, cs_real_t, cs_alloc_mode);
  CS_REALLOC_HD(mq_f->b_sym_flag, n_b_faces_tot, int, cs_alloc_mode);
  CS_REALLOC_HD(mq_f->diipb, n_b_faces_tot, cs_rreal_3_t, cs_alloc_mode);

  cs_field_map_and_init_bcs();
  _field_ibm_reallocate(n_ib_cells);

  /* mq_f points to reallocated fields */
  mq_f->b_face_normal = cs_field_by_name("b_f_face_normal")->val;
  mq_f->b_face_surf = cs_field_by_name("b_f_face_surf")->val;
  mq_f->b_f_face_factor = cs_field_by_name("b_f_face_factor")->val;
  mq_f->b_face_cog = (cs_real_3_t *)cs_field_by_name("b_f_face_cog")->val;
  mq_f->cell_cen = (cs_real_3_t *)cs_field_by_name("cell_f_cen")->val;

  cs_lnum_t *b_face_cells = m->b_face_cells;
  //cs_lnum_t *b_face_vtx_idx = m->b_face_vtx_idx;
  //cs_lnum_t *b_face_vtx_lst = m->b_face_vtx_lst;
  cs_real_t *b_dist = mq_f->b_dist;
  cs_nreal_3_t *b_face_u_normal = mq_f->b_face_u_normal;
  cs_real_3_t *b_face_normal = (cs_real_3_t *)mq_f->b_face_normal;
  cs_real_3_t *b_face_cog = (cs_real_3_t *)mq_f->b_face_cog;
  cs_rreal_3_t *diipb = mq_f->diipb;
  cs_real_t *b_face_surf = mq_f->b_face_surf;

  /* Initialization in the ibm zone for cs_mesh_update_selectors */
  for (cs_lnum_t face_id = n_b_faces_old; face_id < n_b_faces_tot; face_id++) {
    m->b_face_family[face_id] = 1;
  }

  /* Rebuild the boundary zone */
  cs_mesh_free_b_rebuildable(m);

  cs_mesh_init_b_selector();

  CS_REALLOC(m->global_b_face_num, n_b_faces_tot, cs_gnum_t);
  cs_gnum_t *g_b_face_num = m->global_b_face_num;

  fvm_io_num_t *face_io_num = nullptr;
  const cs_gnum_t *face_gnum = nullptr;
  if (cs_glob_n_ranks > 1) {
    /* Order faces by increasing global number */
    face_io_num = fvm_io_num_create_from_scan(n_ib_cells);
    face_gnum = fvm_io_num_get_global_num(face_io_num);
  }

  for (cs_lnum_t i = 0; i < n_ib_cells; i++) {
    const cs_lnum_t c_id = ibcell_cells[i];
    const cs_lnum_t f_id = n_b_faces_old + i;
    const cs_real_t surf = c_w_face_surf[c_id];

    b_face_cells[f_id] = c_id;

    if (cs_glob_n_ranks > 1) {
      g_b_face_num[f_id] = n_g_b_faces_all + face_gnum[i];
    }
    else {
      g_b_face_num[f_id] = f_id+1;
    }

    for (cs_lnum_t j = 0; j < 3; j++) {
      b_face_cog[f_id][j] = c_w_face_cog[c_id][j];
      b_face_normal[f_id][j] = c_w_face_normal[c_id][j];
    }

    /* Unit normal */
    cs_math_3_normalize(c_w_face_normal[c_id], b_face_u_normal[f_id]);

    b_face_surf[f_id] = surf;

    b_dist[f_id] = (c_w_dist_inv[c_id] < DBL_MIN) ?
                    0.:
                    1. / c_w_dist_inv[c_id];

    assert(b_dist[f_id] > 0.);

    // Vector II' for immersed boundaries

    // ---> IF
    cs_real_t ib_vec_if[3] = {c_w_face_cog[c_id][0] - cell_cen[c_id][0],
                              c_w_face_cog[c_id][1] - cell_cen[c_id][1],
                              c_w_face_cog[c_id][2] - cell_cen[c_id][2]};

    // ---> ib_diipb = IF - (IF.NIJ)NIJ
    cs_math_3_orthogonal_projection(b_face_u_normal[f_id], ib_vec_if, diipb[f_id]);

    /*const int n_vtx = 0.5*(w_vtx_idx[c_id+1]-w_vtx_idx[c_id]);
    b_face_vtx_idx[f_id+1] = b_face_vtx_idx[f_id]+n_vtx;

    // vtx order : {n_vertices+1, n_vertices+2,..., n_vertices + n_glob_vtx}
    for (cs_lnum_t j = b_face_vtx_idx[f_id]; j < b_face_vtx_idx[f_id+1]; j++) {
      b_face_vtx_lst[j] = n_vertices+1+(j-b_face_vtx_idx[f_id]);
    }*/

    //n_vertices += n_vtx;
  }

  cs_renumber_b_faces(cs_glob_mesh);

  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);
  assert(m->cell_numbering != nullptr);
  assert(m->i_face_numbering != nullptr);
  assert(m->vtx_numbering != nullptr);
  assert(m->b_face_numbering != nullptr);

  /* Update global n_g_b_faces, and b_cells (from b_face_cells) */
  cs_mesh_update_auxiliary(m);

  cs_gradient_free_quantities();

  cs_mesh_adjacencies_update_cell_b_faces();

  cs_matrix_update_mesh();

  if (cs_glob_n_ranks > 1) {
    fvm_io_num_destroy(face_io_num);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize porous model arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_enable_post(void) {
  if (cs_glob_porous_model_extra_faces->activate_post == false)
    return;

  cs_glob_porous_model_extra_faces->ib_mesh = fvm_nodal_create("ib_mesh", 3);

  const int writer_ids[] = {CS_POST_WRITER_DEFAULT};

  int mesh_id = cs_post_get_free_mesh_id();

  cs_post_define_future_mesh(mesh_id,
                             CS_POST_MESH_BOUNDARY,
                             true,  /* auto variables */
                             1,     /* number of associated writers */
                             writer_ids);

  cs_glob_porous_model_extra_faces->mesh_id = mesh_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-processes the immersed boundary (ib) planes for display
 *         on paraview.
 *
 * \param[in]       n_ib_cells      ib cell number
 * \param[in]       n_glob_vtx      total vertex number
 * \param[in]       ibcell_cells    connectivity ib_cell->cells
 * \param[in]       vtx_ids         vertex ids on both sides of a IB vertex
 *                                  (v0<v1)
 * \param[in]       w_vtx_idx       ib vertex indexes
 * \param[in]       face_vertex_idx vertex indexes of the ib faces
 * \param[in]       w_vtx           ib vertex coordinates
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_post_immmersed_plane(const cs_lnum_t      n_ib_cells,
                                     const cs_lnum_t      n_glob_vtx,
                                     const cs_lnum_t      ibcell_cells[],
                                     const cs_lnum_t      vtx_ids[][2],
                                     const cs_lnum_t      w_vtx_idx[],
                                     const cs_lnum_t      face_vertex_idx[],
                                     const cs_real_t      w_vtx[][3])
{
  if (cs_glob_porous_model_extra_faces->activate_post == false)
    return;

  cs_lnum_t *face_vertex_ids;
  CS_MALLOC(face_vertex_ids, n_glob_vtx, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_glob_vtx; i++)
    face_vertex_ids[i] = i;

  /* Reordering the vertices to obtain a closed contour
     -------------------------------------------------- */

  cs_lnum_t n_w_vtx = 0;
  cs_coord_3_t *edge;
  CS_MALLOC(edge, n_glob_vtx, cs_coord_3_t);

  for (cs_lnum_t c_id_ib = 0; c_id_ib < n_ib_cells; c_id_ib++) {

    const cs_lnum_t c_id = ibcell_cells[c_id_ib];
    const int n_vtx = w_vtx_idx[c_id+1]-w_vtx_idx[c_id];

    /* Edge reference is the first 2 ib vtx by default, then connect other
       ib edges by finding the corresponding vtx ids associated to the
       second ib vertex */

    cs_lnum_t k = 0;
    const cs_lnum_t s_id = w_vtx_idx[c_id];
    const cs_lnum_t e_id = w_vtx_idx[c_id+1];

    const int n_edge = 0.5*n_vtx;
    cs_real_23_t edge_l[n_edge];
    cs_lnum_t vtx_ids_l[n_edge][2][2];

    for (cs_lnum_t j = s_id; j < e_id; j=j+2) {

      /* Store edges and associated vtx_ids of the ib cell */
      for (cs_lnum_t i = 0; i < 3; i++) {
        edge_l[k][0][i] = w_vtx[j][i];
        edge_l[k][1][i] = w_vtx[j+1][i];
      }
      vtx_ids_l[k][0][0] = vtx_ids[j][0];
      vtx_ids_l[k][0][1] = vtx_ids[j][1];

      vtx_ids_l[k][1][0] = vtx_ids[j+1][0];
      vtx_ids_l[k][1][1] = vtx_ids[j+1][1];
      k++;
    }

    /* The reference edge is the first k=0 */
    for (cs_lnum_t i = 0; i < 3; i++) {
      edge[n_w_vtx][i]   = edge_l[0][0][i];
      edge[n_w_vtx+1][i] = edge_l[0][1][i];
    }

    n_w_vtx += 2;

    /* Reference edge id = 0 */
    int edge_id = 0;
    int n_iter = 0;
    int cpt = n_edge - 1;

    /* Fist edge is x0->x1. The vtx to connect is x1. */
    cs_lnum_t vtx_to_connect[2] = {vtx_ids_l[0][1][0], vtx_ids_l[0][1][1]};

    do {

      for (int j = 0; j < n_edge; j++) {

        if (j != edge_id) {

          int new_vtx_id = -1;

          /* x0 of edge j is connected to the vertex */
          if (   vtx_ids_l[j][0][0] == vtx_to_connect[0]
              && vtx_ids_l[j][0][1] == vtx_to_connect[1]) {

            new_vtx_id = 1;
          }
          /* x1 of edge j is connected to the vertex */
          else if (   vtx_ids_l[j][1][0] == vtx_to_connect[0]
                   && vtx_ids_l[j][1][1] == vtx_to_connect[1]) {

            new_vtx_id = 0;
          }
          else {
            continue;
          }

          vtx_to_connect[0] = vtx_ids_l[j][new_vtx_id][0];
          vtx_to_connect[1] = vtx_ids_l[j][new_vtx_id][1];

          for (cs_lnum_t i = 0; i < 3; i++)
            edge[n_w_vtx][i] = edge_l[j][new_vtx_id][i];

          n_w_vtx++;
          edge_id = j;
          cpt -= 1;
          break;
        }
      }

      n_iter++;
      if (n_iter > 100)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem connecting vertices in cell_id = %d"), c_id);

    } while (cpt != 1);

  } /* End loop on ib cells */

  /* Create IBM mesh and give IBM faces and vertices */

  fvm_nodal_t *ib_mesh = cs_glob_porous_model_extra_faces->ib_mesh;
  const cs_lnum_t face_list_shift[2] = {0, n_ib_cells};

  const cs_lnum_t *face_vtx_idx[1] = {face_vertex_idx};
  const cs_lnum_t *face_vtx_ids[1] = {face_vertex_ids};

  fvm_nodal_from_desc_add_faces(ib_mesh,
                                1,
                                n_ib_cells,
                                nullptr,
                                1,
                                face_list_shift,
                                face_vtx_idx,
                                face_vtx_ids,
                                nullptr,
                                nullptr);

  fvm_nodal_transfer_vertices(ib_mesh,
                              reinterpret_cast<cs_coord_t *>(edge));
  edge = nullptr;

  /* Parallel numbering */

  fvm_io_num_t *face_io_num = nullptr;
  fvm_io_num_t *vtx_io_num = nullptr;

  if (cs_glob_n_ranks > 1) {

    /* Order faces by increasing global number */

    face_io_num = fvm_io_num_create_from_scan(n_ib_cells);
    const cs_gnum_t *face_gnum = fvm_io_num_get_global_num(face_io_num);
    fvm_nodal_order_faces(ib_mesh, face_gnum);
    fvm_nodal_init_io_num(ib_mesh, face_gnum, 2);
    fvm_io_num_destroy(face_io_num);

    /* Order vertices by increasing global number */

    vtx_io_num = fvm_io_num_create_from_scan(n_glob_vtx);
    const cs_gnum_t *vertex_gnum = fvm_io_num_get_global_num(vtx_io_num);
    fvm_nodal_order_vertices(ib_mesh, vertex_gnum);
    fvm_nodal_init_io_num(ib_mesh, vertex_gnum, 0);
    fvm_io_num_destroy(vtx_io_num);

  }

  cs_lnum_t *new_parent_id;
  CS_MALLOC(new_parent_id, n_ib_cells, cs_lnum_t);
  cs_lnum_t shift = cs_glob_mesh->n_b_faces;  // Same as n_faces_all here.
  for (cs_lnum_t i = 0; i < n_ib_cells; i++)
    new_parent_id[i] = i + shift;

  fvm_nodal_change_parent_id(ib_mesh, new_parent_id, 2);

  cs_post_assign_existing_mesh(cs_glob_porous_model_extra_faces->mesh_id,
                               ib_mesh,
                               0,
                               true);

  CS_FREE(new_parent_id);
  CS_FREE(face_vertex_ids);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
