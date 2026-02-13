/*============================================================================
 * Adaptive Mesh Refinement (AMR).
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_io_num.h"
#include "fvm/fvm_triangulate.h"

#include "alge/cs_cell_to_vertex.h"
#include "alge/cs_gradient.h"
#include "alge/cs_matrix_default.h"

#include "base/cs_array.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_coarsen.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_refine.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_field_default.h"
#include "base/cs_parall.h"
#include "base/cs_renumber.h"
#include "mesh/cs_redistribute.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_adaptive_refinement.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables.
 *============================================================================*/

static cs_amr_info_t _amr_info = {
  false,
  0,
  {CS_TIME_CONTROL_TIME_STEP, false, false, false, {-1}, {-1}, {-1},
   nullptr, nullptr, false, -1, -1, -1},
  nullptr,
  nullptr,
  nullptr,
  nullptr,
  0
};

cs_amr_info_t *cs_glob_amr_info = &_amr_info;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the cell refinement level array.
 *
 * parameters:
 *   f           <-- pointer to the mesh
 *   c_r_level   --> cell level refinement level
 *----------------------------------------------------------------------------*/

static void
_build_cell_r_level(cs_mesh_t *m,
                    int       *c_r_level)
{
  for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
    c_r_level[i] = 0;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;

  for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {
    for (cs_lnum_t i = 0; i < 2; i++) {
      cs_lnum_t c_id = i_face_cells[f_id][i];
      if (m->i_face_r_gen[f_id] > c_r_level[c_id])
        c_r_level[c_id] = m->i_face_r_gen[f_id];
    }
  }

  if (m->halo != NULL)
    cs_halo_sync_untyped(m->halo,
                         CS_HALO_STANDARD,
                         sizeof(int),
                         c_r_level);
}

/*----------------------------------------------------------------------------
 * In case of refinement during an AMR step, this function perform the
 * memory reallocation for a given field and interpolate its values on
 * the new mesh.
 *
 * parameters:
 *   f           <-- pointer to a field
 *   n_old       <-- old number of elements
 *   o2n_idx     <-- old to new numbering index
 *   new_idx     <-- list of new elements in the new mesh (elements not
 *                   issued from an already existing elements: new internal
 *                   faces created in a refined cell.)
 *   cog         <-- old mesh cells center of gravity
 *   cell_cen    <-- new mesh cells center of gravity
 *   measure     <-- new mesh element measure (volume for cells, surf for faces)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_update_field_refinement(cs_field_t        *f,
                                     const cs_lnum_t    n_old,
                                     const cs_lnum_t    o2n_idx[],
                                     const cs_lnum_t    new_idx[],
                                     const cs_real_3_t  cog[],
                                     const cs_real_3_t  cell_cen[],
                                     const cs_real_t    measure[])
{
  cs_lnum_t n_new = cs_mesh_location_get_n_elts(f->location_id)[0];
  cs_lnum_t n_alloc_new = cs_mesh_location_get_n_elts(f->location_id)[2];
  cs_field_t *f_vel = CS_F_(vel);

  cs_real_t *grad = f->grad;
  const cs_lnum_t dim = f->dim;

  int interpolation_type = _amr_info.interpolation;
  int _field_interp_type = f->get_key_int("amr_interpolation_scheme");

  for (int i = 0; i < f->n_time_vals; i++) {
    // New array for current time vals
    cs_array_2d<cs_real_t> *new_val_i_ptr
      = new cs_array_2d<cs_real_t>(n_alloc_new, f->dim);
    auto new_val_i = new_val_i_ptr->view();
    auto f_val_i = f->view(i);

    // Loop on old elements counts
    for (cs_lnum_t o_id = 0; o_id < n_old; o_id++) {

      // Get list of new elements for current old entity
      cs_lnum_t s_id = o2n_idx[o_id];
      cs_lnum_t e_id = o2n_idx[o_id+1];

      // If the field is extensive, me must weight its new
      // value according to local measure
      cs_real_t o_m = 0.0;
      cs_real_t *_measure_ratio = nullptr;
      cs_real_t measure_ratio_s[64]; // Enough for most cases
      cs_real_t *measure_ratio = measure_ratio_s;
      if (e_id-s_id > 64) {
        CS_MALLOC(_measure_ratio, e_id-s_id, cs_real_t);
        measure_ratio = _measure_ratio;
      }

      for (cs_lnum_t n_id = s_id; n_id < e_id; n_id++) {
        measure_ratio[n_id-s_id] = 1.0;
        o_m += measure[n_id];
      }

      if (f->type == CS_FIELD_EXTENSIVE) {
        for (cs_lnum_t n_id = s_id; n_id < e_id; n_id++)
          measure_ratio[n_id-s_id] = measure[n_id]/o_m;
      }

      /* Some fields (like velocity) are always extrapolated with gradient
       * otherwise specified by user for fields on cells. */
      if (f->location_id == CS_MESH_LOCATION_CELLS &&
          (_field_interp_type == 1 || interpolation_type == 1)) {

        for (cs_lnum_t n_id = s_id; n_id < e_id; n_id++) {
          cs_real_3_t d = {cell_cen[n_id][0] - cog[o_id][0],
                           cell_cen[n_id][1] - cog[o_id][1],
                           cell_cen[n_id][2] - cog[o_id][2]};

          for (cs_lnum_t j = 0; j < dim; j++) {
            cs_lnum_t k = o_id*dim + j;
            new_val_i(n_id, j) =   f_val_i(o_id, j)
                                 + cs_math_3_dot_product(d, grad+3*k);
            new_val_i(n_id, j) *= measure_ratio[n_id-s_id];
          }
        }
      }
      else {
        for (cs_lnum_t n_id = s_id; n_id < e_id; n_id++) {
          for (cs_lnum_t j = 0; j < dim; j++)
            new_val_i(n_id, j) = f_val_i(o_id, j) * measure_ratio[n_id-s_id];
        }
      }

      CS_FREE(_measure_ratio);
    }

    // We delete old ptr, and replace by the new one
    delete f->_vals[i];
    f->_vals[i] = new_val_i_ptr;
  }

  /* For internal faces or vertices fields, new elements are created
   * but do not depend on older elements.
   * For those elements, no interpolation could be performed as such. */

  if (new_idx != nullptr) {
    cs_field_t *f_imf = cs_field_by_name("inner_mass_flux");

    /* When field is defined on vertices, we do not know
     * for now what to impose on newly created vertices.
     * When field is on interior faces, a specific treatment
     * is done for mass fluxes, but other fields are set to
     * 0 for now. */
    if (  f->location_id != CS_MESH_LOCATION_INTERIOR_FACES
        || f->id != f_imf->id) {
      for (int i = 0; i < f->n_time_vals; i++) {
        auto f_view_i = f->view(i);
        for (cs_lnum_t idim = 0; idim < f->dim; idim ++) {
          for (cs_lnum_t n_id = new_idx[0]; n_id < n_new; n_id++) {
            f_view_i(n_id, idim) = 0.;
          }
        }
      }
    }
    else {
      // Specific treatment for interior mass fluxes
      const cs_lnum_2_t *i_face_cells = cs_glob_mesh->i_face_cells;
      cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

      const cs_real_t *weight = mq->weight;
      const cs_real_t *i_face_surf = mq->i_face_surf;
      const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;

      const cs_real_t *crom = CS_F_(rho)->val;
      const cs_real_3_t *vel = (const cs_real_3_t *)(CS_F_(vel)->val);

      for (int i = 0; i < f->n_time_vals; i++) {
        auto f_view_i = f->view(i);

        // Initialisation
        for (cs_lnum_t n_id = new_idx[0]; n_id < n_new; n_id++)
          f_view_i(n_id, 0) = 0;

        for (cs_lnum_t n_id = new_idx[0]; n_id < n_new; n_id++) {
          cs_lnum_t c_id1 = i_face_cells[n_id][0];
          cs_lnum_t c_id2 = i_face_cells[n_id][1];

          cs_real_t w_1 = weight[n_id];
          cs_real_t w_2 = (1. - weight[n_id]);

          cs_real_t q = 0.0;
          for (cs_lnum_t idim = 0; idim < 3; idim ++) {
            cs_real_t qdm1 = crom[c_id1]*vel[c_id1][idim];
            cs_real_t qdm2 = crom[c_id2]*vel[c_id2][idim];
            q +=   (w_1 * qdm1 + w_2 * qdm2)
                 * i_face_u_normal[n_id][idim];
          }

          f_view_i(n_id, 0) += q *  i_face_surf[n_id];
        }
      }
    }
  }

  /* Update public pointers (vals[x], val & val_pre) */
  f->update_public_pointers();
}

/*----------------------------------------------------------------------------
 * In case of refinement during an AMR step, this function perform the
 * memory reallocation for a given field bc_coeffs and interpolate its values
 * on the new mesh.
 *
 * parameters:
 *   f           <-- pointer to a field
 *   n_new       <-- new number of elements
 *   n_old       <-- old number of elements
 *   o2n_idx     <-- old to new numbering index
 *   new_idx     <-- list of new elements in the new mesh (elements not
 *                   issued from an already existing elements: new internal
 *                   faces created in a refined cell.)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_update_bc_coeffs_refinement(const cs_field_t  *f,
                                         const cs_lnum_t    n_old,
                                         const cs_lnum_t    o2n_idx[])
{
  cs_lnum_t n_new = cs_glob_mesh->n_b_faces;

  int n_array = 5;
  cs_real_t *array_list[5] = {f->bc_coeffs->val_f,
                              f->bc_coeffs->val_f_lim,
                              f->bc_coeffs->flux,
                              f->bc_coeffs->flux_lim,
                              f->bc_coeffs->val_f_pre};
  cs_real_t *array_out[5] = {nullptr};

  for (int k = 0; k < n_array; k++) {

    cs_real_t *bc_val   = array_list[k];
    cs_real_t *prev_ptr = nullptr;

    /* In case val_f = val_f_lim or flux_lim = flux, we do not
     * treat those arrays. Simply make them point towards val_f
     * or flux, resp. */
    if (k>0)
      prev_ptr = array_list[k-1];

    if (bc_val != nullptr && bc_val != prev_ptr) {
      cs_real_t *i_vals = nullptr;
      CS_MALLOC(i_vals, n_new * f->dim, cs_real_t);

      /* Fill vals with P0 interpolation */
      for (cs_lnum_t o_id = 0; o_id < n_old; o_id++) {
        for (cs_lnum_t n_id = o2n_idx[o_id];
             n_id < o2n_idx[o_id+1];
             n_id++) {
          for (cs_lnum_t idim = 0; idim < f->dim; idim ++)
            i_vals[f->dim * n_id + idim] = bc_val[f->dim * o_id + idim];
        }
      }

      CS_FREE(bc_val);
      array_out[k] = i_vals;
    }

  }

  /* Make vals member of the current field
     point toward the newly allocated one */
  f->bc_coeffs->val_f = array_out[0];

  if (array_list[1] == array_list[0])
    f->bc_coeffs->val_f_lim = array_out[0];
  else
    f->bc_coeffs->val_f_lim = array_out[1];

  f->bc_coeffs->flux = array_out[2];

  if (array_list[3] == array_list[2])
    f->bc_coeffs->flux_lim = array_out[2];
  else
    f->bc_coeffs->flux_lim = array_out[3];

  f->bc_coeffs->val_f_pre = array_out[4];
}

/*----------------------------------------------------------------------------
 * In case of refinement during an AMR step, this function perform the
 * memory reallocation for a given field and interpolate its values on
 * the new mesh.
 *
 * parameters:
 *   f           <-- pointer to a field
 *   n_old       <-- old number of elements
 *   n2o_idx     <-- new to old numbering index
 *   n2o         <-- new to old values
 *   measure     <-- element measure on old mesh (face surface for
 *                    internal and boundary faces, cell volume for cells)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_update_field_coarsening(cs_field_t      *f,
                                     [[maybe_unused]] const cs_lnum_t n_old,
                                     const cs_lnum_t  n2o_idx[],
                                     const cs_lnum_t  n2o[],
                                     const cs_real_t  measure[])
{
  cs_lnum_t n_new = cs_mesh_location_get_n_elts(f->location_id)[0];

  for (int i = 0; i < f->n_time_vals; i++) {

    /* If current field is the owner of a linked series, it handles the
     * the interpolation and resizing for everyone.
     */
    if (f->has_sub_fields()) {
      /* temporary deep copy of original data */
      auto old_val_i = f->_ns_vals[i]->get_deep_copy();

      /* update_size uses the field's location id to update size */
      f->update_size(i);
      auto new_val_i = f->ns_view(i);

      /* Temporary work array of size n_fields*2 */
      cs_array<cs_real_t> val_sum(f->series_size());
      cs_array<cs_real_t> val_mean(f->series_size());

      for (cs_lnum_t idim = 0; idim < f->dim; idim ++) {

        /* Fill vals with P0 interpolation */
        for (cs_lnum_t n_id = 0; n_id < n_new; n_id++) {
          cs_lnum_t s_id = n2o_idx[n_id], e_id = n2o_idx[n_id+1];

          for (cs_lnum_t tmp_id = 0; tmp_id < f->series_size(); tmp_id++) {
            val_sum[tmp_id] = 0.;
            val_mean[tmp_id] = 0.;
          }

          cs_real_t measure_tot = 0.0;

          /* Compute mean value on old cells */
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t o_id = n2o[j];
            measure_tot += measure[o_id];
            for (cs_lnum_t k = 0; k < f->series_size(); k++) {
              val_mean(k) += old_val_i(k, o_id, idim) * measure[o_id];
              val_sum(k) += old_val_i(k, o_id, idim);
            }
          }
          if (f->type == CS_FIELD_EXTENSIVE)
            for (cs_lnum_t k = 0; k < f->series_size(); k++)
              new_val_i(k, n_id, idim) = val_sum(k);
          else
            for (cs_lnum_t k = 0; k < f->series_size(); k++)
              new_val_i(k, n_id, idim) = val_mean(k) / measure_tot;
        }

      }

    } /* Field is "just" an owner, handle resizing + interpolation */
    else if (f->owner()) {
      /* Copy values to a temporary array and update size of field array. */
      auto old_val_i = f->_vals[i]->get_deep_copy();

      /* update_size uses the field's location id to update size */
      f->update_size(i);
      auto new_val_i = f->view(i);

      for (cs_lnum_t idim = 0; idim < f->dim; idim ++) {

        /* Fill vals with P0 interpolation */
        for (cs_lnum_t n_id = 0; n_id < n_new; n_id++) {
          cs_lnum_t s_id = n2o_idx[n_id], e_id = n2o_idx[n_id+1];

          cs_real_t measure_tot = 0.0;
          cs_real_t val_mean    = 0.0;
          cs_real_t val_sum     = 0.0;

          /* Compute mean value on old cells */
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t o_id = n2o[j];
             measure_tot += measure[o_id];
             val_mean    += old_val_i(o_id,idim) * measure[o_id];
             val_sum     += old_val_i(o_id,idim);
          }
          if (f->type == CS_FIELD_EXTENSIVE)
            new_val_i(n_id, idim) = val_sum;
          else
            new_val_i(n_id, idim) = val_mean/measure_tot;
        }

      }
    }
  }

  /* Update public pointers (vals[x], val & val_pre) */
  f->update_public_pointers();

  /* If current field is a series owner, update data mapping
   * of linked fields. This method updates their cs_array objects
   * and also calls "update_public_pointers" so that val/val_pre/vals are
   * pointing to correct values.
   */
  if (f->is_series_owner())
    f->update_sub_fields_mapping();
}

/*----------------------------------------------------------------------------
 * In case of refinement during an AMR step, this function perform the
 * memory reallocation for a given field bc_coeffs and interpolate its values
 * on the new mesh.
 *
 * parameters:
 *   f           <-- pointer to a field
 *   location_id <-- location id of the field that has to be treated
 *   n_old       <-- old number of elements
 *   n2o_idx     <-- new to old numbering index
 *   n2o         <-- new to old values
 *   measure     <-- element measure on old mesh (face surface for
 *                    internal and boundary faces, cell volume for cells)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_update_bc_coeffs_coarsening(cs_field_t      *f,
                                         [[maybe_unused]] const cs_lnum_t n_old,
                                         const cs_lnum_t  n2o_idx[],
                                         const cs_lnum_t  n2o[],
                                         const cs_real_t  measure[])
{
  cs_lnum_t n_new = cs_glob_mesh->n_b_faces;

  int n_array = 5;
  cs_real_t *array_list[5] = {f->bc_coeffs->val_f,
                              f->bc_coeffs->val_f_lim,
                              f->bc_coeffs->flux,
                              f->bc_coeffs->flux_lim,
                              f->bc_coeffs->val_f_pre};
  cs_real_t *array_out[5] = {nullptr};

  for (int k = 0; k < n_array; k++) {

    cs_real_t *bc_val     = array_list[k];
    cs_real_t *prev_ptr   = nullptr;
    /*In case val_f = val_f_lim or flux_lim = flux, we do not
     * treat those arrays. Simply make them point towards val_f
     * or flux, resp. */
    if (k>0)
      prev_ptr = array_list[k-1];

    if (bc_val != nullptr && bc_val != prev_ptr) {

      cs_real_t *i_vals = nullptr;
      cs_real_t *val_mean = nullptr;

      /* Allocate a new vals array */
      CS_MALLOC(i_vals, n_new * f->dim, cs_real_t);
      CS_MALLOC(val_mean, f->dim, cs_real_t);

      /* Fill vals with P0 interpolation */
      for (cs_lnum_t n_id = 0; n_id < n_new; n_id++) {
        cs_lnum_t s_id = n2o_idx[n_id], e_id = n2o_idx[n_id+1];

        cs_real_t measure_tot = 0.0;
        for (cs_lnum_t idim = 0; idim < f->dim; idim ++)
          val_mean[idim] = 0.0;

        for (cs_lnum_t i = s_id; i < e_id; i++) {
          cs_lnum_t o_id = n2o[i];
          /* Compute mean value */
           measure_tot += measure[o_id];
           for (cs_lnum_t idim = 0; idim < f->dim; idim ++)
             val_mean[idim] += bc_val[f->dim * o_id  + idim]
                                      * measure[o_id];
        }

        for (cs_lnum_t idim = 0; idim < f->dim; idim ++)
          i_vals[f->dim * n_id + idim] = val_mean[idim]/measure_tot;
      }

      CS_FREE(bc_val);
      CS_FREE(val_mean);
      array_out[k] = i_vals;
    }
  }

  /* Make vals member of the current field
     point toward the newly allocated one*/
  f->bc_coeffs->val_f = array_out[0];

  if (array_list[1] == array_list[0])
    f->bc_coeffs->val_f_lim = array_out[0];
  else
    f->bc_coeffs->val_f_lim = array_out[1];

  f->bc_coeffs->flux = array_out[2];

  if (array_list[3] == array_list[2])
    f->bc_coeffs->flux_lim = array_out[2];
  else
    f->bc_coeffs->flux_lim = array_out[3];

  f->bc_coeffs->val_f_pre = array_out[4];
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Sync a standard field
 */
/*--------------------------------------------------------------------------*/
static void
_halo_sync_std_owner
(
  cs_field_t      *f,   /*!<[in] pointer to cs_field_t */
  const cs_halo_t *halo /*!<[in] pointer to cs_halo_t */
)
{
  for (int kk = 0; kk < f->n_time_vals; kk++) {
    cs_halo_sync_untyped(halo,
                         CS_HALO_EXTENDED,
                         f->dim*sizeof(cs_real_t),
                         f->vals[kk]);
    if (f->dim == 3)
      cs_halo_perio_sync_var_vect(halo,
                                  CS_HALO_EXTENDED,
                                  f->vals[kk],
                                  f->dim);
  }
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Sync a field which is a series owner and has subfields
 */
/*--------------------------------------------------------------------------*/
static void
_halo_sync_series_owner
(
  cs_field_t      *f,   /*!<[in] pointer to cs_field_t */
  const cs_halo_t *halo /*!<[in] pointer to cs_halo_t */
)
{
  for (int kk = 0; kk < f->n_time_vals; kk++) {
    auto f_view = f->ns_view(kk);
    /* We can loop on series size since its equal to 1 + n_sub_fields */
    for (int i_s = 0; i_s < f->series_size(); i_s++) {
      cs_real_t *_vals = f_view.sub_array(i_s);

      cs_halo_sync_untyped(halo,
                           CS_HALO_EXTENDED,
                           f->dim*sizeof(cs_real_t),
                           _vals);
      if (f->dim == 3)
        cs_halo_perio_sync_var_vect(halo,
                                    CS_HALO_EXTENDED,
                                    _vals,
                                    f->dim);
    }
  }
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Sync a field
 */
/*--------------------------------------------------------------------------*/
static void
_halo_sync
(
  cs_field_t *f /*!<[in] pointer to cs_field_t */
)
{
  const cs_halo_t *halo = cs_glob_mesh->halo;
  if (halo != nullptr) {
    /* Check if series owner, handle then case for all fields.
     * If not, check if simply owner. Else, nothing to do.
     */
    if (f->has_sub_fields()) {
      _halo_sync_series_owner(f, halo);
    }
    else if (f->owner()) {
      _halo_sync_std_owner(f, halo);
    }
  }
}

/*----------------------------------------------------------------------------
 * In case of refinement during an AMR step, this function perform the
 * memory reallocation for all fields and interpolate fields values on
 * the new mesh.
 *
 * parameters:
 *   mesh        <-- pointer to the mesh
 *   location_id <-- location id of the fields that has to be treated
 *   n_i_elts    <-- old number of elements
 *   o2n_idx     <-- old to new numbering index
 *   new_idx     <-- list of new elements in the new mesh (elements not
 *                   issued from an already existing elements: new internal
 *                   faces created in a refined cell.)
 *   n_i_b_elts  <-- old number of boundary elements (when interpolating field,
 *                   we need to interpolate some bc_coeffs arrays)
 *   o2n_b_idx   <-- old to new numbering index of boundary elements
 *   o_cog       <-- old mesh cells center of gravity
 *   n_cog       <-- new mesh cells center of gravity
 *   measure     <-- new mesh element measure (volume for cells, surf for faces)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_interp_after_refinement([[maybe_unused]] cs_mesh_t  *mesh,
                                     int                location_id,
                                     const cs_lnum_t    n_i_elts,
                                     const cs_lnum_t    o2n_idx[],
                                     const cs_lnum_t    new_idx[],
                                     const cs_real_3_t  o_cog[],
                                     const cs_real_3_t  n_cog[],
                                     const cs_real_t    measure[])
{
  const int n_fields = cs_field_n_fields();

  /* Old elements counts */
  const cs_lnum_t n_old = n_i_elts;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    const cs_mesh_location_type_t field_loc_type
      = cs_mesh_location_get_type(f->location_id);

    /* Handle refinement for series owners only */
    if (field_loc_type == location_id && f->owner()) {
      _realloc_and_update_field_refinement(f, n_old, o2n_idx, new_idx,
                                           o_cog, n_cog, measure);

      /* When dealing with fields on cells, we need to sync halo cells */
      if (field_loc_type == CS_MESH_LOCATION_CELLS)
        _halo_sync(f);
    }

    if (   location_id == CS_MESH_LOCATION_BOUNDARY_FACES
        && field_loc_type == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != nullptr)
      _realloc_and_update_bc_coeffs_refinement(f, n_old, o2n_idx);
  }
}

/*----------------------------------------------------------------------------
 * In case of coarsening during an AMR step, this function perform the
 * memory reallocation for all fields and interpolate fields values on
 * the new mesh.
 *
 * parameters:
 *   mesh        <-- pointer to the mesh
 *   location_id <-- location id of the fields that has to be treated
 *   n_i_elts    <-- old number of elements
 *   n2o_idx     <-- new to old numbering index
 *   n2o         <-- new to old values
 *   measure     <-- element measure on old mesh (face surface for
 *                   internal and boundary faces, cell volume for cells)
 *----------------------------------------------------------------------------*/

static void
_realloc_and_interp_after_coarsening
(
  [[maybe_unused]] cs_mesh_t       *mesh,
  int                               location_id,
  [[maybe_unused]] const cs_lnum_t  n_i_elts,
  const cs_lnum_t                   n2o_idx[],
  const cs_lnum_t                   n2o[],
  cs_real_t                         measure[]
)
{
  const int n_fields = cs_field_n_fields();
  const cs_lnum_t n_old = n_i_elts;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);

    const cs_mesh_location_type_t field_loc_type
      = cs_mesh_location_get_type(f->location_id);

    if ( field_loc_type == location_id ) {
      _realloc_and_update_field_coarsening(f, n_old, n2o_idx, n2o, measure);

      /* When dealing with fields on cells, we need to sync halo cells*/
      if (field_loc_type == CS_MESH_LOCATION_CELLS)
        _halo_sync(f);
    }

    if (   location_id == CS_MESH_LOCATION_BOUNDARY_FACES
        && field_loc_type == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != nullptr)
        _realloc_and_update_bc_coeffs_coarsening(f,
                                                 n_old,
                                                 n2o_idx,
                                                 n2o,
                                                 measure);
  }
}

/*----------------------------------------------------------------------------
 * Refinement operations could introduce jumps in cell refinement level of more
 * than 1. This could be numerically unstable. This subroutine flags the cells
 * that would lead to jump of more than one so that they are also refined.
 *
 * parameters:
 *   mesh        <-- pointer to the mesh
 *   c_r_level   <-- cell refinement level before the refinement process
 *   indic       <-> initial array containing the flagged cells
 *----------------------------------------------------------------------------*/

static void
_propagate_refinement(cs_mesh_t  *m,
                      const int   c_r_level[],
                      int         indic[])
{
  cs_halo_sync_untyped(m->halo, CS_HALO_STANDARD, sizeof(int), indic);

  int reloop = 0;

  do {
    reloop = 0;

    for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id ++) {
      cs_lnum_t c_id1 = m->i_face_cells[f_id][0];
      cs_lnum_t c_id2 = m->i_face_cells[f_id][1];

      int r_level_1 = indic[c_id1] + c_r_level[c_id1];
      int r_level_2 = indic[c_id2] + c_r_level[c_id2];

      if (r_level_1 > r_level_2 + 1) {
        indic[c_id2]++;
        reloop = 1;
      }

      if (r_level_2 > r_level_1 + 1) {
        indic[c_id1]++;
        reloop = 1;
      }
    }

    cs_parall_max(1, CS_INT_TYPE, &reloop);
    if (m->halo != nullptr) {
      cs_halo_sync_untyped(m->halo,
                           CS_HALO_STANDARD,
                           sizeof(int),
                           indic);
    }
  } while (reloop);
}

/*----------------------------------------------------------------------------
 * Write in a CSV file the load imbalance of all ranks.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

static void
_write_load_balance_stats(void)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  int tag = 0;
  int rank = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;
  MPI_Comm comm = cs_glob_mpi_comm;

  if (rank == 0) {
    float *imbalance;
    CS_MALLOC(imbalance, n_ranks, float);

    imbalance[0] = mesh->n_cells;
    for (int i = 1; i < n_ranks; i++) {
      MPI_Recv(&imbalance[i], 1, MPI_FLOAT, i, tag, comm, MPI_STATUS_IGNORE);
    }

    float balance = mesh->n_g_cells / cs_glob_n_ranks;

    for (int i = 0; i < n_ranks; i++)
      imbalance[i] = fabsf(imbalance[i] - balance) / balance * 100;

    char fname[256] = "load_balance_stats.csv";
    FILE *fh;
    static int balance_iter = 0;
    balance_iter++;
    if (balance_iter == 1) {
      fh = fopen(fname, "w");
      fprintf(fh, "time ");
      for (int i = 0; i < n_ranks; i++)
        fprintf(fh, "P%d ", i);
      fprintf(fh, "\n");
    }
    else {
      fh = fopen(fname, "a");
    }

    fprintf(fh, "%d ", balance_iter);

    for (int i = 0; i < n_ranks; i++)
      fprintf(fh, "%f ", imbalance[i]);
    fprintf(fh, "\n");

    fclose(fh);

  }
  else {
    float to_send = (float)mesh->n_cells;
    MPI_Send(&to_send, 1, MPI_FLOAT, 0, tag, comm);
  }
}

/*----------------------------------------------------------------------------
 * Redistribute the mesh so as to ensure a balanced load of cells between
 * all the procs.
 *
 * For now, a cell distribution map is computed based on a Morton space-filling
 * curve.
 *
 * TODO: redistribute using the algorithm initially used to read-in the mesh
 * in cs_partition.
 *
 * parameters:
 *   write_stats <-- switch to write load imbalance stats
 *----------------------------------------------------------------------------*/

static void
_load_balance(bool write_stats)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_coord_t *cell_center = (cs_coord_t *)mq->cell_cen;
  int dim = 3;
  fvm_io_num_sfc_t sfc_type = FVM_IO_NUM_SFC_MORTON_CUBE;

  fvm_io_num_t *cell_io_num = fvm_io_num_create_from_sfc(cell_center,
                                                         dim,
                                                         mesh->n_cells,
                                                         sfc_type);

  const cs_gnum_t *cell_num = fvm_io_num_get_global_num(cell_io_num);

  int n_ranks = cs_glob_n_ranks;

  int *cell_rank = nullptr;
  CS_MALLOC(cell_rank, mesh->n_cells, int);

  // Non-uniform block size.

  cs_gnum_t cells_per_rank = mesh->n_g_cells / n_ranks;
  cs_lnum_t rmdr = mesh->n_g_cells - cells_per_rank * (cs_gnum_t)n_ranks;
  if (rmdr == 0) {
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      cell_rank[i] = (cell_num[i] - 1) / cells_per_rank;
  }
  else {
    cs_gnum_t n_ranks_rmdr = n_ranks - rmdr;
    cs_gnum_t n_ranks_cells_per_rank = n_ranks_rmdr * cells_per_rank;
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      if ((cell_num[i] - 1) < n_ranks_cells_per_rank)
        cell_rank[i] = (cell_num[i] - 1) / cells_per_rank;
      else
        cell_rank[i] = (cell_num[i] + n_ranks_rmdr - 1) / (cells_per_rank + 1);
    }
  }

  fvm_io_num_destroy(cell_io_num);

  cs_redistribute(cell_rank, nullptr);

  CS_FREE(cell_rank);

  if (write_stats)
    _write_load_balance_stats();
}

#endif // defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Refinement step function.
 *----------------------------------------------------------------------------*/

static void
_refine_step(void)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  const int c_r_level_max = _amr_info.n_layers;
  cs_lnum_t   n_selected_cells = 0;
  cs_lnum_t  *selected_cells = nullptr;
  int        *s = nullptr;
  int        *indic = nullptr;
  int        *c_r_level = nullptr;

  CS_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);
  CS_MALLOC(s, mesh->n_cells, int);
  CS_MALLOC(indic, mesh->n_cells_with_ghosts, int);
  CS_MALLOC(c_r_level, mesh->n_cells_with_ghosts, int);

  /* Determine current cell refinement level */

  _build_cell_r_level(mesh, c_r_level);

  /* Compute the user criteria */

  memset(s, 0, mesh->n_cells*sizeof(int));
  memset(indic, 0, mesh->n_cells_with_ghosts*sizeof(int));

  _amr_info.indic_func(_amr_info.indic_input, s);

  /* Flag cells and eventually propagate refinement
   * to avoid jump in refinement levels*/

  int n_threads = cs_parall_n_threads(mesh->n_cells, CS_THR_MIN);
# pragma omp parallel for  num_threads(n_threads)
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id ++) {
    if (s[c_id] > 0 && c_r_level[c_id] < c_r_level_max)
      indic[c_id] = 1;
    else
      indic[c_id] = 0;
  }

  _propagate_refinement(mesh, c_r_level, indic);

  /* Refinement */

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id ++) {
    if (indic[c_id] > 0) {
      selected_cells[n_selected_cells] = c_id;
      n_selected_cells++;
    }
  }

  cs_mesh_refine_simple_selected(mesh,
                                 false,
                                 n_selected_cells,
                                 selected_cells);

  CS_FREE(selected_cells);
  CS_FREE(s);
  CS_FREE(indic);
  CS_FREE(c_r_level);

  /* Renumber mesh based on code options */

  cs_renumber_mesh(mesh, nullptr, nullptr, nullptr, nullptr);

  /* Free and re-compute mesh related quantities */

  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);
  cs_mesh_quantities_compute(mesh, cs_glob_mesh_quantities);

  /* Reset boundary conditions */

  cs_boundary_conditions_realloc();
  cs_field_map_and_init_bcs();
}

/*----------------------------------------------------------------------------
 * Coarsening step function.
 *----------------------------------------------------------------------------*/

static void
_coarsen_step(void)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  cs_lnum_t   n_selected_cells = 0;
  cs_lnum_t  *selected_cells = nullptr;
  int        *s = nullptr;
  int        *indic = nullptr;

  CS_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);
  CS_MALLOC(s, mesh->n_cells, int);
  CS_MALLOC(indic, mesh->n_cells_with_ghosts, int);

  /* Compute the user criteria */

  memset(s, 0, mesh->n_cells*sizeof(int));
  memset(indic, 0, mesh->n_cells_with_ghosts*sizeof(int));

  _amr_info.indic_func(_amr_info.indic_input, s);

  /* First we tag cell according to indicator.
   * Eventually, block coarsening operations
   * to avoid jump in cell refinement levels */

  int n_threads = cs_parall_n_threads(mesh->n_cells, CS_THR_MIN);
# pragma omp parallel for  num_threads(n_threads)
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id ++) {
    if (s[c_id] > 0)
      indic[c_id] = 0;
    else
      indic[c_id] = 1;
  }

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id ++) {
    if (indic[c_id] > 0) {
      selected_cells[n_selected_cells] = c_id;
      n_selected_cells++;
    }
  }

  cs_mesh_coarsen_simple_selected(mesh,
                                  n_selected_cells,
                                  selected_cells);

  CS_FREE(selected_cells);
  CS_FREE(s);
  CS_FREE(indic);

  /* Renumber mesh based on code options */

  cs_renumber_mesh(mesh, nullptr, nullptr, nullptr, nullptr);

  /* Free and compute mesh related quantities */

  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);
  cs_mesh_quantities_compute(cs_glob_mesh,cs_glob_mesh_quantities);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_update_selectors(cs_glob_mesh);
  cs_mesh_location_build(cs_glob_mesh, -1);
  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);

  /* Update for boundary conditions */

  cs_boundary_conditions_realloc();
  cs_field_map_and_init_bcs();
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activation and Initializing af the global AMR strcture
 *
 * \param[in]  n_layers     number of layers of refinement around user criteria
 * \param[in]  nt_interval  time_step interval between adaptation steps
 * \param[in]  indic_func   user function marking cells that should be refined
 * \param[in]  indic_input  input for indic_func or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_adaptive_refinement_define(int                  n_layers,
                              int                  nt_interval,
                              cs_amr_indicator_t  *indic_func,
                              const void          *indic_input,
                              int                  interpolation)
{
  _amr_info.is_set = true;
  _amr_info.n_layers = n_layers;
  cs_time_control_init_by_time_step(&_amr_info.time_control,
                                    -1, -1, nt_interval,
                                    false, false);

  _amr_info.indic_func = indic_func;
  _amr_info.indic_input = indic_input;
  _amr_info.fields_interp_refinement_func = _realloc_and_interp_after_refinement;
  _amr_info.fields_interp_coarsening_func = _realloc_and_interp_after_coarsening;
  _amr_info.interpolation = interpolation;

  cs_glob_mesh->time_dep = CS_MESH_TRANSIENT_CONNECT;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute field gradients when necessary (for fields interpolation)
 */
/*----------------------------------------------------------------------------*/

void
cs_adaptive_refinement_update_gradients(void)
{
  const int n_fields = cs_field_n_fields();
  int interpolation_type = _amr_info.interpolation;
  cs_field_t *fv = cs_field_by_name("velocity");

  /* Velocity field is always interpolated using P1
   * to reconstruct mass flux at newly created internal
   * faces */
  cs_field_allocate_gradient(fv);

  cs_real_33_t *gradu = (cs_real_33_t *)fv->grad;
  cs_field_gradient_vector(fv,
                           false,  /* use_previous_t */
                           1,     /* inc */
                           gradu);

  /* In cas of P1 interpolation requested by user, we compute
   * all fields gradients for cell based fields.
   * This is also done to update the gradient in case it was
   * set by user. */

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);

    const cs_mesh_location_type_t field_loc_type
      = cs_mesh_location_get_type(f->location_id);

    if (   field_loc_type == CS_MESH_LOCATION_CELLS
        && f->id != fv->id
        && (interpolation_type == 1 || f->grad != nullptr)) {
      cs_field_allocate_gradient(f);

      switch (f->dim) {
      case 1:
        {
          cs_real_3_t *grad = (cs_real_3_t *)f->grad;
          cs_field_gradient_scalar(f, false, 1, grad);
        }
        break;
      case 3:
        {
          cs_real_33_t *gradv = (cs_real_33_t *)f->grad;
          cs_field_gradient_vector(f, false, 1, gradv);
        }
        break;
      case 6:
        {
          cs_real_63_t *gradt = (cs_real_63_t *)f->grad;
          cs_field_gradient_tensor(f, false, 1, gradt);
        }
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  _("Field '%s' has dimension %d, which is not compatible\n"
                    "with P1 interpolation for AMR for now..\n"),
                  f->name, f->dim);
      }

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free field gradients.
 */
/*----------------------------------------------------------------------------*/

void
cs_adaptive_refinement_free_gradients(void)
{
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->grad != nullptr)
      CS_FREE(f->grad);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a refinement / coarsening step
 */
/*----------------------------------------------------------------------------*/

void
cs_adaptive_refinement_step(void)
{
  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "Starting mesh adaptation (AMR)\n"
                  "------------------------------\n\n"));

  _refine_step();

  _coarsen_step();

  /* Post-adaptation updates */

  cs_mesh_adjacencies_update_mesh();
  cs_gradient_free_quantities();
  cs_matrix_update_mesh();

  /* Free gradients
   * WARNING : for now stored gradients are not necessary
   * but one should remain cautious here if gradient are
   * needed afterwards */

  cs_adaptive_refinement_free_gradients();

  /* Perform load balancing */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    bool write_stats = true;
    _load_balance(write_stats);
  }
#endif

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "Completed mesh adaptation\n"
                  "-------------------------\n\n"));
}

/*----------------------------------------------------------------------------*/
