/*============================================================================
 * Gradient reconstruction.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <chrono>

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "alge/cs_bad_cells_regularisation.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary_conditions.h"
#include "alge/cs_cell_to_vertex.h"
#include "base/cs_dispatch.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "alge/cs_gradient_boundary.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_porous_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_timer.h"
#include "base/cs_timer_stats.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_gradient.h"
#include "alge/cs_gradient_priv.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_gradient.cpp
 * \brief Gradient reconstruction.
 *
 * Please refer to the
 * <a href="../../theory.pdf#gradreco"><b>gradient reconstruction</b></a>
 * section of the theory guide for more informations.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* SIMD unit size to ensure SIMD alignement (2 to 8 required on most
 * current architectures, so 16 should be enough on most architectures) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure associated to gradient quantities management */

typedef struct {

  cs_real_33_t  *cocg_it;          /* Interleaved cocg matrix
                                      for iterative gradients */

  cs_cocg_6_t   *cocgb_s_lsq;      /* coupling of gradient components for
                                      least-square reconstruction at boundary */
  cs_cocg_6_t   *cocg_lsq;         /* Interleaved cocg matrix
                                      for least square gradients */

  cs_cocg_6_t   *cocgb_s_lsq_ext;  /* coupling of gradient components for
                                      least-square reconstruction at boundary */
  cs_cocg_6_t   *cocg_lsq_ext;     /* Interleaved cocg matrix for least
                                      squares gradients with ext. neighbors */

} cs_gradient_quantities_t;

/* Basic per gradient computation options and logging */
/*----------------------------------------------------*/

typedef struct _cs_gradient_info_t {

  char                *name;               /* System name */
  cs_gradient_type_t   type;               /* Gradient type */

  unsigned             n_calls;            /* Number of times system solved */

  int                  n_iter_min;         /* Minimum number of iterations */
  int                  n_iter_max;         /* Minimum number of iterations */
  unsigned long        n_iter_tot;         /* Total number of iterations */

  cs_timer_counter_t   t_tot;              /* Total time used */

} cs_gradient_info_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int _gradient_n_systems = 0;      /* Current number of systems */
static int _gradient_n_max_systems = 0;  /* Max. number of systems for
                                            _gradient_systems. */

/* System info array */
static cs_gradient_info_t **_gradient_systems = nullptr;

/* Short names for gradient computation types */

const char *cs_gradient_type_name[]
  = {N_("Green-Gauss, iterative handling of non-orthogonalities"),
     N_("Least-squares"),
     N_("Green-Gauss, least-squares gradient face values"),
     N_("Green-Gauss, vertex-based face interpolation"),
     N_("Green-Gauss, multiplied by renormalization")};

/* Timer statistics */

static cs_timer_counter_t   _gradient_t_tot;     /* Total time in gradients */
static int _gradient_stat_id = -1;

/* Gradient quantities */

static int                        _n_gradient_quantities = 0;
static cs_gradient_quantities_t  *_gradient_quantities = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 symmetric matrix (with symmetric storage)
 *         using Cramer's rule
 *
 * \param[in, out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
_math_6_inv_cramer_sym(cs_cocg_t  a[6],
                       cs_cocg_t  cocg[6])
{
  cs_real_t a00 = a[1]*a[2] - a[4]*a[4];
  cs_real_t a01 = a[4]*a[5] - a[3]*a[2];
  cs_real_t a02 = a[3]*a[4] - a[1]*a[5];
  cs_real_t a11 = a[0]*a[2] - a[5]*a[5];
  cs_real_t a12 = a[3]*a[5] - a[0]*a[4];
  cs_real_t a22 = a[0]*a[1] - a[3]*a[3];

  double det_inv = 1. / (a[0]*a00 + a[3]*a01 + a[5]*a02);

  cocg[0] = a00 * det_inv;
  cocg[1] = a11 * det_inv;
  cocg[2] = a22 * det_inv;
  cocg[3] = a01 * det_inv;
  cocg[4] = a12 * det_inv;
  cocg[5] = a02 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 symmetric matrix (with symmetric storage)
 *         in place, using Cramer's rule
 *
 * \param[in, out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
_math_6_inv_cramer_sym_in_place(cs_cocg_t  a[6])
{
  cs_real_t a00 = a[1]*a[2] - a[4]*a[4];
  cs_real_t a01 = a[4]*a[5] - a[3]*a[2];
  cs_real_t a02 = a[3]*a[4] - a[1]*a[5];
  cs_real_t a11 = a[0]*a[2] - a[5]*a[5];
  cs_real_t a12 = a[3]*a[5] - a[0]*a[4];
  cs_real_t a22 = a[0]*a[1] - a[3]*a[3];

  double det_inv = 1. / (a[0]*a00 + a[3]*a01 + a[5]*a02);

  a[0] = a00 * det_inv;
  a[1] = a11 * det_inv;
  a[2] = a22 * det_inv;
  a[3] = a01 * det_inv;
  a[4] = a12 * det_inv;
  a[5] = a02 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a gradient quantities structure, adding one if needed
 *
 * \param[in]  id  id of structure to return

 * \return  pointer to gradient quantities structure
 */
/*----------------------------------------------------------------------------*/

static cs_gradient_quantities_t  *
_gradient_quantities_get(int  id)
{
  assert(id >= 0);

  if (id >= _n_gradient_quantities) {

    CS_REALLOC(_gradient_quantities, id+1, cs_gradient_quantities_t);

    for (int i = _n_gradient_quantities; i < id+1; i++) {
      cs_gradient_quantities_t  *gq = _gradient_quantities + i;

      gq->cocg_it = nullptr;
      gq->cocgb_s_lsq = nullptr;
      gq->cocg_lsq = nullptr;
      gq->cocgb_s_lsq_ext = nullptr;
      gq->cocg_lsq_ext = nullptr;
    }

    _n_gradient_quantities = id+1;

  }

  return _gradient_quantities + id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_quantities_destroy(void)
{
  for (int i = 0; i < _n_gradient_quantities; i++) {

    cs_gradient_quantities_t  *gq = _gradient_quantities + i;

    CS_FREE(gq->cocg_it);
    CS_FREE(gq->cocgb_s_lsq);
    CS_FREE(gq->cocg_lsq);
    CS_FREE(gq->cocgb_s_lsq_ext);
    CS_FREE(gq->cocg_lsq_ext);

  }

  CS_FREE(_gradient_quantities);
  _n_gradient_quantities = 0;
}

/*----------------------------------------------------------------------------
 * Return pointer to new gradient computation info structure.
 *
 * parameters:
 *   name <-- system name
 *   type <-- resolution method
 *
 * returns:
 *   pointer to newly created linear system info structure
 *----------------------------------------------------------------------------*/

static cs_gradient_info_t *
_gradient_info_create(const char          *name,
                      cs_gradient_type_t   type)
{
  cs_gradient_info_t *new_info = nullptr;

  CS_MALLOC(new_info, 1, cs_gradient_info_t);
  CS_MALLOC(new_info->name, strlen(name) + 1, char);

  strcpy(new_info->name, name);
  new_info->type = type;

  new_info->n_calls = 0;
  new_info->n_iter_min = 0;
  new_info->n_iter_max = 0;
  new_info->n_iter_tot = 0;

  CS_TIMER_COUNTER_INIT(new_info->t_tot);

  return new_info;
}

/*----------------------------------------------------------------------------
 * Destroy gradient computation info structure.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure pointer
 *----------------------------------------------------------------------------*/

static void
_gradient_info_destroy(cs_gradient_info_t  **this_info)
{
  if (*this_info != nullptr) {
    CS_FREE((*this_info)->name);
    CS_FREE(*this_info);
  }
}

/*----------------------------------------------------------------------------
 * Update the number of sweeps for gradient information
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure
 *   n_iter    <-- number of iterations
 *----------------------------------------------------------------------------*/

static void
_gradient_info_update_iter(cs_gradient_info_t  *this_info,
                           int                  n_iter)
{
  if (n_iter > this_info->n_iter_max) {
    this_info->n_iter_max = n_iter;
    /* for first pass: */
    if (this_info->n_calls == 0)
      this_info->n_iter_min = n_iter;
  }
  else if (n_iter < this_info->n_iter_min)
    this_info->n_iter_min = n_iter;

  this_info->n_iter_tot += n_iter;
}

/*----------------------------------------------------------------------------
 * Output information regarding gradient computation.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure
 *----------------------------------------------------------------------------*/

static void
_gradient_info_dump(cs_gradient_info_t *this_info)
{
  int n_calls = this_info->n_calls;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Summary of gradient computations for \"%s\":\n\n"
                  "  Reconstruction type:   %s\n"
                  "  Number of calls:       %d\n"),
                this_info->name, cs_gradient_type_name[this_info->type],
                n_calls);
  if (this_info->n_iter_tot > 0)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Number of iterations:  %d mean, %d min., %d max.\n"),
                  (int)(this_info->n_iter_tot / (unsigned long)n_calls),
                  this_info->n_iter_min,
                  this_info->n_iter_max);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  Total elapsed time:    %.3f\n"),
                this_info->t_tot.nsec*1e-9);
}

/*----------------------------------------------------------------------------
 * Return pointer to gradient computation info.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name --> system name
 *   type --> resolution method
 *----------------------------------------------------------------------------*/

static cs_gradient_info_t *
_find_or_add_system(const char          *name,
                    cs_gradient_type_t   type)
{
  int ii, start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find system */

  start_id = 0;
  end_id = _gradient_n_systems - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((_gradient_systems[mid_id])->name, name);
    if (cmp_ret == 0)
      cmp_ret = (_gradient_systems[mid_id])->type - type;
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If found, return */

  if (cmp_ret == 0)
    return _gradient_systems[mid_id];

  /* Reallocate global array if necessary */

  if (_gradient_n_systems >= _gradient_n_max_systems) {

    if (_gradient_n_max_systems == 0)
      _gradient_n_max_systems = 10;
    else
      _gradient_n_max_systems *= 2;

    CS_REALLOC(_gradient_systems,
               _gradient_n_max_systems,
               cs_gradient_info_t *);

  }

  /* Insert in sorted list */

  for (ii = _gradient_n_systems; ii > mid_id; ii--)
    _gradient_systems[ii] = _gradient_systems[ii - 1];

  _gradient_systems[mid_id] = _gradient_info_create(name, type);
  _gradient_n_systems += 1;

  return _gradient_systems[mid_id];
}

/*----------------------------------------------------------------------------
 * Compute L2 norm.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- array of 3-vectors
 *----------------------------------------------------------------------------*/

static double
_l2_norm_1(cs_lnum_t            n_elts,
           cs_real_t  *restrict x)
{
  double s = cs_dot(n_elts, x, x);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double _s;
    MPI_Allreduce(&s, &_s, 1, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    s = _s;
  }

#endif /* defined(HAVE_MPI) */

  return (sqrt(s));
}

/*----------------------------------------------------------------------------
 * Synchronize strided gradient ghost cell values.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   halo_type      <-- halo type (extended or not)
 *   grad           --> gradient of a variable
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sync_strided_gradient_halo(const cs_mesh_t         *m,
                            cs_halo_type_t           halo_type,
                            cs_real_t (*restrict grad)[stride][3])
{
  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, stride*3);
    if (m->have_rotation_perio) {
      if (stride == 1)
        cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)grad, 3);
      else if (stride == 3)
        cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
      else if (stride == 6)
        cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                             halo_type,
                                             (cs_real_t *)grad);
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute face-based gradient clipping factor based on per-cell factor.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   ma             <-- mesh adjacencies
 *   halo_type      <-- halo type (extended or not)
 *   factor         <-- per cell input factor
 *   clip_factor    <-> final clip factor (initialized to 1)
 *----------------------------------------------------------------------------*/

static void
_gradient_update_face_clip_factor(const cs_mesh_t              *m,
                                  const cs_mesh_adjacencies_t  *ma,
                                  cs_halo_type_t                halo_type,
                                  cs_real_t           *restrict factor,
                                  cs_real_t           *restrict clip_factor)
{
  const cs_lnum_t n_cells = m->n_cells;
  const int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

  const size_t block_size = 128;
  const size_t n_blocks = cs_parall_block_count(n_cells, block_size);

  if (m->halo != nullptr) {
    cs_halo_sync_var(m->halo, halo_type, factor);
  }

# pragma omp parallel for  if  (n_cells > CS_THR_MIN)
  for (size_t b_id = 0; b_id < n_blocks; b_id++) {

    cs_lnum_t s_id = b_id*block_size, e_id = (b_id+1)*block_size;
    if (e_id > n_cells) e_id = n_cells;

    /* Remark:
       denum: maximum l2 norm of the variation of the gradient squared
       denom: maximum l2 norm of the variation of the variable squared */

    for (int adj_id = 0; adj_id < n_adj; adj_id++) {

      const cs_lnum_t *restrict cell_cells_idx;
      const cs_lnum_t *restrict cell_cells;

      if (adj_id == 0) {
        cell_cells_idx = ma->cell_cells_idx;
        cell_cells = ma->cell_cells;
      }
      else if (ma->cell_cells_e_idx != nullptr) {
        cell_cells_idx = ma->cell_cells_e_idx;
        cell_cells = ma->cell_cells_e;
      }
      else
        break;

      for (cs_lnum_t c_id1 = s_id; c_id1 < e_id; c_id1++) {

        cs_real_t l_min_factor = factor[c_id1];

        for (cs_lnum_t cidx = cell_cells_idx[c_id1];
             cidx < cell_cells_idx[c_id1+1];
             cidx++) {

          cs_lnum_t c_id2 = cell_cells[cidx];

          l_min_factor = std::min(l_min_factor, factor[c_id2]);

        }

        clip_factor[c_id1] = std::min(clip_factor[c_id1], l_min_factor);

      }  /* End of loop on block elements */

    } /* End of loop on adjacency type */

  } /* End of (parallel) loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return pointer to clipping factor field values if present
 *
 * The field, if present, is expected to be called
 * "algo:grad_clip_factor_<variable_name>".
 *
 * \param[in]  var_name  variable name
 *
 * \return  point to iteration count values, or nullptr
 */
/*----------------------------------------------------------------------------*/

static cs_real_t *
_get_clip_factor_try(const char  *var_name)
{
  cs_real_t *c_iter = nullptr;

  if (var_name != nullptr) {
    cs_field_t *f_c_iter
      = cs_field_by_composite_name_try("algo:grad_clip_factor",
                                       var_name);
    if (f_c_iter != nullptr)
      c_iter = f_c_iter->val;
  }

  return c_iter;
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a scalar if necessary.
 * This function deals with the standard or extended neighborhood.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   ma             <-- mesh adjacencies
 *   halo_type      <-- halo type (extended or not)
 *   clip_mode      <-- type of clipping for the computation of the gradient
 *   verbosity      <-- output level
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   pvar           <-- variable
 *   grad           <-> gradient of pvar (du/dx_j : grad[][j])
 *----------------------------------------------------------------------------*/

static void
_scalar_gradient_clipping(const cs_mesh_t              *m,
                          const cs_mesh_quantities_t   *fvq,
                          const cs_mesh_adjacencies_t  *ma,
                          cs_halo_type_t                halo_type,
                          int                           clip_mode,
                          int                           verbosity,
                          cs_real_t                     climgp,
                          const char                   *var_name,
                          const cs_real_t              *restrict  pvar,
                          cs_real_t                   (*restrict grad)[3])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;

  const cs_halo_t *halo = m->halo;

  if (clip_mode < 0 || climgp < 0)
    return;

  cs_dispatch_context ctx;

  /* The gradient and the variable must be already synchronized */

  cs_real_t *restrict clip_factor = _get_clip_factor_try(var_name);
  cs_real_t *_clip_factor = nullptr;

  if (clip_factor == nullptr) {
    CS_MALLOC(_clip_factor, n_cells_ext, cs_real_t);
    clip_factor = _clip_factor;
  }

  const int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

  const size_t block_size = 128;
  const size_t n_blocks = cs_parall_block_count(n_cells, block_size);

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clip_mode == CS_GRADIENT_LIMIT_CELL) {

#   pragma omp parallel for  if  (n_cells > CS_THR_MIN)
    for (size_t b_id = 0; b_id < n_blocks; b_id++) {

      cs_lnum_t s_id = b_id*block_size, e_id = (b_id+1)*block_size;
      if (e_id > n_cells) e_id = n_cells;

      /* Remark:
         denum: maximum l2 norm of the variation of the gradient squared
         denom: maximum l2 norm of the variation of the variable squared */

      cs_real_t denum[block_size];
      cs_real_t denom[block_size];

      cs_lnum_t b_e_id = e_id - s_id;

      for (cs_lnum_t i = 0; i < b_e_id; i++) {
        denum[i] = 0;
        denom[i] = 0;
      }

      for (int adj_id = 0; adj_id < n_adj; adj_id++) {

        const cs_lnum_t *restrict cell_cells_idx;
        const cs_lnum_t *restrict cell_cells;

        if (adj_id == 0) {
          cell_cells_idx = ma->cell_cells_idx;
          cell_cells = ma->cell_cells;
        }
        else if (ma->cell_cells_e_idx != nullptr) {
          cell_cells_idx = ma->cell_cells_e_idx;
          cell_cells = ma->cell_cells_e;
        }
        else
          break;

        for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */

          cs_lnum_t c_id1 = i + s_id;

          for (cs_lnum_t cidx = cell_cells_idx[c_id1];
               cidx < cell_cells_idx[c_id1+1];
               cidx++) {

            cs_lnum_t c_id2 = cell_cells[cidx];

            cs_real_t dist[3];

            for (cs_lnum_t k = 0; k < 3; k++)
              dist[k] = cell_cen[c_id1][k] - cell_cen[c_id2][k];

            cs_real_t dist1 = std::abs(  dist[0]*grad[c_id1][0]
                                       + dist[1]*grad[c_id1][1]
                                       + dist[2]*grad[c_id1][2]);
            cs_real_t dvar = std::abs(pvar[c_id1] - pvar[c_id2]);

            denum[i] = std::max(denum[i], dist1);
            denom[i] = std::max(denom[i], dvar);

          }

        }  /* End of loop on block elements */

      } /* End of loop on adjacency type */

      for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */
        cs_real_t factor1 = 1.;
        if (denum[i] > climgp * denom[i])
          factor1 = climgp * denom[i]/denum[i];

        clip_factor[s_id + i] = factor1;
      }

    } /* End of (parallel) loop on blocks */

  } /* End for clip_mode == CS_GRADIENT_LIMIT_CELL */

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clip_mode == CS_GRADIENT_LIMIT_FACE) {

    cs_real_t *factor;
    CS_MALLOC(factor, n_cells_ext, cs_real_t);
    cs_array_real_set_scalar(n_cells_ext, DBL_MAX, factor);

#   pragma omp parallel for  if  (n_cells > CS_THR_MIN)
    for (size_t b_id = 0; b_id < n_blocks; b_id++) {

      cs_lnum_t s_id = b_id*block_size, e_id = (b_id+1)*block_size;
      if (e_id > n_cells) e_id = n_cells;

      /* Remark:
         denum: maximum l2 norm of the variation of the gradient squared
         denom: maximum l2 norm of the variation of the variable squared */

      cs_real_t denum[block_size];
      cs_real_t denom[block_size];

      cs_lnum_t b_e_id = e_id - s_id;

      for (cs_lnum_t i = 0; i < b_e_id; i++) {
        denum[i] = 0;
        denom[i] = 0;
        clip_factor[s_id + i] = 1;
      }

      for (int adj_id = 0; adj_id < n_adj; adj_id++) {

        const cs_lnum_t *restrict cell_cells_idx;
        const cs_lnum_t *restrict cell_cells;

        if (adj_id == 0) {
          cell_cells_idx = ma->cell_cells_idx;
          cell_cells = ma->cell_cells;
        }
        else if (ma->cell_cells_e_idx != nullptr) {
          cell_cells_idx = ma->cell_cells_e_idx;
          cell_cells = ma->cell_cells_e;
        }
        else
          break;

        for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */

          cs_lnum_t c_id1 = i + s_id;

          for (cs_lnum_t cidx = cell_cells_idx[c_id1];
               cidx < cell_cells_idx[c_id1+1];
               cidx++) {

            cs_lnum_t c_id2 = cell_cells[cidx];

            cs_real_t dist[3], gradm[3];

            for (cs_lnum_t k = 0; k < 3; k++)
              dist[k] = cell_cen[c_id1][k] - cell_cen[c_id2][k];

            for (cs_lnum_t k = 0; k < 3; k++)
              gradm[k] = grad[c_id1][k] + grad[c_id2][k];

            cs_real_t dist1 = 0.5 * std::abs(cs_math_3_dot_product(dist, gradm));
            cs_real_t dvar = std::abs(pvar[c_id1] - pvar[c_id2]);

            denum[i] = std::max(denum[i], dist1);
            denom[i] = std::max(denom[i], dvar);

          }

        }  /* End of loop on block elements */

      } /* End of loop on adjacency type */

      for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */
        cs_real_t factor1 = 1.;
        if (denum[i] > climgp * denom[i])
          factor1 = climgp * denom[i]/denum[i];

        factor[s_id + i] = factor1;
      }

    } /* End of (parallel) loop on blocks */

    /* Now compute clip factor (kernel common to scalar and strided clases */

    _gradient_update_face_clip_factor(m, ma, halo_type, factor, clip_factor);

    CS_FREE(factor);

  } /* End for clip_mode == CS_GRADIENT_LIMIT_FACE */

  /* Synchronize variable */

  if (halo != nullptr) {
    cs_halo_sync_var(m->halo, halo_type, clip_factor);
  }

  /* Apply clip factor to gradient
     ----------------------------- */

  if (verbosity > 1) {

    cs_gnum_t n_clip = 0;
    cs_real_t min_factor = 0, max_factor = 0, mean_factor = 0;

    cs_array_reduce_simple_stats_l(ctx, n_cells, 1, nullptr, clip_factor,
                                   &min_factor,
                                   &max_factor,
                                   &mean_factor);

#   pragma omp parallel for reduction(+:n_clip) if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t factor1 = clip_factor[c_id];
      if (factor1 < 1) {
        for (cs_lnum_t j = 0; j < 3; j++)
          grad[c_id][j] *= factor1;
        n_clip += 1;
      }
    }

    cs_real_t buf[2] = {-min_factor, max_factor};
    cs_parall_max(2, CS_REAL_TYPE, buf);
    min_factor = -buf[0];
    max_factor =  buf[1];

    cs_parall_sum(1, CS_REAL_TYPE, &mean_factor);
    mean_factor /= (cs_real_t)(m->n_g_cells);

    cs_parall_counter(&n_clip, 1);

    bft_printf
      (_(" Variable: %s; gradient limitation in %llu cells\n"
         "   minimum factor = %g; maximum factor = %g; mean factor = %g\n"),
       var_name,
       (unsigned long long)n_clip, min_factor, max_factor, mean_factor);
  }

  else { /* Avoid unneeded reduction if no logging */

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t factor1 = clip_factor[c_id];
      if (factor1 < 1) {
        for (cs_lnum_t j = 0; j < 3; j++)
          grad[c_id][j] *= factor1;
      }

    } /* End of loop on cells */

  }

  /* Synchronize the updated Gradient */

  cs_halo_sync_r(m->halo, halo_type, false, grad);

  CS_FREE(_clip_factor);
}

/*----------------------------------------------------------------------------
 * Initialize gradient and right-hand side for scalar gradient reconstruction.
 *
 * A non-reconstructed gradient is computed at this stage.
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   w_stride       <-- stride for weighting coefficient
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_initialize_scalar_gradient(const cs_mesh_t                *m,
                            const cs_mesh_quantities_t     *fvq,
                            int                             w_stride,
                            int                             hyd_p_flag,
                            cs_real_t                       inc,
                            const cs_real_3_t               f_ext[],
                            const cs_field_bc_coeffs_t     *bc_coeffs,
                            const cs_real_t                 pvar[],
                            const cs_real_t                 c_weight[],
                            cs_real_3_t           *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  const cs_real_t *c_weight_s = nullptr;
  const cs_real_6_t *c_weight_t = nullptr;

  if (c_weight != nullptr) {
    if (w_stride == 1)
      c_weight_s = c_weight;
    else if (w_stride == 6)
      c_weight_t = (const cs_real_6_t *)c_weight;
  }

  /* Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  cs_lnum_t is_porous = 0;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  /* Initialize gradient */
  /*---------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      grad[cell_id][j] = 0.0;
  }

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  if (hyd_p_flag == 1) {

    /* Contribution from interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t ktpond = weight[f_id]; /* no cell weighting */
          /* if cell weighting is active */
          if (c_weight_s != nullptr) {
            ktpond = weight[f_id] * c_weight_s[ii]
                      / (       weight[f_id] * c_weight_s[ii]
                         + (1.0-weight[f_id])* c_weight_s[jj]);
          }
          else if (c_weight_t != nullptr) {
            cs_real_t sum[6];
            cs_real_t inv_sum[6];

            for (cs_lnum_t kk = 0; kk < 6; kk++)
              sum[kk] = weight[f_id]*c_weight_t[ii][kk]
                 +(1.0-weight[f_id])*c_weight_t[jj][kk];

            cs_math_sym_33_inv_cramer(sum, inv_sum);

            ktpond = weight[f_id] / 3.0
                     * (  inv_sum[0]*c_weight_t[ii][0]
                        + inv_sum[1]*c_weight_t[ii][1]
                        + inv_sum[2]*c_weight_t[ii][2]
                        + 2.0 * (  inv_sum[3]*c_weight_t[ii][3]
                                 + inv_sum[4]*c_weight_t[ii][4]
                                 + inv_sum[5]*c_weight_t[ii][5]));
          }

          cs_real_2_t poro = {
            i_poro_duq_0[is_porous*f_id],
            i_poro_duq_1[is_porous*f_id]
          };

          /*
             Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                      + (1-\alpha_\ij) \varia_\cellj\f$
                     but for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                     and for the cell \f$ \cellj \f$ we remove
                     \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
          */

          /* Reconstruction part */
          cs_real_t pfaci
            =  ktpond
                 * (  (i_face_cog[f_id][0] - cell_cen[ii][0])*f_ext[ii][0]
                    + (i_face_cog[f_id][1] - cell_cen[ii][1])*f_ext[ii][1]
                    + (i_face_cog[f_id][2] - cell_cen[ii][2])*f_ext[ii][2]
                    + poro[0])
            +  (1.0 - ktpond)
                 * (  (i_face_cog[f_id][0] - cell_cen[jj][0])*f_ext[jj][0]
                    + (i_face_cog[f_id][1] - cell_cen[jj][1])*f_ext[jj][1]
                    + (i_face_cog[f_id][2] - cell_cen[jj][2])*f_ext[jj][2]
                    + poro[1]);

          cs_real_t pfacj = pfaci;

          pfaci += (1.0-ktpond) * (pvar[jj] - pvar[ii]);
          pfacj -=      ktpond  * (pvar[jj] - pvar[ii]);

          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[f_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[f_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t f_id = b_group_index[t_id*2];
           f_id < b_group_index[t_id*2 + 1];
           f_id++) {

        cs_lnum_t ii = b_face_cells[f_id];

        cs_real_t poro = b_poro_duq[is_porous*f_id];

        /*
          Remark: for the cell \f$ \celli \f$ we remove
          \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
        */

        /* Reconstruction part */
        cs_real_t pfac
          = coefap[f_id] * inc
            + coefbp[f_id]
              * (  (b_face_cog[f_id][0] - cell_cen[ii][0])*f_ext[ii][0]
                 + (b_face_cog[f_id][1] - cell_cen[ii][1])*f_ext[ii][1]
                 + (b_face_cog[f_id][2] - cell_cen[ii][2])*f_ext[ii][2]
                 + poro);

        pfac += (coefbp[f_id] - 1.0) * pvar[ii];

        for (cs_lnum_t j = 0; j < 3; j++)
          grad[ii][j] += pfac * b_f_face_normal[f_id][j];

      } /* loop on faces */

    } /* loop on threads */

  } /* End of test on hydrostatic pressure */

  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  else {

    /* Contribution from interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t ktpond = weight[f_id]; /* no cell weighting */
          /* if cell weighting is active */
          if (w_stride == 1 && c_weight != nullptr) {
            ktpond = weight[f_id] * c_weight_s[ii]
                      / (      weight[f_id] * c_weight_s[ii]
                         + (1.0-weight[f_id])* c_weight_s[jj]);
          }
          else if (w_stride == 6 && c_weight != nullptr) {
            cs_real_t sum[6];
            cs_real_t inv_sum[6];

            for (cs_lnum_t kk = 0; kk < 6; kk++)
              sum[kk] = weight[f_id]*c_weight_t[ii][kk]
                 +(1.0-weight[f_id])*c_weight_t[jj][kk];

            cs_math_sym_33_inv_cramer(sum, inv_sum);

            ktpond =   weight[f_id] / 3.0
                     * (  inv_sum[0]*c_weight_t[ii][0]
                        + inv_sum[1]*c_weight_t[ii][1]
                        + inv_sum[2]*c_weight_t[ii][2]
                        + 2.0*(  inv_sum[3]*c_weight_t[ii][3]
                               + inv_sum[4]*c_weight_t[ii][4]
                               + inv_sum[5]*c_weight_t[ii][5]));
          }

          /*
             Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                      + (1-\alpha_\ij) \varia_\cellj\f$
                     but for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                     and for the cell \f$ \cellj \f$ we remove
                     \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
          */
          cs_real_t pfaci = (1.0-ktpond) * (pvar[jj] - pvar[ii]);
          cs_real_t pfacj =     -ktpond  * (pvar[jj] - pvar[ii]);

          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[f_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[f_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t f_id = b_group_index[t_id*2];
           f_id < b_group_index[t_id*2 + 1];
           f_id++) {

        cs_lnum_t ii = b_face_cells[f_id];

        /*
          Remark: for the cell \f$ \celli \f$ we remove
                  \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
        */

        cs_real_t pfac =   inc*coefap[f_id]
                         + (coefbp[f_id]-1.0)*pvar[ii];

        for (cs_lnum_t j = 0; j < 3; j++)
          grad[ii][j] += pfac * b_f_face_normal[f_id][j];

      } /* loop on faces */

    } /* loop on threads */

  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_vol[cell_id];
    else
      dvol = 0.;

    for (cs_lnum_t j = 0; j < 3; j++)
      grad[cell_id][j] *= dvol;
  }

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_EXTENDED, false, grad);
}

/*----------------------------------------------------------------------------
 * Renomalize scalar gradient by multiplying with a renormalization matrix
 * (so this gradient is not conservative on non Cartesian grids)
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_renormalize_scalar_gradient(const cs_mesh_t                *m,
                             const cs_mesh_quantities_t     *fvq,
                             int                             hyd_p_flag,
                             cs_real_3_t           *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;
  const int *bc_type = cs_glob_bc_type;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict cell_f_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict cell_cen = mq_g->cell_cen;
  const cs_real_3_t *restrict cell_f_cen = fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = mq_g->b_face_cog;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mq_g->b_face_normal;

  /* Correction matrix */
  cs_real_33_t *cor_mat;
  CS_MALLOC(cor_mat, n_cells_ext, cs_real_33_t);

  /* InitializatioC */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cor_mat[c_id][i][j] = 0.;
      }
    }
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           f_id++) {
        cs_lnum_t ii = i_face_cells[f_id][0];
        cs_lnum_t jj = i_face_cells[f_id][1];
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            cor_mat[ii][i][j] +=   (  i_face_cog[f_id][i]
                                    - cell_f_cen[ii][i])
                                 * i_f_face_normal[f_id][j];
            cor_mat[jj][i][j] -=   (  i_face_cog[f_id][i]
                                    - cell_f_cen[jj][i])
                                 * i_f_face_normal[f_id][j];
          }
        }
      }
    }
  }

  /* Contribution from boundary faces */
  for (int t_id = 0; t_id < n_b_threads; t_id++) {
    for (cs_lnum_t face_id = b_group_index[t_id*2];
        face_id < b_group_index[t_id*2 + 1];
        face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      if (bc_type[face_id] == CS_OUTLET || bc_type[face_id] == CS_SYMMETRY) {
        cs_real_t xip[3], xij[3];
        for (int k = 0; k < 3; k++) {
          xip[k] = cell_f_cen[ii][k];
          xij[k] = b_face_cog[face_id][k];
        }

        /* Special treatment for outlets */
        cs_real_t xi[3], nn[3];
        for (int k = 0; k < 3; k++) {
          xi[k] = cell_cen[ii][k];
          nn[k] = b_face_normal[face_id][k];
        }
        cs_real_t psca1 = cs_math_3_distance_dot_product(xip, xij, nn);
        cs_real_t psca2 = cs_math_3_distance_dot_product(xi,  xij, nn);
        cs_real_t lambda = psca1 / psca2;
        for (int k = 0; k < 3; k++) {
          xij[k] = xip[k] + lambda * (xij[k]-xi[k]);
        }
        /* End special treatment for outlets */

        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            cor_mat[ii][k][l] += (xij[k] - xip[k]) * b_f_face_normal[face_id][l];
      }
    }
  }

  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)cor_mat, 9);
    if (cs_glob_mesh->have_rotation_perio)
      cs_halo_perio_sync_var_tens(m->halo, CS_HALO_EXTENDED,
                                  (cs_real_t *)cor_mat);
  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
        cor_mat[cell_id][i][j] *= dvol;

    cs_math_33_inv_cramer_in_place(cor_mat[cell_id]);

    cs_real_t grd[3];
    for (cs_lnum_t j = 0; j < 3; j++)
      grd[j] = grad[cell_id][j];

    if (hyd_p_flag == 1) {
      cs_math_33_3_product(cor_mat[cell_id], grd, grad[cell_id]);
    }
  }

  CS_FREE(cor_mat);

  /* Synchronize halos */
  cs_halo_sync_r(m->halo, CS_HALO_EXTENDED, false, grad);
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the iterative algorithm
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <--  mesh quantities
 *   gq   <->  gradient quantities
 *
 * returns:
 *   pointer to cocg matrix (handled in main or coupled mesh quantities)
 *----------------------------------------------------------------------------*/

static cs_real_33_t *
_compute_cell_cocg_it(const cs_mesh_t               *m,
                      const cs_mesh_quantities_t    *fvq,
                      cs_gradient_quantities_t      *gq)
{
  /* Local variables */

  const int n_cells = m->n_cells;
  const int n_cells_with_ghosts = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;

  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  cs_real_33_t *restrict cocg = gq->cocg_it;

  cs_lnum_t  cell_id, f_id, i, j;
  cs_real_t  pfac, vecfac;
  cs_real_t  dvol1, dvol2;

  if (cocg == nullptr) {
    CS_MALLOC_HD(cocg, n_cells_with_ghosts, cs_real_33_t, cs_alloc_mode);
    gq->cocg_it = cocg;
  }

  /* compute the dimensionless matrix COCG for each cell*/

  for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    cocg[cell_id][0][0]= 1.0;
    cocg[cell_id][0][1]= 0.0;
    cocg[cell_id][0][2]= 0.0;
    cocg[cell_id][1][0]= 0.0;
    cocg[cell_id][1][1]= 1.0;
    cocg[cell_id][1][2]= 0.0;
    cocg[cell_id][2][0]= 0.0;
    cocg[cell_id][2][1]= 0.0;
    cocg[cell_id][2][2]= 1.0;
  }

  /* Interior face treatment */

  for (f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t cell_id1 = i_face_cells[f_id][0];
    cs_lnum_t cell_id2 = i_face_cells[f_id][1];

    if (cell_vol[cell_id1] > 0)
      dvol1 = 1. / cell_vol[cell_id1];
    else
      dvol1 = 0.;
    if (cell_vol[cell_id2] > 0)
      dvol2 = 1. / cell_vol[cell_id2];
    else
      dvol2 = 0.;

    for (i = 0; i < 3; i++) {

      pfac = -0.5*dofij[f_id][i];

      for (j = 0; j < 3; j++) {
        vecfac = pfac*i_face_normal[f_id][j];
        cocg[cell_id1][i][j] += vecfac * dvol1;
        cocg[cell_id2][i][j] -= vecfac * dvol2;
      }
    }
  }

  /* 3x3 Matrix inversion */

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells; cell_id++)
    cs_math_33_inv_cramer_in_place(cocg[cell_id]);

  return cocg;
}

/*----------------------------------------------------------------------------
 * Compute cell gradient using iterative reconstruction for non-orthogonal
 * meshes (nswrgp > 1).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * parameters:
 *   m               <-- pointer to associated mesh structure
 *   fvq             <-- pointer to associated finite volume quantities
 *   w_stride        <-- stride for weighting coefficient
 *   var_name        <-- variable name
 *   gradient_info   <-- pointer to performance logging structure, or nullptr
 *   nswrgp          <-- number of sweeps for gradient reconstruction
 *   hyd_p_flag      <-- flag for hydrostatic pressure
 *   verbosity       <-- verbosity level
 *   inc             <-- if 0, solve on increment; 1 otherwise
 *   epsrgp          <-- relative precision for gradient reconstruction
 *   f_ext           <-- exterior force generating pressure
 *   bc_coeffs       <-- B.C. structure for boundary face normals
 *   pvar            <-- variable
 *   c_weight        <-- weighted gradient coefficient variable
 *   grad            <-> gradient of pvar (halo prepared for periodicity
 *                       of rotation)
 *----------------------------------------------------------------------------*/

static void
_iterative_scalar_gradient(const cs_mesh_t                *m,
                           const cs_mesh_quantities_t     *fvq,
                           int                             w_stride,
                           const char                     *var_name,
                           cs_gradient_info_t             *gradient_info,
                           int                             nswrgp,
                           int                             hyd_p_flag,
                           int                             verbosity,
                           cs_real_t                       inc,
                           cs_real_t                       epsrgp,
                           const cs_real_3_t               f_ext[],
                           const cs_field_bc_coeffs_t     *bc_coeffs,
                           const cs_real_t                 pvar[],
                           const cs_real_t                 c_weight[],
                           cs_real_3_t           *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  cs_real_3_t *rhs;

  int n_sweeps = 0;
  cs_real_t l2_residual = 0.;

  if (nswrgp < 1) {
    if (gradient_info != nullptr)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  cs_gradient_quantities_t  *gq = _gradient_quantities_get(0);

  cs_real_33_t *restrict cocg = gq->cocg_it;
  if (cocg == nullptr)
    cocg = _compute_cell_cocg_it(m, fvq, gq);

  const cs_real_t *c_weight_s = nullptr;
  const cs_real_6_t *c_weight_t = nullptr;

  if (c_weight != nullptr) {
    if (w_stride == 1)
      c_weight_s = c_weight;
    else if (w_stride == 6)
      c_weight_t = (const cs_real_6_t *)c_weight;
  }

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  int is_porous = 0;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  /* Reconstruct gradients for non-orthogonal meshes */
  /*-------------------------------------------------*/

  /* Semi-implicit resolution on the whole mesh  */

  /* Compute normalization residual */

  cs_real_t  rnorm = _l2_norm_1(3*n_cells, (cs_real_t *)grad);

  if (rnorm <= cs_math_epzero) {
    if (gradient_info != nullptr)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  CS_MALLOC(rhs, n_cells_ext, cs_real_3_t);

  /* Vector OijFij is computed in CLDijP */

  /* Start iterations */
  /*------------------*/

  for (n_sweeps = 1; n_sweeps < nswrgp; n_sweeps++) {

    /* Case with hydrostatic pressure */
    /*--------------------------------*/

    if (hyd_p_flag == 1) {

    /* Compute right hand side */

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        rhs[c_id][0] = -(grad[c_id][0] - f_ext[c_id][0]) * cell_vol[c_id];
        rhs[c_id][1] = -(grad[c_id][1] - f_ext[c_id][1]) * cell_vol[c_id];
        rhs[c_id][2] = -(grad[c_id][2] - f_ext[c_id][2]) * cell_vol[c_id];
      }

      /* Contribution from interior faces */

      for (int g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_i_threads; t_id++) {

          for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               f_id++) {

            cs_lnum_t c_id1 = i_face_cells[f_id][0];
            cs_lnum_t c_id2 = i_face_cells[f_id][1];

            cs_real_t ktpond = weight[f_id]; /* no cell weighting */
            /* if cell weighting is active */
            if (c_weight_s != nullptr) {
              ktpond = weight[f_id] * c_weight_s[c_id1]
                        / (       weight[f_id] * c_weight_s[c_id1]
                           + (1.0-weight[f_id])* c_weight_s[c_id2]);
            }
            else if (c_weight_t != nullptr) {
              cs_real_t sum[6];
              cs_real_t inv_sum[6];

              for (cs_lnum_t ii = 0; ii < 6; ii++)
                sum[ii] =        weight[f_id] *c_weight_t[c_id1][ii]
                          + (1.0-weight[f_id])*c_weight_t[c_id2][ii];

              cs_math_sym_33_inv_cramer(sum, inv_sum);

              ktpond =   weight[f_id] / 3.0
                       * (  inv_sum[0]*c_weight_t[c_id1][0]
                          + inv_sum[1]*c_weight_t[c_id1][1]
                          + inv_sum[2]*c_weight_t[c_id1][2]
                          + 2.0 * (  inv_sum[3]*c_weight_t[c_id1][3]
                                   + inv_sum[4]*c_weight_t[c_id1][4]
                                   + inv_sum[5]*c_weight_t[c_id1][5]));
            }

            cs_real_2_t poro = {
              i_poro_duq_0[is_porous*f_id],
              i_poro_duq_1[is_porous*f_id]
            };

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                but for the cell \f$ \celli \f$ we remove
                \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                and for the cell \f$ \cellj \f$ we remove
                \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
                We also use the property
                \f$ \sum_\face \centf \otimes \vect{S}_\face = \vol \tens{I} \f$
                to remove the contributions
                \f$ ( \centi - \centf ) \cdot \vect{f}_\celli \f$
                \f$ ( \centj - \centf ) \cdot \vect{f}_\cellj \f$
            */

            /* Reconstruction part */
            cs_real_t dpfaci = pvar[c_id1]
                   + cs_math_3_distance_dot_product(cell_cen[c_id1],
                                                    i_face_cog[f_id],
                                                    f_ext[c_id1]);

            cs_real_t dpfacj = pvar[c_id2]
                   + cs_math_3_distance_dot_product(cell_cen[c_id2],
                                                    i_face_cog[f_id],
                                                    f_ext[c_id2]);

            cs_real_t pfaci = ktpond*poro[0] + (1.0-ktpond)*poro[1]
                 + 0.5*( cs_math_3_distance_dot_product(f_ext[c_id1],
                                                        grad[c_id1],
                                                        dofij[f_id])
                       + cs_math_3_distance_dot_product(f_ext[c_id2],
                                                        grad[c_id2],
                                                        dofij[f_id]));

            cs_real_t pfacj = pfaci;

            pfaci += (1.0-ktpond) * (dpfacj - dpfaci);
            pfacj -= ktpond * (dpfacj - dpfaci);

            for (cs_lnum_t j = 0; j < 3; j++) {
              rhs[c_id1][j] += pfaci * i_f_face_normal[f_id][j];
              rhs[c_id2][j] -= pfacj * i_f_face_normal[f_id][j];
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from boundary faces */

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[t_id*2];
             f_id < b_group_index[t_id*2 + 1];
             f_id++) {

          cs_lnum_t c_id = b_face_cells[f_id];

          cs_real_t poro = b_poro_duq[is_porous*f_id];

          /*
            Remark: for the cell \f$ \celli \f$ we remove
                    \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
          */

          /* Reconstruction part */
          cs_real_t pfac
            =   coefap[f_id] * inc
              + coefbp[f_id]
                 * (  cs_math_3_distance_dot_product(f_ext[c_id],
                                                     grad[c_id],
                                                     diipb[f_id])
                    + poro);

          cs_real_t dpfac = pvar[c_id]
                  + cs_math_3_distance_dot_product(cell_cen[c_id],
                                                   b_face_cog[f_id],
                                                   f_ext[c_id]);

          pfac += (coefbp[f_id] - 1.0) * dpfac;

          rhs[c_id][0] += pfac * b_f_face_normal[f_id][0];
          rhs[c_id][1] += pfac * b_f_face_normal[f_id][1];
          rhs[c_id][2] += pfac * b_f_face_normal[f_id][2];

        } /* loop on faces */

      } /* loop on threads */

    } /* End of test on hydrostatic pressure */

    /* Standard case, without hydrostatic pressure */
    /*---------------------------------------------*/

    else {

      /* Compute right hand side */

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        rhs[c_id][0] = -grad[c_id][0] * cell_vol[c_id];
        rhs[c_id][1] = -grad[c_id][1] * cell_vol[c_id];
        rhs[c_id][2] = -grad[c_id][2] * cell_vol[c_id];
      }

      /* Contribution from interior faces */

      for (int g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_i_threads; t_id++) {

          for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               f_id++) {

            cs_lnum_t c_id1 = i_face_cells[f_id][0];
            cs_lnum_t c_id2 = i_face_cells[f_id][1];

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            /* Reconstruction part */
            cs_real_t pfaci
              = 0.5 * (  dofij[f_id][0] * (grad[c_id1][0]+grad[c_id2][0])
                       + dofij[f_id][1] * (grad[c_id1][1]+grad[c_id2][1])
                       + dofij[f_id][2] * (grad[c_id1][2]+grad[c_id2][2]));
            cs_real_t pfacj = pfaci;

            cs_real_t ktpond = weight[f_id]; /* no cell weighting */
            /* if cell weighting is active */
            if (c_weight_s != nullptr) {
              ktpond =   weight[f_id] * c_weight_s[c_id1]
                       / (       weight[f_id] * c_weight_s[c_id1]
                          + (1.0-weight[f_id])* c_weight_s[c_id2]);
            }
            else if (c_weight_t != nullptr) {
              cs_real_t sum[6];
              cs_real_t inv_sum[6];

              for (cs_lnum_t ii = 0; ii < 6; ii++)
                sum[ii] = weight[f_id]*c_weight_t[c_id1][ii]
                   +(1.0-weight[f_id])*c_weight_t[c_id2][ii];

              cs_math_sym_33_inv_cramer(sum, inv_sum);

              ktpond =   weight[f_id] / 3.0
                       * (  inv_sum[0]*c_weight_t[c_id1][0]
                          + inv_sum[1]*c_weight_t[c_id1][1]
                          + inv_sum[2]*c_weight_t[c_id1][2]
                          + 2.0 * (  inv_sum[3]*c_weight_t[c_id1][3]
                                   + inv_sum[4]*c_weight_t[c_id1][4]
                                   + inv_sum[5]*c_weight_t[c_id1][5]));
            }

            pfaci += (1.0-ktpond) * (pvar[c_id2] - pvar[c_id1]);
            pfacj -=      ktpond  * (pvar[c_id2] - pvar[c_id1]);

            for (cs_lnum_t j = 0; j < 3; j++) {
              rhs[c_id1][j] += pfaci * i_f_face_normal[f_id][j];
              rhs[c_id2][j] -= pfacj * i_f_face_normal[f_id][j];
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from boundary faces */

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[t_id*2];
             f_id < b_group_index[t_id*2 + 1];
             f_id++) {

          cs_lnum_t c_id = b_face_cells[f_id];

          /*
            Remark: for the cell \f$ \celli \f$ we remove
            \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
          */

          /* Reconstruction part */
          cs_real_t pfac
            =      coefap[f_id] * inc
                   + coefbp[f_id]
                     * (  diipb[f_id][0] * grad[c_id][0]
                        + diipb[f_id][1] * grad[c_id][1]
                        + diipb[f_id][2] * grad[c_id][2]);

          pfac += (coefbp[f_id] -1.0) * pvar[c_id];

          rhs[c_id][0] += pfac * b_f_face_normal[f_id][0];
          rhs[c_id][1] += pfac * b_f_face_normal[f_id][1];
          rhs[c_id][2] += pfac * b_f_face_normal[f_id][2];

        } /* loop on faces */

      } /* loop on threads */

    }

    /* Increment gradient */
    /*--------------------*/

#   pragma omp parallel for
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t dvol;
      /* Is the cell disabled (for solid or porous)? Not the case if coupled */
      if (has_dc * c_disable_flag[has_dc * c_id] == 0)
        dvol = 1. / cell_vol[c_id];
      else
        dvol = 0.;

      rhs[c_id][0] *= dvol;
      rhs[c_id][1] *= dvol;
      rhs[c_id][2] *= dvol;

      grad[c_id][0] +=   rhs[c_id][0] * cocg[c_id][0][0]
                       + rhs[c_id][1] * cocg[c_id][1][0]
                       + rhs[c_id][2] * cocg[c_id][2][0];
      grad[c_id][1] +=   rhs[c_id][0] * cocg[c_id][0][1]
                       + rhs[c_id][1] * cocg[c_id][1][1]
                       + rhs[c_id][2] * cocg[c_id][2][1];
      grad[c_id][2] +=   rhs[c_id][0] * cocg[c_id][0][2]
                       + rhs[c_id][1] * cocg[c_id][1][2]
                       + rhs[c_id][2] * cocg[c_id][2][2];
    }

    /* Synchronize halos */

    cs_halo_sync_r(m->halo, CS_HALO_STANDARD, false, grad);

    /* Convergence test */

    l2_residual = _l2_norm_1(3*n_cells, (cs_real_t *)rhs);

    if (l2_residual < epsrgp*rnorm) {
      if (verbosity >= 2)
        bft_printf(_(" %s; variable: %s; converged in %d sweeps\n"
                     " %*s  normed residual: %11.4e; norm: %11.4e\n"),
                   __func__, var_name, n_sweeps,
                   (int)(strlen(__func__)), " ", l2_residual/rnorm, rnorm);
      break;
    }

  } /* Loop on sweeps */

  if (l2_residual >= epsrgp*rnorm && verbosity > -1) {
    bft_printf(_(" Warning:\n"
                 " --------\n"
                 "   %s; variable: %s; sweeps: %d\n"
                 "   %*s  normed residual: %11.4e; norm: %11.4e\n"),
               __func__, var_name, n_sweeps,
               (int)(strlen(__func__)), " ", l2_residual/rnorm, rnorm);
  }

  if (gradient_info != nullptr)
    _gradient_info_update_iter(gradient_info, n_sweeps);

  CS_FREE(rhs);
}

/*----------------------------------------------------------------------------
 * Add compute 3x3 cocg for least squares algorithm contribution from hidden
 * faces (as pure homogeneous Neumann BC's) to a single cell
 *
 * parameters:
 *   c_id               <-- cell id
 *   cell_hb_faces_idx  <-- cells -> hidden boundary faces index
 *   cell_hb_faces      <-- cells -> hidden boundary faces adjacency
 *   b_face_u_normal    <-- boundary faces unit normals
 *   cocg               <-> cocg covariance matrix for given cell
 *----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
_add_hb_faces_cocg_lsq_cell(cs_lnum_t         c_id,
                            const cs_lnum_t   cell_hb_faces_idx[],
                            const cs_lnum_t   cell_hb_faces[],
                            const cs_nreal_t  b_face_u_normal[][3],
                            cs_cocg_t         cocg[6])

{
  cs_lnum_t s_id = cell_hb_faces_idx[c_id];
  cs_lnum_t e_id = cell_hb_faces_idx[c_id+1];

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    cs_lnum_t f_id = cell_hb_faces[i];
    const cs_nreal_t *dddij = b_face_u_normal[f_id];

    cocg[0] += dddij[0]*dddij[0];
    cocg[1] += dddij[1]*dddij[1];
    cocg[2] += dddij[2]*dddij[2];
    cocg[3] += dddij[0]*dddij[1];
    cocg[4] += dddij[1]*dddij[2];
    cocg[5] += dddij[0]*dddij[2];

  }
}

/*----------------------------------------------------------------------------
 * Add compute 3x3 cocg for least squares algorithm contribution from hidden
 * faces (as pure homogeneous Neumann BC's)
 *
 * parameters:
 *   m    <-- mesh
 *   fvq  <-- mesh quantities
 *   ma   <-- mesh adjacencies
 *   ctx  <-- Reference to dispatch context
 *   cocg <-> cocg covariance matrix
 *----------------------------------------------------------------------------*/

static void
_add_hb_faces_cell_cocg_lsq(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            const cs_mesh_adjacencies_t  *ma,
                            cs_dispatch_context          &ctx,
                            cs_cocg_6_t                  *cocg)
{
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;

  const cs_lnum_t  *restrict cell_hb_faces_idx = ma->cell_hb_faces_idx;
  const cs_lnum_t  *restrict cell_hb_faces = ma->cell_hb_faces;

  ctx.parallel_for(m->n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    _add_hb_faces_cocg_lsq_cell(c_id,
                                cell_hb_faces_idx,
                                cell_hb_faces,
                                b_face_u_normal,
                                cocg[c_id]);
  });
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient least squares algorithm
 *
 * parameters:
 *   m             <--  mesh
 *   extended      <--  true if extended neighborhood used
 *   fvq           <--  mesh quantities
 *   gq            <->  gradient quantities
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cocg_lsq(const cs_mesh_t               *m,
                       bool                           extended,
                       const cs_mesh_quantities_t    *fvq,
                       cs_gradient_quantities_t      *gq)
{
  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict c2b_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict c2b = ma->cell_b_faces;

  const cs_lnum_t *restrict b_cells = m->b_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;

  cs_cocg_6_t  *restrict cocgb = nullptr, *restrict cocg = nullptr;

   cs_dispatch_context ctx;

  /* Map cocg/cocgb to correct structure, reallocate if needed */

  if (extended) {
    cocg = gq->cocg_lsq_ext;
    cocgb = gq->cocgb_s_lsq_ext;
  }
  else {
    cocg = gq->cocg_lsq;
    cocgb =gq->cocgb_s_lsq;
  }

  if (cocg == nullptr) {

    assert(cocgb == nullptr);

    CS_MALLOC_HD(cocg, n_cells_ext, cs_cocg_6_t, cs_alloc_mode);
    CS_MALLOC_HD(cocgb, m->n_b_cells, cs_cocg_6_t, cs_alloc_mode);

    if (extended) {
      gq->cocg_lsq_ext = cocg;
      gq->cocgb_s_lsq_ext = cocgb;
    }
    else {
      gq->cocg_lsq = cocg;
      gq->cocgb_s_lsq = cocgb;
    }

  }

  /* Initialize cocg in ghost cells */

  cs_lnum_t n_halo_cells = n_cells_ext - n_cells;
  if (n_halo_cells > 0) {
    ctx.parallel_for(n_halo_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t idx) {
      auto _cocg = cocg[n_cells + idx];
      _cocg[0] = 0; _cocg[1] = 0; _cocg[2] = 0;
      _cocg[3] = 0; _cocg[4] = 0; _cocg[5] = 0;
    });
  }

  /* Add contributions from neighbor cells (standard and extended)
     ------------------------------------------------------------- */

  int n_steps = (extended) ? 2 : 1;

  for (int step_id = 0; step_id < n_steps; step_id++) {

    const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
    const cs_lnum_t *c2c = ma->cell_cells;
    if (step_id == 1) {
      c2c_idx = ma->cell_cells_e_idx;
      c2c = ma->cell_cells_e;
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      auto _cocg = cocg[c_id];

      if (step_id == 0) {
        _cocg[0] = 0; _cocg[1] = 0; _cocg[2] = 0;
        _cocg[3] = 0; _cocg[4] = 0; _cocg[5] = 0;
      }

      cs_lnum_t s_id = c2c_idx[c_id];
      cs_lnum_t e_id = c2c_idx[c_id + 1];

      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t c_id1 = c2c[i];

        cs_real_t dc[3] = {cell_cen[c_id1][0] - cell_cen[c_id][0],
                           cell_cen[c_id1][1] - cell_cen[c_id][1],
                           cell_cen[c_id1][2] - cell_cen[c_id][2]};

        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        _cocg[0] += dc[0] * dc[0] * ddc;
        _cocg[1] += dc[1] * dc[1] * ddc;
        _cocg[2] += dc[2] * dc[2] * ddc;
        _cocg[3] += dc[0] * dc[1] * ddc;
        _cocg[4] += dc[1] * dc[2] * ddc;
        _cocg[5] += dc[0] * dc[2] * ddc;
      }

    });

  } /* Loop on standard/extended neighborhood */

  /* Contribution from hidden boundary faces, if present */

  if (m->n_b_faces_all > m->n_b_faces)
    _add_hb_faces_cell_cocg_lsq(m, fvq, ma, ctx, cocg);

  /* Contribution from boundary faces
     --------------------------------

     Assume symmetry everywhere so as to avoid obtaining
     a non-invertible matrix in 2D cases. */

  ctx.parallel_for(m->n_b_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t idx) {

    cs_lnum_t c_id = b_cells[idx];

    auto _cocg = cocg[c_id];

    /* Save partial cocg at interior faces of boundary cells */
    for (cs_lnum_t ll = 0; ll < 6; ll++) {
      cocgb[idx][ll] = _cocg[ll];
    }

    cs_lnum_t s_id = c2b_idx[c_id];
    cs_lnum_t e_id = c2b_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) { /* loop on boundary faces */
      cs_lnum_t f_id = c2b[i];
      const cs_nreal_t *u_normal = b_face_u_normal[f_id];

      _cocg[0] += u_normal[0] * u_normal[0];
      _cocg[1] += u_normal[1] * u_normal[1];
      _cocg[2] += u_normal[2] * u_normal[2];
      _cocg[3] += u_normal[0] * u_normal[1];
      _cocg[4] += u_normal[1] * u_normal[2];
      _cocg[5] += u_normal[0] * u_normal[2];
    }

  });

  /* Invert for all cells
     -------------------- */

  /* The cocg term for interior cells only changes if the mesh does */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    _math_6_inv_cramer_sym_in_place(cocg[c_id]);
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Return current symmetric 3x3 matrix cocg for least squares algorithm
 *
 * parameters:
 *   m          <--  mesh
 *   halo_type  <--  halo type
 *   accel      <--  use accelerator device (if true, cocg and cocgb
 *                   pointers returned are device pointers)
 *   fvq        <--  mesh quantities
 *   cocg       -->  coupling coefficients (covariance matrices)
 *   cocgb      -->  partial boundary coupling coeffients, or nullptr
 *----------------------------------------------------------------------------*/

static void
_get_cell_cocg_lsq(const cs_mesh_t               *m,
                   cs_halo_type_t                 halo_type,
                   bool                           accel,
                   const cs_mesh_quantities_t    *fvq,
                   cs_cocg_6_t                   *restrict *cocg,
                   cs_cocg_6_t                   *restrict *cocgb)
{
  cs_gradient_quantities_t  *gq = _gradient_quantities_get(0);

  cs_cocg_6_t *_cocg = nullptr;

  bool extended = (   halo_type == CS_HALO_EXTENDED
                   && m->cell_cells_idx) ? true : false;

  if (extended) {
    _cocg = gq->cocg_lsq_ext;
  }
  else {
    _cocg = gq->cocg_lsq;
  }

  /* Compute if not present yet.
   *
   * TODO: when using accelerators, this implies a first computation will be
   *       run on the host. This will usually be amortized, but could be
   *       further improved. */

  if (_cocg == nullptr)
    _compute_cell_cocg_lsq(m, extended, fvq, gq);

  /* If used on accelerator, ensure arrays are available there */

#if defined(HAVE_ACCEL)

  if (accel) {

    cs_alloc_mode_t alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;

    if (extended) {
      cs_set_alloc_mode_r(gq->cocg_lsq_ext, alloc_mode);
      cs_set_alloc_mode_r(gq->cocgb_s_lsq_ext, alloc_mode);
    }
    else {
      cs_set_alloc_mode_r(gq->cocg_lsq, alloc_mode);
      cs_set_alloc_mode_r(gq->cocgb_s_lsq, alloc_mode);
    }

  }

#endif

  /* Set pointers */

  if (extended)
    *cocg = gq->cocg_lsq_ext;
  else
    *cocg = gq->cocg_lsq;

  if (cocgb != nullptr) {
    if (extended)
      *cocgb = gq->cocgb_s_lsq_ext;
    else
      *cocgb = gq->cocgb_s_lsq;
  }

  /* If used on accelerator, copy/prefetch values and switch to
     device pointers */

  if (accel) {
    cs_sync_h2d(*cocg);
    *cocg = (cs_cocg_6_t *)cs_get_device_ptr(*cocg);

    if (cocgb != nullptr) {
      cs_sync_h2d(*cocgb);
      *cocgb = (cs_cocg_6_t *)cs_get_device_ptr(*cocgb);
    }
  }
}

/*----------------------------------------------------------------------------
 * Recompute scalar least-squares cocg at boundaries, using saved cocgb.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   ctx            <-- Reference to dispatch context
 *   cocgb          <-- saved B.C. coefficients for boundary cells
 *   cocgb          <-> B.C. coefficients, updated at boundary cells
 *----------------------------------------------------------------------------*/

static void
_recompute_lsq_scalar_cocg(const cs_mesh_t                *m,
                           const cs_mesh_quantities_t     *fvq,
                           const cs_field_bc_coeffs_t     *bc_coeffs,
                           cs_dispatch_context            &ctx,
                           const cs_cocg_t               (*restrict cocgb)[6],
                           cs_cocg_t                     (*restrict cocg)[6])
{
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_b_faces = ma->cell_b_faces;
  const cs_lnum_t *restrict b_cells = m->b_cells;

  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /* Recompute cocg at boundaries, using saved cocgb */

  ctx.parallel_for(m->n_b_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    cs_lnum_t c_id = b_cells[ii];

    cs_cocg_t _cocg[6];
    for (cs_lnum_t ll = 0; ll < 6; ll++)
      _cocg[ll] = cocgb[ii][ll];

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) { /* loop on boundary faces */

      cs_lnum_t f_id = cell_b_faces[i];

      cs_real_t umcbdd = (1. - coefbp[f_id]) / b_dist[f_id];

      cs_real_t dddij[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dddij[ll] =   b_face_u_normal[f_id][ll]
                    + umcbdd * diipb[f_id][ll];

      _cocg[0] += dddij[0]*dddij[0];
      _cocg[1] += dddij[1]*dddij[1];
      _cocg[2] += dddij[2]*dddij[2];
      _cocg[3] += dddij[0]*dddij[1];
      _cocg[4] += dddij[1]*dddij[2];
      _cocg[5] += dddij[0]*dddij[2];

    } /* loop on boundary faces */

    _math_6_inv_cramer_sym(_cocg, cocg[c_id]);

  }); /* loop on boundary cells */
}

/*----------------------------------------------------------------------------
 * Compute cell gradient using least-squares reconstruction.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   recompute_cocg <-- flag to recompute cocg
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or nullptr
 *   grad           --> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_gradient(const cs_mesh_t                *m,
                     const cs_mesh_quantities_t     *fvq,
                     cs_halo_type_t                  halo_type,
                     bool                            recompute_cocg,
                     cs_real_t                       inc,
                     const cs_field_bc_coeffs_t     *bc_coeffs,
                     const cs_real_t                 pvar[],
                     const cs_real_t       *restrict c_weight,
                     cs_real_3_t           *restrict grad)
{
  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx = m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst = m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_t *restrict weight = fvq->weight;

  cs_cocg_6_t  *restrict cocgb = nullptr;
  cs_cocg_6_t  *restrict cocg = nullptr;

  /* Weighting requires face weights, so cannot be applied consistently
     with an extended neighborhhood. */
  if (c_weight != nullptr)
    halo_type = CS_HALO_STANDARD;

#if defined(HAVE_CUDA)
  bool accel = (cs_get_device_id() > -1) ? true : false;
#else
  bool accel = false;
#endif

  _get_cell_cocg_lsq(m, halo_type, accel, fvq, &cocg, &cocgb);

#if defined(HAVE_CUDA)

  if (accel) {
    cs_gradient_scalar_lsq_cuda(m,
                                fvq,
                                halo_type,
                                recompute_cocg,
                                inc,
                                bc_coeffs,
                                pvar,
                                c_weight,
                                cocg,
                                cocgb,
                                grad);

    return;
  }

#endif

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  /* Compute cocg and save contribution at boundaries */

  if (recompute_cocg) {
    cs_dispatch_context ctx;
    ctx.set_use_gpu(false);
    _recompute_lsq_scalar_cocg(m,
                               fvq,
                               bc_coeffs,
                               ctx,
                               cocgb,
                               cocg);
    ctx.wait();
  }

  /* Compute Right-Hand Side */
  /*-------------------------*/

  cs_real_4_t  *restrict rhsv;
  CS_MALLOC(rhsv, n_cells_ext, cs_real_4_t);

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    rhsv[c_id][0] = 0.0;
    rhsv[c_id][1] = 0.0;
    rhsv[c_id][2] = 0.0;
    rhsv[c_id][3] = pvar[c_id];
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           f_id++) {

        cs_lnum_t ii = i_face_cells[f_id][0];
        cs_lnum_t jj = i_face_cells[f_id][1];

        cs_real_t pond = weight[f_id];

        cs_real_t pfac, dc[3], fctb[4];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        if (c_weight != nullptr) {
          /* (P_j - P_i) / ||d||^2 */
          pfac =   (rhsv[jj][3] - rhsv[ii][3])
                 / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            fctb[ll] = dc[ll] * pfac;

          cs_real_t denom = 1. / (  pond       *c_weight[ii]
                                  + (1. - pond)*c_weight[jj]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ii][ll] +=  c_weight[jj] * denom * fctb[ll];

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[jj][ll] +=  c_weight[ii] * denom * fctb[ll];
        }
        else {
          /* (P_j - P_i) / ||d||^2 */
          pfac =   (rhsv[jj][3] - rhsv[ii][3])
                 / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            fctb[ll] = dc[ll] * pfac;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ii][ll] += fctb[ll];

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[jj][ll] += fctb[ll];
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from extended neighborhood */

  if (halo_type == CS_HALO_EXTENDED && cell_cells_idx != nullptr) {

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
      for (cs_lnum_t cidx = cell_cells_idx[ii];
           cidx < cell_cells_idx[ii+1];
           cidx++) {

        cs_lnum_t jj = cell_cells_lst[cidx];

        cs_real_t pfac, dc[3], fctb[4];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        pfac =   (rhsv[jj][3] - rhsv[ii][3])
               / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          fctb[ll] = dc[ll] * pfac;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ii][ll] += fctb[ll];

      }
    }

  } /* End for extended neighborhood */

  /* Contribution from boundary faces */

# pragma omp parallel for
  for (int t_id = 0; t_id < n_b_threads; t_id++) {

    for (cs_lnum_t f_id = b_group_index[t_id*2];
         f_id < b_group_index[t_id*2 + 1];
         f_id++) {

      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t unddij = 1. / b_dist[f_id];
      cs_real_t umcbdd = (1. - coefbp[f_id]) * unddij;

      cs_real_t dsij[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dsij[ll] =   b_face_u_normal[f_id][ll]
                   + umcbdd*diipb[f_id][ll];

      cs_real_t pfac =   (coefap[f_id]*inc + (coefbp[f_id] -1.)
                       * rhsv[ii][3]) * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ii][ll] += dsij[ll] * pfac;

    } /* loop on faces */

  } /* loop on threads */

  /* Compute gradient */
  /*------------------*/

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    grad[c_id][0] =   cocg[c_id][0] *rhsv[c_id][0]
                    + cocg[c_id][3] *rhsv[c_id][1]
                    + cocg[c_id][5] *rhsv[c_id][2];
    grad[c_id][1] =   cocg[c_id][3] *rhsv[c_id][0]
                    + cocg[c_id][1] *rhsv[c_id][1]
                    + cocg[c_id][4] *rhsv[c_id][2];
    grad[c_id][2] =   cocg[c_id][5] *rhsv[c_id][0]
                    + cocg[c_id][4] *rhsv[c_id][1]
                    + cocg[c_id][2] *rhsv[c_id][2];
  }

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_STANDARD, false, grad);

  CS_FREE(rhsv);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient by least-squares reconstruction with a volume force
 * generating a hydrostatic pressure component.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   recompute_cocg <-- flag to recompute cocg
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight_s     <-- weighted gradient coefficient variable,
 *                      or nullptr
 *   grad           --> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_gradient_hyd_p(const cs_mesh_t                *m,
                           const cs_mesh_quantities_t     *fvq,
                           cs_halo_type_t                  halo_type,
                           bool                            recompute_cocg,
                           cs_real_t                       inc,
                           const cs_real_3_t               f_ext[],
                           const cs_field_bc_coeffs_t     *bc_coeffs,
                           const cs_real_t                 pvar[],
                           const cs_real_t       *restrict c_weight_s,
                           cs_real_3_t           *restrict grad)
{
  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx = m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst = m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_t *restrict weight = fvq->weight;

  std::chrono::high_resolution_clock::time_point t_start, t_cocg, t_init, \
    t_i_faces, t_ext_cells, t_b_faces, t_grad, t_stop;

  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  /* Weighting requires face weights, so cannot be applied consistently
     with an extended neighborhhod.  */
  if (c_weight_s != nullptr)
    halo_type = CS_HALO_STANDARD;

  cs_dispatch_context ctx;
  cs_dispatch_context ctx_b;

  ctx_b.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_b.use_gpu())
    ctx_b.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  cs_cocg_6_t  *restrict cocgb = nullptr;
  cs_cocg_6_t  *restrict cocg = nullptr;

  /* Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0 = nullptr;
  cs_real_t *i_poro_duq_1 = nullptr;
  cs_real_t *b_poro_duq = nullptr;

  int is_porous = 0;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  }

  _get_cell_cocg_lsq(m,
                     halo_type,
                     ctx.use_gpu(),
                     fvq,
                     &cocg,
                     &cocgb);

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  /* Compute cocg and save contribution at boundaries */

  if (recompute_cocg)
    _recompute_lsq_scalar_cocg(m,
                               fvq,
                               bc_coeffs,
                               ctx,
                               cocgb,
                               cocg);

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();
    t_cocg = std::chrono::high_resolution_clock::now();
  }

  /* Compute Right-Hand Side */
  /*-------------------------*/

  cs_real_4_t  *restrict rhsv;
  CS_MALLOC_HD(rhsv, n_cells_ext, cs_real_4_t, cs_alloc_mode);

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    rhsv[c_id][0] = 0.0;
    rhsv[c_id][1] = 0.0;
    rhsv[c_id][2] = 0.0;
    rhsv[c_id][3] = pvar[c_id];
  });

  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_init = std::chrono::high_resolution_clock::now();

  /* Contribution from interior faces
     -------------------------------- */

  if (c_weight_s != nullptr) {  /* With cell weighting */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t dvarij, pfac, dc[3];

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

      /* In this case, use rhsv, as access patterns lead to better
         caching behavior (or did when last tested) */
      dvarij = rhsv[jj][3] - rhsv[ii][3];

      if (is_porous) {
        pfac =   dvarij
               + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[ii],
                                                f_ext[ii])
               + i_poro_duq_0[f_id]
               - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[jj],
                                                f_ext[jj])
               - i_poro_duq_1[f_id];
      }
      else {
        pfac =   dvarij
               + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[ii],
                                                f_ext[ii])
               - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[jj],
                                                f_ext[jj]);
      }
      pfac /= cs_math_3_square_norm(dc);

      cs_real_t pond = weight[f_id];

      pfac /= (  pond       *c_weight_s[ii]
               + (1. - pond)*c_weight_s[jj]);

      cs_real_t pfac_ii = pfac * c_weight_s[jj];
      cs_real_t pfac_jj = pfac * c_weight_s[ii];

      cs_real_t rhsv_ii[3], rhsv_jj[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++) {
        rhsv_ii[ll] = dc[ll] * pfac_ii;
        rhsv_jj[ll] = dc[ll] * pfac_jj;
      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhsv[ii], rhsv_ii, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhsv[jj], rhsv_jj, i_sum_type);

    });

  }
  else { /* Without cell weights */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t dvarij, pfac, dc[3];

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

      /* In this case, use rhsv, as access patterns lead to better
         caching behavior (or did when last tested) */
      dvarij = rhsv[jj][3] - rhsv[ii][3];

      if (is_porous) {
        pfac =   dvarij
               + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[ii],
                                                f_ext[ii])
               + i_poro_duq_0[f_id]
               - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[jj],
                                                f_ext[jj])
               - i_poro_duq_1[is_porous*f_id];
      }
      else {
        pfac =   dvarij
               + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[ii],
                                                f_ext[ii])
               - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                cell_cen[jj],
                                                f_ext[jj]);
      }
      pfac /= cs_math_3_square_norm(dc);

      cs_real_t fctb[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++) {
        fctb[ll] = dc[ll] * pfac;
      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhsv[ii], fctb, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhsv[jj], fctb, i_sum_type);

    }); /* End of loop on faces */

  }  /* End of test on weighting */

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();  // Using an event would be preferred here
    t_i_faces = std::chrono::high_resolution_clock::now();
  }

  /* Contribution from extended neighborhood;
     We assume that the middle of the segment joining cell centers
     may replace the center of gravity of a fictitious face. */

  if (halo_type == CS_HALO_EXTENDED && cell_cells_idx != nullptr) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  ii) {
      for (cs_lnum_t cidx = cell_cells_idx[ii];
           cidx < cell_cells_idx[ii+1];
           cidx++) {

        cs_lnum_t jj = cell_cells_lst[cidx];

        /* Note: replaced the expressions:
         *  a) ptmid = 0.5 * (cell_cen[jj] - cell_cen[ii])
         *  b)   (cell_cen[ii] - ptmid) * f_ext[ii]
         *  c) - (cell_cen[jj] - ptmid) * f_ext[jj]
         * with:
         *  a) dc = cell_cen[jj] - cell_cen[ii]
         *  b) - 0.5 * dc * f_ext[ii]
         *  c) - 0.5 * dc * f_ext[jj]
         */

        cs_real_t pfac, dc[3], fctb[4];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        pfac =   (  rhsv[jj][3] - rhsv[ii][3]
                  - 0.5 * cs_math_3_dot_product(dc, f_ext[ii])
                  - 0.5 * cs_math_3_dot_product(dc, f_ext[jj]))
                / cs_math_3_square_norm(dc);

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          fctb[ll] = dc[ll] * pfac;

        cs_dispatch_sum<3>(rhsv[ii], fctb, i_sum_type);

      }
    });

  } /* End for extended neighborhood */

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();  // Using an event would be preferred here
    t_ext_cells = std::chrono::high_resolution_clock::now();
  }

  /* Contribution from boundary faces */

  ctx_b.set_use_gpu(false);
  ctx_b.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = b_face_cells[f_id];

    cs_real_t poro = (is_porous) ? b_poro_duq[f_id] : 0.;

    cs_real_t unddij = 1. / b_dist[f_id];
    cs_real_t umcbdd = (1. - coefbp[f_id]) * unddij;

    cs_real_t dsij[3];
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      dsij[ll] =   b_face_u_normal[f_id][ll]
                 + umcbdd*diipb[f_id][ll];

    cs_real_t pfac
      =   (coefap[f_id]*inc
          + (  (coefbp[f_id] -1.)
             * (  rhsv[ii][3]
                 /* (b_face_cog - cell_cen).f_ext, or IF.F_i */
                + cs_math_3_distance_dot_product(cell_cen[ii],
                                                 b_face_cog[f_id],
                                                 f_ext[ii])
                + poro)))
          * unddij;

    cs_real_t fctb[3];
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      fctb[ll] = dsij[ll] * pfac;

    cs_dispatch_sum<3>(rhsv[ii], fctb, b_sum_type);

  }); /* loop on faces */

  ctx_b.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_b_faces = std::chrono::high_resolution_clock::now();

  /* Compute gradient */
  /*------------------*/

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
    grad[c_id][0] =   cocg[c_id][0] *rhsv[c_id][0]
                    + cocg[c_id][3] *rhsv[c_id][1]
                    + cocg[c_id][5] *rhsv[c_id][2]
                    + f_ext[c_id][0];
    grad[c_id][1] =   cocg[c_id][3] *rhsv[c_id][0]
                    + cocg[c_id][1] *rhsv[c_id][1]
                    + cocg[c_id][4] *rhsv[c_id][2]
                    + f_ext[c_id][1];
    grad[c_id][2] =   cocg[c_id][5] *rhsv[c_id][0]
                    + cocg[c_id][4] *rhsv[c_id][1]
                    + cocg[c_id][2] *rhsv[c_id][2]
                    + f_ext[c_id][2];
  });

  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_grad = std::chrono::high_resolution_clock::now();

  CS_FREE(rhsv);

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_STANDARD, ctx.use_gpu(), grad);

  if (cs_glob_timer_kernels_flag > 0) {
    t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_cocg - t_start);
    printf(", cocg = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_init - t_cocg);
    printf(", init = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_i_faces - t_init);
    printf(", i_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_ext_cells - t_i_faces);
    printf(", ext_cells = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_b_faces - t_ext_cells);
    printf(", b_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_grad - t_b_faces);
    printf(", gradient = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_grad);
    printf(", halo = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Compute cell gradient by least-squares reconstruction with a volume force
 * generating a hydrostatic pressure component, with a gather algorithm.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight_s     <-- weighted gradient coefficient variable,
 *                      or nullptr
 *   grad           --> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_gradient_hyd_p_gather(const cs_mesh_t                *m,
                                  const cs_mesh_quantities_t     *fvq,
                                  cs_halo_type_t                  halo_type,
                                  bool                            recompute_cocg,
                                  cs_real_t                       inc,
                                  const cs_real_3_t               f_ext[],
                                  const cs_field_bc_coeffs_t     *bc_coeffs,
                                  const cs_real_t                 pvar[],
                                  const cs_real_t       *restrict c_weight_s,
                                  cs_real_3_t           *restrict grad)
{
  std::chrono::high_resolution_clock::time_point t_start, t_cocg, t_cells, \
    t_stop;

  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  /* Weighting requires face weights, so cannot be applied consistently
     with an extended neighborhhood. */
  if (c_weight_s != nullptr)
    halo_type = CS_HALO_STANDARD;

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_t *restrict weight = fvq->weight;

  cs_dispatch_context ctx;

  /* Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0 = nullptr;
  cs_real_t *i_poro_duq_1 = nullptr;
  cs_real_t *b_poro_duq = nullptr;

  bool is_porous = false;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = true;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  }

  cs_cocg_6_t  *restrict cocgb = nullptr;
  cs_cocg_6_t  *restrict cocg = nullptr;

  _get_cell_cocg_lsq(m,
                     halo_type,
                     ctx.use_gpu(),
                     fvq,
                     &cocg,
                     &cocgb);

  /* Compute cocg and save contribution at boundaries */

  if (recompute_cocg) {
    _recompute_lsq_scalar_cocg(m,
                               fvq,
                               bc_coeffs,
                               ctx,
                               cocgb,
                               cocg);
    ctx.wait();
  }

  if (cs_glob_timer_kernels_flag > 0)
    t_cocg = std::chrono::high_resolution_clock::now();

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2c = ma->cell_cells;
  const cs_lnum_t *c2f = ma->cell_i_faces;
  short int *c2f_sgn = ma->cell_i_faces_sgn;
  if (c2f == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
    c2f = ma->cell_i_faces;
    c2f_sgn = ma->cell_i_faces_sgn;
  }
  const cs_lnum_t *c2c_e_idx = nullptr;
  const cs_lnum_t *c2c_e = nullptr;
  if (halo_type == CS_HALO_EXTENDED) {
    c2c_e_idx = ma->cell_cells_e_idx;
    c2c_e = ma->cell_cells_e;
  }
  const cs_lnum_t *restrict c2b_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict c2b = ma->cell_b_faces;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {

    cs_real_t rhsv[3] = {0, 0, 0};

    /* Contribution from interior faces
       -------------------------------- */

    cs_lnum_t s_id = c2c_idx[ii];
    cs_lnum_t e_id = c2c_idx[ii+1];

    if (c_weight_s != nullptr) {  /* With cell weighting */

      const cs_real_t w_ii = c_weight_s[ii];

      for (cs_lnum_t i = s_id; i < e_id; i++) {
        const cs_lnum_t jj = c2c[i];
        const cs_lnum_t f_id = c2f[i];
        const cs_real_t w_jj = c_weight_s[jj];

        cs_real_t poro[2] = {0, 0};
        if (is_porous) {
          poro[0] = i_poro_duq_0[f_id];
          poro[1] = i_poro_duq_1[f_id];
        }

        cs_real_t dc[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

        cs_real_t pfac =   (  pvar[jj] - pvar[ii]
                            + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                             cell_cen[ii],
                                                             f_ext[ii])
                            + poro[0]
                            - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                             cell_cen[jj],
                                                             f_ext[jj])
                            - poro[1])
                         * ddc;

        cs_real_t pond = (c2f_sgn[i] > 0) ? weight[f_id] : 1. - weight[f_id];

        pfac /= (  pond       *w_ii
                 + (1. - pond)*w_jj);

        cs_real_t fctb[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          fctb[ll] = dc[ll] * pfac;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ll] += w_jj * fctb[ll];
      }

    }
    else {  /* Without cell weighting */

      for (cs_lnum_t i = s_id; i < e_id; i++) {
        const cs_lnum_t jj = c2c[i];
        const cs_lnum_t f_id = c2f[i];

        cs_real_t poro[2] = {0, 0};
        if (is_porous) {
          poro[0] = i_poro_duq_0[f_id];
          poro[1] = i_poro_duq_1[f_id];
        }

        cs_real_t dc[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

        cs_real_t pfac =   (  pvar[jj] - pvar[ii]
                            + cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                             cell_cen[ii],
                                                             f_ext[ii])
                            + poro[0]
                            - cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                             cell_cen[jj],
                                                             f_ext[jj])
                            - poro[1])
                         * ddc;

        cs_real_t fctb[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          fctb[ll] = dc[ll] * pfac;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ll] += fctb[ll];
      }

    }  /* End for face-adjacent cells */

    /* Contribution from extended neighborhood;
       We assume that the middle of the segment joining cell centers
       may replace the center of gravity of a fictitious face. */

    if (c2c_e_idx != nullptr) {

      s_id = c2c_e_idx[ii];
      e_id = c2c_e_idx[ii+1];

      for (cs_lnum_t cidx = s_id; cidx < e_id; cidx++) {

        cs_lnum_t jj = c2c_e[cidx];

        /* Note: replaced the expressions:
         *  a) ptmid = 0.5 * (cell_cen[jj] - cell_cen[ii])
         *  b)   (cell_cen[ii] - ptmid) * f_ext[ii]
         *  c) - (cell_cen[jj] - ptmid) * f_ext[jj]
         * with:
         *  a) dc = cell_cen[jj] - cell_cen[ii]
         *  b) - 0.5 * dc * f_ext[ii]
         *  c) - 0.5 * dc * f_ext[jj]
         */

        cs_real_t pfac, dc[3], fctb[4];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

        pfac =   (  pvar[jj] - pvar[ii]
                  - 0.5 * cs_math_3_dot_product(dc, f_ext[ii])
                  - 0.5 * cs_math_3_dot_product(dc, f_ext[jj]))
                * ddc;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          fctb[ll] = dc[ll] * pfac;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ll] += fctb[ll];

      }

    } /* End for extended neighborhood */

    /* Contribution from boundary faces */

    s_id = c2b_idx[ii];
    e_id = c2b_idx[ii+1];

    for (cs_lnum_t fidx = s_id; fidx < e_id; fidx++) {
      const cs_lnum_t f_id = c2b[fidx];

      cs_real_t poro = (is_porous) ? b_poro_duq[f_id] : 0.;

      cs_real_t unddij = 1. / b_dist[f_id];
      cs_real_t umcbdd = (1. - coefbp[f_id]) * unddij;

      cs_real_t dddij[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dddij[ll] =   b_face_u_normal[f_id][ll]
                    + umcbdd * diipb[f_id][ll];

      cs_real_t pfac
        =   (coefap[f_id]*inc
            + (  (coefbp[f_id] -1.)
               * (  pvar[ii]
                   /* (b_face_cog - cell_cen).f_ext, or IF.F_i */
                  + cs_math_3_distance_dot_product(cell_cen[ii],
                                                   b_face_cog[f_id],
                                                   f_ext[ii])
                  + poro)))
            * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ll] += dddij[ll] * pfac;

    } /* loop on faces */

    /* Compute gradient */
    /*------------------*/

    grad[ii][0] =   cocg[ii][0] *rhsv[0]
                  + cocg[ii][3] *rhsv[1]
                  + cocg[ii][5] *rhsv[2]
                  + f_ext[ii][0];
    grad[ii][1] =   cocg[ii][3] *rhsv[0]
                  + cocg[ii][1] *rhsv[1]
                  + cocg[ii][4] *rhsv[2]
                  + f_ext[ii][1];
    grad[ii][2] =   cocg[ii][5] *rhsv[0]
                  + cocg[ii][4] *rhsv[1]
                  + cocg[ii][2] *rhsv[2]
                  + f_ext[ii][2];

  }); /* loop on cells */
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_cells = std::chrono::high_resolution_clock::now();

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_STANDARD, ctx.use_gpu(), grad);

  if (cs_glob_timer_kernels_flag > 0) {
    t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_cocg - t_start);
    printf(", cocgb = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_cells - t_cocg);
    printf(", cells= %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_cells);
    printf(", halo = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Compute cell gradient using least-squares reconstruction for non-orthogonal
 * meshes in the anisotropic case.
 *
 * cocg is computed to account for variable B.C.'s (flux).
 *
 * parameters:
 *   m            <-- pointer to associated mesh structure
 *   fvq          <-- pointer to associated finite volume quantities
 *   inc          <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs    <-- B.C. structure for boundary face normals
 *   pvar         <-- variable
 *   c_weight     <-- weighted gradient coefficient variable
 *   grad         <-> gradient of pvar (halo prepared for periodicity
 *                    of rotation)
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_gradient_ani(const cs_mesh_t               *m,
                         const cs_mesh_quantities_t    *fvq,
                         cs_real_t                      inc,
                         const cs_field_bc_coeffs_t    *bc_coeffs,
                         const cs_real_t                pvar[],
                         const cs_real_t              (*restrict c_weight)[6],
                         cs_real_t                    (*restrict grad)[3])
{
  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_t *restrict weight = fvq->weight;

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2c = ma->cell_cells;
  const cs_lnum_t *c2f = ma->cell_i_faces;
  short int *c2f_sgn = ma->cell_i_faces_sgn;
  if (c2f == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
    c2f = ma->cell_i_faces;
    c2f_sgn = ma->cell_i_faces_sgn;
  }
  const cs_lnum_t *restrict c2b_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict c2b = ma->cell_b_faces;

  cs_dispatch_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {

    cs_real_t rhsv[3] = {0, 0, 0};
    cs_real_t cocg[6] = {0, 0, 0, 0, 0, 0};

    /* Contribution from interior faces
       -------------------------------- */

    cs_lnum_t s_id = c2c_idx[ii];
    cs_lnum_t e_id = c2c_idx[ii+1];

    const cs_real_t *cell_cen_ii = cell_cen[ii];
    const cs_real_t *wi = c_weight[ii];

    for (cs_lnum_t i = s_id; i < e_id; i++) {
      const cs_lnum_t jj = c2c[i];
      const cs_lnum_t f_id = c2f[i];
      const cs_real_t *cell_cen_jj = cell_cen[jj];
      const cs_real_t *wj = c_weight[jj];

      cs_real_t fw = (c2f_sgn[i] > 0) ? weight[f_id] : 1. - weight[f_id];

      cs_real_t dc[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dc[ll] = cell_cen_jj[ll] - cell_cen_ii[ll];

      cs_real_t sum[6];
      for (cs_lnum_t ll = 0; ll < 6; ll++)
        sum[ll] = fw*wi[ll] + (1. - fw)*wj[ii];

      /* cocg contribution (inverse of the face viscosity tensor and
         anisotropic vector taking into account the weight coefficients) */

      /* Note: K_i.K_f^-1 = SUM.K_j^-1
       *       K_j.K_f^-1 = SUM.K_i^-1
       * So: K_i d = SUM.K_j^-1.IJ */

      cs_real_t inv_wj[6], _d[3], ki_d[3];
      cs_math_sym_33_inv_cramer(wj, inv_wj);
      cs_math_sym_33_3_product(inv_wj, dc,  _d);
      cs_math_sym_33_3_product(sum, _d, ki_d);

      /* 1 / ||Ki. K_f^-1. IJ||^2 */

      cs_real_t i_dci = 1. / cs_math_3_square_norm(ki_d);

      cocg[0] += ki_d[0] * ki_d[0] * i_dci;
      cocg[1] += ki_d[1] * ki_d[1] * i_dci;
      cocg[2] += ki_d[2] * ki_d[2] * i_dci;
      cocg[3] += ki_d[0] * ki_d[1] * i_dci;
      cocg[4] += ki_d[1] * ki_d[2] * i_dci;
      cocg[5] += ki_d[0] * ki_d[2] * i_dci;

      /* RHS contribution */

      /* (P_j - P_i)*/
      cs_real_t p_diff = (pvar[jj] - pvar[ii]);

      for (cs_lnum_t ll = 0; ll < 3; ll++) {
        rhsv[ll] += p_diff * ki_d[ll] * i_dci;
      }

    }

    /* Contribution from boundary faces */

    s_id = c2b_idx[ii];
    e_id = c2b_idx[ii+1];

    for (cs_lnum_t fidx = s_id; fidx < e_id; fidx++) {
      const cs_lnum_t f_id = c2b[fidx];

      cs_real_t umcbdd = (1. - coefbp[f_id]) / b_dist[f_id];
      cs_real_t unddij = 1. / b_dist[f_id];

      cs_real_t dsij[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dsij[ll] =   b_face_u_normal[f_id][ll]
                   + umcbdd*diipb[f_id][ll];

      /* cocg contribution */

      cocg[0] += dsij[0]*dsij[0];
      cocg[1] += dsij[1]*dsij[1];
      cocg[2] += dsij[2]*dsij[2];
      cocg[3] += dsij[0]*dsij[1];
      cocg[4] += dsij[1]*dsij[2];
      cocg[5] += dsij[0]*dsij[2];

      /* RHS contribution */

      cs_real_t pfac =   (coefap[f_id]*inc + (coefbp[f_id] -1.)*pvar[ii])
                       * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ll] += dsij[ll] * pfac;

    } /* loop on boundary faces */

    /* Invert cocg
       ----------- */

    _math_6_inv_cramer_sym_in_place(cocg);

    /* Compute gradient
       ---------------- */

    grad[ii][0] =   cocg[0] * rhsv[0]
                  + cocg[3] * rhsv[1]
                  + cocg[5] * rhsv[2];
    grad[ii][1] =   cocg[3] * rhsv[0]
                  + cocg[1] * rhsv[1]
                  + cocg[4] * rhsv[2];
    grad[ii][2] =   cocg[5] * rhsv[0]
                  + cocg[4] * rhsv[1]
                  + cocg[2] * rhsv[2];

  }); // End of parallel dispatch on cells

  /* Synchronize halos */

  ctx.wait();
  cs_halo_sync_r(m->halo, CS_HALO_STANDARD, false, grad);
}

/*----------------------------------------------------------------------------
 * Reconstruct the gradient of a scalar using a given gradient of
 * this scalar (typically lsq).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   w_stride       <-- stride for weighting coefficient
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. strucutre for boundary face normals
 *   c_weight       <-- weighted gradient coefficient variable
 *   c_var          <-- variable
 *   r_grad         <-- gradient used for reconstruction
 *   grad           <-> gradient of c_var (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_reconstruct_scalar_gradient(const cs_mesh_t                 *m,
                             const cs_mesh_quantities_t      *fvq,
                             int                              w_stride,
                             int                              hyd_p_flag,
                             cs_real_t                        inc,
                             const cs_real_t                  f_ext[][3],
                             const cs_field_bc_coeffs_t      *bc_coeffs,
                             const cs_real_t                  c_weight[],
                             const cs_real_t                  c_var[],
                             cs_real_3_t            *restrict r_grad,
                             cs_real_3_t            *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  const cs_real_3_t *restrict dofij = fvq->dofij;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const cs_real_33_t *restrict corr_grad_lin = fvq->corr_grad_lin;

  const cs_real_t *c_weight_s = nullptr;
  const cs_real_6_t *c_weight_t = nullptr;

  if (c_weight != nullptr) {
    if (w_stride == 1)
      c_weight_s = c_weight;
    else if (w_stride == 6)
      c_weight_t = (const cs_real_6_t *)c_weight;
  }

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0 = nullptr;
  cs_real_t *i_poro_duq_1 = nullptr;
  cs_real_t *b_poro_duq = nullptr;

  cs_lnum_t is_porous = 0;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  }

  bool warped_correction = (  cs_glob_mesh_quantities_flag
                            & CS_BAD_CELLS_WARPED_CORRECTION) ? true : false;

#if defined(HAVE_CUDA)
  bool accel = (cs_get_device_id() > -1 && hyd_p_flag == 0) ? true : false;

  if (accel) {
    const cs_mesh_adjacencies_t *madj = cs_glob_mesh_adjacencies;

    cs_gradient_strided_gg_r_cuda(m,
                                  madj,
                                  fvq,
                                  CS_HALO_EXTENDED,
                                  inc,
                                  warped_correction,
                                  (const cs_real_t (*)[1])coefap,
                                  (const cs_real_t (*)[1][1])coefbp,
                                  (const cs_real_t (*)[1])c_var,
                                  c_weight,
                                  (const cs_real_t (*)[1][3])r_grad,
                                  (cs_real_t (*)[1][3])grad);
    return;
  }
#endif

  std::chrono::high_resolution_clock::time_point t_start, t_init, t_i_faces, \
    t_b_faces, t_rescale, t_stop;

  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_dispatch_context ctx;
  cs_dispatch_context ctx_b;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx_b.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_b.use_gpu())
    ctx_b.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  /* Initialize gradient */
  /*---------------------*/

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    for (cs_lnum_t j = 0; j < 3; j++)
      grad[cell_id][j] = 0.0;
  });
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_init = std::chrono::high_resolution_clock::now();

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  if (hyd_p_flag == 1) {

    /* Contribution from interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t c_id1 = i_face_cells[f_id][0];
      cs_lnum_t c_id2 = i_face_cells[f_id][1];

      cs_real_t ktpond = weight[f_id]; /* no cell weighting */
      /* if cell weighting is active */
      if (c_weight_s != nullptr) {
        ktpond =   weight[f_id] * c_weight_s[c_id1]
                 / (       weight[f_id] * c_weight_s[c_id1]
                    + (1.0-weight[f_id])* c_weight_s[c_id2]);
      }
      else if (c_weight_t != nullptr) {
        cs_real_t sum[6], inv_sum[6];

        for (cs_lnum_t ii = 0; ii < 6; ii++)
          sum[ii] =        weight[f_id] *c_weight_t[c_id1][ii]
                    + (1.0-weight[f_id])*c_weight_t[c_id2][ii];

        cs_math_sym_33_inv_cramer(sum, inv_sum);

        ktpond =   weight[f_id] / 3.0
                 * (  inv_sum[0]*c_weight_t[c_id1][0]
                    + inv_sum[1]*c_weight_t[c_id1][1]
                    + inv_sum[2]*c_weight_t[c_id1][2]
                    + 2.0 * (  inv_sum[3]*c_weight_t[c_id1][3]
                             + inv_sum[4]*c_weight_t[c_id1][4]
                             + inv_sum[5]*c_weight_t[c_id1][5]));
      }

      cs_real_t poro[2] = {0., 0.};
      if (is_porous) {
        poro[0] = i_poro_duq_0[f_id];
        poro[1] = i_poro_duq_1[f_id];
      };

      cs_real_t  fexd[3];
      fexd[0] = 0.5 * (f_ext[c_id1][0] + f_ext[c_id2][0]);
      fexd[1] = 0.5 * (f_ext[c_id1][1] + f_ext[c_id2][1]);
      fexd[2] = 0.5 * (f_ext[c_id1][2] + f_ext[c_id2][2]);

      cs_real_t d_oi[3], d_oj[3];
      d_oi[0] = i_face_cog[f_id][0] - cell_cen[c_id1][0];
      d_oi[1] = i_face_cog[f_id][1] - cell_cen[c_id1][1];
      d_oi[2] = i_face_cog[f_id][2] - cell_cen[c_id1][2];
      d_oj[0] = i_face_cog[f_id][0] - cell_cen[c_id2][0];
      d_oj[1] = i_face_cog[f_id][1] - cell_cen[c_id2][1];
      d_oj[2] = i_face_cog[f_id][2] - cell_cen[c_id2][2];

      /*
        Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                  + (1-\alpha_\ij) \varia_\cellj\f$
                 but for the cell \f$ \celli \f$ we remove
                 \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                 and for the cell \f$ \cellj \f$ we remove
                 \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
      */

      cs_real_t pfaci
        =  ktpond
             * (  d_oi[0]*f_ext[c_id1][0]
                + d_oi[1]*f_ext[c_id1][1]
                + d_oi[2]*f_ext[c_id1][2]
                + poro[0])
        +  (1.0 - ktpond)
             * (  d_oj[0]*f_ext[c_id2][0]
                + d_oj[1]*f_ext[c_id2][1]
                + d_oj[2]*f_ext[c_id2][2]
                + poro[1]);

      cs_real_t pfacj = pfaci;
      cs_real_t d_var = c_var[c_id2] - c_var[c_id1];

      pfaci += (1.0-ktpond) * d_var;
      pfacj -=      ktpond  * d_var;

      /* Reconstruction part */
      cs_real_t rfac =
               - weight[f_id]      * cs_math_3_dot_product(d_oi, fexd)
           -  (1.0 - weight[f_id]) * cs_math_3_dot_product(d_oj, fexd)
           + (  dofij[f_id][0] * (r_grad[c_id1][0]+r_grad[c_id2][0])
              + dofij[f_id][1] * (r_grad[c_id1][1]+r_grad[c_id2][1])
              + dofij[f_id][2] * (r_grad[c_id1][2]+r_grad[c_id2][2])) * 0.5;

      cs_real_t face_normal[3], rhsv1[3], rhsv2[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        face_normal[j] = i_f_face_normal[f_id][j];
        rhsv1[j] =   (pfaci + rfac) * face_normal[j];
        rhsv2[j] = - (pfacj + rfac) * face_normal[j];
      }

      cs_dispatch_sum<3>(grad[c_id1], rhsv1, i_sum_type);
      cs_dispatch_sum<3>(grad[c_id2], rhsv2, i_sum_type);

    }); /* loop on faces */

    if (cs_glob_timer_kernels_flag > 0) {
      ctx.wait();  // Using an event would be preferred here
      t_i_faces = std::chrono::high_resolution_clock::now();
    }

    /* Contribution from boundary faces */

    ctx_b.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t c_id = b_face_cells[f_id];

      cs_real_t poro = (is_porous) ? b_poro_duq[f_id] : 0.;

      /*
        Remark: for the cell \f$ \celli \f$ we remove
                 \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
      */

      cs_real_t pfac
        = coefap[f_id] * inc
          + coefbp[f_id]
              /* (b_face_cog - cell_cen).f_ext, or IF.F_i */
            * (  cs_math_3_distance_dot_product(cell_cen[c_id],
                                                b_face_cog[f_id],
                                                f_ext[c_id])
               + poro);

      pfac += (coefbp[f_id] - 1.0) * c_var[c_id];

      /* Reconstruction part */
      cs_real_t
        rfac = coefbp[f_id]
               * (  diipb[f_id][0] * (r_grad[c_id][0] - f_ext[c_id][0])
                  + diipb[f_id][1] * (r_grad[c_id][1] - f_ext[c_id][1])
                  + diipb[f_id][2] * (r_grad[c_id][2] - f_ext[c_id][2]));

      cs_real_t rhsv[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        rhsv[j] = (pfac + rfac) * b_f_face_normal[f_id][j];
      }
      cs_dispatch_sum<3>(grad[c_id], rhsv, b_sum_type);

    }); /* loop on faces */

  } /* End of test on hydrostatic pressure */

  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  else {

    /* Contribution from interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t c_id1 = i_face_cells[f_id][0];
      cs_lnum_t c_id2 = i_face_cells[f_id][1];

      cs_real_t ktpond = weight[f_id]; /* no cell weighting */
      /* if cell weighting is active */
      if (c_weight_s != nullptr) {
        ktpond =   weight[f_id] * c_weight_s[c_id1]
                 / (       weight[f_id] *c_weight_s[c_id1]
                    + (1.0-weight[f_id])*c_weight_s[c_id2]);
      }
      else if (c_weight_t != nullptr) {
        cs_real_t sum[6], inv_sum[6];

        for (cs_lnum_t ii = 0; ii < 6; ii++)
          sum[ii] =        weight[f_id] *c_weight_t[c_id1][ii]
                    + (1.0-weight[f_id])*c_weight_t[c_id2][ii];

        cs_math_sym_33_inv_cramer(sum, inv_sum);

        ktpond =   weight[f_id] / 3.0
                 * (  inv_sum[0]*c_weight_t[c_id1][0]
                    + inv_sum[1]*c_weight_t[c_id1][1]
                    + inv_sum[2]*c_weight_t[c_id1][2]
                    + 2.0 * (  inv_sum[3]*c_weight_t[c_id1][3]
                             + inv_sum[4]*c_weight_t[c_id1][4]
                             + inv_sum[5]*c_weight_t[c_id1][5]));
      }

      /*
         Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                    + (1-\alpha_\ij) \varia_\cellj\f$
                 but for the cell \f$ \celli \f$ we remove
                 \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                 and for the cell \f$ \cellj \f$ we remove
                 \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
      */

      cs_real_t d_var = c_var[c_id2] - c_var[c_id1];
      cs_real_t pfaci = (1.0-ktpond) * d_var;
      cs_real_t pfacj =     -ktpond  * d_var;

      /* Reconstruction part */
      cs_real_t rfac = 0.5 *
                ( dofij[f_id][0]*(r_grad[c_id1][0]+r_grad[c_id2][0])
                 +dofij[f_id][1]*(r_grad[c_id1][1]+r_grad[c_id2][1])
                 +dofij[f_id][2]*(r_grad[c_id1][2]+r_grad[c_id2][2]));

      cs_real_t face_normal[3], rhsv1[3], rhsv2[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        face_normal[j] = i_f_face_normal[f_id][j];
        rhsv1[j] =   (pfaci + rfac) * face_normal[j];
        rhsv2[j] = - (pfacj + rfac) * face_normal[j];
      }

      cs_dispatch_sum<3>(grad[c_id1], rhsv1, i_sum_type);
      cs_dispatch_sum<3>(grad[c_id2], rhsv2, i_sum_type);

    }); /* loop on faces */

    if (cs_glob_timer_kernels_flag > 0) {
      ctx.wait();  // Using an event would be preferred here
      t_i_faces = std::chrono::high_resolution_clock::now();
    }

    /* Contribution from boundary faces */

    ctx_b.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t c_id = b_face_cells[f_id];

      /*
        Remark: for the cell \f$ \celli \f$ we remove
                \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
      */

      cs_real_t pfac =   inc*coefap[f_id]
                       + (coefbp[f_id]-1.0)*c_var[c_id];

      /* Reconstruction part */
      cs_real_t
        rfac =   coefbp[f_id]
               * (  diipb[f_id][0] * r_grad[c_id][0]
                  + diipb[f_id][1] * r_grad[c_id][1]
                  + diipb[f_id][2] * r_grad[c_id][2]);

      cs_real_t rhsv[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        rhsv[j] = (pfac + rfac) * b_f_face_normal[f_id][j];
      }
      cs_dispatch_sum<3>(grad[c_id], rhsv, b_sum_type);

    }); /* loop on faces */

  }

  ctx_b.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_b_faces = std::chrono::high_resolution_clock::now();

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    grad[c_id][0] *= dvol;
    grad[c_id][1] *= dvol;
    grad[c_id][2] *= dvol;

    if (warped_correction) {
      cs_real_t gradpa[3];
      for (cs_lnum_t i = 0; i < 3; i++) {
        gradpa[i] = grad[c_id][i];
        grad[c_id][i] = 0.;
      }

      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          grad[c_id][i] += corr_grad_lin[c_id][i][j] * gradpa[j];
      }
    }
  });

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();  // Using an event would be preferred here
    t_rescale = std::chrono::high_resolution_clock::now();
  }

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_EXTENDED, ctx.use_gpu(), grad);

  if (cs_glob_timer_kernels_flag > 0) {
    t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s (hyd_p %d)", cs_glob_rank_id, __func__, hyd_p_flag);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_init - t_start);
    printf(", init = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_i_faces - t_init);
    printf(", i_faces (hyd_p %d) = %ld", hyd_p_flag, elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_b_faces - t_i_faces);
    printf(", b_faces (hyd_p %d) = %ld", hyd_p_flag, elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_rescale - t_b_faces);
    printf(", rescale = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_rescale);
    printf(", halo = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Compute boundary face scalar values using least-squares reconstruction
 * for non-orthogonal meshes.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   nswrgp         <-- number of sweeps for gradient reconstruction
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs      <-- B.C. strucutre for boundary face normals
 *   c_var          <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or nullptr
 *   b_f_var        --> boundary face value.
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_b_face_val(const cs_mesh_t             *m,
                       const cs_mesh_quantities_t  *fvq,
                       cs_halo_type_t               halo_type,
                       cs_real_t                    inc,
                       const cs_field_bc_coeffs_t  *bc_coeffs,
                       const cs_real_t              c_var[],
                       const cs_real_t              c_weight[],
                       cs_real_t          *restrict b_f_var)
{
  const cs_lnum_t n_b_cells = m->n_b_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_b_faces = ma->cell_b_faces;

  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  cs_field_bc_coeffs_t *bc_coeffs_loc = nullptr;

  if (inc < 1) {
    CS_MALLOC(bc_coeffs_loc, 1, cs_field_bc_coeffs_t);
    cs_field_bc_coeffs_shallow_copy(bc_coeffs, bc_coeffs_loc);

    CS_MALLOC(bc_coeffs_loc->a, m->n_b_faces, cs_real_t);
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      bc_coeffs_loc->a[i] = 0;
  }

  const cs_field_bc_coeffs_t *_bc_coeffs
    = (bc_coeffs_loc != nullptr) ?
      (const cs_field_bc_coeffs_t *)bc_coeffs_loc :
      (const cs_field_bc_coeffs_t *)bc_coeffs;

  /* Reconstruct gradients using least squares for non-orthogonal meshes */

# pragma omp parallel for if (n_b_cells > CS_THR_MIN)
  for (cs_lnum_t ci = 0; ci < n_b_cells; ci++) {

    cs_real_t grad[3];
    cs_lnum_t c_id = m->b_cells[ci];

    cs_gradient_scalar_cell(m,
                            fvq,
                            c_id,
                            halo_type,
                            _bc_coeffs,
                            c_var,
                            c_weight,
                            grad);

    /* Update boundary face value */

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t f_id = cell_b_faces[i];

      cs_real_t pip =   c_var[c_id]
                      + cs_math_3_dot_product(diipb[f_id], grad);
      b_f_var[f_id] = _bc_coeffs->a[f_id]*inc + pip*_bc_coeffs->b[f_id];

    }
  }

  if (bc_coeffs_loc != nullptr) {
    cs_field_bc_coeffs_free_copy(bc_coeffs, bc_coeffs_loc);
    CS_FREE(bc_coeffs_loc);
  }
}

/*----------------------------------------------------------------------------
 * Compute boundary face scalar values using least-squares reconstruction
 * for non-orthogonal meshes in the presence of a volume force generating
 * a hydrostatic pressure component.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   nswrgp         <-- number of sweeps for gradient reconstruction
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   c_var          <-- cell variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or nullptr
 *   b_f_var        --> boundary face value.
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_b_face_val_phyd(const cs_mesh_t             *m,
                            const cs_mesh_quantities_t  *fvq,
                            cs_halo_type_t               halo_type,
                            cs_real_t                    inc,
                            const cs_real_t              f_ext[][3],
                            const cs_field_bc_coeffs_t  *bc_coeffs,
                            const cs_real_t              c_var[],
                            const cs_real_t              c_weight[],
                            cs_real_t          *restrict b_f_var)
{
  const cs_real_t *bc_coeff_a = bc_coeffs->a;
  const cs_real_t *bc_coeff_b = bc_coeffs->b;

  const cs_lnum_t n_b_cells = m->n_b_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict cell_cells_idx = ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_e_idx = ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells = ma->cell_cells;
  const cs_lnum_t *restrict cell_cells_e = ma->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces = ma->cell_b_faces;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /*Additional terms due to porosity */

  cs_field_t *f_b_poro_duq_0 = cs_field_by_name_try("b_poro_duq_0");

  cs_real_t *b_poro_duq = nullptr;
  cs_lnum_t is_porous = false;
  if (f_b_poro_duq_0 != nullptr) {
    is_porous = 1;
    b_poro_duq = f_b_poro_duq_0->val;
  }

  /* Reconstruct gradients using least squares for non-orthogonal meshes */

# pragma omp parallel for if (n_b_cells > CS_THR_MIN)
  for (cs_lnum_t ci = 0; ci < n_b_cells; ci++) {

    cs_lnum_t c_id = m->b_cells[ci];

    cs_real_t cocg[6] = {0., 0., 0., 0., 0., 0.};
    cs_real_t rhsv[3] = {0., 0., 0.};

    int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

    for (int adj_id = 0; adj_id < n_adj; adj_id++) {

      const cs_lnum_t *restrict cell_cells_p;
      cs_lnum_t s_id, e_id;

      if (adj_id == 0) {
        s_id = cell_cells_idx[c_id];
        e_id = cell_cells_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *)(cell_cells);
      }
      else if (cell_cells_e_idx != nullptr) {
        s_id = cell_cells_e_idx[c_id];
        e_id = cell_cells_e_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *)(cell_cells_e);
      }
      else
        break;

      if (c_weight == nullptr) {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[c_id1][ll] - cell_cen[c_id][ll];

          cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          cs_real_t pfac = (  c_var[c_id1] - c_var[c_id]
                            - 0.5 * dc[0] * f_ext[c_id][0]
                            - 0.5 * dc[1] * f_ext[c_id][1]
                            - 0.5 * dc[2] * f_ext[c_id][2]
                            - 0.5 * dc[0] * f_ext[c_id1][0]
                            - 0.5 * dc[1] * f_ext[c_id1][1]
                            - 0.5 * dc[2] * f_ext[c_id1][2]) * ddc;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ll] += dc[ll] * pfac;

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

        }

      }
      else {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[c_id1][ll] - cell_cen[c_id][ll];

          cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          cs_real_t pfac = (  c_var[c_id1] - c_var[c_id]
                            - 0.5 * dc[0] * f_ext[c_id][0]
                            - 0.5 * dc[1] * f_ext[c_id][1]
                            - 0.5 * dc[2] * f_ext[c_id][2]
                            - 0.5 * dc[0] * f_ext[c_id1][0]
                            - 0.5 * dc[1] * f_ext[c_id1][1]
                            - 0.5 * dc[2] * f_ext[c_id1][2]) * ddc;

          cs_real_t _weight =   2. * c_weight[c_id1]
                              / (c_weight[c_id] + c_weight[c_id1]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ll] += dc[ll] * pfac* _weight;

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

        }

      }

    } /* end of contribution from interior and extended cells */

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    /* Contribution from boundary faces */

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_real_t  dsij[3];

      cs_lnum_t f_id = cell_b_faces[i];

      cs_real_t poro = (is_porous) ? b_poro_duq[f_id] : 0;

      cs_real_t unddij = 1. / b_dist[f_id];
      cs_real_t umcbdd = (1. -bc_coeff_b[f_id]) * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dsij[ll] = b_face_u_normal[f_id][ll] + umcbdd*diipb[f_id][ll];

      /* (b_face_cog - cell_cen).f_ext, or IF.F_i */
      cs_real_t c_f_ext
        = cs_math_3_distance_dot_product(cell_cen[c_id],
                                         b_face_cog[f_id],
                                         f_ext[c_id]);

      cs_real_t pfac =  (  bc_coeff_a[f_id]*inc + (bc_coeff_b[f_id] -1.)
                         * (c_var[c_id] + c_f_ext + poro))
                       * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ll] += dsij[ll] * pfac;

      cocg[0] += dsij[0]*dsij[0];
      cocg[1] += dsij[1]*dsij[1];
      cocg[2] += dsij[2]*dsij[2];
      cocg[3] += dsij[0]*dsij[1];
      cocg[4] += dsij[1]*dsij[2];
      cocg[5] += dsij[0]*dsij[2];

    } // end of contribution from boundary cells

    /* Invert */

    cs_real_t a11 = cocg[1]*cocg[2] - cocg[4]*cocg[4];
    cs_real_t a12 = cocg[4]*cocg[5] - cocg[3]*cocg[2];
    cs_real_t a13 = cocg[3]*cocg[4] - cocg[1]*cocg[5];
    cs_real_t a22 = cocg[0]*cocg[2] - cocg[5]*cocg[5];
    cs_real_t a23 = cocg[3]*cocg[5] - cocg[0]*cocg[4];
    cs_real_t a33 = cocg[0]*cocg[1] - cocg[3]*cocg[3];

    cs_real_t det_inv = 1. / (cocg[0]*a11 + cocg[3]*a12 + cocg[5]*a13);

    cs_real_t grad[3];

    grad[0] =  (  a11 * rhsv[0]
                + a12 * rhsv[1]
                + a13 * rhsv[2]) * det_inv
                + f_ext[c_id][0];
    grad[1] =  (  a12 * rhsv[0]
                + a22 * rhsv[1]
                + a23 * rhsv[2]) * det_inv
                + f_ext[c_id][1];
    grad[2] =  (  a13 * rhsv[0]
                + a23 * rhsv[1]
                + a33 * rhsv[2]) * det_inv
                + f_ext[c_id][2];

    /* Update boundary face value */

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t f_id = cell_b_faces[i];

      cs_real_t pip =   c_var[c_id]
                      + cs_math_3_dot_product(diipb[f_id], grad);
      b_f_var[f_id] = bc_coeff_a[f_id]*inc + pip*bc_coeff_b[f_id];

    }
  }
}

/*----------------------------------------------------------------------------
 * Compute gradient using vertex-based face values for scalar gradient
 * reconstruction.
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   w_stride       <-- stride for weighting coefficient
 *   halo_type      <-- halo type (extended or not)
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   c_var          <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_fv_vtx_based_scalar_gradient(const cs_mesh_t                *m,
                              const cs_mesh_quantities_t     *fvq,
                              int                             w_stride,
                              int                             hyd_p_flag,
                              cs_real_t                       inc,
                              const cs_real_3_t               f_ext[],
                              const cs_field_bc_coeffs_t     *bc_coeffs,
                              const cs_real_t                 c_var[],
                              const cs_real_t                 c_weight[],
                              cs_real_3_t           *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_t *bc_coeff_b = (const cs_real_t *)bc_coeffs->b;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  const cs_real_t *c_weight_s = nullptr;
  const cs_real_6_t *c_weight_t = nullptr;

  if (c_weight != nullptr) {
    if (w_stride == 1)
      c_weight_s = c_weight;
    else if (w_stride == 6)
      c_weight_t = (const cs_real_6_t *)c_weight;
  }

  /* Additional terms due to porosity */

  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  cs_lnum_t is_porous = 0;
  if (f_i_poro_duq_0 != nullptr) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name_try("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name_try("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  /* Initialize gradient
     ------------------- */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      grad[c_id][j] = 0.0;
  }

  /* Pre-compute values at boundary using least squares */

  cs_real_t *b_f_var;
  CS_MALLOC(b_f_var, m->n_b_faces, cs_real_t);

  if (hyd_p_flag == 1)
    _lsq_scalar_b_face_val_phyd(m,
                                fvq,
                                CS_HALO_STANDARD,
                                inc,
                                f_ext,
                                bc_coeffs,
                                c_var,
                                c_weight,
                                b_f_var);
  else
    _lsq_scalar_b_face_val(m,
                           fvq,
                           CS_HALO_STANDARD,
                           inc,
                           bc_coeffs,
                           c_var,
                           c_weight,
                           b_f_var);

  /* Compute vertex-based values
     --------------------------- */

  cs_real_t *v_var;
  CS_MALLOC(v_var, m->n_vertices, cs_real_t);

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0, /* verbosity */
                    1, /* var_dim */
                    0, /* tr_dim */
                    c_weight,
                    c_var,
                    b_f_var,
                    v_var);

  /* Interpolate to face-based values
     -------------------------------- */

  cs_real_t *i_f_var;
  CS_MALLOC(i_f_var, m->n_i_faces, cs_real_t);

  for (int f_t = 0; f_t < 2; f_t++) {

    const cs_lnum_t n_faces = (f_t == 0) ? m->n_i_faces : m->n_b_faces;
    const cs_lnum_t *f2v_idx= nullptr, *f2v_ids = nullptr;
    cs_real_t *f_var = nullptr;

    if (f_t == 0) {
      f2v_idx = m->i_face_vtx_idx;
      f2v_ids = m->i_face_vtx_lst;
      f_var = i_f_var;
    }
    else {
      f2v_idx = m->b_face_vtx_idx;
      f2v_ids = m->b_face_vtx_lst;
      f_var = b_f_var;
    }

#   pragma omp parallel for if (n_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
      cs_lnum_t s_id = f2v_idx[f_id];
      cs_lnum_t e_id = f2v_idx[f_id+1];
      cs_real_t s = 0;
      for (cs_lnum_t i = s_id; i < e_id; i++)
        s += v_var[f2v_ids[i]];
      f_var[f_id] = s / (e_id-s_id);
    }

  }

  /* Vertex values are not needed after this stage */

  CS_FREE(v_var);

  /* Case with hydrostatic pressure
     ------------------------------ */

  if (hyd_p_flag == 1) {

    /* Contribution from interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t ktpond = weight[f_id]; /* no cell weighting */
          /* if cell weighting is active */
          if (c_weight_s != nullptr) {
            ktpond =   weight[f_id] * c_weight_s[ii]
                     / (       weight[f_id] * c_weight_s[ii]
                        + (1.0-weight[f_id])* c_weight_s[jj]);
          }
          else if (c_weight_t != nullptr) {
            cs_real_t sum[6], inv_sum[6];

            for (cs_lnum_t kk = 0; kk < 6; kk++)
              sum[kk] =        weight[f_id]*c_weight_t[ii][kk]
                        + (1.0-weight[f_id])*c_weight_t[jj][kk];

            cs_math_sym_33_inv_cramer(sum, inv_sum);

            ktpond =   weight[f_id] / 3.0
                     * (  inv_sum[0]*c_weight_t[ii][0]
                        + inv_sum[1]*c_weight_t[ii][1]
                        + inv_sum[2]*c_weight_t[ii][2]
                        + 2.0 * (  inv_sum[3]*c_weight_t[ii][3]
                                 + inv_sum[4]*c_weight_t[ii][4]
                                 + inv_sum[5]*c_weight_t[ii][5]));
          }

          cs_real_2_t poro = {
            i_poro_duq_0[is_porous*f_id],
            i_poro_duq_1[is_porous*f_id]
          };

          cs_real_t pfaci
            =  ktpond
                 * (  (i_face_cog[f_id][0] - cell_cen[ii][0])*f_ext[ii][0]
                    + (i_face_cog[f_id][1] - cell_cen[ii][1])*f_ext[ii][1]
                    + (i_face_cog[f_id][2] - cell_cen[ii][2])*f_ext[ii][2]
                    + poro[0])
            +  (1.0 - ktpond)
                 * (  (i_face_cog[f_id][0] - cell_cen[jj][0])*f_ext[jj][0]
                    + (i_face_cog[f_id][1] - cell_cen[jj][1])*f_ext[jj][1]
                    + (i_face_cog[f_id][2] - cell_cen[jj][2])*f_ext[jj][2]
                    + poro[1]);
          cs_real_t pfacj = pfaci;

          pfaci += i_f_var[f_id] - c_var[ii];
          pfacj += i_f_var[f_id] - c_var[jj];

          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[f_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[f_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             f_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t c_id = b_face_cells[f_id];

          cs_real_t poro = b_poro_duq[is_porous*f_id];

          /*
            Remark: for the cell \f$ \celli \f$ we remove
                    \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
          */

          cs_real_t pfac = b_f_var[f_id] - c_var[c_id];

          pfac +=  bc_coeff_b[f_id]
                  * (  cs_math_3_distance_dot_product(cell_cen[c_id],
                                                      b_face_cog[f_id],
                                                      f_ext[c_id])
                     + poro);

          for (cs_lnum_t j = 0; j < 3; j++)
            grad[c_id][j] += pfac * b_f_face_normal[f_id][j];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  } /* End of test on hydrostatic pressure */

  /* Standard case, without hydrostatic pressure
     ------------------------------------------- */

  else {

    /* Contribution from interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t pfaci = i_f_var[f_id] - c_var[ii];
          cs_real_t pfacj = i_f_var[f_id] - c_var[jj];

          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[f_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[f_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             f_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = b_face_cells[f_id];

          /*
            Remark: for the cell \f$ \celli \f$ we remove
                    \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
          */

          cs_real_t pfac = b_f_var[f_id] - c_var[ii];

          for (cs_lnum_t j = 0; j < 3; j++)
            grad[ii][j] += pfac * b_f_face_normal[f_id][j];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  }

  CS_FREE(i_f_var);
  CS_FREE(b_f_var);

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    for (cs_lnum_t j = 0; j < 3; j++)
      grad[c_id][j] *= dvol;
  }

  /* Synchronize halos */

  cs_halo_sync_r(m->halo, CS_HALO_EXTENDED, false, grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square of the Frobenius norm of a symmetric tensor
 *
 * \param[in]  t
 *
 * \return the square of the norm
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_tensor_norm_2(const cs_real_t  t[6])
{
  cs_real_t retval =     t[0]*t[0] +   t[1]*t[1] +   t[2]*t[2]
                     + 2*t[3]*t[3] + 2*t[4]*t[4] + 2*t[5]*t[5];
  return retval;
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a vector or tensor if necessary.
 * This function deals with the standard or extended neighborhood.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   ma             <-- mesh adjacencies
 *   halo_type      <-- halo type (extended or not)
 *   clip_mode      <-- type of clipping for the computation of the gradient
 *   verbosity      <-- output level
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   var_name       <-- name of the current variable
 *   pvar           <-- variable
 *   grad           <-> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_strided_gradient_clipping(const cs_mesh_t              *m,
                           const cs_mesh_quantities_t   *fvq,
                           const cs_mesh_adjacencies_t  *ma,
                           cs_halo_type_t                halo_type,
                           int                           clip_mode,
                           int                           verbosity,
                           cs_real_t                     climgp,
                           const char                   *var_name,
                           const cs_real_t    (*restrict pvar)[stride],
                           cs_real_t          (*restrict grad)[stride][3])
{
  const cs_real_t  clipp_coef_sq = climgp*climgp;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;

  const cs_halo_t *halo = m->halo;

  if (clip_mode < 0 || climgp < 0)
    return;

  cs_dispatch_context ctx;

  /* The gradient and the variable must be already synchronized */

  cs_real_t *restrict clip_factor = _get_clip_factor_try(var_name);
  cs_real_t *_clip_factor = nullptr;

  if (clip_factor == nullptr) {
    CS_MALLOC(_clip_factor, n_cells_ext, cs_real_t);
    clip_factor = _clip_factor;
  }

  const int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

  const size_t block_size = 128;
  const size_t n_blocks = cs_parall_block_count(n_cells, block_size);

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clip_mode == CS_GRADIENT_LIMIT_CELL) {

#   pragma omp parallel for  if  (n_cells > CS_THR_MIN)
    for (size_t b_id = 0; b_id < n_blocks; b_id++) {

      cs_lnum_t s_id = b_id*block_size, e_id = (b_id+1)*block_size;
      if (e_id > n_cells) e_id = n_cells;

      /* Remark:
         denum: maximum l2 norm of the variation of the gradient squared
         denom: maximum l2 norm of the variation of the variable squared */

      cs_real_t denum[block_size];
      cs_real_t denom[block_size];

      cs_lnum_t b_e_id = e_id - s_id;

      for (cs_lnum_t i = 0; i < b_e_id; i++) {
        denum[i] = 0;
        denom[i] = 0;
      }

      for (int adj_id = 0; adj_id < n_adj; adj_id++) {

        const cs_lnum_t *restrict cell_cells_idx;
        const cs_lnum_t *restrict cell_cells;

        if (adj_id == 0) {
          cell_cells_idx = ma->cell_cells_idx;
          cell_cells = ma->cell_cells;
        }
        else if (ma->cell_cells_e_idx != nullptr) {
          cell_cells_idx = ma->cell_cells_e_idx;
          cell_cells = ma->cell_cells_e;
        }
        else
          break;

        for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */

          cs_lnum_t c_id1 = i + s_id;

          for (cs_lnum_t cidx = cell_cells_idx[c_id1];
               cidx < cell_cells_idx[c_id1+1];
               cidx++) {

            cs_lnum_t c_id2 = cell_cells[cidx];

            cs_real_t dist[3], grad_dist1[stride], var_dist[stride];

            for (cs_lnum_t k = 0; k < 3; k++)
              dist[k] = cell_cen[c_id1][k] - cell_cen[c_id2][k];

            for (cs_lnum_t k = 0; k < stride; k++) {

              grad_dist1[k] =   grad[c_id1][k][0] * dist[0]
                              + grad[c_id1][k][1] * dist[1]
                              + grad[c_id1][k][2] * dist[2];

              var_dist[k] = pvar[c_id1][k] - pvar[c_id2][k];

            }

            cs_real_t  dvar_sq, dist_sq1;

            if (stride == 3) {
              dist_sq1 = cs_math_3_square_norm(grad_dist1);
              dvar_sq = cs_math_3_square_norm(var_dist);
            }
            else if (stride == 6) {
              dist_sq1 = _tensor_norm_2(grad_dist1);
              dvar_sq = _tensor_norm_2(var_dist);
            }

            denum[i] = std::max(denum[i], dist_sq1);
            denom[i] = std::max(denom[i], dvar_sq);

          }

        }  /* End of loop on block elements */

      } /* End of loop on adjacency type */

      for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */
        cs_real_t factor1 = 1.;
        if (denum[i] > clipp_coef_sq * denom[i])
          factor1 = sqrt(clipp_coef_sq * denom[i]/denum[i]);

        clip_factor[s_id + i] = factor1;
      }

    } /* End of (parallel) loop on blocks */

  } /* End for clip_mode == CS_GRADIENT_LIMIT_CELL */

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clip_mode == CS_GRADIENT_LIMIT_FACE) {

    cs_real_t *factor;
    CS_MALLOC(factor, n_cells_ext, cs_real_t);
    cs_array_real_set_scalar(n_cells_ext, DBL_MAX, factor);

#   pragma omp parallel for  if  (n_cells > CS_THR_MIN)
    for (size_t b_id = 0; b_id < n_blocks; b_id++) {

      cs_lnum_t s_id = b_id*block_size, e_id = (b_id+1)*block_size;
      if (e_id > n_cells) e_id = n_cells;

      /* Remark:
         denum: maximum l2 norm of the variation of the gradient squared
         denom: maximum l2 norm of the variation of the variable squared */

      cs_real_t denum[block_size];
      cs_real_t denom[block_size];

      cs_lnum_t b_e_id = e_id - s_id;

      for (cs_lnum_t i = 0; i < b_e_id; i++) {
        denum[i] = 0;
        denom[i] = 0;
        clip_factor[s_id + i] = 1;
      }

      for (int adj_id = 0; adj_id < n_adj; adj_id++) {

        const cs_lnum_t *restrict cell_cells_idx;
        const cs_lnum_t *restrict cell_cells;

        if (adj_id == 0) {
          cell_cells_idx = ma->cell_cells_idx;
          cell_cells = ma->cell_cells;
        }
        else if (ma->cell_cells_e_idx != nullptr) {
          cell_cells_idx = ma->cell_cells_e_idx;
          cell_cells = ma->cell_cells_e;
        }
        else
          break;

        for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */

          cs_lnum_t c_id1 = i + s_id;

          for (cs_lnum_t cidx = cell_cells_idx[c_id1];
               cidx < cell_cells_idx[c_id1+1];
               cidx++) {

            cs_lnum_t c_id2 = cell_cells[cidx];

            cs_real_t dist[3], grad_dist1[stride], var_dist[stride];

            for (cs_lnum_t k = 0; k < 3; k++)
              dist[k] = cell_cen[c_id1][k] - cell_cen[c_id2][k];

            for (cs_lnum_t k = 0; k < stride; k++) {
              grad_dist1[k]
                = 0.5 * (  (grad[c_id1][k][0] + grad[c_id2][k][0]) * dist[0]
                         + (grad[c_id1][k][1] + grad[c_id2][k][1]) * dist[1]
                         + (grad[c_id1][k][2] + grad[c_id2][k][2]) * dist[2]);
              var_dist[k] = pvar[c_id1][k] - pvar[c_id2][k];
            }

            cs_real_t dist_sq1, dvar_sq;

            if (stride == 3) {
              dist_sq1 = cs_math_3_square_norm(grad_dist1);
              dvar_sq = cs_math_3_square_norm(var_dist);
            }
            else if (stride == 6) {
              dist_sq1 = _tensor_norm_2(grad_dist1);
              dvar_sq = _tensor_norm_2(var_dist);
            }

            denum[i] = std::max(denum[i], dist_sq1);
            denom[i] = std::max(denom[i], dvar_sq);

          }

        }  /* End of loop on block elements */

      } /* End of loop on adjacency type */

      for (cs_lnum_t i = 0; i < b_e_id; i++) { /* Loop on block elements */
        cs_real_t factor1 = 1.;
        if (denum[i] > clipp_coef_sq * denom[i])
          factor1 = sqrt(clipp_coef_sq * denom[i]/denum[i]);

        factor[s_id + i] = factor1;
      }

    } /* End of (parallel) loop on blocks */

    /* Now compute clip factor (kernel common to scalar and strided clases */

    _gradient_update_face_clip_factor(m, ma, halo_type, factor, clip_factor);

    CS_FREE(factor);

  } /* End for clip_mode == CS_GRADIENT_LIMIT_FACE */

  /* Synchronize variable */

  if (halo != nullptr) {
    cs_halo_sync_var(m->halo, halo_type, clip_factor);
  }

  /* Apply clip factor to gradient
     ----------------------------- */

  if (verbosity > 1) {

    cs_gnum_t n_clip = 0;
    cs_real_t min_factor = 0, max_factor = 0, mean_factor = 0;

    cs_array_reduce_simple_stats_l(ctx, 1, n_cells, nullptr, clip_factor,
                                   &min_factor,
                                   &max_factor,
                                   &mean_factor);

#   pragma omp parallel for reduction(+:n_clip) if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t factor1 = clip_factor[c_id];
      if (factor1 < 1) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            grad[c_id][i][j] *= factor1;
        }
        n_clip += 1;
      }
    }

    cs_real_t buf[2] = {-min_factor, max_factor};
    cs_parall_max(2, CS_REAL_TYPE, buf);
    min_factor = -buf[0];
    max_factor =  buf[1];

    cs_parall_sum(1, CS_REAL_TYPE, &mean_factor);
    mean_factor /= (cs_real_t)(m->n_g_cells);

    cs_parall_counter(&n_clip, 1);

    bft_printf
      (_(" Variable: %s; gradient limitation in %llu cells\n"
         "   minimum factor = %g; maximum factor = %g; mean factor = %g\n"),
       var_name,
       (unsigned long long)n_clip, min_factor, max_factor, mean_factor);
  }

  else { /* Avoid unneeded reduction if no logging */

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t factor1 = clip_factor[c_id];
      if (factor1 < 1) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            grad[c_id][i][j] *= factor1;
        }
      }

    } /* End of loop on cells */

  }

  /* Synchronize the updated Gradient */

  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, stride*3);

    if (cs_glob_mesh->have_rotation_perio) {
      if (stride == 3)
        cs_halo_perio_sync_var_tens(m->halo,
                                    halo_type,
                                    (cs_real_t *)grad);
      else if (stride == 6)
        cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                             halo_type,
                                             (cs_real_t *)grad);
    }
  }

  CS_FREE(_clip_factor);
}

/*----------------------------------------------------------------------------
 * Initialize the gradient of a vector for gradient reconstruction.
 *
 * A non-reconstructed gradient is computed at this stage.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs_v    <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_initialize_vector_gradient(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            int                           inc,
                            const cs_field_bc_coeffs_t   *bc_coeffs_v,
                            const cs_real_3_t   *restrict pvar,
                            const cs_real_t     *restrict c_weight,
                            cs_real_33_t        *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_3_t  *restrict coefav
    = (const cs_real_3_t *)bc_coeffs_v->a;
  const cs_real_33_t *restrict coefbv
    = (const cs_real_33_t *)bc_coeffs_v->b;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;

  /* Computation without reconstruction */
  /*------------------------------------*/

  /* Initialization */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] = 0.0;
    }
  }

  /* Interior faces contribution */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           f_id++) {

        cs_lnum_t c_id1 = i_face_cells[f_id][0];
        cs_lnum_t c_id2 = i_face_cells[f_id][1];

        cs_real_t pond = weight[f_id];

        cs_real_t ktpond = (c_weight == nullptr) ?
          pond :                    // no cell weighting
          pond * c_weight[c_id1] // cell weighting active
            / (      pond * c_weight[c_id1]
              + (1.0-pond)* c_weight[c_id2]);

        /*
           Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                    + (1-\alpha_\ij) \varia_\cellj\f$
                   but for the cell \f$ \celli \f$ we remove
                   \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                   and for the cell \f$ \cellj \f$ we remove
                   \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
        */
        for (cs_lnum_t i = 0; i < 3; i++) {
          cs_real_t pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
          cs_real_t pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[c_id1][i][j] += pfaci * i_f_face_normal[f_id][j];
            grad[c_id2][i][j] -= pfacj * i_f_face_normal[f_id][j];
          }
        }

      } /* End of loop on faces */

    } /* End of loop on threads */

  } /* End of loop on thread groups */

  /* Boundary face treatment */

# pragma omp parallel for
  for (int t_id = 0; t_id < n_b_threads; t_id++) {

    for (cs_lnum_t f_id = b_group_index[t_id*2];
         f_id < b_group_index[t_id*2 + 1];
         f_id++) {

      cs_lnum_t c_id = b_face_cells[f_id];

      for (cs_lnum_t i = 0; i < 3; i++) {
        cs_real_t pfac = inc*coefav[f_id][i];

        for (cs_lnum_t k = 0; k < 3; k++) {
          if (i == k)
            pfac += (coefbv[f_id][i][k] - 1.0) * pvar[c_id][k];
          else
            pfac += coefbv[f_id][i][k] * pvar[c_id][k];
        }

        for (cs_lnum_t j = 0; j < 3; j++)
          grad[c_id][i][j] += pfac * b_f_face_normal[f_id][j];
      }

    } /* loop on faces */

  } /* loop on threads */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] *= dvol;
    }
  }

  /* Periodicity and parallelism treatment */

  _sync_strided_gradient_halo<3>(m, halo_type, grad);
}

/*----------------------------------------------------------------------------
 * Green-Gauss reconstruction of the gradient of a vector or tensor using
 * an initial gradient of this quantity (typically lsq).
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   pvar           <-- variable
 *   val_f          <-- face value for gradient
 *   c_weight       <-- weighted gradient coefficient variable
 *   r_grad         <-- gradient used for reconstruction
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_reconstruct_strided_gradient
(
  const cs_mesh_t                               *m,
  [[maybe_unused]] const cs_mesh_adjacencies_t  *madj,
  const cs_mesh_quantities_t                    *fvq,
  cs_halo_type_t                                 halo_type,
  const cs_real_t                     (*restrict pvar)[stride],
  const cs_real_t                     (*restrict val_f)[stride],
  const cs_real_t                      *restrict c_weight,
  cs_real_t                           (*restrict r_grad)[stride][3],
  cs_real_t                           (*restrict grad)[stride][3]
)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  bool warped_correction = (  cs_glob_mesh_quantities_flag
                            & CS_BAD_CELLS_WARPED_CORRECTION) ? true : false;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *restrict dofij = fvq->dofij;
  const cs_real_33_t *restrict corr_grad_lin = fvq->corr_grad_lin;

#if defined(HAVE_CUDA)
  bool accel = (cs_get_device_id() > -1) ? true : false;

  if (accel) {
    cs_gradient_strided_gg_r_cuda(m,
                                  madj,
                                  fvq,
                                  halo_type,
                                  warped_correction,
                                  val_f,
                                  pvar,
                                  c_weight,
                                  r_grad,
                                  grad);

    return;
  }
#endif

  cs_dispatch_context ctx;
  cs_dispatch_context ctx_b;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx_b.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_b.use_gpu())
    ctx_b.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  /* Initialize gradient */
  /*---------------------*/

  std::chrono::high_resolution_clock::time_point t_start, t_init, t_i_faces, \
    t_b_faces, t_rescale, t_stop;

  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  /* Initialization */

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] = 0.0;
    }
  });
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_init = std::chrono::high_resolution_clock::now();

  /* Interior faces contribution */

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    cs_real_t pond = weight[f_id];

    cs_real_t ktpond = (c_weight == nullptr) ?
      pond :                    // no cell weighting
      pond * c_weight[c_id1] // cell weighting active
           / (      pond * c_weight[c_id1]
              + (1.0-pond)* c_weight[c_id2]);

    /*
      Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                + (1-\alpha_\ij) \varia_\cellj\f$
               but for the cell \f$ \celli \f$ we remove
               \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
               and for the cell \f$ \cellj \f$ we remove
               \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
    */

    for (cs_lnum_t i = 0; i < stride; i++) {

      cs_real_t pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
      cs_real_t pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

      /* Reconstruction part */
      cs_real_t rfac = 0.5 * (  dofij[f_id][0]*(  r_grad[c_id1][i][0]
                                                + r_grad[c_id2][i][0])
                              + dofij[f_id][1]*(  r_grad[c_id1][i][1]
                                                + r_grad[c_id2][i][1])
                              + dofij[f_id][2]*(  r_grad[c_id1][i][2]
                                                + r_grad[c_id2][i][2]));

      cs_real_t ctb1[3], ctb2[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        ctb1[j] =   (pfaci + rfac) * i_f_face_normal[f_id][j];
        ctb2[j] = - (pfacj + rfac) * i_f_face_normal[f_id][j];
      }

      if (c_id1 < n_cells)
        cs_dispatch_sum<3>(grad[c_id1][i], ctb1, i_sum_type);
      if (c_id2 < n_cells)
        cs_dispatch_sum<3>(grad[c_id2][i], ctb2, i_sum_type);

    }

  }); /* End of loop on faces */

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();  // Using an event would be preferred here
    t_i_faces = std::chrono::high_resolution_clock::now();
  }

  /* Boundary face treatment */

  ctx_b.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t c_id = b_face_cells[f_id];

    /*
      Remark: for the cell \f$ \celli \f$ we remove
              \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
    */

    for (cs_lnum_t i = 0; i < stride; i++) {

      cs_real_t ctb[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        ctb[j] =  (val_f[f_id][i] - pvar[c_id][i]) * b_f_face_normal[f_id][j];
      }
      cs_dispatch_sum<3>(grad[c_id][i], ctb, b_sum_type);
    }

  }); /* loop on faces */
  ctx_b.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_b_faces = std::chrono::high_resolution_clock::now();

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] *= dvol;
    }

    if (warped_correction) {
      cs_real_t gradpa[3];
      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          gradpa[j] = grad[c_id][i][j];
          grad[c_id][i][j] = 0.;
        }

        for (cs_lnum_t j = 0; j < 3; j++)
          for (cs_lnum_t k = 0; k < 3; k++)
            grad[c_id][i][j] += corr_grad_lin[c_id][j][k] * gradpa[k];
      }
    }
  });
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0)
    t_rescale = std::chrono::high_resolution_clock::now();

  /* Periodicity and parallelism treatment */

  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, stride*3);
    if (cs_glob_mesh->have_rotation_perio) {
      if (stride == 3)
        cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
      else if (stride == 6)
        cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                             halo_type,
                                             (cs_real_t *)grad);
    }
  }

  if (cs_glob_timer_kernels_flag > 0) {
    t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s<%d>", cs_glob_rank_id, __func__, stride);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_init - t_start);
    printf(", init = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_i_faces - t_init);
    printf(", i_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_b_faces - t_i_faces);
    printf(", b_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_rescale - t_b_faces);
    printf(", rescale = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_rescale);
    printf(", halo = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Compute the gradient of a vector with an iterative technique in order to
 * handle non-orthoganalities (n_r_sweeps > 1).
 *
 * We do not take into account any volumic force here.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-> pointer to associated finite volume quantities
 *   var_name       <-- variable's name
 *   gradient_info  <-- pointer to performance logging structure, or nullptr
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_r_sweeps     --> >1: with reconstruction
 *   verbosity      --> verbosity level
 *   epsrgp         --> precision for iterative gradient calculation
 *   coefav         <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           <-> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_iterative_vector_gradient(const cs_mesh_t               *m,
                           const cs_mesh_quantities_t    *fvq,
                           const char                    *var_name,
                           cs_gradient_info_t            *gradient_info,
                           cs_halo_type_t                 halo_type,
                           int                            inc,
                           int                            n_r_sweeps,
                           int                            verbosity,
                           cs_real_t                      epsrgp,
                           const cs_field_bc_coeffs_t    *bc_coeffs_v,
                           const cs_real_3_t   *restrict  pvar,
                           const cs_real_t               *c_weight,
                           cs_real_33_t         *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_3_t  *restrict coefav
    = (const cs_real_3_t *)bc_coeffs_v->a;
  const cs_real_33_t *restrict coefbv
    = (const cs_real_33_t *)bc_coeffs_v->b;

  int isweep = 0;

  cs_real_33_t *rhs;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  cs_gradient_quantities_t  *gq = _gradient_quantities_get(0);

  cs_real_33_t *restrict cocg = gq->cocg_it;
  if (cocg == nullptr)
    cocg = _compute_cell_cocg_it(m, fvq, gq);

  CS_MALLOC(rhs, n_cells_ext, cs_real_33_t);

  /* Gradient reconstruction to handle non-orthogonal meshes */
  /*---------------------------------------------------------*/

  /* L2 norm */

  cs_real_t l2_norm = _l2_norm_1(9*n_cells, (cs_real_t *)grad);
  cs_real_t l2_residual = l2_norm;

  if (l2_norm > cs_math_epzero) {

    /* Iterative process */
    /*-------------------*/

    for (isweep = 1;
         isweep < n_r_sweeps && l2_residual > epsrgp*l2_norm;
         isweep++) {

      /* Computation of the Right Hand Side*/

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            rhs[c_id][i][j] = -grad[c_id][i][j] * cell_vol[c_id];
        }
      }

      /* Interior face treatment */

      for (int g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_i_threads; t_id++) {

          for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               f_id++) {

            cs_lnum_t c_id1 = i_face_cells[f_id][0];
            cs_lnum_t c_id2 = i_face_cells[f_id][1];
            cs_real_t pond = weight[f_id];

            cs_real_t ktpond = (c_weight == nullptr) ?
              pond :                     // no cell weighting
              pond  * c_weight[c_id1] // cell weighting active
                / (      pond  * c_weight[c_id1]
                  + (1.0-pond) * c_weight[c_id2]);

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            for (cs_lnum_t i = 0; i < 3; i++) {

              /* Reconstruction part */
              cs_real_t
                pfaci = 0.5 * (    (grad[c_id1][i][0] + grad[c_id2][i][0])
                                 * dofij[f_id][0]
                               +   (grad[c_id1][i][1] + grad[c_id2][i][1])
                                 * dofij[f_id][1]
                               +   (grad[c_id1][i][2] + grad[c_id2][i][2])
                                 * dofij[f_id][2]);
              cs_real_t pfacj = pfaci;

              pfaci += (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
              pfacj -=      ktpond  * (pvar[c_id2][i] - pvar[c_id1][i]);

              for (cs_lnum_t j = 0; j < 3; j++) {
                rhs[c_id1][i][j] += pfaci * i_f_face_normal[f_id][j];
                rhs[c_id2][i][j] -= pfacj * i_f_face_normal[f_id][j];
              }
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Boundary face treatment */

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[t_id*2];
             f_id < b_group_index[t_id*2 + 1];
             f_id++) {

          cs_lnum_t c_id = b_face_cells[f_id];

          for (cs_lnum_t i = 0; i < 3; i++) {

            /*
              Remark: for the cell \f$ \celli \f$ we remove
              \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            cs_real_t pfac = inc*coefav[f_id][i];

            for (cs_lnum_t k = 0; k < 3; k++) {
              /* Reconstruction part */
              cs_real_t vecfac =   grad[c_id][k][0] * diipb[f_id][0]
                                 + grad[c_id][k][1] * diipb[f_id][1]
                                 + grad[c_id][k][2] * diipb[f_id][2];
              pfac += coefbv[f_id][i][k] * vecfac;

              if (i == k)
                pfac += (coefbv[f_id][i][k] - 1.0) * pvar[c_id][k];
              else
                pfac += coefbv[f_id][i][k] * pvar[c_id][k];
            }

            for (cs_lnum_t j = 0; j < 3; j++)
              rhs[c_id][i][j] += pfac * b_f_face_normal[f_id][j];

          }

        } /* loop on faces */

      } /* loop on threads */

      /* Increment of the gradient */

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t dvol;
        /* Is the cell disabled (for solid or porous)? Not the case if coupled */
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_vol[c_id];
        else
          dvol = 0.;

        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            rhs[c_id][i][j] *= dvol;
        }

        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            for (cs_lnum_t k = 0; k < 3; k++)
              grad[c_id][i][j] += rhs[c_id][i][k] * cocg[c_id][k][j];
          }
        }
      }

      /* Periodicity and parallelism treatment */

      if (m->halo != nullptr) {
        cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 9);
        if (cs_glob_mesh->have_rotation_perio)
          cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
      }

      /* Convergence test (L2 norm) */

      l2_residual = _l2_norm_1(9*n_cells, (cs_real_t *)rhs);

    } /* End of the iterative process */

    /* Printing */

    if (l2_residual < epsrgp*l2_norm) {
      if (verbosity >= 2) {
        bft_printf
          (_(" %s: isweep = %d, normed residual: %e, norm: %e, var: %s\n"),
           __func__, isweep, l2_residual/l2_norm, l2_norm, var_name);
      }
    }
    else if (isweep >= n_r_sweeps) {
      if (verbosity >= 0) {
        bft_printf(_(" Warning:\n"
                     " --------\n"
                     "   %s; variable: %s; sweeps: %d\n"
                     "   %*s  normed residual: %11.4e; norm: %11.4e\n"),
                   __func__, var_name, isweep,
                   (int)(strlen(__func__)), " ", l2_residual/l2_norm, l2_norm);
      }
    }

  }

  if (gradient_info != nullptr)
    _gradient_info_update_iter(gradient_info, isweep);

  CS_FREE(rhs);
}

/*----------------------------------------------------------------------------
 * Compute the gradient of a vector with an iterative technique in order to
 * handle non-orthoganalities (n_r_sweeps > 1).
 *
 * We do not take into account any volumic force here.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-> pointer to associated finite volume quantities
 *   var_name       <-- variable's name
 *   gradient_info  <-- pointer to performance logging structure, or nullptr
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_r_sweeps     --> >1: with reconstruction
 *   verbosity      --> verbosity level
 *   epsrgp         --> precision for iterative gradient calculation
 *   bc_coeffs_ts   <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   grad          <-> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_iterative_tensor_gradient(const cs_mesh_t              *m,
                           const cs_mesh_quantities_t   *fvq,
                           const char                   *var_name,
                           cs_gradient_info_t           *gradient_info,
                           cs_halo_type_t                halo_type,
                           int                           inc,
                           int                           n_r_sweeps,
                           int                           verbosity,
                           cs_real_t                     epsrgp,
                           const cs_field_bc_coeffs_t   *bc_coeffs_ts,
                           const cs_real_6_t   *restrict pvar,
                           cs_real_63_t        *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_6_t  *restrict coefat
    = (const cs_real_6_t *)bc_coeffs_ts->a;
  const cs_real_66_t *restrict coefbt
    = (const cs_real_66_t *)bc_coeffs_ts->b;

  int isweep = 0;

  cs_real_63_t *rhs;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  cs_gradient_quantities_t  *gq = _gradient_quantities_get(0);

  cs_real_33_t *restrict cocg = gq->cocg_it;
  if (cocg == nullptr)
    cocg = _compute_cell_cocg_it(m, fvq, gq);

  CS_MALLOC(rhs, n_cells_ext, cs_real_63_t);

  /* Gradient reconstruction to handle non-orthogonal meshes */
  /*---------------------------------------------------------*/

  /* L2 norm */

  cs_real_t l2_norm = _l2_norm_1(18*n_cells, (cs_real_t *)grad);
  cs_real_t l2_residual = l2_norm ;

  if (l2_norm > cs_math_epzero) {

    /* Iterative process */
    /*-------------------*/

    for (isweep = 1;
         isweep < n_r_sweeps && l2_residual > epsrgp*l2_norm;
         isweep++) {

      /* Computation of the Right Hand Side*/

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        for (cs_lnum_t i = 0; i < 6; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            rhs[c_id][i][j] = - cell_vol[c_id] * grad[c_id][i][j];
        }
      }

      /* Interior face treatment */

      for (int g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_i_threads; t_id++) {

          for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               f_id++) {

            cs_lnum_t c_id1 = i_face_cells[f_id][0];
            cs_lnum_t c_id2 = i_face_cells[f_id][1];
            cs_real_t pond = weight[f_id];

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            for (cs_lnum_t i = 0; i < 6; i++) {

              /* Reconstruction part */
              cs_real_t
                pfaci = 0.5 * (    (grad[c_id1][i][0] + grad[c_id2][i][0])
                                 * dofij[f_id][0]
                               +   (grad[c_id1][i][1] + grad[c_id2][i][1])
                                 * dofij[f_id][1]
                               +   (grad[c_id1][i][2] + grad[c_id2][i][2])
                                 * dofij[f_id][2]);
              cs_real_t pfacj = pfaci;

              pfaci += (1.0-pond) * (pvar[c_id2][i] - pvar[c_id1][i]);
              pfacj -=       pond * (pvar[c_id2][i] - pvar[c_id1][i]);
              for (cs_lnum_t j = 0; j < 3; j++) {
                rhs[c_id1][i][j] += pfaci * i_f_face_normal[f_id][j];
                rhs[c_id2][i][j] -= pfacj * i_f_face_normal[f_id][j];
              }
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Boundary face treatment */

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t f_id = b_group_index[t_id*2];
             f_id < b_group_index[t_id*2 + 1];
             f_id++) {

          cs_lnum_t c_id = b_face_cells[f_id];

          for (cs_lnum_t i = 0; i < 6; i++) {

            /*
              Remark: for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            cs_real_t pfac = inc*coefat[f_id][i];

            for (cs_lnum_t k = 0; k < 6; k++) {
              /* Reconstruction part */
              cs_real_t vecfac =   grad[c_id][k][0] * diipb[f_id][0]
                                 + grad[c_id][k][1] * diipb[f_id][1]
                                 + grad[c_id][k][2] * diipb[f_id][2];
              pfac += coefbt[f_id][i][k] * vecfac;

              if (i == k)
                pfac += (coefbt[f_id][i][k] - 1.0) * pvar[c_id][k];
              else
                pfac += coefbt[f_id][i][k] * pvar[c_id][k];
            }

            for (cs_lnum_t j = 0; j < 3; j++)
              rhs[c_id][i][j] += pfac * b_f_face_normal[f_id][j];

          }

        } /* loop on faces */

      } /* loop on threads */

      /* Increment of the gradient */

#     pragma omp parallel for
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t dvol;
        /* Is the cell disabled (for solid or porous)? Not the case if coupled */
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_vol[c_id];
        else
          dvol = 0.;

        for (cs_lnum_t i = 0; i < 6; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            rhs[c_id][i][j] *= dvol;
        }

        for (cs_lnum_t i = 0; i < 6; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            for (cs_lnum_t k = 0; k < 3; k++)
              grad[c_id][i][j] += rhs[c_id][i][k] * cocg[c_id][k][j];
          }
        }
      }

      /* Periodicity and parallelism treatment */

      if (m->halo != nullptr) {
        cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 18);
        if (cs_glob_mesh->have_rotation_perio)
          cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                               halo_type,
                                               (cs_real_t *)grad);
      }

      /* Convergence test (L2 norm) */

      ///FIXME
      l2_residual = _l2_norm_1(18*n_cells, (cs_real_t *)rhs);

    } /* End of the iterative process */

    /* Printing */

    if (l2_residual < epsrgp*l2_norm) {
      if (verbosity >= 2) {
        bft_printf
          (_(" %s: isweep = %d, normed residual: %e, norm: %e, var: %s\n"),
           __func__, isweep, l2_residual/l2_norm, l2_norm, var_name);
      }
    }
    else if (isweep >= n_r_sweeps) {
      if (verbosity >= 0) {
        bft_printf(_(" Warning:\n"
                     " --------\n"
                     "   %s; variable: %s; sweeps: %d\n"
                     "   %*s  normed residual: %11.4e; norm: %11.4e\n"),
                   __func__, var_name, isweep,
                   (int)(strlen(__func__)), " ", l2_residual/l2_norm, l2_norm);
      }
    }

  }

  if (gradient_info != nullptr)
    _gradient_info_update_iter(gradient_info, isweep);

  CS_FREE(rhs);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of a vector or tensor using least-squares
 * reconstruction for non-orthogonal meshes.
 *
 * template parameters:
 *   e2n           type of assembly algorithm used
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   pvar           <-- variable
 *   val_f          <-- face value for gradient
 *   c_weight       <-- weighted gradient coefficient variable, or nullptr
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <const cs_e2n_sum_t e2n, cs_lnum_t stride>
static void
_lsq_strided_gradient(const cs_mesh_t             *m,
                      const cs_mesh_adjacencies_t *madj,
                      const cs_mesh_quantities_t  *fvq,
                      cs_halo_type_t               halo_type,
                      const cs_real_t (*restrict pvar)[stride],
                      const cs_real_t (*restrict val_f)[stride],
                      const cs_real_t *restrict c_weight,
                      cs_real_t (*restrict grad)[stride][3])
{
  using grad_t = cs_real_t[stride][3];
  const cs_lnum_t n_cells     = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_cells = m->n_b_cells;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;

  const cs_lnum_t *restrict cell_cells_e_idx = madj->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx = madj->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells_e = madj->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces = madj->cell_b_faces;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  std::chrono::high_resolution_clock::time_point t_start, t_init, t_i_faces, \
    t_ext_n, t_b_faces, t_gradient, t_b_correction, t_halo, t_stop;

  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

#if defined(HAVE_CUDA)
  bool accel = (cs_get_device_id() > -1) ? true : false;
#else
  bool accel = false;
#endif

  cs_cocg_6_t *restrict cocgb = nullptr;
  cs_cocg_6_t *restrict cocg = nullptr;

  _get_cell_cocg_lsq(m, halo_type, accel, fvq, &cocg, &cocgb);

#if defined(HAVE_CUDA)

  if (accel) {
    cs_gradient_strided_lsq_cuda(m,
                                 madj,
                                 fvq,
                                 halo_type,
                                 val_f,
                                 pvar,
                                 c_weight,
                                 cocgb,
                                 cocg,
                                 grad);

    return;
  }

#endif

  /* Initialize RHS
     (note this could be done locally in the gather variants */

  grad_t *rhs;
  CS_MALLOC_HD(rhs, n_cells_ext, grad_t, cs_alloc_mode);

  if (cs_glob_timer_kernels_flag > 0)
    t_init = std::chrono::high_resolution_clock::now();

  /* Contribution from interior faces
     -------------------------------- */

  if (e2n != CS_E2N_SUM_GATHER) {

    cs_array_real_fill_zero(n_cells_ext*stride*3, (cs_real_t *)rhs);

    cs_lnum_t n_i_groups, n_i_threads;
    cs_mesh_i_faces_thread_block_count(m, e2n, 0, &n_i_groups, &n_i_threads);

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        cs_lnum_t s_id, e_id;
        cs_mesh_i_faces_thread_block_range(m, e2n, g_id, t_id, n_i_threads, 0,
                                           &s_id, &e_id);

        for (cs_lnum_t f_id = s_id; f_id < e_id; f_id++) {

          cs_lnum_t c_id1 = i_face_cells[f_id][0];
          cs_lnum_t c_id2 = i_face_cells[f_id][1];

          cs_real_t  dc[3];
          for (cs_lnum_t i = 0; i < 3; i++)
            dc[i] = cell_cen[c_id2][i] - cell_cen[c_id1][i];

          cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          if (c_weight != nullptr) {
            cs_real_t pond = weight[f_id];
            cs_real_t denom = 1. / (  pond       *c_weight[c_id1]
                                    + (1. - pond)*c_weight[c_id2]);

            for (cs_lnum_t i = 0; i < stride; i++) {
              cs_real_t pfac = denom * (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;

              if (e2n == CS_E2N_SUM_SCATTER) {
                cs_real_t fctb[3];
                for (cs_lnum_t j = 0; j < 3; j++) {
                  fctb[j] = dc[j] * pfac;
                  rhs[c_id1][i][j] += c_weight[c_id2] * fctb[j];
                  rhs[c_id2][i][j] += c_weight[c_id1] * fctb[j];
                }
              }
              else if (e2n == CS_E2N_SUM_SCATTER_ATOMIC) {
                cs_real_t fctb[3];
                for (cs_lnum_t j = 0; j < 3; j++) {
                  fctb[j] = dc[j] * pfac;
                  #pragma omp atomic
                  rhs[c_id1][i][j] += c_weight[c_id2] * fctb[j];
                  #pragma omp atomic
                  rhs[c_id2][i][j] += c_weight[c_id1] * fctb[j];
                }
              }

            }
          }

          else { /* if (c_weight == nullptr) */
            for (cs_lnum_t i = 0; i < stride; i++) {
              cs_real_t pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;

              if (e2n == CS_E2N_SUM_SCATTER) {
                cs_real_t fctb[3];
                for (cs_lnum_t j = 0; j < 3; j++) {
                  fctb[j] = dc[j] * pfac;
                  rhs[c_id1][i][j] += fctb[j];
                  rhs[c_id2][i][j] += fctb[j];
                }
              }
              else if (e2n == CS_E2N_SUM_SCATTER_ATOMIC) {
                cs_real_t fctb[3];
                for (cs_lnum_t j = 0; j < 3; j++) {
                  fctb[j] = dc[j] * pfac;
                  #pragma omp atomic
                  rhs[c_id1][i][j] += fctb[j];
                  #pragma omp atomic
                  rhs[c_id2][i][j] += fctb[j];
                }
              }
            }
          } /* c_weight */

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Assembly for 2-stage algorithm */

  } /* Test on e2n (template argument for algorithm type) */

  else if (e2n == CS_E2N_SUM_GATHER) {

    const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
    const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
    const cs_lnum_t *c2c = ma->cell_cells;
    const cs_lnum_t *c2f = ma->cell_i_faces;
    short int *c2f_sgn = ma->cell_i_faces_sgn;
    if (c2f == nullptr) {
      cs_mesh_adjacencies_update_cell_i_faces();
      c2f = ma->cell_i_faces;
      c2f_sgn = ma->cell_i_faces_sgn;
    }

    if (c_weight != nullptr) {  /* With cell weighting */

#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {

        const cs_lnum_t s_id = c2c_idx[ii];
        const cs_lnum_t e_id = c2c_idx[ii+1];

        const cs_real_t w_ii = c_weight[ii];

        for (cs_lnum_t k = 0; k < stride; k++) {
          for (cs_lnum_t l = 0; l < 3; l++)
            rhs[ii][k][l] = 0.;
        }

        for (cs_lnum_t i = s_id; i < e_id; i++) {
          const cs_lnum_t jj = c2c[i];
          const cs_lnum_t f_id = c2f[i];
          const cs_real_t w_jj = c_weight[jj];

          cs_real_t  dc[3];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

          cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          cs_real_t pond = (c2f_sgn[i] > 0) ? weight[f_id] : 1. - weight[f_id];

          cs_real_t denom = 1. / (  pond       *w_ii
                                  + (1. - pond)*w_jj);

          for (cs_lnum_t k = 0; k < stride; k++) {
            cs_real_t pfac = denom * (pvar[jj][k] - pvar[ii][k]) * ddc;

            cs_real_t fctb[3];
            for (cs_lnum_t l = 0; l < 3; l++) {
              fctb[l] = dc[l] * pfac;
              rhs[ii][k][l] += c_weight[jj] * fctb[l];
            }
          }
        }

      }  /* loop on cells */

    }

    else { /* if (c_weight == nullptr) */

#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {

        const cs_lnum_t s_id = c2c_idx[ii];
        const cs_lnum_t e_id = c2c_idx[ii+1];

        for (cs_lnum_t k = 0; k < stride; k++) {
          for (cs_lnum_t l = 0; l < 3; l++)
            rhs[ii][k][l] = 0.;
        }

        for (cs_lnum_t i = s_id; i < e_id; i++) {
          const cs_lnum_t jj = c2c[i];

          cs_real_t  dc[3];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

          cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          for (cs_lnum_t k = 0; k < stride; k++) {
            cs_real_t pfac = (pvar[jj][k] - pvar[ii][k]) * ddc;

            for (cs_lnum_t l = 0; l < 3; l++) {
              rhs[ii][k][l] += dc[l] * pfac;
            }
          }
        }

      }  /* loop on cells */

    }

  } /* (e2n == CS_E2N_SUM_GATHER) */

  if (cs_glob_timer_kernels_flag > 0)
    t_i_faces = std::chrono::high_resolution_clock::now();

  /* Contribution from extended neighborhood */

  if (cell_cells_e_idx != nullptr && halo_type == CS_HALO_EXTENDED) {

#   pragma omp parallel for
    for (cs_lnum_t c_id1 = 0; c_id1 < n_cells; c_id1++) {
      for (cs_lnum_t cidx = cell_cells_e_idx[c_id1];
           cidx < cell_cells_e_idx[c_id1+1];
           cidx++) {

        cs_lnum_t c_id2 = cell_cells_e[cidx];

        cs_real_t dc[3];

        for (cs_lnum_t i = 0; i < 3; i++)
          dc[i] = cell_cen[c_id2][i] - cell_cen[c_id1][i];

        cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (cs_lnum_t i = 0; i < stride; i++) {

          cs_real_t pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;

          for (cs_lnum_t j = 0; j < 3; j++) {
            rhs[c_id1][i][j] += dc[j] * pfac;
          }
        }
      }
    }

  } /* End for extended neighborhood */

  if (cs_glob_timer_kernels_flag > 0)
    t_ext_n = std::chrono::high_resolution_clock::now();

  /* Contribution from boundary faces */

# pragma omp parallel for  if (n_b_cells > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_b_cells; ii++) {

    cs_lnum_t c_id = m->b_cells[ii];

    for (cs_lnum_t ll = 0; ll < 6; ll++)
      cocg[c_id][ll] = cocgb[ii][ll];

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) { /* loop on boundary faces */

      cs_lnum_t f_id = cell_b_faces[i];

      cs_real_t dif[3];
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dif[ll] = b_face_cog[f_id][ll] - cell_cen[c_id][ll];

      cs_real_t ddif = 1. / cs_math_3_square_norm(dif);

      cocg[c_id][0] += dif[0]*dif[0]*ddif;
      cocg[c_id][1] += dif[1]*dif[1]*ddif;
      cocg[c_id][2] += dif[2]*dif[2]*ddif;
      cocg[c_id][3] += dif[0]*dif[1]*ddif;
      cocg[c_id][4] += dif[1]*dif[2]*ddif;
      cocg[c_id][5] += dif[0]*dif[2]*ddif;

      for (cs_lnum_t kk = 0; kk < stride; kk++) {
        cs_real_t pfac = (val_f[f_id][kk] - pvar[c_id][kk]) * ddif;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhs[c_id][kk][ll] += dif[ll] * pfac;
      }
    } /* loop on faces */

    _math_6_inv_cramer_sym_in_place(cocg[c_id]);

  } /* loop on boundary cells */

  if (cs_glob_timer_kernels_flag > 0)
    t_b_faces = std::chrono::high_resolution_clock::now();

  /* Compute gradient */
  /*------------------*/

  #pragma omp parallel for if(n_cells >= CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      grad[c_id][i][0] =   rhs[c_id][i][0] * cocg[c_id][0]
                         + rhs[c_id][i][1] * cocg[c_id][3]
                         + rhs[c_id][i][2] * cocg[c_id][5];

      grad[c_id][i][1] =   rhs[c_id][i][0] * cocg[c_id][3]
                         + rhs[c_id][i][1] * cocg[c_id][1]
                         + rhs[c_id][i][2] * cocg[c_id][4];

      grad[c_id][i][2] =   rhs[c_id][i][0] * cocg[c_id][5]
                         + rhs[c_id][i][1] * cocg[c_id][4]
                         + rhs[c_id][i][2] * cocg[c_id][2];

    }
  }

  if (cs_glob_timer_kernels_flag > 0)
    t_gradient = std::chrono::high_resolution_clock::now();

  /* Synchronize halos */

  _sync_strided_gradient_halo<stride>(m, halo_type, grad);

  if (cs_glob_timer_kernels_flag > 0)
    t_halo = std::chrono::high_resolution_clock::now();

  CS_FREE(rhs);

  if (cs_glob_timer_kernels_flag > 0) {
    t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s<%d>", cs_glob_rank_id, __func__, stride);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_init - t_start);
    printf(", init = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_i_faces - t_init);
    printf(", i_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_ext_n - t_i_faces);
    printf(", ext_n = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_b_faces - t_ext_n);
    printf(", b_faces = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_gradient - t_b_faces);
    printf(", gradient = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_halo - t_b_correction);
    printf(", halo = %ld", elapsed.count());

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of a scalar using least-squares reconstruction for
 * non-orthogonal meshes (n_r_sweeps > 1).
 * Calls the strided version with the correct types.
 *
 * template parameters:
 *   e2n           type of assembly algorithm used
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   pvar           <-- variable
 *   val_f          <-- face value for gradient
 *   c_weight       <-- weighted gradient coefficient variable, or nullptr
 *   gradv          --> gradient of pvar (du_i/dx_j : gradv[][j])
 *----------------------------------------------------------------------------*/

template <const cs_e2n_sum_t e2n>
static void
_lsq_strided_gradient(const cs_mesh_t             *m,
                      const cs_mesh_adjacencies_t *madj,
                      const cs_mesh_quantities_t  *fvq,
                      cs_halo_type_t               halo_type,
                      const cs_real_t    *restrict pvar,
                      const cs_real_t    *restrict val_f,
                      const cs_real_t    *restrict c_weight,
                      cs_real_3_t        *restrict gradv)
{
  _lsq_strided_gradient<e2n>(
    m,
    madj,
    fvq,
    halo_type,
    reinterpret_cast<const cs_real_t(*)[1]>(pvar),
    reinterpret_cast<const cs_real_t(*)[1]>(val_f),
    c_weight,
    reinterpret_cast<const cs_real_t(*)[1][3]>(gradv)
  );
}

template <size_t stride>
using lsq_strided_gradient_t = void(const cs_mesh_t *,
                                    const cs_mesh_adjacencies_t *,
                                    const cs_mesh_quantities_t *,
                                    cs_halo_type_t,
                                    const cs_real_t (*)[stride],
                                    const cs_real_t (*)[stride],
                                    const cs_real_t *restrict,
                                    cs_real_t (*)[stride][3]);

/*----------------------------------------------------------------------------
 * Compute gradient using vertex-based face values for scalar gradient
 * reconstruction.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   c_var          <-- variable
 *   val_f          <-- face value for gradient
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           --> gradient of c_var (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_fv_vtx_based_strided_gradient(const cs_mesh_t               *m,
                               const cs_mesh_quantities_t    *fvq,
                               const cs_real_t (*restrict c_var)[stride],
                               const cs_real_t (*restrict val_f)[stride],
                               const cs_real_t                c_weight[],
                               cs_real_t (*restrict grad)[stride][3])
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  using var_t = cs_real_t[stride];

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;

  /* Initialize gradient
     ------------------- */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[cell_id][i][j] = 0.0;
    }
  }

  /* Pre-compute gradient at boundary using least squares */

  var_t *b_f_var;
  CS_MALLOC(b_f_var, m->n_b_faces, var_t);
  cs_array_copy<cs_real_t>(stride*m->n_b_faces,
                           (const cs_real_t *)val_f,
                           (cs_real_t *)b_f_var);

  /* Compute vertex-based values
     --------------------------- */

  var_t *v_var;
  CS_MALLOC(v_var, m->n_vertices, var_t);

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0, /* verbosity */
                    stride, /* var_dim */
                    0,
                    c_weight,
                    (const cs_real_t *)c_var,
                    (const cs_real_t *)val_f,
                    (cs_real_t *)v_var);

  /* Interpolate to face-based values
     -------------------------------- */

  var_t *i_f_var;
  CS_MALLOC(i_f_var, m->n_i_faces, var_t);

  for (int f_t = 0; f_t < 2; f_t++) {

    const cs_lnum_t n_faces = (f_t == 0) ? m->n_i_faces : m->n_b_faces;
    const cs_lnum_t *f2v_idx= nullptr, *f2v_ids = nullptr;
    var_t *f_var = nullptr;

    if (f_t == 0) {
      f2v_idx = m->i_face_vtx_idx;
      f2v_ids = m->i_face_vtx_lst;
      f_var = i_f_var;
    }
    else {
      f2v_idx = m->b_face_vtx_idx;
      f2v_ids = m->b_face_vtx_lst;
      f_var = b_f_var;
    }

#   pragma omp parallel for if (n_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
      cs_lnum_t s_id = f2v_idx[f_id];
      cs_lnum_t e_id = f2v_idx[f_id+1];
      cs_real_t s[stride];
      for (cs_lnum_t k = 0; k < stride; k++)
        s[k] = 0;
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        for (cs_lnum_t k = 0; k < stride; k++)
          s[k] += v_var[f2v_ids[i]][k];
      }
      for (cs_lnum_t k = 0; k < stride; k++)
        f_var[f_id][k] = s[k] / (e_id-s_id);
    }

  }

  /* Vertex values are not needed after this stage */

  CS_FREE(v_var);

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           f_id++) {

        cs_lnum_t ii = i_face_cells[f_id][0];
        cs_lnum_t jj = i_face_cells[f_id][1];

        for (cs_lnum_t k = 0; k < stride; k++) {

          cs_real_t pfaci = i_f_var[f_id][k] - c_var[ii][k];
          cs_real_t pfacj = i_f_var[f_id][k] - c_var[jj][k];

          for (cs_lnum_t l = 0; l < 3; l++) {
            grad[ii][k][l] += pfaci * i_f_face_normal[f_id][l];
            grad[jj][k][l] -= pfacj * i_f_face_normal[f_id][l];
          }

        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from boundary faces */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t f_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           f_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           f_id++) {

        cs_lnum_t ii = b_face_cells[f_id];

        /*
          Remark: for the cell \f$ \celli \f$ we remove
          \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
        */

        for (cs_lnum_t k = 0; k < stride; k++) {

          cs_real_t pfac = b_f_var[f_id][k] - c_var[ii][k];

          for (cs_lnum_t l = 0; l < 3; l++)
            grad[ii][k][l] += pfac * b_f_face_normal[f_id][l];

        }

      } /* loop on faces */

    }

  }

  CS_FREE(i_f_var);
  CS_FREE(b_f_var);

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] *= dvol;
    }
  }

  /* Synchronize halos */

  _sync_strided_gradient_halo<stride>(m, CS_HALO_EXTENDED, grad);
}

/*----------------------------------------------------------------------------
 * Initialize the gradient of a tensor for gradient reconstruction.
 *
 * A non-reconstructed gradient is computed at this stage.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs_ts   <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   grad          --> gradient of pvar (dts_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_initialize_tensor_gradient(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            int                           inc,
                            const cs_field_bc_coeffs_t   *bc_coeffs_ts,
                            const cs_real_6_t   *restrict pvar,
                            cs_real_63_t        *restrict grad)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_real_6_t  *restrict coefat
    = (const cs_real_6_t *)bc_coeffs_ts->a;
  const cs_real_66_t *restrict coefbt
    = (const cs_real_66_t *)bc_coeffs_ts->b;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_vol = mq_g->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)fvq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;

  /* Computation without reconstruction */
  /*------------------------------------*/

  /* Initialization */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t i = 0; i < 6; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] = 0.0;
    }
  }

  /* Interior faces contribution */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           f_id++) {

        cs_lnum_t c_id1 = i_face_cells[f_id][0];
        cs_lnum_t c_id2 = i_face_cells[f_id][1];

        cs_real_t pond = weight[f_id];

        /*
           Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                    + (1-\alpha_\ij) \varia_\cellj\f$
                   but for the cell \f$ \celli \f$ we remove
                   \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                   and for the cell \f$ \cellj \f$ we remove
                   \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
        */
        for (cs_lnum_t i = 0; i < 6; i++) {
          cs_real_t pfaci = (1.0-pond) * (pvar[c_id2][i] - pvar[c_id1][i]);
          cs_real_t pfacj =     - pond * (pvar[c_id2][i] - pvar[c_id1][i]);
          for (cs_lnum_t j = 0; j < 3; j++) {
            grad[c_id1][i][j] += pfaci * i_f_face_normal[f_id][j];
            grad[c_id2][i][j] -= pfacj * i_f_face_normal[f_id][j];
          }
        }

      } /* End of loop on faces */

    } /* End of loop on threads */

  } /* End of loop on thread groups */

  /* Boundary face treatment */

# pragma omp parallel for
  for (int t_id = 0; t_id < n_b_threads; t_id++) {

    for (cs_lnum_t f_id = b_group_index[t_id*2];
         f_id < b_group_index[t_id*2 + 1];
         f_id++) {

      cs_lnum_t c_id = b_face_cells[f_id];

      /*
        Remark: for the cell \f$ \celli \f$ we remove
                 \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
      */
      for (cs_lnum_t i = 0; i < 6; i++) {
        cs_real_t pfac = inc*coefat[f_id][i];

        for (cs_lnum_t k = 0; k < 6; k++) {
          if (i == k)
            pfac += (coefbt[f_id][i][k] - 1.0) * pvar[c_id][k];
          else
            pfac += coefbt[f_id][i][k] * pvar[c_id][k] ;

        }

        for (cs_lnum_t j = 0; j < 3; j++)
          grad[c_id][i][j] += pfac * b_f_face_normal[f_id][j];
      }

    } /* loop on faces */

  } /* loop on threads */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * c_id] == 0)
      dvol = 1. / cell_vol[c_id];
    else
      dvol = 0.;

    for (cs_lnum_t i = 0; i < 6; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[c_id][i][j] *= dvol;
    }
  }

  /* Periodicity and parallelism treatment */

  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 18);
    if (cs_glob_mesh->have_rotation_perio)
      cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                           halo_type,
                                           (cs_real_t *)grad);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * This variant of the \ref cs_gradient_scalar function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]     var_name         variable name
 * \param[in]     gradient_info    performance logging structure, or nullptr
 * \param[in]     gradient_type    gradient type
 * \param[in]     halo_type        halo type
 * \param[in]     inc              if 0, solve on increment; 1 otherwise
 * \param[in]     check_recompute_cocg  should boundary COCG be recomputed ?
 * \param[in]     n_r_sweeps       if > 1, number of reconstruction sweeps
 * \param[in]     hyd_p_flag       flag for hydrostatic pressure
 * \param[in]     w_stride         stride for weighting coefficient
 * \param[in]     verbosity        verbosity level
 * \param[in]     clip_mode        clipping mode
 * \param[in]     epsilon          precision for iterative gradient calculation
 * \param[in]     clip_coeff       clipping coefficient
 * \param[in]     f_ext            exterior force generating
 *                                 the hydrostatic pressure
 * \param[in]     bc_coeffs        boundary condition structure of the variable
 * \param[in]     var              gradient's base variable
 * \param[in]     c_weight         weighted gradient coefficient variable,
 *                                 or nullptr
 * \param[in]     cpl              structure associated with internal coupling,
 *                                 or nullptr
 * \param[out]    grad             gradient
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_scalar(const char                    *var_name,
                 cs_gradient_info_t            *gradient_info,
                 cs_gradient_type_t             gradient_type,
                 cs_halo_type_t                 halo_type,
                 int                            inc,
                 bool                           check_recompute_cocg,
                 int                            n_r_sweeps,
                 int                            hyd_p_flag,
                 int                            w_stride,
                 int                            verbosity,
                 cs_gradient_limit_t            clip_mode,
                 double                         epsilon,
                 double                         clip_coeff,
                 const cs_real_3_t             *f_ext,
                 const cs_field_bc_coeffs_t    *bc_coeffs,
                 const cs_real_t                var[],
                 const cs_real_t                c_weight[],
                 const cs_internal_coupling_t  *cpl,
                 cs_real_t           (*restrict grad)[3])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t n_b_faces = mesh->n_b_faces;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  static int last_fvm_count = 0;
  static char var_name_prev[96] = "";

  bool recompute_cocg = true;

  if (check_recompute_cocg) {
    /* We may reuse the boundary COCG values
       if we last computed a gradient for this same field or array,
       and if we are solving in increment. */

    if (strncmp(var_name_prev, var_name, 95) == 0 && inc == 0)
      recompute_cocg = false;

    int prev_fvq_count = last_fvm_count;
    last_fvm_count = cs_mesh_quantities_compute_count();
    if (last_fvm_count != prev_fvq_count)
      recompute_cocg = true;
  }
  strncpy(var_name_prev, var_name, 95);

  /* For internal coupling, find field BC Coefficients
     matching the current variable.
     FIXME: this should also work with the iterative gradient,
     but needs extra checking. */

  /* Use Neumann BC's as default if not provided */

  cs_field_bc_coeffs_t *bc_coeffs_loc = nullptr;

  if (bc_coeffs == nullptr) {
    CS_MALLOC(bc_coeffs_loc, 1, cs_field_bc_coeffs_t);
    cs_field_bc_coeffs_init(bc_coeffs_loc);

    CS_MALLOC(bc_coeffs_loc->a, n_b_faces, cs_real_t);
    CS_MALLOC(bc_coeffs_loc->b, n_b_faces, cs_real_t);

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      bc_coeffs_loc->a[i] = 0.;
      bc_coeffs_loc->b[i] = 1.;
    }

    bc_coeffs = bc_coeffs_loc;

  }

  if (bc_coeffs != nullptr) {
    if (bc_coeffs->b == nullptr) {
      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        bc_coeffs->b[i] = 1.;
      }
    }
  }

  /* Update of local BC. coefficients for internal coupling */

  if (cpl != nullptr) {

    if (bc_coeffs_loc == nullptr) {
      cs_real_t *bc_coeff_a = bc_coeffs->a;
      cs_real_t *bc_coeff_b = bc_coeffs->b;

      CS_MALLOC(bc_coeffs_loc, 1, cs_field_bc_coeffs_t);
      cs_field_bc_coeffs_shallow_copy(bc_coeffs, bc_coeffs_loc);

      CS_MALLOC(bc_coeffs_loc->a, n_b_faces, cs_real_t);
      CS_MALLOC(bc_coeffs_loc->b, n_b_faces, cs_real_t);

      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        bc_coeffs_loc->a[i] = inc * bc_coeff_a[i];
        bc_coeffs_loc->b[i] = bc_coeff_b[i];
      }

      bc_coeffs = bc_coeffs_loc;
    }

    inc = 1;  /* Local _bc_coeff_a already multiplied by inc = 0 above for
                 uncoupled faces, and bc_coeff_a used for coupled faces. */

    cs_dispatch_context ctx;

    cs_real_t _clip_coeff = (clip_mode >= 0) ? clip_coeff : -1;
    cs_internal_coupling_update_bc_coeffs_s(ctx,
                                            bc_coeffs,
                                            cpl,
                                            halo_type,
                                            w_stride,
                                            _clip_coeff,
                                            var,
                                            c_weight);

    cpl = nullptr;  /* Coupling not needed in lower functions in this case.
                     * TODO check for reconstruction case */

  }

  /* Allocate work arrays */

  /* Compute gradient */

  switch (gradient_type) {

  case CS_GRADIENT_GREEN_ITER:
    _initialize_scalar_gradient(mesh,
                                fvq,
                                w_stride,
                                hyd_p_flag,
                                inc,
                                f_ext,
                                bc_coeffs,
                                var,
                                c_weight,
                                grad);

    _iterative_scalar_gradient(mesh,
                               fvq,
                               w_stride,
                               var_name,
                               gradient_info,
                               n_r_sweeps,
                               hyd_p_flag,
                               verbosity,
                               inc,
                               epsilon,
                               f_ext,
                               bc_coeffs,
                               var,
                               c_weight,
                               grad);
    break;

  case CS_GRADIENT_GREEN_R:
    _initialize_scalar_gradient(mesh,
                                fvq,
                                w_stride,
                                hyd_p_flag,
                                inc,
                                f_ext,
                                bc_coeffs,
                                var,
                                c_weight,
                                grad);

    _renormalize_scalar_gradient(mesh,
                                 fvq,
                                 hyd_p_flag,
                                 grad);

    break;

  case CS_GRADIENT_LSQ:
    [[fallthrough]];
  case CS_GRADIENT_GREEN_LSQ:
    {
      cs_real_3_t  *restrict r_grad;
      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        cs_alloc_mode_t amode = CS_ALLOC_HOST;
#if defined(HAVE_ACCEL)
        if (cs_get_device_id() > -1)
          amode = CS_ALLOC_DEVICE;
#endif
        CS_MALLOC_HD(r_grad, n_cells_ext, cs_real_3_t, amode);
      }
      else
        r_grad = grad;

      if (w_stride == 6 && c_weight != nullptr)
        _lsq_scalar_gradient_ani(mesh,
                                 fvq,
                                 inc,
                                 bc_coeffs,
                                 var,
                                 (const cs_real_6_t *)c_weight,
                                 r_grad);

      else if (hyd_p_flag) {
        cs_e2n_sum_t e2n_sum_type = cs_glob_e2n_sum_type;
#if defined(HAVE_ACCEL)
        if (cs_check_device_ptr(var) > CS_ALLOC_HOST)
          e2n_sum_type = CS_E2N_SUM_GATHER;
#endif
        if (e2n_sum_type == CS_E2N_SUM_GATHER)
          _lsq_scalar_gradient_hyd_p_gather
            (mesh,
             fvq,
             halo_type,
             recompute_cocg,
             inc,
             f_ext,
             bc_coeffs,
             var,
             c_weight,
             r_grad);
        else
          _lsq_scalar_gradient_hyd_p
            (mesh,
             fvq,
             halo_type,
             recompute_cocg,
             inc,
             f_ext,
             bc_coeffs,
             var,
             c_weight,
             r_grad);
      }

      else
        _lsq_scalar_gradient(mesh,
                             fvq,
                             halo_type,
                             recompute_cocg,
                             inc,
                             bc_coeffs,
                             var,
                             c_weight,
                             r_grad);

      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        _reconstruct_scalar_gradient(mesh,
                                     fvq,
                                     w_stride,
                                     hyd_p_flag,
                                     inc,
                                     (const cs_real_3_t *)f_ext,
                                     bc_coeffs,
                                     c_weight,
                                     var,
                                     r_grad,
                                     grad);

        CS_FREE_HD(r_grad);
      }
    }
    break;

    case CS_GRADIENT_GREEN_VTX:
      _fv_vtx_based_scalar_gradient(mesh,
                                    fvq,
                                    w_stride,
                                    hyd_p_flag,
                                    inc,
                                    (const cs_real_3_t *)f_ext,
                                    bc_coeffs,
                                    var,
                                    c_weight,
                                    grad);
   break;

  }

  if (bc_coeffs_loc != nullptr) {
    CS_FREE(bc_coeffs_loc->a);
    CS_FREE(bc_coeffs_loc->b);
    CS_FREE(bc_coeffs_loc);
  }

  _scalar_gradient_clipping(mesh,
                            fvq,
                            cs_glob_mesh_adjacencies,
                            halo_type,
                            clip_mode,
                            verbosity,
                            clip_coeff,
                            var_name,
                            var,
                            grad);

  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION)
    cs_bad_cells_regularisation_vector(grad, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_info   performance logging structure, or nullptr
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_v     boundary condition structure
 * \param[in]       var             gradient's base variable
 * \param[in]       val_f           face value for gradient computation
 * \param[in]       c_weight        weighted gradient coefficient variable,
 *                                  or nullptr
 * \param[out]      grad            gradient
                                    (\f$ \der{u_i}{x_j} \f$ is grad[][i][j])
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_vector(const char                     *var_name,
                 cs_gradient_info_t             *gradient_info,
                 cs_gradient_type_t              gradient_type,
                 cs_halo_type_t                  halo_type,
                 int                             inc,
                 int                             n_r_sweeps,
                 int                             verbosity,
                 int                             clip_mode,
                 double                          epsilon,
                 double                          clip_coeff,
                 const cs_field_bc_coeffs_t     *bc_coeffs_v,
                 const cs_real_3_t     *restrict var,
                 const cs_real_t                 val_f[][3],
                 const cs_real_t       *restrict c_weight,
                 cs_real_33_t          *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_adjacencies_t *madj = cs_glob_mesh_adjacencies;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  /* Compute gradient */

  switch (gradient_type) {

  case CS_GRADIENT_GREEN_ITER:

    _initialize_vector_gradient(mesh,
                                fvq,
                                halo_type,
                                inc,
                                bc_coeffs_v,
                                var,
                                c_weight,
                                grad);

    /* If reconstructions are required */

    if (n_r_sweeps > 1)
      _iterative_vector_gradient(mesh,
                                 fvq,
                                 var_name,
                                 gradient_info,
                                 halo_type,
                                 inc,
                                 n_r_sweeps,
                                 verbosity,
                                 epsilon,
                                 bc_coeffs_v,
                                 var,
                                 c_weight,
                                 grad);

    break;

  case CS_GRADIENT_LSQ:
    [[fallthrough]];
  case CS_GRADIENT_GREEN_LSQ:
    {
      cs_real_33_t *restrict r_grad;
      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        cs_alloc_mode_t amode = CS_ALLOC_HOST;
#if defined(HAVE_CUDA)
        if (cs_get_device_id() > -1)
          amode = CS_ALLOC_DEVICE;
#endif
        CS_MALLOC_HD(r_grad, n_cells_ext, cs_real_33_t, amode);
      }
      else
        r_grad = grad;

      lsq_strided_gradient_t<3> *gradient_f;
      if (cs_glob_e2n_sum_type == CS_E2N_SUM_SCATTER)
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_SCATTER>;
      else if (cs_glob_e2n_sum_type == CS_E2N_SUM_SCATTER_ATOMIC)
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_SCATTER_ATOMIC>;
      else
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_GATHER>;
      gradient_f(mesh,
                 madj,
                 fvq,
                 halo_type,
                 var,
                 val_f,
                 c_weight,
                 r_grad);

      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        _reconstruct_strided_gradient<3>(mesh,
                                         madj,
                                         fvq,
                                         halo_type,
                                         var,
                                         val_f,
                                         c_weight,
                                         r_grad,
                                         grad);

        CS_FREE_HD(r_grad);
      }
    }
    break;

  case CS_GRADIENT_GREEN_VTX:
    _fv_vtx_based_strided_gradient<3>(mesh,
                                      fvq,
                                      var,
                                      val_f,
                                      c_weight,
                                      grad);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: gradient type \"%s\" not handled."),
              __func__, cs_gradient_type_name[gradient_type]);
  }

  _strided_gradient_clipping(mesh,
                             fvq,
                             madj,
                             halo_type,
                             clip_mode,
                             verbosity,
                             clip_coeff,
                             var_name,
                             var,
                             grad);

  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION)
    cs_bad_cells_regularisation_tensor((cs_real_9_t *)grad, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_info   performance logging structure, or nullptr
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_ts    boundary condition structure
 * \param[in]       var             gradient's base variable
 * \param[in]       val_f           face value for gradient
 * \param[out]      grad            gradient
                                      (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_tensor(const char                 *var_name,
                 cs_gradient_info_t         *gradient_info,
                 cs_gradient_type_t          gradient_type,
                 cs_halo_type_t              halo_type,
                 int                         inc,
                 int                         n_r_sweeps,
                 int                         verbosity,
                 int                         clip_mode,
                 double                      epsilon,
                 double                      clip_coeff,
                 const cs_field_bc_coeffs_t *bc_coeffs_ts,
                 const cs_real_6_t          *restrict var,
                 const cs_real_6_t          *restrict val_f,
                 cs_real_63_t               *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_adjacencies_t *madj = cs_glob_mesh_adjacencies;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  /* Compute gradient */

  switch (gradient_type) {

  case CS_GRADIENT_GREEN_ITER:

    _initialize_tensor_gradient(mesh,
                                fvq,
                                halo_type,
                                inc,
                                bc_coeffs_ts,
                                var,
                                grad);

    /* If reconstructions are required */

    if (n_r_sweeps > 1)
      _iterative_tensor_gradient(mesh,
                                 fvq,
                                 var_name,
                                 gradient_info,
                                 halo_type,
                                 inc,
                                 n_r_sweeps,
                                 verbosity,
                                 epsilon,
                                 bc_coeffs_ts,
                                 var,
                                 grad);

    break;

  case CS_GRADIENT_LSQ:
    [[fallthrough]];
  case CS_GRADIENT_GREEN_LSQ:
    {
      cs_real_63_t *restrict r_grad;
      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        cs_alloc_mode_t amode = CS_ALLOC_HOST;
#if defined(HAVE_CUDA)
        if (cs_get_device_id() > -1)
          amode = CS_ALLOC_DEVICE;
#endif
        CS_MALLOC_HD(r_grad, mesh->n_cells_with_ghosts, cs_real_63_t, amode);
      }
      else
        r_grad = grad;

      lsq_strided_gradient_t<6> *gradient_f;
      if (cs_glob_e2n_sum_type == CS_E2N_SUM_SCATTER)
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_SCATTER>;
      else if (cs_glob_e2n_sum_type == CS_E2N_SUM_SCATTER_ATOMIC)
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_SCATTER_ATOMIC>;
      else
        gradient_f = _lsq_strided_gradient<CS_E2N_SUM_GATHER>;
      gradient_f(mesh,
                 madj,
                 fvq,
                 halo_type,
                 var,
                 val_f,
                 nullptr, /* c_weight */
                 r_grad);

      if (gradient_type == CS_GRADIENT_GREEN_LSQ) {
        _reconstruct_strided_gradient<6>(mesh,
                                         madj,
                                         fvq,
                                         halo_type,
                                         var,
                                         val_f,
                                         nullptr, /* c_weight */
                                         r_grad,
                                         grad);

        CS_FREE_HD(r_grad);
      }
    }
    break;

  case CS_GRADIENT_GREEN_VTX:
    _fv_vtx_based_strided_gradient<6>(mesh,
                                      fvq,
                                      var,
                                      val_f,
                                      nullptr, /* c_weight */
                                      grad);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: gradient type \"%s\" not handled."),
              __func__, cs_gradient_type_name[gradient_type]);
  }

  _strided_gradient_clipping(mesh,
                             fvq,
                             madj,
                             halo_type,
                             clip_mode,
                             verbosity,
                             clip_coeff,
                             var_name,
                             (const cs_real_6_t *)var,
                             grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a strided field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs_v     boundary condition structure, or nullptr
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  c_grad          cell gradient
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_gradient_strided_cell(const cs_mesh_t             *m,
                       const cs_mesh_quantities_t  *fvq,
                       cs_lnum_t                    c_id,
                       cs_halo_type_t               halo_type,
                       const cs_field_bc_coeffs_t  *bc_coeffs_v,
                       const cs_real_t              var[][stride],
                       const cs_real_t              c_weight[],
                       cs_real_t                    c_grad[stride][3])
{
  CS_UNUSED(m);

  using a_t = cs_real_t[stride];
  using b_t = cs_real_t[stride][stride];

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

  const cs_lnum_t *restrict cell_cells_idx = ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_e_idx = ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells = ma->cell_cells;
  const cs_lnum_t *restrict cell_cells_e = ma->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces = ma->cell_b_faces;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /* Compute covariance matrix and Right-Hand Side */

  cs_real_t cocg[6] = {0., 0., 0., 0., 0., 0.};
  cs_real_t rhs[stride][3];

  for (int i = 0; i < stride; i++) {
    for (int j = 0; j < 3; j++)
      rhs[i][j] = 0;
  }

  /* Contribution from interior and extended cells */

  int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

  for (int adj_id = 0; adj_id < n_adj; adj_id++) {

    const cs_lnum_t *restrict cell_cells_p;
    cs_lnum_t s_id, e_id;

    if (adj_id == 0) {
      s_id = cell_cells_idx[c_id];
      e_id = cell_cells_idx[c_id+1];
      cell_cells_p = (const cs_lnum_t *)(cell_cells);
    }
    else if (cell_cells_e_idx != nullptr) {
      s_id = cell_cells_e_idx[c_id];
      e_id = cell_cells_e_idx[c_id+1];
      cell_cells_p = (const cs_lnum_t *)(cell_cells_e);
    }
    else
      break;

    if (c_weight == nullptr) {

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_real_t dc[3];
        cs_lnum_t c_id1 = cell_cells_p[i];
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        cocg[0] += dc[0]*dc[0]*ddc;
        cocg[1] += dc[1]*dc[1]*ddc;
        cocg[2] += dc[2]*dc[2]*ddc;
        cocg[3] += dc[0]*dc[1]*ddc;
        cocg[4] += dc[1]*dc[2]*ddc;
        cocg[5] += dc[0]*dc[2]*ddc;

        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          cs_real_t pfac = (var[c_id1][kk] - var[c_id][kk]) * ddc;
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[kk][ll] += dc[ll] * pfac;
        }

      }

    }
    else {

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_real_t dc[3];
        cs_lnum_t c_id1 = cell_cells_p[i];
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        cs_real_t _weight =   2. * c_weight[c_id1]
                            / (c_weight[c_id] + c_weight[c_id1]);

        cocg[0] += dc[0]*dc[0]*ddc;
        cocg[1] += dc[1]*dc[1]*ddc;
        cocg[2] += dc[2]*dc[2]*ddc;
        cocg[3] += dc[0]*dc[1]*ddc;
        cocg[4] += dc[1]*dc[2]*ddc;
        cocg[5] += dc[0]*dc[2]*ddc;

        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          cs_real_t pfac = (var[c_id1][kk] - var[c_id][kk]) * ddc;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[kk][ll] += dc[ll] * pfac * _weight;
        }

      }
    }

  } /* end of contribution from interior and extended cells */

  /* Contribution from hidden boundary faces, if present */

  if (ma->cell_hb_faces_idx != nullptr) {
    cs_lnum_t s_id = ma->cell_hb_faces_idx[c_id];
    cs_lnum_t e_id = ma->cell_hb_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t f_id = ma->cell_hb_faces[i];

      const cs_nreal_t *dddij = fvq->b_face_u_normal[f_id];

      cocg[0] += dddij[0]*dddij[0];
      cocg[1] += dddij[1]*dddij[1];
      cocg[2] += dddij[2]*dddij[2];
      cocg[3] += dddij[0]*dddij[1];
      cocg[4] += dddij[1]*dddij[2];
      cocg[5] += dddij[0]*dddij[2];

    }
  }

  /* Contribution from boundary conditions */

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

  const a_t *bc_coeff_a = nullptr;
  const b_t *bc_coeff_b = nullptr;

  if (bc_coeffs_v != nullptr) {
    bc_coeff_a = (const a_t *)bc_coeffs_v->a;
    bc_coeff_b = (const b_t *)bc_coeffs_v->b;
  }

  if (e_id > s_id) {

    /* Case with known BC's */

    if (bc_coeff_a != nullptr) {

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_lnum_t f_id = cell_b_faces[i];

        cs_real_t dif[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dif[ll] = b_face_cog[f_id][ll] - cell_cen[c_id][ll];

        cs_real_t ddif = 1. / cs_math_3_square_norm(dif);

        cocg[0] += dif[0]*dif[0]*ddif;
        cocg[1] += dif[1]*dif[1]*ddif;
        cocg[2] += dif[2]*dif[2]*ddif;
        cocg[3] += dif[0]*dif[1]*ddif;
        cocg[4] += dif[1]*dif[2]*ddif;
        cocg[5] += dif[0]*dif[2]*ddif;

        cs_real_t var_f[stride];
        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          var_f[kk] = bc_coeff_a[f_id][kk];
        }

        /* Initial prediction using non-reconstructed value for Neumann */
        if (bc_coeff_b != nullptr) {
          for (cs_lnum_t kk = 0; kk < stride; kk++) {
            for (cs_lnum_t ll = 0; ll < 3; ll++) {
              var_f[kk] += bc_coeff_b[f_id][ll][kk] * var[c_id][ll];
            }
          }
        }

        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          cs_real_t pfac = (var_f[kk] - var[c_id][kk]) * ddif;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[kk][ll] += dif[ll] * pfac;
        }

      }

    }

    /* Case with no boundary conditions or known face values */

    else {

      /* Use homogeneous Neumann BC; as above, but
         pfac and RHS contribution become zero */

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_lnum_t f_id = cell_b_faces[i];

        cs_real_t dif[3];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dif[ll] = b_face_cog[f_id][ll] - cell_cen[c_id][ll];

        cs_real_t ddif = 1. / cs_math_3_square_norm(dif);

        cocg[0] += dif[0]*dif[0]*ddif;
        cocg[1] += dif[1]*dif[1]*ddif;
        cocg[2] += dif[2]*dif[2]*ddif;
        cocg[3] += dif[0]*dif[1]*ddif;
        cocg[4] += dif[1]*dif[2]*ddif;
        cocg[5] += dif[0]*dif[2]*ddif;

      }

    }

  }

  /* Invert */

  _math_6_inv_cramer_sym_in_place(cocg);

  for (cs_lnum_t i = 0; i < stride; i++) {
    c_grad[i][0] =   rhs[i][0] * cocg[0]
                   + rhs[i][1] * cocg[3]
                   + rhs[i][2] * cocg[5];

    c_grad[i][1] =   rhs[i][0] * cocg[3]
                   + rhs[i][1] * cocg[1]
                   + rhs[i][2] * cocg[4];

    c_grad[i][2] =   rhs[i][0] * cocg[5]
                   + rhs[i][1] * cocg[4]
                   + rhs[i][2] * cocg[2];
  }

  /* Correct gradient in case of Neumann BC's */

  if (e_id > s_id && bc_coeff_b != nullptr) {

    cs_real_t grad_0[stride][3], grad_i[stride][3];

    memcpy(grad_0, c_grad, sizeof(cs_real_t)*stride*3);
    memcpy(grad_i, c_grad, sizeof(cs_real_t)*stride*3);

    /* Compute norm for convergence testing. */

    cs_real_t ref_norm = 0;
    for (cs_lnum_t kk = 0; kk < stride; kk++) {
      for (cs_lnum_t ll = 0; ll < 3; ll++)
        ref_norm += std::abs(c_grad[kk][ll]);
    }

    /* Iterate over boundary condition contributions. */

    const int n_c_iter_max = 30;
    const cs_real_t c_eps = 1e-4, eps_dvg = 1e-2;

    cs_real_t c_norm = 0;

    int n_c_it;
    for (n_c_it = 0; n_c_it < n_c_iter_max; n_c_it++) {

      cs_real_t rhs_c[stride][3];

      for (cs_lnum_t ll = 0; ll < stride; ll++) {
        rhs_c[ll][0] = 0;
        rhs_c[ll][1] = 0;
        rhs_c[ll][2] = 0;
      }

      /* Loop on boundary faces */

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_lnum_t f_id = cell_b_faces[i];

        /* Remark: we could avoid recomputing dif and ddif in most cases by
           saving at least a few values from the previous step, in
           fixed-size buffers indexed by i_rel = i - s_id. */

        cs_real_t dif[3];
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dif[ii] = b_face_cog[f_id][ii] - cell_cen[c_id][ii];

        cs_real_t ddif = 1. / cs_math_3_square_norm(dif);

        /* Note that the contribution to the right-hand side from
           bc_coeff_a[c_f_id] + (bc_coeff_b -1).var[c_id] has already been
           counted, so it does not appear in the following terms.
           We should perhaps check whether the numerical sensitivity
           is lower when proceeding thus or when computing the full
           face value at each step. */

        cs_real_t var_ip_f[stride];

        for (cs_lnum_t ll = 0; ll < stride; ll++) {
          var_ip_f[ll] = cs_math_3_dot_product(c_grad[ll], diipb[f_id]);
        }

        const cs_real_t *b =   ((const cs_real_t *)bc_coeff_b)
                             + (f_id*stride*stride);

        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          cs_real_t pfac = 0;
          for (cs_lnum_t ll = 0; ll < stride; ll++) {
            pfac += b[kk*stride + ll] * var_ip_f[ll] * ddif;
          }

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs_c[kk][ll] += dif[ll] * pfac;
        }

      } /* End of loop on boundary faces */

      /* Compute gradient correction */

      cs_real_t grad_c[stride][3];

      for (cs_lnum_t ii = 0; ii < stride; ii++) {

        grad_c[ii][0] = (  cocg[0] * rhs_c[ii][0]
                         + cocg[3] * rhs_c[ii][1]
                         + cocg[5] * rhs_c[ii][2]);
        grad_c[ii][1] = (  cocg[3] * rhs_c[ii][0]
                         + cocg[1] * rhs_c[ii][1]
                         + cocg[4] * rhs_c[ii][2]);
        grad_c[ii][2] = (  cocg[5] * rhs_c[ii][0]
                         + cocg[4] * rhs_c[ii][1]
                         + cocg[2] * rhs_c[ii][2]);

      }

      /* Update gradient values */

      c_norm = 0;
      for (cs_lnum_t ii = 0; ii < stride; ii++) {
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          c_grad[ii][jj] = grad_0[ii][jj] + grad_c[ii][jj];
          c_norm += std::abs(c_grad[ii][jj] - grad_i[ii][jj]);
          grad_i[ii][jj] = c_grad[ii][jj];
        }
      }

#if 0
      printf("grads %d: it %d, ref_norm %g, it_norm %g\n",
             c_id, n_c_it, ref_norm, c_norm);
#endif

      /* Use of Frobenius norm is not rotation-independent,
         but is cheaper to compute and assumed "good enough".
         Note that the comparison to cs_math_epzero is done for
         the gradient correction itself, so is already
         independent of the cell size. */
      if (c_norm < ref_norm * c_eps || c_norm < cs_math_epzero)
        break;

    } /* End of loop on iterations */

    /* If the last correction was too large, we suspect
       the the algorithm did not converge at all/diverged,
       so we simply restore the non-reconstructed value
       (additional precaution, not encountered in testing). */

    if (c_norm > eps_dvg * ref_norm) {
      memcpy(c_grad, grad_0, sizeof(cs_real_t)*stride*3);
#if 0
      printf("%s: non-convergence for boundary cell %ld/%ld\n"
             "  use non-recontruced value\n", __func__,
             (long)c_idx, (long)c_id);
#endif
      n_c_it *= -1;
    }

  }  /* End of correction in case of Neumann BC's */
}

BEGIN_C_DECLS

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize gradient computation API.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_initialize(void)
{
  assert(cs_glob_mesh != nullptr);

  CS_TIMER_COUNTER_INIT(_gradient_t_tot);

  int stats_root = cs_timer_stats_id_by_name("operations");
  if (stats_root > -1) {
    _gradient_stat_id = cs_timer_stats_create("operations",
                                              "gradients",
                                              "gradients reconstruction");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize gradient computation API.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_finalize(void)
{
  _gradient_quantities_destroy();

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Total elapsed time for all gradient computations:  %.3f s\n"),
                _gradient_t_tot.nsec*1e-9);

  /* Free system info */

  for (int ii = 0; ii < _gradient_n_systems; ii++) {
    _gradient_info_dump(_gradient_systems[ii]);
    _gradient_info_destroy(&(_gradient_systems[ii]));
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  CS_FREE(_gradient_systems);

  _gradient_n_systems = 0;
  _gradient_n_max_systems = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free saved gradient quantities.
 *
 * This is required when the mesh changes, so that the on-demand computation
 * will be updated.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_free_quantities(void)
{
  for (int i = 0; i < _n_gradient_quantities; i++) {

    cs_gradient_quantities_t  *gq = _gradient_quantities + i;

    CS_FREE(gq->cocg_it);
    CS_FREE(gq->cocgb_s_lsq);
    CS_FREE(gq->cocg_lsq);
    CS_FREE(gq->cocgb_s_lsq_ext);
    CS_FREE(gq->cocg_lsq_ext);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * \param[in]       var_name       variable name
 * \param[in]       gradient_type  gradient type
 * \param[in]       halo_type      halo type
 * \param[in]       inc            if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps     if > 1, number of reconstruction sweeps
 *                                 (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       hyd_p_flag     flag for hydrostatic pressure
 * \param[in]       w_stride       stride for weighting coefficient
 * \param[in]       verbosity      verbosity level
 * \param[in]       clip_mode      clipping mode
 * \param[in]       epsilon        precision for iterative gradient calculation
 * \param[in]       clip_coeff     clipping coefficient
 * \param[in]       f_ext          exterior force generating the
 *                                 hydrostatic pressure
 * \param[in]       bc_coeffs      boundary condition structure
 * \param[in, out]  var            gradient's base variable
 * \param[in, out]  c_weight       cell variable weight, or nullptr
 * \param[in]       cpl            associated internal coupling, or nullptr
 * \param[out]      grad           gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar(const char                    *var_name,
                   cs_gradient_type_t             gradient_type,
                   cs_halo_type_t                 halo_type,
                   int                            inc,
                   int                            n_r_sweeps,
                   int                            hyd_p_flag,
                   int                            w_stride,
                   int                            verbosity,
                   cs_gradient_limit_t            clip_mode,
                   double                         epsilon,
                   double                         clip_coeff,
                   cs_real_3_t                    f_ext[],
                   const cs_field_bc_coeffs_t    *bc_coeffs,
                   cs_real_t                      var[],
                   cs_real_t                     *c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t           (*restrict grad)[3])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  /* Synchronize variable */

  if (mesh->halo != nullptr) {

    bool on_device = cs_mem_is_device_ptr(var);

    cs_halo_sync(mesh->halo, halo_type, on_device, var);

    if (c_weight != nullptr) {
      if (w_stride == 6) {
        cs_real_6_t *c_weight_t = (cs_real_6_t *)c_weight;
        cs_halo_sync_r(mesh->halo, halo_type, on_device, c_weight_t);
      }
      else
        cs_halo_sync(mesh->halo, halo_type, on_device, c_weight);
    }

    if (hyd_p_flag == 1) {
      cs_halo_sync_r(mesh->halo, halo_type, on_device, f_ext);
    }

  }

  _gradient_scalar(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   false, /* Do not use previous cocg at boundary, recompute */
                   n_r_sweeps,
                   hyd_p_flag,
                   w_stride,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   f_ext,
                   bc_coeffs,
                   var,
                   c_weight,
                   cpl,
                   grad);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_v     boundary condition structure
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        cell variable weight, or nullptr
 * \param[in]       cpl             associated internal coupling, or nullptr
 * \param[out]      gradv           gradient
                                    (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector(const char                    *var_name,
                   cs_gradient_type_t             gradient_type,
                   cs_halo_type_t                 halo_type,
                   int                            inc,
                   int                            n_r_sweeps,
                   int                            verbosity,
                   cs_gradient_limit_t            clip_mode,
                   double                         epsilon,
                   double                         clip_coeff,
                   const cs_field_bc_coeffs_t    *bc_coeffs_v,
                   cs_real_t                      var[][3],
                   cs_real_t        *restrict     c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t                      gradv[][3][3])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true) {
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (mesh->halo != nullptr) {
    bool on_device = cs_mem_is_device_ptr(var);

    cs_halo_sync_r(mesh->halo, halo_type, on_device, var);
    if (c_weight != nullptr)
      cs_halo_sync(mesh->halo, halo_type, on_device, c_weight);
  }

  /* Compute face value for gradient
     ------------------------------- */

  cs_real_3_t *val_ip = nullptr, *val_f = nullptr;
  cs_real_3_t *val_f_hmg = nullptr, *val_f_wrk = nullptr;

  /* For internal coupling, find field BC Coefficients
     matching the current variable.
     FIXME: this should also work with the iterative gradient,
     but needs extra checking. */

  /* Update of local BC. coefficients for internal coupling */

  cs_real_3_t *bc_coeff_loc_cpl_a = nullptr;

  if (cpl != nullptr) {

    if (bc_coeffs_v != nullptr) {

      cs_real_3_t *bc_coeff_a = (cs_real_3_t *)bc_coeffs_v->a;
      cs_field_bc_coeffs_t bc_coeffs_loc_cpl;
      cs_field_bc_coeffs_shallow_copy(bc_coeffs_v, &bc_coeffs_loc_cpl);
      CS_MALLOC_HD(bc_coeffs_loc_cpl.a, 3*n_b_faces, cs_real_t, cs_alloc_mode);

      bc_coeff_loc_cpl_a = (cs_real_3_t *)bc_coeffs_loc_cpl.a;

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          bc_coeff_loc_cpl_a[face_id][i] = inc * bc_coeff_a[face_id][i];
        }
      }

      bc_coeffs_v = &bc_coeffs_loc_cpl;

      inc = 1;  /* Local _bc_coeff_a already multiplied by inc = 0 above for
                   uncoupled faces, and bc_coeff_a used for coupled faces. */

      cs_dispatch_context ctx;

      cs_real_t _clip_coeff = (clip_mode >= 0) ? clip_coeff : -1;
      cs_internal_coupling_update_bc_coeffs_strided<3>(ctx,
                                                       bc_coeffs_v,
                                                       cpl,
                                                       halo_type,
                                                       _clip_coeff,
                                                       nullptr,
                                                       var,
                                                       c_weight);

      cpl = nullptr;  /* Coupling not needed in lower functions in this case. */
    }
  }

  cs_real_3_t  *bc_coeff_loc_a = nullptr;
  cs_real_33_t *bc_coeff_loc_b = nullptr;

  // Use Neumann BC's as default if not provided
  if (bc_coeffs_v == nullptr) {

    cs_field_bc_coeffs_t bc_coeffs_v_loc;
    cs_field_bc_coeffs_init(&bc_coeffs_v_loc);

    CS_MALLOC(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t);
    CS_MALLOC(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t);

    bc_coeff_loc_a = (cs_real_3_t  *)bc_coeffs_v_loc.a;
    bc_coeff_loc_b = (cs_real_33_t *)bc_coeffs_v_loc.b;

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        bc_coeff_loc_a[i][j] = 0;
        for (cs_lnum_t k = 0; k < 3; k++)
          bc_coeff_loc_b[i][j][k] = 0;

        bc_coeff_loc_b[i][j][j] = 1;
      }
    }

    bc_coeffs_v = &bc_coeffs_v_loc;

    if (gradient_type != CS_GRADIENT_GREEN_ITER) {
      // else only above standard coefa&b are used

      CS_MALLOC_HD(val_f_hmg, n_b_faces, cs_real_3_t, cs_alloc_mode);

      cs_dispatch_context ctx;

      /* Compute var_iprime (val_f = var_iprime for hmg Neumann) */
      cs_gradient_boundary_iprime_lsq_strided<3>(ctx,
                                                 mesh,
                                                 fvq,
                                                 n_b_faces,
                                                 nullptr,
                                                 halo_type,
                                                 -1,
                                                 nullptr,
                                                 bc_coeffs_v,
                                                 c_weight,
                                                 var,
                                                 val_f_hmg,
                                                 nullptr);

      val_f = val_f_hmg;
    }
  }
  else { // bc_coeffs_v != nullptr, enter here for cpl

    if (bc_coeffs_v->val_f != nullptr) {
      val_f = (cs_real_3_t *)bc_coeffs_v->val_f;
    }
    else { // work array (momemtum, ...)

      /* Compute face value for gradient exept for iterative_gradient
         which compute val_f with iterative process via coeffa&b
         (copied by cs_field_bc_coeffs_shallow_copy)

         Remark: We cannot copy val_f in cs_field_bc_coeffs_shallow_copy.
                 For Momentum for example, coefa is multiplied by rho so
                 val_f must be not copied but updated with the new coefa.
                 It is what we do bellow.
                 TODO : give val_f in argument of this function */

      if (gradient_type != CS_GRADIENT_GREEN_ITER) {
        // else standard coefa&b are used

        CS_MALLOC_HD(val_ip, n_b_faces, cs_real_3_t, cs_alloc_mode);
        CS_MALLOC_HD(val_f_wrk, n_b_faces, cs_real_3_t, cs_alloc_mode);

        cs_dispatch_context ctx;

        cs_gradient_boundary_iprime_lsq_strided<3>(ctx,
                                                   mesh,
                                                   fvq,
                                                   n_b_faces,
                                                   nullptr,
                                                   halo_type,
                                                   -1,
                                                   nullptr,
                                                   bc_coeffs_v,
                                                   c_weight,
                                                   var,
                                                   val_ip,
                                                   nullptr);

        cs_real_3_t  *coefa = (cs_real_3_t  *)bc_coeffs_v->a;
        cs_real_33_t *coefb = (cs_real_33_t *)bc_coeffs_v->b;

        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

          for (cs_lnum_t i = 0; i < 3; i++) {
            val_f_wrk[face_id][i] = inc*coefa[face_id][i];
            for (cs_lnum_t j = 0; j < 3; j++) {
              val_f_wrk[face_id][i] += coefb[face_id][j][i]*val_ip[face_id][j];
            }
          }
        } /* End loop on boundary faces */

        val_f = val_f_wrk;
      }
    }
  }

  /* Compute gradient */

  _gradient_vector(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   n_r_sweeps,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   bc_coeffs_v,
                   (const cs_real_3_t *)var,
                   (const cs_real_3_t *)val_f,
                   (const cs_real_t *)c_weight,
                   gradv);

  CS_FREE_HD(val_ip);
  CS_FREE_HD(val_f_hmg);
  CS_FREE_HD(val_f_wrk);
  CS_FREE_HD(bc_coeff_loc_a);
  CS_FREE_HD(bc_coeff_loc_b);
  CS_FREE_HD(bc_coeff_loc_cpl_a);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_ts    boundary condition structure
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor(const char                  *var_name,
                   cs_gradient_type_t           gradient_type,
                   cs_halo_type_t               halo_type,
                   int                          inc,
                   int                          n_r_sweeps,
                   int                          verbosity,
                   cs_gradient_limit_t          clip_mode,
                   double                       epsilon,
                   double                       clip_coeff,
                   const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                   cs_real_6_t        *restrict var,
                   cs_real_63_t       *restrict grad)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true) {
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (mesh->halo != nullptr) {
    bool on_device = cs_mem_is_device_ptr(var);
    cs_halo_sync_r(mesh->halo, halo_type, on_device, var);
  }

  /* Use Neumann BC's as default if not provided */

  cs_real_6_t  *bc_coeff_loc_a = nullptr;
  cs_real_66_t *bc_coeff_loc_b = nullptr;
  cs_real_6_t *val_f = nullptr, *val_ip = nullptr;
  cs_real_6_t *val_f_hmg = nullptr, *val_f_wrk = nullptr;

  if (bc_coeffs_ts == nullptr) {

    cs_field_bc_coeffs_t bc_coeffs_ts_loc;
    cs_field_bc_coeffs_init(&bc_coeffs_ts_loc);

    CS_MALLOC(bc_coeffs_ts_loc.a, 6*n_b_faces, cs_real_t);
    CS_MALLOC(bc_coeffs_ts_loc.b, 36*n_b_faces, cs_real_t);

    bc_coeff_loc_a = (cs_real_6_t  *)bc_coeffs_ts_loc.a;
    bc_coeff_loc_b = (cs_real_66_t *)bc_coeffs_ts_loc.b;

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 6; j++) {
        bc_coeff_loc_a[i][j] = 0;
        for (cs_lnum_t k = 0; k < 6; k++)
          bc_coeff_loc_b[i][j][k] = 0;

        bc_coeff_loc_b[i][j][j] = 1;
      }
    }

    bc_coeffs_ts = &bc_coeffs_ts_loc;

    if (gradient_type != CS_GRADIENT_GREEN_ITER) {
      // else only above standard coefa&b are used

      CS_MALLOC_HD(val_f_hmg, n_b_faces, cs_real_6_t, cs_alloc_mode);

      cs_dispatch_context ctx;

      /* Compute var_iprime (val_f = var_iprime for hmg Neumann) */
      cs_gradient_boundary_iprime_lsq_strided<6>(ctx,
                                                 mesh,
                                                 fvq,
                                                 n_b_faces,
                                                 nullptr,
                                                 halo_type,
                                                 -1,
                                                 nullptr,
                                                 bc_coeffs_ts,
                                                 nullptr, // c_weight,
                                                 var,
                                                 val_f_hmg,
                                                 nullptr);

      val_f = val_f_hmg;
    }
  }
  else {

    if (bc_coeffs_ts->val_f != nullptr)
      val_f = (cs_real_6_t *)bc_coeffs_ts->val_f;
    else { // work array

    /* Compute face value for gradient exept for iterative_gradient
       which compute val_f with iterative process */

      if (gradient_type != CS_GRADIENT_GREEN_ITER) {
        // else standard coefa&b are used

        CS_MALLOC_HD(val_ip, n_b_faces, cs_real_6_t, cs_alloc_mode);
        CS_MALLOC_HD(val_f_wrk, n_b_faces, cs_real_6_t, cs_alloc_mode);

        cs_real_6_t  *coefav = (cs_real_6_t  *)bc_coeffs_ts->a;
        cs_real_66_t *coefbv = (cs_real_66_t *)bc_coeffs_ts->b;

        cs_dispatch_context ctx;

        /* Compute var_iprime */
        cs_gradient_boundary_iprime_lsq_strided<6>(ctx,
                                                   mesh,
                                                   fvq,
                                                   n_b_faces,
                                                   nullptr,
                                                   halo_type,
                                                   -1,
                                                   nullptr,
                                                   bc_coeffs_ts,
                                                   nullptr, // c_weight
                                                   var,
                                                   val_ip,
                                                   nullptr);

        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
          for (cs_lnum_t i = 0; i < 6; i++) {
            val_f_wrk[face_id][i] = inc*coefav[face_id][i];
            for (cs_lnum_t j = 0; j < 6; j++) {
              val_f_wrk[face_id][i] += coefbv[face_id][j][i]*val_ip[face_id][j];
            }
          }
        }

        val_f = val_f_wrk;

      }
    }
  }

  /* Compute gradient */

  _gradient_tensor(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   n_r_sweeps,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   bc_coeffs_ts,
                   (const cs_real_6_t *)var,
                   (const cs_real_6_t *)val_f,
                   grad);

  CS_FREE_HD(val_ip);
  CS_FREE_HD(val_f_hmg);
  CS_FREE_HD(val_f_wrk);
  CS_FREE_HD(bc_coeff_loc_a);
  CS_FREE_HD(bc_coeff_loc_b);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * This variant of the \ref cs_gradient_scalar function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]   var_name        variable name
 * \param[in]   gradient_type   gradient type
 * \param[in]   halo_type       halo type
 * \param[in]   inc             if 0, solve on increment; 1 otherwise
 * \param[in]   n_r_sweeps      if > 1, number of reconstruction sweeps
 *                              (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]   hyd_p_flag      flag for hydrostatic pressure
 * \param[in]   w_stride        stride for weighting coefficient
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   f_ext           exterior force generating the
 *                              hydrostatic pressure
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[in]   cpl             associated internal coupling, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_synced_input(const char                 *var_name,
                                cs_gradient_type_t          gradient_type,
                                cs_halo_type_t              halo_type,
                                int                         inc,
                                int                         n_r_sweeps,
                                int                         hyd_p_flag,
                                int                         w_stride,
                                int                         verbosity,
                                cs_gradient_limit_t         clip_mode,
                                double                      epsilon,
                                double                      clip_coeff,
                                cs_real_t                   f_ext[][3],
                                const cs_field_bc_coeffs_t *bc_coeffs,
                                const cs_real_t             var[],
                                const cs_real_t             c_weight[],
                                const cs_internal_coupling_t  *cpl,
                                cs_real_t                   grad[][3])
{
  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  if (hyd_p_flag == 2)
    hyd_p_flag = 0;

  if (hyd_p_flag == 1) {
    bool on_device = cs_mem_is_device_ptr(f_ext);
    const cs_halo_t *halo = cs_glob_mesh->halo;
    cs_halo_sync_r(halo, halo_type, on_device, f_ext);
  }

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  _gradient_scalar(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   true,  /* Check recompute of cocg */
                   n_r_sweeps,
                   hyd_p_flag,
                   w_stride,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   f_ext,
                   bc_coeffs,
                   var,
                   c_weight,
                   cpl,
                   grad);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * This variant of the \ref cs_gradient_vector function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]   var_name        variable name
 * \param[in]   gradient_type   gradient type
 * \param[in]   halo_type       halo type
 * \param[in]   inc             if 0, solve on increment; 1 otherwise
 * \param[in]   n_r_sweeps      if > 1, number of reconstruction sweeps
 *                              (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   bc_coeffs_v     boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[in]   cpl             associated internal coupling, or nullptr
 * \param[out]  grad            gradient
                                (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_synced_input(const char                 *var_name,
                                cs_gradient_type_t          gradient_type,
                                cs_halo_type_t              halo_type,
                                int                         inc,
                                int                         n_r_sweeps,
                                int                         verbosity,
                                cs_gradient_limit_t         clip_mode,
                                double                      epsilon,
                                double                      clip_coeff,
                                const cs_field_bc_coeffs_t *bc_coeffs_v,
                                const cs_real_t             var[][3],
                                const cs_real_t             val_f[][3],
                                const cs_real_t             c_weight[],
                                cs_real_t                   grad[][3][3])
{
  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  /* Compute gradient */

  _gradient_vector(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   n_r_sweeps,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   bc_coeffs_v,
                   var,
                   val_f,
                   c_weight,
                   grad);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * This variant of the \ref cs_gradient_tensor function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_ts    boundary condition structure
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_synced_input(const char                  *var_name,
                                cs_gradient_type_t           gradient_type,
                                cs_halo_type_t               halo_type,
                                int                          inc,
                                int                          n_r_sweeps,
                                int                          verbosity,
                                cs_gradient_limit_t          clip_mode,
                                double                       epsilon,
                                double                       clip_coeff,
                                const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                                const cs_real_t              var[][6],
                                const cs_real_t              val_f[][6],
                                cs_real_63_t                *grad)
{
  cs_gradient_info_t *gradient_info = nullptr;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  /* Compute gradient */

  _gradient_tensor(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   n_r_sweeps,
                   verbosity,
                   clip_mode,
                   epsilon,
                   clip_coeff,
                   bc_coeffs_ts,
                   var,
                   val_f,
                   grad);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&_gradient_t_tot, &t0, &t1);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  if (_gradient_stat_id > -1)
    cs_timer_stats_add_diff(_gradient_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a scalar field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs,
                        const cs_real_t              var[],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3])
{
  CS_UNUSED(m);

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict cell_cells_idx = ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_e_idx = ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells = ma->cell_cells;
  const cs_lnum_t *restrict cell_cells_e = ma->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces = ma->cell_b_faces;

  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;

  /* Reconstruct gradients using least squares for non-orthogonal meshes */

  cs_real_t cocg[6] = {0., 0., 0., 0., 0., 0.};
  cs_real_t rhsv[3] = {0., 0., 0.};

  int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

  for (int adj_id = 0; adj_id < n_adj; adj_id++) {

    const cs_lnum_t *restrict cell_cells_p;
    cs_lnum_t s_id, e_id;

    if (adj_id == 0) {
      s_id = cell_cells_idx[c_id];
      e_id = cell_cells_idx[c_id+1];
      cell_cells_p = cell_cells;
    }
    else if (cell_cells_e_idx != nullptr) {
      s_id = cell_cells_e_idx[c_id];
      e_id = cell_cells_e_idx[c_id+1];
      cell_cells_p = cell_cells_e;
    }
    else
      break;

    if (c_weight == nullptr) {

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_real_t dc[3];
        cs_lnum_t c_id1 = cell_cells_p[i];
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        cs_real_t pfac = (var[c_id1] - var[c_id]) * ddc;

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ll] += dc[ll] * pfac;

        cocg[0] += dc[0]*dc[0]*ddc;
        cocg[1] += dc[1]*dc[1]*ddc;
        cocg[2] += dc[2]*dc[2]*ddc;
        cocg[3] += dc[0]*dc[1]*ddc;
        cocg[4] += dc[1]*dc[2]*ddc;
        cocg[5] += dc[0]*dc[2]*ddc;

      }

    }
    else {

      for (cs_lnum_t i = s_id; i < e_id; i++) {

        cs_real_t dc[3];
        cs_lnum_t c_id1 = cell_cells_p[i];
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        cs_real_t pfac = (var[c_id1] - var[c_id]) * ddc;

        cs_real_t _weight =   2. * c_weight[c_id1]
                            / (c_weight[c_id] + c_weight[c_id1]);

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhsv[ll] += dc[ll] * pfac * _weight;

        cocg[0] += dc[0]*dc[0]*ddc;
        cocg[1] += dc[1]*dc[1]*ddc;
        cocg[2] += dc[2]*dc[2]*ddc;
        cocg[3] += dc[0]*dc[1]*ddc;
        cocg[4] += dc[1]*dc[2]*ddc;
        cocg[5] += dc[0]*dc[2]*ddc;

      }
    }

  } /* end of contribution from interior and extended cells */

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

  /* Contribution from hidden boundary faces, if present */

  if (ma->cell_hb_faces_idx != nullptr)
    _add_hb_faces_cocg_lsq_cell(c_id,
                                ma->cell_hb_faces_idx,
                                ma->cell_hb_faces,
                                fvq->b_face_u_normal,
                                cocg);

  /* Contribution from boundary conditions */

  const cs_real_t *bc_coeff_a = nullptr;
  const cs_real_t *bc_coeff_b = nullptr;

  if (bc_coeffs != nullptr) {
    bc_coeff_a = bc_coeffs->a;
    bc_coeff_b = bc_coeffs->b;
  }

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
    const cs_real_t *restrict b_dist = fvq->b_dist;
    const cs_rreal_3_t *restrict diipb = fvq->diipb;

    cs_real_t  dsij[3];

    cs_lnum_t f_id = cell_b_faces[i];

    for (cs_lnum_t ll = 0; ll < 3; ll++)
      dsij[ll] = b_face_u_normal[f_id][ll];

    if (bc_coeff_a != nullptr && bc_coeff_b != nullptr) { /* Known face BC's */

      cs_real_t unddij = 1. / b_dist[f_id];
      cs_real_t umcbdd = (1. -bc_coeff_b[f_id]) * unddij;

      cs_real_t pfac =  (  bc_coeff_a[f_id] + (bc_coeff_b[f_id] -1.)
                         * var[c_id]) * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        dsij[ll] += umcbdd*diipb[f_id][ll];

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ll] += dsij[ll] * pfac;

      cocg[0] += dsij[0] * dsij[0];
      cocg[1] += dsij[1] * dsij[1];
      cocg[2] += dsij[2] * dsij[2];
      cocg[3] += dsij[0] * dsij[1];
      cocg[4] += dsij[1] * dsij[2];
      cocg[5] += dsij[0] * dsij[2];

    }
    else if (bc_coeff_a != nullptr) { /* Known face values */

      const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

      cs_real_t dc[3];
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        dc[ii] = b_face_cog[f_id][ii] - cell_cen[c_id][ii];

      cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

      cs_real_t pfac = (bc_coeff_a[f_id] - var[c_id]) * ddc;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhsv[ll] += dc[ll] * pfac;

      cocg[0] += dc[0]*dc[0]*ddc;
      cocg[1] += dc[1]*dc[1]*ddc;
      cocg[2] += dc[2]*dc[2]*ddc;
      cocg[3] += dc[0]*dc[1]*ddc;
      cocg[4] += dc[1]*dc[2]*ddc;
      cocg[5] += dc[0]*dc[2]*ddc;

    }
    else { /* Assign cell values as face values (homogeneous Neumann);
              as above, pfac cancels out, so does contribution to RHS */

      const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

      cs_real_t dc[3];
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        dc[ii] = b_face_cog[f_id][ii] - cell_cen[c_id][ii];

      cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

      cocg[0] += dc[0]*dc[0]*ddc;
      cocg[1] += dc[1]*dc[1]*ddc;
      cocg[2] += dc[2]*dc[2]*ddc;
      cocg[3] += dc[0]*dc[1]*ddc;
      cocg[4] += dc[1]*dc[2]*ddc;
      cocg[5] += dc[0]*dc[2]*ddc;

    }

  } // end of contribution from boundary cells

  /* Invert */

  cs_real_t a00 = cocg[1]*cocg[2] - cocg[4]*cocg[4];
  cs_real_t a01 = cocg[4]*cocg[5] - cocg[3]*cocg[2];
  cs_real_t a02 = cocg[3]*cocg[4] - cocg[1]*cocg[5];
  cs_real_t a11 = cocg[0]*cocg[2] - cocg[5]*cocg[5];
  cs_real_t a12 = cocg[3]*cocg[5] - cocg[0]*cocg[4];
  cs_real_t a22 = cocg[0]*cocg[1] - cocg[3]*cocg[3];

  cs_real_t det_inv = 1. / (cocg[0]*a00 + cocg[3]*a01 + cocg[5]*a02);

  grad[0] = (  a00 * rhsv[0]
             + a01 * rhsv[1]
             + a02 * rhsv[2]) * det_inv;
  grad[1] = (  a01 * rhsv[0]
             + a11 * rhsv[1]
             + a12 * rhsv[2]) * det_inv;
  grad[2] = (  a02 * rhsv[0]
             + a12 * rhsv[1]
             + a22 * rhsv[2]) * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a vector field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs_v,
                        const cs_real_t              var[][3],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3][3])
{
  _gradient_strided_cell<3>(m,
                            fvq,
                            c_id,
                            halo_type,
                            bc_coeffs_v,
                            var,
                            c_weight,
                            grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a tensor field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs_ts    boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                        const cs_real_t              var[][6],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[6][3])
{
  _gradient_strided_cell<6>(m,
                            fvq,
                            c_id,
                            halo_type,
                            bc_coeffs_ts,
                            var,
                            c_weight,
                            grad);
}

/*----------------------------------------------------------------------------
 * Determine gradient type by Fortran "imrgra" value
 *
 * parameters:
 *   imrgra         <-- Fortran gradient option
 *   gradient_type  --> gradient type
 *   halo_type      --> halo type
 *----------------------------------------------------------------------------*/

void
cs_gradient_type_by_imrgra(int                  imrgra,
                           cs_gradient_type_t  *gradient_type,
                           cs_halo_type_t      *halo_type)
{
  *halo_type = CS_HALO_STANDARD;
  *gradient_type = CS_GRADIENT_GREEN_ITER;

  switch (imrgra) {
  case 0:
    *gradient_type = CS_GRADIENT_GREEN_ITER;
    break;
  case 1:
    *gradient_type = CS_GRADIENT_LSQ;
    break;
  case 2:
  case 3:
    *gradient_type = CS_GRADIENT_LSQ;
    *halo_type = CS_HALO_EXTENDED;
    break;
  case 4:
    *gradient_type = CS_GRADIENT_GREEN_LSQ;
    break;
  case 5:
  case 6:
    *gradient_type = CS_GRADIENT_GREEN_LSQ;
    *halo_type = CS_HALO_EXTENDED;
    break;
  case 7:
    *gradient_type = CS_GRADIENT_GREEN_VTX;
    break;
  case 8:
    *gradient_type = CS_GRADIENT_GREEN_R;
    break;
  default:
    *gradient_type = CS_GRADIENT_GREEN_ITER;
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute the steady balance due to porous modelling for the pressure
 *        gradient.
 *
 * \param[in]  inc  if 0, solve on increment; 1 otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_porosity_balance(int inc)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;
  const cs_halo_t  *halo = m->halo;

  const cs_real_t *restrict cell_vol = mq->cell_vol;
  cs_real_2_t *i_f_face_factor = mq->i_f_face_factor;
  cs_real_t *b_f_face_factor = mq->b_f_face_factor;
  cs_real_t *i_massflux = cs_field_by_name("inner_mass_flux")->val;
  cs_real_t *b_massflux = cs_field_by_name("boundary_mass_flux")->val;
  const cs_nreal_3_t *restrict i_face_u_normal = mq_g->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = mq_g->b_face_u_normal;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const int *restrict c_disable_flag = mq->c_disable_flag;
  cs_lnum_t has_dc = mq->has_disable_flag; /* Has cells disabled? */

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  if (f_i_poro_duq_0 == nullptr)
    return;

  cs_real_t *i_poro_duq_0 = f_i_poro_duq_0->val;
  cs_real_t *i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
  cs_real_t *b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  cs_real_3_t *c_poro_div_duq
    = (cs_real_3_t *)cs_field_by_name("poro_div_duq")->val;

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      c_poro_div_duq[c_id][i] = 0.;
  }

  if (inc == 1) {

    /* Inner faces corrections */
    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          const cs_nreal_t *normal = i_face_u_normal[f_id];

          cs_real_t *vel_i = &(CS_F_(vel)->val_pre[3*ii]);
          cs_real_t *vel_j = &(CS_F_(vel)->val_pre[3*jj]);

          cs_real_t veli_dot_n =    (1. - i_f_face_factor[f_id][0])
                                  * cs_math_3_dot_product(vel_i, normal);
          cs_real_t velj_dot_n =    (1. - i_f_face_factor[f_id][1])
                                  * cs_math_3_dot_product(vel_j, normal);

          cs_real_t d_f_surf = 0.;
          /* Is the cell disabled (for solid or porous)?
             Not the case if coupled */
          if (   has_dc * c_disable_flag[has_dc * ii] == 0
              && has_dc * c_disable_flag[has_dc * jj] == 0)
            d_f_surf = 1.;

          i_poro_duq_0[f_id] = veli_dot_n * i_massflux[f_id] * d_f_surf;
          i_poro_duq_1[f_id] = velj_dot_n * i_massflux[f_id] * d_f_surf;

          for (cs_lnum_t i = 0; i < 3; i++) {
            c_poro_div_duq[ii][i] +=   i_poro_duq_0[f_id]
                                     * i_face_u_normal[f_id][i];
            c_poro_div_duq[jj][i] -=   i_poro_duq_1[f_id]
                                     * i_face_u_normal[f_id][i];
          }
        }
      }

    }

    /* Boundary faces corrections */

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t f_id = b_group_index[t_id*2];
           f_id < b_group_index[t_id*2 + 1];
           f_id++) {

        cs_lnum_t ii = b_face_cells[f_id];

        const cs_nreal_t *normal = b_face_u_normal[f_id];

        cs_real_t *vel_i = &(CS_F_(vel)->val_pre[3*ii]);

        cs_real_t veli_dot_n =   (1. - b_f_face_factor[f_id])
                               * cs_math_3_dot_product(vel_i, normal);

        cs_real_t d_f_surf = 0.;
        /* Is the cell disabled (for solid or porous)?
           Not the case if coupled */
        if (has_dc * c_disable_flag[has_dc * ii] == 0)
          d_f_surf = 1.;

        b_poro_duq[f_id] = veli_dot_n * b_massflux[f_id] * d_f_surf;

        for (cs_lnum_t i = 0; i < 3; i++)
          c_poro_div_duq[ii][i] +=   b_poro_duq[f_id]
                                   * b_face_u_normal[f_id][i];
      }

      /* Finalization of cell terms */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        /* Is the cell disabled (for solid or porous)?
           Not the case if coupled */
        cs_real_t dvol = 0.;
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_vol[c_id];

        for (cs_lnum_t i = 0; i < 3; i++)
          c_poro_div_duq[c_id][i] *= dvol;
      }
    }

    /* Handle parallelism and periodicity */
    if (halo != nullptr)
      cs_halo_sync_var_strided(halo,
                               CS_HALO_STANDARD,
                               (cs_real_t *)c_poro_div_duq,
                               3);

  }
  else {
#   pragma omp parallel for
    for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {
      i_poro_duq_0[f_id] = 0.;
      i_poro_duq_1[f_id] = 0.;
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
