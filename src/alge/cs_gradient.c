/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_internal_coupling.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_internal_coupling.h"
#include "cs_bad_cells_regularisation.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_gradient.c
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

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Basic per gradient computation options and logging */
/*---------------------------------------------------*/

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

static int cs_glob_gradient_n_systems = 0;      /* Current number of systems */
static int cs_glob_gradient_n_max_systems = 0;  /* Max. number of systems for
                                                   cs_glob_gradient_systems. */

/* System info array */
static cs_gradient_info_t **cs_glob_gradient_systems = NULL;

/* Short names for gradient computation types */

const char *cs_gradient_type_name[]
  = {N_("Iterative reconstruction"),
     N_("Least-squares"),
     N_("conservative, reconstruction with Least-squares"),
     N_("Iterative (old)")};

/* Timer statistics */

static cs_timer_counter_t   _gradient_t_tot;     /* Total time in gradients */
static int _gradient_stat_id = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start and past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static void
_thread_range(cs_lnum_t   n,
              cs_lnum_t  *s_id,
              cs_lnum_t  *e_id)
{
#if defined(HAVE_OPENMP)
  int t_id = omp_get_thread_num();
  int n_t = omp_get_num_threads();
  cs_lnum_t t_n = (n + n_t - 1) / n_t;
  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------
 * Factorize dense p*p symmetric matrices.
 * Only the lower triangular part is stored and the factorization is performed
 * in place (original coefficients are replaced).
 * Crout Factorization is performed (A = L D t(L)).
 *
 * parameters:
 *   d_size   <--  matrix size (p)
 *   ad       <--> symmetric matrix to be factorized
 *----------------------------------------------------------------------------*/

inline static void
_fact_crout_pp(const int   d_size,
               cs_real_t  *ad)
{
  cs_real_t aux[d_size];
  for (int kk = 0; kk < d_size - 1; kk++) {
    int kk_d_size = kk*(kk + 1)/2;
    for (int ii = kk + 1; ii < d_size; ii++) {
      int ii_d_size = ii*(ii + 1)/2;
      aux[ii] = ad[ii_d_size + kk];
      ad[ii_d_size + kk] =   ad[ii_d_size + kk]
                           / ad[kk_d_size + kk];
      for (int jj = kk + 1; jj < ii + 1; jj++) {
        ad[ii_d_size + jj] = ad[ii_d_size + jj] - ad[ii_d_size + kk]*aux[jj];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Solve forward and backward linear systems of the form L D t(L) x = b.
 * Matrix L D t(L) should be given by a Crout factorization.
 *
 * parameters:
 *   mat      <--  symmetric factorized (Crout) matrix
 *   d_size   <--  matrix size (p)
 *   x         --> linear system vector solution
 *   b        <--  linear system vector right hand side
 *----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_ldtl_pp(const cs_real_t mat[],
                   const int       d_size,
                         cs_real_t x[],
                   const cs_real_t b[])
{
  cs_real_t  aux[d_size];

  /* forward (stricly lower + identity) */
  for (int ii = 0; ii < d_size; ii++) {
    int ii_d_size = ii*(ii + 1)/2;
    aux[ii] = b[ii];
    for (int jj = 0; jj < ii; jj++) {
      aux[ii] -= aux[jj]*mat[ii_d_size + jj];
    }
  }

  /* diagonal */
  for (int ii = 0; ii < d_size; ii++) {
    int ii_d_size = ii*(ii + 1)/2;
    aux[ii] /= mat[ii_d_size + ii];
  }

  /* backward (transposed of strictly lower + identity) */
  for (int ii = d_size - 1; ii >= 0; ii--) {
    x[ii] = aux[ii];
    for (int jj = d_size - 1; jj > ii; jj--) {
      int jj_d_size = jj*(jj + 1)/2;
      x[ii] -= x[jj]*mat[jj_d_size + ii];
    }
  }
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
  cs_gradient_info_t *new_info = NULL;

  BFT_MALLOC(new_info, 1, cs_gradient_info_t);
  BFT_MALLOC(new_info->name, strlen(name) + 1, char);

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
  if (*this_info != NULL) {
    BFT_FREE((*this_info)->name);
    BFT_FREE(*this_info);
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
                  "Summary of gradient computations for \"%s\" (%s):\n\n"
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
                this_info->t_tot.wall_nsec*1e-9);
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
  end_id = cs_glob_gradient_n_systems - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((cs_glob_gradient_systems[mid_id])->name, name);
    if (cmp_ret == 0)
      cmp_ret = (cs_glob_gradient_systems[mid_id])->type - type;
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
    return cs_glob_gradient_systems[mid_id];

  /* Reallocate global array if necessary */

  if (cs_glob_gradient_n_systems >= cs_glob_gradient_n_max_systems) {

    if (cs_glob_gradient_n_max_systems == 0)
      cs_glob_gradient_n_max_systems = 10;
    else
      cs_glob_gradient_n_max_systems *= 2;

    BFT_REALLOC(cs_glob_gradient_systems,
                cs_glob_gradient_n_max_systems,
                cs_gradient_info_t *);

  }

  /* Insert in sorted list */

  for (ii = cs_glob_gradient_n_systems; ii > mid_id; ii--)
    cs_glob_gradient_systems[ii] = cs_glob_gradient_systems[ii - 1];

  cs_glob_gradient_systems[mid_id] = _gradient_info_create(name,
                                                           type);
  cs_glob_gradient_n_systems += 1;

  return cs_glob_gradient_systems[mid_id];
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
 * Compute triple L2 norm, summing result over axes.
 *
 * The input array is assumed to be interleaved with block of 4 values,
 * of which the first 3 are used.
 *
 * A superblock algorithm is used for better precision.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- array of 3-vectors
 *----------------------------------------------------------------------------*/

static double
_l2_norm_3(cs_lnum_t              n_elts,
           cs_real_4_t  *restrict x)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t n_blocks = n_elts / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double s[3];
  double s1 = 0.0, s2 = 0.0, s3 = 0.0;

# pragma omp parallel for reduction(+:s1, s2, s3)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot1 = 0.0;
    double sdot2 = 0.0;
    double sdot3 = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot1 = 0.0, cdot2 = 0.0, cdot3 = 0.0;
      for (cs_lnum_t ii = start_id; ii < end_id; ii++) {
        cdot1 += x[ii][0] * x[ii][0];
        cdot2 += x[ii][1] * x[ii][1];
        cdot3 += x[ii][2] * x[ii][2];
      }
      sdot1 += cdot1;
      sdot2 += cdot2;
      sdot3 += cdot3;
    }

    s1 += sdot1;
    s2 += sdot2;
    s3 += sdot3;

  }

  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n_elts;
  double cdot1 = 0.0, cdot2 = 0.0, cdot3 = 0.0;
  for (cs_lnum_t ii = start_id; ii < end_id; ii++) {
    cdot1 += x[ii][0] * x[ii][0];
    cdot2 += x[ii][1] * x[ii][1];
    cdot3 += x[ii][2] * x[ii][2];
  }

  s[0] = s1 + cdot1;
  s[1] = s2 + cdot2;
  s[2] = s3 + cdot3;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double _s[3];
    MPI_Allreduce(s, _s, 3, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    s[0] = _s[0];
    s[1] = _s[1];
    s[2] = _s[2];
  }

#endif /* defined(HAVE_MPI) */

  return (sqrt(s[0]) + sqrt(s[1]) + sqrt(s[2]));
}

/*----------------------------------------------------------------------------
 * Update R.H.S. for lsq gradient taking into account the weight coefficients.
 *
 * parameters:
 *   wi     <-- Weight coefficient of cell i
 *   wj     <-- Weight coefficient of cell j
 *   p_diff <-- R.H.S.
 *   d      <-- R.H.S.
 *   a      <-- geometric weight J'F/I'J'
 *   resi   --> Updated R.H.S. for cell i
 *   resj   --> Updated R.H.S. for cell j
 *----------------------------------------------------------------------------*/

static inline void
_compute_ani_weighting(const cs_real_t  wi[],
                       const cs_real_t  wj[],
                       const cs_real_t  p_diff,
                       const cs_real_t  d[],
                       const cs_real_t  a,
                       cs_real_t        resi[],
                       cs_real_t        resj[])
{
  int ii;
  cs_real_t ki_d[3] = {0., 0., 0.};
  cs_real_t kj_d[3] = {0., 0., 0.};

  cs_real_6_t sum;
  cs_real_6_t inv_wi;
  cs_real_6_t inv_wj;
  cs_real_t _d[3];

  for (ii = 0; ii < 6; ii++)
    sum[ii] = a*wi[ii] + (1. - a)*wj[ii];

  cs_math_sym_33_inv_cramer(wi, inv_wi);

  cs_math_sym_33_inv_cramer(wj, inv_wj);

  cs_math_sym_33_3_product(inv_wj, d,  _d);
  cs_math_sym_33_3_product(sum, _d, ki_d);
  cs_math_sym_33_3_product(inv_wi, d, _d);
  cs_math_sym_33_3_product(sum, _d, kj_d);

  /* 1 / ||Ki. K_f^-1. IJ||^2 */
  cs_real_t normi = 1. / cs_math_3_dot_product(ki_d, ki_d);
  /* 1 / ||Kj. K_f^-1. IJ||^2 */
  cs_real_t normj = 1. / cs_math_3_dot_product(kj_d, kj_d);

  for (ii = 0; ii < 3; ii++) {
    resi[ii] += p_diff * ki_d[ii] * normi;
    resj[ii] += p_diff * kj_d[ii] * normj;
  }
}

/*----------------------------------------------------------------------------
 * Compute the inverse of the face viscosity tensor and anisotropic vector
 * taking into account the weight coefficients to update cocg for lsq gradient.
 *
 * parameters:
 *   wi     <-- Weight coefficient of cell i
 *   wj     <-- Weight coefficient of cell j
 *   d      <-- IJ direction
 *   a      <-- geometric weight J'F/I'J'
 *   ki_d   --> Updated vector for cell i
 *   kj_d   --> Updated vector for cell j
 *----------------------------------------------------------------------------*/

static inline void
_compute_ani_weighting_cocg(const cs_real_t  wi[],
                            const cs_real_t  wj[],
                            const cs_real_t  d[],
                            const cs_real_t  a,
                            cs_real_t        ki_d[],
                            cs_real_t        kj_d[])
{
  int ii;
  cs_real_6_t sum;
  cs_real_6_t inv_wi;
  cs_real_6_t inv_wj;
  cs_real_t _d[3];

  for (ii = 0; ii < 6; ii++)
    sum[ii] = a*wi[ii] + (1. - a)*wj[ii];

  cs_math_sym_33_inv_cramer(wi,
                            inv_wi);
  cs_math_sym_33_inv_cramer(wj,
                            inv_wj);

  /* Note: K_i.K_f^-1 = SUM.K_j^-1
   *       K_j.K_f^-1 = SUM.K_i^-1
   * So: K_i d = SUM.K_j^-1.IJ */

  cs_math_sym_33_3_product(inv_wj,
                           d,
                           _d);
  cs_math_sym_33_3_product(sum,
                           _d,
                           ki_d);
  cs_math_sym_33_3_product(inv_wi,
                           d,
                           _d);
  cs_math_sym_33_3_product(sum,
                           _d,
                           kj_d);
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient least squares algorithm with
 * weighting coefficients.
 * The cocg are recomputed at this step to have the latest weighting
 * coefficients. Performance optimisation may be done.
 *
 * parameters:
 *   m         <--  mesh
 *   c_weight  <--  weighting coefficients
 *   fvq       <--  mesh quantities
 *   cpl       <--  pointer to coupling entity
 *   cocgb     -->  weighted cocgb
 *   cocg      -->  weighted cocg
 *----------------------------------------------------------------------------*/

static void
_compute_weighted_cell_cocg_s_lsq(const cs_mesh_t              *m,
                                  const cs_real_t              *c_weight,
                                  const cs_mesh_quantities_t   *fvq,
                                  const cs_internal_coupling_t *cpl,
                                  cs_real_33_t                 *cocgb,
                                  cs_real_33_t                 *cocg)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;

  const cs_real_t *restrict weight = fvq->weight;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Initialization */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t ll = 0; ll < 3; ll++) {
      for (cs_lnum_t mm = 0; mm < 3; mm++)
        cocg[cell_id][ll][mm] = 0.0;
    }
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        cs_real_t pond = weight[face_id];

        cs_real_t dc_i[3];
        cs_real_t dc_j[3];
        cs_real_t dc[3];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

        _compute_ani_weighting_cocg(&c_weight[ii*6],
                                    &c_weight[jj*6],
                                    dc,
                                    pond,
                                    dc_i,
                                    dc_j);

        cs_real_t i_dci = 1. / (dc_i[0]*dc_i[0] + dc_i[1]*dc_i[1] + dc_i[2]*dc_i[2]);
        cs_real_t i_dcj = 1. / (dc_j[0]*dc_j[0] + dc_j[1]*dc_j[1] + dc_j[2]*dc_j[2]);

        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[ii][ll][mm] += dc_i[mm] * dc_i[ll] * i_dci;
        }
        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[jj][ll][mm] += dc_j[mm] * dc_j[ll] * i_dcj;
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from extended neighborhood */

  if (m->halo_type == CS_HALO_EXTENDED) {

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
      for (cs_lnum_t cidx = cell_cells_idx[ii];
           cidx < cell_cells_idx[ii+1];
           cidx++) {

        cs_lnum_t jj = cell_cells_lst[cidx];

        cs_real_t dc[3];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];
        cs_real_t uddij2 = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[ii][ll][mm] += dc[ll] * dc[mm] * uddij2;
        }

      }
    }

  } /* End for extended neighborhood */

  /* Contribution from coupled faces */
  if (cpl != NULL)
    cs_internal_coupling_lsq_cocg_weighted(cpl,
                                           c_weight,
                                           cocg);

  /* Save partial cocg at interior faces of boundary cells */

# pragma omp parallel for
  for (cs_lnum_t ii = 0; ii < m->n_b_cells; ii++) {
    cs_lnum_t cell_id = m->b_cells[ii];
    for (cs_lnum_t ll = 0; ll < 3; ll++) {
      for (cs_lnum_t mm = 0; mm < 3; mm++)
        cocgb[ii][ll][mm] = cocg[cell_id][ll][mm];
    }
  }

  /* Contribution from boundary faces, assuming symmetry everywhere
     so as to avoid obtaining a non-invertible matrix in 2D cases. */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        if (cpl == NULL || !coupled_faces[face_id]) {

          cs_lnum_t ii = b_face_cells[face_id];

          cs_real_3_t normal;
          /* Normal is vector 0 if the b_face_normal norm is too small */
          cs_math_3_normalise(b_face_normal[face_id], normal);

          for (cs_lnum_t ll = 0; ll < 3; ll++) {
            for (cs_lnum_t mm = 0; mm < 3; mm++)
              cocg[ii][ll][mm] += normal[ll] * normal[mm];
          }

        } /* face without internal coupling */

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Invert for all cells. */
  /*-----------------------*/

  /* The cocg term for interior cells only changes if the mesh does */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    cs_math_33_inv_cramer_in_place(cocg[cell_id]);
}

/*----------------------------------------------------------------------------
 * Synchronize halos for scalar gradients.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   halo_type      <-- halo type (extended or not)
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_sync_scalar_gradient_halo(const cs_mesh_t  *m,
                           cs_halo_type_t    halo_type,
                           int               idimtr,
                           cs_real_3_t       grad[])
{
  if (m->halo != NULL) {
    if (idimtr == 0) {
      cs_halo_sync_var_strided
        (m->halo, halo_type, (cs_real_t *)grad, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect
          (m->halo, halo_type, (cs_real_t *)grad, 3);
    }
    else
      cs_halo_sync_components_strided(m->halo,
                                      halo_type,
                                      CS_HALO_ROTATION_IGNORE,
                                      (cs_real_t *)grad,
                                      3);
  }
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a scalar if necessary. This function deals with
 * the standard or extended neighborhood.
 *
 * parameters:
 *   halo_type      <-- halo type (extended or not)
 *   clip_mode      <-- type of clipping for the computation of the gradient
 *   iwarnp         <-- output level
 *   idimtr         <-- 0 for scalars or without rotational periodicity,
 *                      1 or 2 for vectors or tensors in case of rotational
 *                      periodicity
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   var            <-- variable
 *   grad           --> components of the pressure gradient
 *----------------------------------------------------------------------------*/

static void
_scalar_gradient_clipping(cs_halo_type_t         halo_type,
                          int                    clip_mode,
                          int                    verbosity,
                          int                    idimtr,
                          cs_real_t              climgp,
                          const cs_real_t        var[],
                          cs_real_3_t  *restrict grad)
{
  cs_gnum_t  t_n_clip;
  cs_real_t  global_min_factor, global_max_factor;
  cs_real_t  t_min_factor, t_max_factor;

  cs_gnum_t  n_clip = 0, n_g_clip = 0;
  cs_real_t  min_factor = 1, max_factor = 0;
  cs_real_t  *buf = NULL, *restrict clip_factor = NULL;
  cs_real_t  *restrict denom = NULL, *restrict denum = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const int n_i_groups = mesh->i_face_numbering->n_groups;
  const int n_i_threads = mesh->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = mesh->i_face_numbering->group_index;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;

  const cs_halo_t *halo = mesh->halo;

  if (clip_mode < 0)
    return;

  /* Synchronize variable */

  if (halo != NULL) {

    /* Exchange for the gradients. Not useful for working array */

    if (clip_mode == 1) {

      if (idimtr > 0)
        cs_halo_sync_components_strided(halo,
                                        halo_type,
                                        CS_HALO_ROTATION_IGNORE,
                                        (cs_real_t *restrict)grad,
                                        3);
      else {
        cs_halo_sync_var_strided(halo,
                                 halo_type,
                                 (cs_real_t *restrict)grad,
                                 3);
        cs_halo_perio_sync_var_vect(halo,
                                    halo_type,
                                    (cs_real_t *restrict)grad,
                                    3);
      }

    } /* End if clip_mode == 1 */

  } /* End if halo */

  /* Allocate and initialize working buffers */

  if (clip_mode == 1)
    BFT_MALLOC(buf, 3*n_cells_ext, cs_real_t);
  else
    BFT_MALLOC(buf, 2*n_cells_ext, cs_real_t);

  denum = buf;
  denom = buf + n_cells_ext;

  if (clip_mode == 1)
    clip_factor = buf + 2*n_cells_ext;

# pragma omp parallel for
  for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
    denum[ii] = 0;
    denom[ii] = 0;
  }

  /* First computations:
      denum holds the maximum variation of the gradient
      denom holds the maximum variation of the variable */

  if (clip_mode == 0) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t dist[3];
          for (int ll = 0; ll < 3; ll++)
            dist[ll] = cell_cen[ii][ll] - cell_cen[jj][ll];

          cs_real_t dist1 = CS_ABS(  dist[0]*grad[ii][0]
                                   + dist[1]*grad[ii][1]
                                   + dist[2]*grad[ii][2]);
          cs_real_t dist2 = CS_ABS(  dist[0]*grad[jj][0]
                                   + dist[1]*grad[jj][1]
                                   + dist[2]*grad[jj][2]);

          cs_real_t dvar = CS_ABS(var[ii] - var[jj]);

          denum[ii] = CS_MAX(denum[ii], dist1);
          denum[jj] = CS_MAX(denum[jj], dist2);
          denom[ii] = CS_MAX(denom[ii], dvar);
          denom[jj] = CS_MAX(denom[jj], dvar);

        } /* End of loop on faces */

      } /* End of loop on threads */

    } /* End of loop on thread groups */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
        for (cs_lnum_t cidx = cell_cells_idx[ii];
             cidx < cell_cells_idx[ii+1];
             cidx++) {

          cs_lnum_t jj = cell_cells_lst[cidx];

          cs_real_t dist[3];
          for (int ll = 0; ll < 3; ll++)
            dist[ll] = cell_cen[ii][ll] - cell_cen[jj][ll];

          cs_real_t dist1 = CS_ABS(  dist[0]*grad[ii][0]
                                   + dist[1]*grad[ii][1]
                                   + dist[2]*grad[ii][2]);
          cs_real_t dvar = CS_ABS(var[ii] - var[jj]);

          denum[ii] = CS_MAX(denum[ii], dist1);
          denom[ii] = CS_MAX(denom[ii], dvar);

        }
      }

    } /* End for extended halo */

  }
  else if (clip_mode == 1) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t dist[3];
          for (int ll = 0; ll < 3; ll++)
            dist[ll] = cell_cen[ii][ll] - cell_cen[jj][ll];

          cs_real_t dpdxf = 0.5 * (grad[ii][0] + grad[jj][0]);
          cs_real_t dpdyf = 0.5 * (grad[ii][1] + grad[jj][1]);
          cs_real_t dpdzf = 0.5 * (grad[ii][2] + grad[jj][2]);

          cs_real_t dist1 = CS_ABS(  dist[0]*dpdxf
                                   + dist[1]*dpdyf
                                   + dist[2]*dpdzf);
          cs_real_t dvar = CS_ABS(var[ii] - var[jj]);

          denum[ii] = CS_MAX(denum[ii], dist1);
          denum[jj] = CS_MAX(denum[jj], dist1);
          denom[ii] = CS_MAX(denom[ii], dvar);
          denom[jj] = CS_MAX(denom[jj], dvar);

        } /* End of loop on faces */

      } /* End of loop on threads */

    } /* End of loop on thread groups */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
        for (cs_lnum_t cidx = cell_cells_idx[ii];
             cidx < cell_cells_idx[ii+1];
             cidx++) {

          cs_lnum_t jj = cell_cells_lst[cidx];

          cs_real_t dist[3];
          for (int ll = 0; ll < 3; ll++)
            dist[ll] = cell_cen[ii][ll] - cell_cen[jj][ll];

          cs_real_t dpdxf = 0.5 * (grad[ii][0] + grad[jj][0]);
          cs_real_t dpdyf = 0.5 * (grad[ii][1] + grad[jj][1]);
          cs_real_t dpdzf = 0.5 * (grad[ii][2] + grad[jj][2]);

          cs_real_t dist1 = CS_ABS(  dist[0]*dpdxf
                                   + dist[1]*dpdyf
                                   + dist[2]*dpdzf);
          cs_real_t dvar = CS_ABS(var[ii] - var[jj]);

          denum[ii] = CS_MAX(denum[ii], dist1);
          denom[ii] = CS_MAX(denom[ii], dvar);

        }
      }

    } /* End for extended neighborhood */

  } /* End if clip_mode == 1 */

  /* Clipping of the gradient if denum/denom > climgp */

  if (clip_mode == 0) {

    t_min_factor = min_factor;
    t_max_factor = max_factor;

#   pragma omp parallel private(t_min_factor, t_max_factor, t_n_clip)
    {
      t_n_clip = 0;
      t_min_factor = min_factor; t_max_factor = max_factor;

#     pragma omp for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {

        if (denum[ii] > climgp * denom[ii]) {

          cs_real_t factor1 = climgp * denom[ii]/denum[ii];
          grad[ii][0] *= factor1;
          grad[ii][1] *= factor1;
          grad[ii][2] *= factor1;

          t_min_factor = CS_MIN(factor1, t_min_factor);
          t_max_factor = CS_MAX(factor1, t_max_factor);
          t_n_clip++;

        } /* If clipping */

      } /* End of loop on cells */

#     pragma omp critical
      {
        min_factor = CS_MIN(min_factor, t_min_factor);
        max_factor = CS_MAX(max_factor, t_max_factor);
        n_clip += t_n_clip;
      }
    } /* End of omp parallel construct */

  }
  else if (clip_mode == 1) {

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++)
      clip_factor[ii] = (cs_real_t)DBL_MAX;

    /* Synchronize variable */

    if (halo != NULL) {
      if (idimtr > 0) {
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, denom);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, denum);
      }
      else {
        cs_halo_sync_var(halo, halo_type, denom);
        cs_halo_sync_var(halo, halo_type, denum);
      }
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(min_factor)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t factor1 = 1.0;
          if (denum[ii] > climgp * denom[ii])
            factor1 = climgp * denom[ii]/denum[ii];

          cs_real_t factor2 = 1.0;
          if (denum[jj] > climgp * denom[jj])
            factor2 = climgp * denom[jj]/denum[jj];

          min_factor = CS_MIN(factor1, factor2);

          clip_factor[ii] = CS_MIN(clip_factor[ii], min_factor);
          clip_factor[jj] = CS_MIN(clip_factor[jj], min_factor);

        } /* End of loop on faces */

      } /* End of loop on threads */

    } /* End of loop on thread groups */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {

        cs_real_t factor1 = 1.0;

        for (cs_lnum_t cidx = cell_cells_idx[ii];
             cidx < cell_cells_idx[ii+1];
             cidx++) {

          cs_lnum_t jj = cell_cells_lst[cidx];

          cs_real_t factor2 = 1.0;

          if (denum[jj] > climgp * denom[jj])
            factor2 = climgp * denom[jj]/denum[jj];

          factor1 = CS_MIN(factor1, factor2);

        }

        clip_factor[ii] = CS_MIN(clip_factor[ii], factor1);

      } /* End of loop on cells */

    } /* End for extended neighborhood */

#   pragma omp parallel private(t_min_factor, t_max_factor, t_n_clip)
    {
      t_n_clip = 0;
      t_min_factor = min_factor; t_max_factor = max_factor;

#     pragma omp for
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {

        for (int ll = 0; ll < 3; ll++)
          grad[ii][ll] *= clip_factor[ii];

        if (clip_factor[ii] < 0.99) {
          t_max_factor = CS_MAX(t_max_factor, clip_factor[ii]);
          t_min_factor = CS_MIN(t_min_factor, clip_factor[ii]);
          t_n_clip++;
        }

      } /* End of loop on cells */

#     pragma omp critical
      {
        min_factor = CS_MIN(min_factor, t_min_factor);
        max_factor = CS_MAX(max_factor, t_max_factor);
        n_clip += t_n_clip;
      }
    } /* End of omp parallel construct */

  } /* End if clip_mode == 1 */

  /* Update min/max and n_clip in case of parallelism */

#if defined(HAVE_MPI)

  if (mesh->n_domains > 1) {

    assert(sizeof(cs_real_t) == sizeof(double));

    /* Global Max */

    MPI_Allreduce(&max_factor, &global_max_factor, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    max_factor = global_max_factor;

    /* Global min */

    MPI_Allreduce(&min_factor, &global_min_factor, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    min_factor = global_min_factor;

    /* Sum number of clippings */

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_GNUM,
                  MPI_SUM, cs_glob_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(HAVE_MPI) */

  /* Output warning if necessary */

  if (verbosity > 1)
    bft_printf(_(" Gradient limitation in %llu cells\n"
                 "   minimum factor = %14.5e; maximum factor = %14.5e\n"),
               (unsigned long long)n_clip, min_factor, max_factor);

  /* Synchronize grad */

  if (halo != NULL) {

    if (idimtr > 0) {

      /* If the gradient is not treated as a "true" vector */

      cs_halo_sync_components_strided(halo,
                                      halo_type,
                                      CS_HALO_ROTATION_IGNORE,
                                      (cs_real_t *restrict)grad,
                                      3);

    }
    else {

      cs_halo_sync_var_strided(halo,
                               halo_type,
                               (cs_real_t *restrict)grad,
                               3);

      cs_halo_perio_sync_var_vect(halo,
                                  halo_type,
                                  (cs_real_t *restrict)grad,
                                  3);

    }

  }

  BFT_FREE(buf);
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
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or NULL
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *   rhsv           <-> interleaved array for gradient RHS components
 *                      (0, 1, 2) and variable copy (3)
 *----------------------------------------------------------------------------*/

static void
_initialize_scalar_gradient_old(const cs_mesh_t             *m,
                                cs_mesh_quantities_t        *fvq,
                                int                          idimtr,
                                int                          hyd_p_flag,
                                cs_real_t                    inc,
                                const cs_real_3_t            f_ext[],
                                const cs_real_t              coefap[],
                                const cs_real_t              coefbp[],
                                const cs_real_t              pvar[],
                                const cs_real_t              c_weight[],
                                cs_real_3_t        *restrict grad,
                                cs_real_4_t        *restrict rhsv)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  cs_lnum_t  cell_id, face_id, ii, jj;
  int        g_id, t_id;
  cs_real_t  pfac;
  cs_real_t  ktpond;
  cs_real_4_t  fctb;

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  int is_porous = 0;
  if (f_i_poro_duq_0 != NULL) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name_try("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name_try("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  /* Initialize gradient */
  /*---------------------*/

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    rhsv[cell_id][0] = 0.0;
    rhsv[cell_id][1] = 0.0;
    rhsv[cell_id][2] = 0.0;
    rhsv[cell_id][3] = pvar[cell_id];
  }

  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  if (hyd_p_flag == 0 || hyd_p_flag == 2) {

    /* Pressure gradient with weighting activated */

    if (c_weight != NULL) {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, jj, ktpond, pfac, fctb)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            ktpond =    weight[face_id] * c_weight[ii]
                     / (       weight[face_id] * c_weight[ii]
                        + (1.0-weight[face_id])* c_weight[jj]);
            pfac  =        ktpond  * rhsv[ii][3]
                    + (1.0-ktpond) * rhsv[jj][3];
            fctb[0] = pfac * i_f_face_normal[face_id][0];
            fctb[1] = pfac * i_f_face_normal[face_id][1];
            fctb[2] = pfac * i_f_face_normal[face_id][2];
            rhsv[ii][0] += fctb[0];
            rhsv[ii][1] += fctb[1];
            rhsv[ii][2] += fctb[2];
            rhsv[jj][0] -= fctb[0];
            rhsv[jj][1] -= fctb[1];
            rhsv[jj][2] -= fctb[2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    }
    else {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, jj, pfac, fctb)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            pfac  =        weight[face_id]  * rhsv[ii][3]
                    + (1.0-weight[face_id]) * rhsv[jj][3];
            fctb[0] = pfac * i_f_face_normal[face_id][0];
            fctb[1] = pfac * i_f_face_normal[face_id][1];
            fctb[2] = pfac * i_f_face_normal[face_id][2];
            rhsv[ii][0] += fctb[0];
            rhsv[ii][1] += fctb[1];
            rhsv[ii][2] += fctb[2];
            rhsv[jj][0] -= fctb[0];
            rhsv[jj][1] -= fctb[1];
            rhsv[jj][2] -= fctb[2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    } /* loop on contribution for interior faces without weighting */

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(face_id, ii, pfac)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          pfac = inc*coefap[face_id] + coefbp[face_id]*rhsv[ii][3];
          rhsv[ii][0] += pfac * b_f_face_normal[face_id][0];
          rhsv[ii][1] += pfac * b_f_face_normal[face_id][1];
          rhsv[ii][2] += pfac * b_f_face_normal[face_id][2];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  }

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  else {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(face_id, ii, jj, pfac, fctb)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          cs_real_2_t poro = {
            i_poro_duq_0[is_porous*face_id],
            i_poro_duq_1[is_porous*face_id]
          };

          pfac
            =   (  weight[face_id]
                 * (  rhsv[ii][3]
                    - (cell_cen[ii][0] - i_face_cog[face_id][0])*f_ext[ii][0]
                    - (cell_cen[ii][1] - i_face_cog[face_id][1])*f_ext[ii][1]
                    - (cell_cen[ii][2] - i_face_cog[face_id][2])*f_ext[ii][2]
                    + poro[0]
                    ))
              + ( (1.0 - weight[face_id])
                 * (  rhsv[jj][3]
                    - (cell_cen[jj][0] - i_face_cog[face_id][0])*f_ext[jj][0]
                    - (cell_cen[jj][1] - i_face_cog[face_id][1])*f_ext[jj][1]
                    - (cell_cen[jj][2] - i_face_cog[face_id][2])*f_ext[jj][2]
                    - poro[1]
                    ));

          fctb[0] = pfac * i_f_face_normal[face_id][0];
          fctb[1] = pfac * i_f_face_normal[face_id][1];
          fctb[2] = pfac * i_f_face_normal[face_id][2];
          rhsv[ii][0] += fctb[0];
          rhsv[ii][1] += fctb[1];
          rhsv[ii][2] += fctb[2];
          rhsv[jj][0] -= fctb[0];
          rhsv[jj][1] -= fctb[1];
          rhsv[jj][2] -= fctb[2];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(face_id, ii, pfac)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          cs_real_t poro = b_poro_duq[is_porous*face_id];

          pfac
            =      coefap[face_id] * inc
              + (  coefbp[face_id]
                 * (  rhsv[ii][3]
                    - (cell_cen[ii][0] - b_face_cog[face_id][0])*f_ext[ii][0]
                    - (cell_cen[ii][1] - b_face_cog[face_id][1])*f_ext[ii][1]
                    - (cell_cen[ii][2] - b_face_cog[face_id][2])*f_ext[ii][2]
                    + poro));

          rhsv[ii][0] += pfac * b_f_face_normal[face_id][0];
          rhsv[ii][1] += pfac * b_f_face_normal[face_id][1];
          rhsv[ii][2] += pfac * b_f_face_normal[face_id][2];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  } /* End of test on hydrostatic pressure */

  /* Compute gradient */
  /*------------------*/

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell fully solid? */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    grad[cell_id][0] = rhsv[cell_id][0] * dvol;
    grad[cell_id][1] = rhsv[cell_id][1] * dvol;
    grad[cell_id][2] = rhsv[cell_id][2] * dvol;
  }

  /* Synchronize halos */

  _sync_scalar_gradient_halo(m, CS_HALO_EXTENDED, idimtr, grad);
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
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   f_ext          <-- exterior force generating pressure
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_initialize_scalar_gradient(const cs_mesh_t                *m,
                            cs_mesh_quantities_t           *fvq,
                            const cs_internal_coupling_t   *cpl,
                            int                             idimtr,
                            int                             hyd_p_flag,
                            cs_real_t                       inc,
                            const cs_real_3_t               f_ext[],
                            const cs_real_t                 coefap[],
                            const cs_real_t                 coefbp[],
                            const cs_real_t                 pvar[],
                            const cs_real_t                 c_weight[],
                            cs_real_3_t           *restrict grad)
{
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  cs_lnum_t  ii, jj;
  int        g_id, t_id;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  int is_porous = 0;
  if (f_i_poro_duq_0 != NULL) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name_try("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name_try("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  /* Initialize gradient */
  /*---------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int j = 0; j < 3; j++)
      grad[cell_id][j] = 0.0;
  }

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  if (hyd_p_flag == 1) {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(ii, jj)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          cs_real_t ktpond = (c_weight == NULL) ?
             weight[face_id] :              /* no cell weightening */
             weight[face_id] * c_weight[ii] /* cell weightening active */
               / (      weight[face_id] * c_weight[ii]
                 + (1.0-weight[face_id])* c_weight[jj]);

          cs_real_2_t poro = {
            i_poro_duq_0[is_porous*face_id],
            i_poro_duq_1[is_porous*face_id]
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
                 * (  (i_face_cog[face_id][0] - cell_cen[ii][0])*f_ext[ii][0]
                    + (i_face_cog[face_id][1] - cell_cen[ii][1])*f_ext[ii][1]
                    + (i_face_cog[face_id][2] - cell_cen[ii][2])*f_ext[ii][2]
                    + poro[0])
            +  (1.0 - ktpond)
                 * (  (i_face_cog[face_id][0] - cell_cen[jj][0])*f_ext[jj][0]
                    + (i_face_cog[face_id][1] - cell_cen[jj][1])*f_ext[jj][1]
                    + (i_face_cog[face_id][2] - cell_cen[jj][2])*f_ext[jj][2]
                    - poro[1]);

          cs_real_t pfacj = pfaci;

          pfaci += (1.0-ktpond) * (pvar[jj] - pvar[ii]);
          pfacj -=      ktpond  * (pvar[jj] - pvar[ii]);

          for (int j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[face_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[face_id][j];
          }


        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(ii)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          cs_real_t poro = b_poro_duq[is_porous*face_id];

          /*
             Remark: for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
           */

          /* Reconstruction part */
          cs_real_t pfac
            = coefap[face_id] * inc
            + coefbp[face_id]
              * ( (b_face_cog[face_id][0] - cell_cen[ii][0])*f_ext[ii][0]
                + (b_face_cog[face_id][1] - cell_cen[ii][1])*f_ext[ii][1]
                + (b_face_cog[face_id][2] - cell_cen[ii][2])*f_ext[ii][2]
                + poro);

          pfac += (coefbp[face_id] - 1.0) * pvar[ii];

          for (int j = 0; j < 3; j++)
            grad[ii][j] += pfac * b_f_face_normal[face_id][j];

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  } /* End of test on hydrostatic pressure */


  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  else {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(ii, jj)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          cs_real_t ktpond = (c_weight == NULL) ?
             weight[face_id] :              /* no cell weightening */
             weight[face_id] * c_weight[ii] /* cell weightening active */
               / (      weight[face_id] * c_weight[ii]
                 + (1.0-weight[face_id])* c_weight[jj]);

          /*
             Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                      + (1-\alpha_\ij) \varia_\cellj\f$
                     but for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                     and for the cell \f$ \cellj \f$ we remove
                     \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
          */
          cs_real_t pfaci = (1.0-ktpond) * (pvar[jj] - pvar[ii]);
          cs_real_t pfacj = - ktpond * (pvar[jj] - pvar[ii]);

          for (int j = 0; j < 3; j++) {
            grad[ii][j] += pfaci * i_f_face_normal[face_id][j];
            grad[jj][j] -= pfacj * i_f_face_normal[face_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from coupled faces */
    if (cpl != NULL)
      cs_internal_coupling_initialize_scalar_gradient
        (cpl, c_weight, pvar, grad);

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(ii)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          if (cpl == NULL || !coupled_faces[face_id]) {

            ii = b_face_cells[face_id];

            /*
               Remark: for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
             */

            cs_real_t pfac = inc*coefap[face_id] + (coefbp[face_id]-1.0)*pvar[ii];

            for (int j = 0; j < 3; j++)
              grad[ii][j] += pfac * b_f_face_normal[face_id][j];

          } /* face without internal coupling */

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    for (int j = 0; j < 3; j++)
      grad[cell_id][j] *= dvol;
  }

  /* Synchronize halos */

  _sync_scalar_gradient_halo(m, CS_HALO_EXTENDED, idimtr, grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient using iterative reconstruction for non-orthogonal
 * meshes (nswrgp > 1).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * cocg is computed to account for variable B.C.'s (flux).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   var_name       <-- variable name
 *   gradient_info  <-- pointer to performance logging structure, or NULL
 *   recompute_cocg <-- flag to recompute cocg
 *   nswrgp         <-- number of sweeps for gradient reconstruction
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   verbosity      <-- verbosity level
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   epsrgp         <-- relative precision for gradient reconstruction
 *   extrap         <-- gradient extrapolation coefficient
 *   f_ext          <-- exterior force generating pressure
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *   rhsv           <-> interleaved array for gradient RHS components
 *                      (0, 1, 2) and variable copy (3)
 *----------------------------------------------------------------------------*/

static void
_iterative_scalar_gradient_old(const cs_mesh_t             *m,
                               cs_mesh_quantities_t        *fvq,
                               const char                  *var_name,
                               cs_gradient_info_t          *gradient_info,
                               bool                         recompute_cocg,
                               int                          nswrgp,
                               int                          idimtr,
                               int                          hyd_p_flag,
                               int                          verbosity,
                               cs_real_t                    inc,
                               cs_real_t                    epsrgp,
                               cs_real_t                    extrap,
                               const cs_real_3_t            f_ext[],
                               const cs_real_t              coefap[],
                               const cs_real_t              coefbp[],
                               cs_real_3_t        *restrict grad,
                               cs_real_4_t        *restrict rhsv)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  cs_real_33_t   *restrict cocgb = fvq->cocgb_s_it;
  cs_real_33_t   *restrict cocg = fvq->cocg_s_it;

  cs_lnum_t  cell_id, face_id, ii, jj, ll, mm;
  int        g_id, t_id;
  cs_real_t  rnorm;
  cs_real_t  pfac, pfac0, pfac1, pip;
  cs_real_3_t  fexd;
  cs_real_4_t  fctb;

  int nswmax = nswrgp;
  int n_sweeps = 0;
  cs_real_t residue = 0.;

  if (nswrgp <  1) {
    if (gradient_info != NULL)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  /* Reconstruct gradients for non-orthogonal meshes */
  /*-------------------------------------------------*/

  /* Semi-implicit resolution on the whole mesh  */

  /* If cocg must be recomputed, only do it for boundary cells,
     with saved cocgb */

  if (recompute_cocg) {

#   pragma omp parallel for private(cell_id, ll, mm)
    for (ii = 0; ii < m->n_b_cells; ii++) {
      cell_id = m->b_cells[ii];
      for (ll = 0; ll < 3; ll++) {
        for (mm = 0; mm < 3; mm++)
          cocg[cell_id][ll][mm] = cocgb[ii][ll][mm];
      }
    }

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(face_id, ii, ll, mm)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          for (ll = 0; ll < 3; ll++) {
            for (mm = 0; mm < 3; mm++)
              cocg[ii][ll][mm] -= (  coefbp[face_id]*diipb[face_id][mm]
                                   * b_f_face_normal[face_id][ll]);
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

#   pragma omp parallel for private(cell_id)
    for (ii = 0; ii < m->n_b_cells; ii++) {
      cell_id = m->b_cells[ii];
      cs_math_33_inv_cramer_in_place(cocg[cell_id]);
    }

  } /* End of test on ale, mobile mesh, or call counter */

  /* Compute normalization residue */

  rnorm = _l2_norm_3(n_cells, rhsv);

  if (fvq->max_vol > 1)
    rnorm /= fvq->max_vol;

  if (rnorm <= cs_math_epzero) {
    if (gradient_info != NULL)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  /* Vector OijFij is computed in CLDijP */

  /* Start iterations */
  /*------------------*/

  for (n_sweeps = 1; n_sweeps < nswmax; n_sweeps++) {

    /* Compute right hand side */

#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      rhsv[cell_id][0] = -grad[cell_id][0] * cell_f_vol[cell_id];
      rhsv[cell_id][1] = -grad[cell_id][1] * cell_f_vol[cell_id];
      rhsv[cell_id][2] = -grad[cell_id][2] * cell_f_vol[cell_id];
    }

    /* Standard case, without hydrostatic pressure */
    /*---------------------------------------------*/

    if (hyd_p_flag == 0 || hyd_p_flag == 2) {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, jj, pfac, fctb)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            pfac  =        weight[face_id]  * rhsv[ii][3]
                    + (1.0-weight[face_id]) * rhsv[jj][3]
                    + ( dofij[face_id][0] * (grad[ii][0]+grad[jj][0])
                    +   dofij[face_id][1] * (grad[ii][1]+grad[jj][1])
                    +   dofij[face_id][2] * (grad[ii][2]+grad[jj][2])) * 0.5;
            fctb[0] = pfac * i_f_face_normal[face_id][0];
            fctb[1] = pfac * i_f_face_normal[face_id][1];
            fctb[2] = pfac * i_f_face_normal[face_id][2];
            rhsv[ii][0] += fctb[0];
            rhsv[ii][1] += fctb[1];
            rhsv[ii][2] += fctb[2];
            rhsv[jj][0] -= fctb[0];
            rhsv[jj][1] -= fctb[1];
            rhsv[jj][2] -= fctb[2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from boundary faces */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, pip, pfac0, pfac1, pfac)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pip =   rhsv[ii][3]
                  + diipb[face_id][0] * grad[ii][0]
                  + diipb[face_id][1] * grad[ii][1]
                  + diipb[face_id][2] * grad[ii][2];

            pfac0 =   coefap[face_id] * inc
                    + coefbp[face_id] * pip;

            pfac1 =   rhsv[ii][3]
                    + (b_face_cog[face_id][0]-cell_cen[ii][0]) * grad[ii][0]
                    + (b_face_cog[face_id][1]-cell_cen[ii][1]) * grad[ii][1]
                    + (b_face_cog[face_id][2]-cell_cen[ii][2]) * grad[ii][2];

            pfac =          coefbp[face_id]  *(extrap*pfac1 + (1.0-extrap)*pfac0)
                   + (1.0 - coefbp[face_id]) * pfac0;

            rhsv[ii][0] = rhsv[ii][0] + pfac * b_f_face_normal[face_id][0];
            rhsv[ii][1] = rhsv[ii][1] + pfac * b_f_face_normal[face_id][1];
            rhsv[ii][2] = rhsv[ii][2] + pfac * b_f_face_normal[face_id][2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    }

    /* Case with hydrostatic pressure */
    /*--------------------------------*/

    else {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, jj, pfac, fexd, fctb)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            fexd[0] = 0.5 * (f_ext[ii][0] - f_ext[jj][0]);
            fexd[1] = 0.5 * (f_ext[ii][1] - f_ext[jj][1]);
            fexd[2] = 0.5 * (f_ext[ii][2] - f_ext[jj][2]);

            /* Note: changed expression from:
             *   fmean = 0.5 * (f_ext[ii] + f_ext[jj])
             *   fii = f_ext[ii] - fmean
             *   fjj = f_ext[jj] - fmean
             * to:
             *   fexd = 0.5 * (f_ext[ii] - f_ext[jj])
             *   fii =  fexd
             *   fjj = -fexd
             */

            pfac
              =   (  weight[face_id]
                   * (  rhsv[ii][3]
                      - (cell_cen[ii][0]-i_face_cog[face_id][0])*fexd[0]
                      - (cell_cen[ii][1]-i_face_cog[face_id][1])*fexd[1]
                      - (cell_cen[ii][2]-i_face_cog[face_id][2])*fexd[2]))
              +   (  (1.0 - weight[face_id])
                   * (  rhsv[jj][3]
                      + (cell_cen[jj][0]-i_face_cog[face_id][0])*fexd[0]
                      + (cell_cen[jj][1]-i_face_cog[face_id][1])*fexd[1]
                      + (cell_cen[jj][2]-i_face_cog[face_id][2])*fexd[2]))
              +   (  dofij[face_id][0] * (grad[ii][0]+grad[jj][0])
                   + dofij[face_id][1] * (grad[ii][1]+grad[jj][1])
                   + dofij[face_id][2] * (grad[ii][2]+grad[jj][2]))*0.5;

            fctb[0] = pfac * i_f_face_normal[face_id][0];
            fctb[1] = pfac * i_f_face_normal[face_id][1];
            fctb[2] = pfac * i_f_face_normal[face_id][2];

            rhsv[ii][0] += fctb[0];
            rhsv[ii][1] += fctb[1];
            rhsv[ii][2] += fctb[2];
            rhsv[jj][0] -= fctb[0];
            rhsv[jj][1] -= fctb[1];
            rhsv[jj][2] -= fctb[2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from boundary faces */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id, ii, pip, pfac0, pfac1, pfac)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pip =   rhsv[ii][3]
                  + diipb[face_id][0] * grad[ii][0]
                  + diipb[face_id][1] * grad[ii][1]
                  + diipb[face_id][2] * grad[ii][2];

            pfac0 =      coefap[face_id] * inc
                    +    coefbp[face_id]
                       * (  pip
                          - (  cell_cen[ii][0]
                             - b_face_cog[face_id][0]
                             + diipb[face_id][0]) * f_ext[ii][0]
                          - (  cell_cen[ii][1]
                             - b_face_cog[face_id][1]
                             + diipb[face_id][1]) * f_ext[ii][1]
                          - (  cell_cen[ii][2]
                             - b_face_cog[face_id][2]
                             + diipb[face_id][2]) * f_ext[ii][2]);

            pfac1 =   rhsv[ii][3]
                    + (b_face_cog[face_id][0]-cell_cen[ii][0]) * grad[ii][0]
                    + (b_face_cog[face_id][1]-cell_cen[ii][1]) * grad[ii][1]
                    + (b_face_cog[face_id][2]-cell_cen[ii][2]) * grad[ii][2];

            pfac =          coefbp[face_id]  *(extrap*pfac1 + (1.0-extrap)*pfac0)
                   + (1.0 - coefbp[face_id]) * pfac0;

            rhsv[ii][0] += pfac * b_f_face_normal[face_id][0];
            rhsv[ii][1] += pfac * b_f_face_normal[face_id][1];
            rhsv[ii][2] += pfac * b_f_face_normal[face_id][2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    } /* End of test on hydrostatic pressure */

    /* Increment gradient */
    /*--------------------*/

#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      grad[cell_id][0] +=   cocg[cell_id][0][0] * rhsv[cell_id][0]
                          + cocg[cell_id][0][1] * rhsv[cell_id][1]
                          + cocg[cell_id][0][2] * rhsv[cell_id][2];
      grad[cell_id][1] +=   cocg[cell_id][1][0] * rhsv[cell_id][0]
                          + cocg[cell_id][1][1] * rhsv[cell_id][1]
                          + cocg[cell_id][1][2] * rhsv[cell_id][2];
      grad[cell_id][2] +=   cocg[cell_id][2][0] * rhsv[cell_id][0]
                          + cocg[cell_id][2][1] * rhsv[cell_id][1]
                          + cocg[cell_id][2][2] * rhsv[cell_id][2];
    }

    /* Synchronize halos */

    _sync_scalar_gradient_halo(m, CS_HALO_STANDARD, idimtr, grad);

    /* Convergence test */

    residue = _l2_norm_3(n_cells, rhsv);

    if (fvq->max_vol > 1)
      residue /= fvq->max_vol;

    if (residue < epsrgp*rnorm) {
      if (verbosity > 1)
        bft_printf(_(" %s; variable: %s; converged in %d sweeps\n"
                     " %*s  normed residual: %11.4e; norm: %11.4e\n"),
                   __func__, var_name, n_sweeps,
                   (int)(strlen(__func__)), " ", residue/rnorm, rnorm);
      break;
    }

  } /* Loop on sweeps */

  if (residue >= epsrgp*rnorm && verbosity > -1) {
    bft_printf(_(" Warning:\n"
                 " --------\n"
                 "   %s; variable: %s; sweeps: %d\n"
                 "   %*s  normed residual: %11.4e; norm: %11.4e\n"),
               __func__, var_name, n_sweeps,
               (int)(strlen(__func__)), " ", residue/rnorm, rnorm);
  }

  if (gradient_info != NULL)
    _gradient_info_update_iter(gradient_info, n_sweeps);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient using least-squares reconstruction for non-orthogonal
 * meshes (nswrgp > 1).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * cocg is computed to account for variable B.C.'s (flux).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   halo_type      <-- halo type (extended or not)
 *   recompute_cocg <-- flag to recompute cocg
 *   nswrgp         <-- number of sweeps for gradient reconstruction
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   w_stride       <-- stride for weighting coefficient
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   extrap         <-- gradient extrapolation coefficient
 *   fextx          <-- x component of exterior force generating pressure
 *   fexty          <-- y component of exterior force generating pressure
 *   fextz          <-- z component of exterior force generating pressure
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or NULL
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_lsq_scalar_gradient(const cs_mesh_t                *m,
                     cs_mesh_quantities_t           *fvq,
                     const cs_internal_coupling_t   *cpl,
                     cs_halo_type_t                  halo_type,
                     bool                            recompute_cocg,
                     int                             nswrgp,
                     int                             idimtr,
                     int                             hyd_p_flag,
                     int                             w_stride,
                     cs_real_t                       inc,
                     cs_real_t                       extrap,
                     const cs_real_3_t               f_ext[],
                     const cs_real_t                 coefap[],
                     const cs_real_t                 coefbp[],
                     const cs_real_t                 pvar[],
                     const cs_real_t       *restrict c_weight,
                     cs_real_3_t           *restrict grad,
                     cs_real_4_t           *restrict rhsv)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)fvq->b_dist;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_int_t *isympa = fvq->b_sym_flag;
  const cs_real_t *restrict weight = fvq->weight;

  cs_real_33_t   *restrict cocgb = (cpl == NULL) ?
    fvq->cocgb_s_lsq :
    cpl->cocgb_s_lsq;
  cs_real_33_t   *restrict cocg = fvq->cocg_lsq;
  cs_real_33_t   *restrict _cocg = NULL;
  cs_real_33_t   *restrict _cocgb = NULL;

  int        g_id, t_id;
  cs_real_t  pfac;
  cs_real_t  extrab, unddij, umcbdd, udbfs;
  cs_real_3_t  dc, dddij, dsij;
  cs_real_4_t  fctb;

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t _f_ext = 0.;

  int is_porous = 0;
  if (f_i_poro_duq_0 != NULL) {
    is_porous = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name_try("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name_try("b_poro_duq")->val;
  } else {
    i_poro_duq_0 = &_f_ext;
    i_poro_duq_1 = &_f_ext;
    b_poro_duq = &_f_ext;
  }

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Remark:

     for 2D calculations, if we extrapolate the pressure gradient,
     we obtain a non-invertible cocg matrix, because of the third
     direction.

     To avoid this, we multiply extrap by isympa which is zero for
     symmetries: the gradient is thus not extrapolated on those faces. */

   if (c_weight != NULL) {
     if (w_stride == 6) {
       BFT_MALLOC(_cocgb, m->n_b_cells, cs_real_33_t);
       BFT_MALLOC(_cocg, n_cells_ext, cs_real_33_t);
       _compute_weighted_cell_cocg_s_lsq(cs_glob_mesh,
                                         c_weight,
                                         cs_glob_mesh_quantities,
                                         cpl,
                                         _cocgb,
                                         _cocg);
       cocg = _cocg;
       cocgb = _cocgb;
       recompute_cocg = true;
     }
   }

  /* Initialize gradient */
  /*---------------------*/

  if (nswrgp <= 1) {

    _initialize_scalar_gradient(m,
                                fvq,
                                cpl,
                                idimtr,
                                hyd_p_flag,
                                inc,
                                f_ext,
                                coefap,
                                coefbp,
                                pvar,
                                c_weight,
                                grad);

    return;

  }

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  /* Compute cocg and save contribution at boundaries */

  if (recompute_cocg) {

    /* Recompute cocg at boundaries, using saved cocgb */

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < m->n_b_cells; ii++) {
      cs_lnum_t cell_id = m->b_cells[ii];
      for (cs_lnum_t ll = 0; ll < 3; ll++) {
        for (cs_lnum_t mm = 0; mm < 3; mm++)
          cocg[cell_id][ll][mm] = cocgb[ii][ll][mm];
      }
    }

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(extrab, umcbdd, udbfs, dddij)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          if (cpl == NULL || !coupled_faces[face_id]) {

            cs_lnum_t ii = b_face_cells[face_id];

            extrab = 1. - isympa[face_id]*extrap*coefbp[face_id];

            umcbdd = extrab * (1. - coefbp[face_id]) / b_dist[face_id];
            udbfs = extrab / b_face_surf[face_id];

            for (cs_lnum_t ll = 0; ll < 3; ll++)
              dddij[ll] =   udbfs * b_face_normal[face_id][ll]
                          + umcbdd * diipb[face_id][ll];

            for (cs_lnum_t ll = 0; ll < 3; ll++) {
              for (cs_lnum_t mm = 0; mm < 3; mm++)
                cocg[ii][ll][mm] += dddij[ll]*dddij[mm];
            }

          }  /* face without internal coupling */

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < m->n_b_cells; ii++) {
      cs_lnum_t cell_id = m->b_cells[ii];
      cs_math_33_inv_cramer_sym_in_place(cocg[cell_id]);
    }

  } /* End of recompute_cocg */

  /* Compute Right-Hand Side */
  /*-------------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    rhsv[cell_id][0] = 0.0;
    rhsv[cell_id][1] = 0.0;
    rhsv[cell_id][2] = 0.0;
    rhsv[cell_id][3] = pvar[cell_id];
  }

  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  if (hyd_p_flag == 0 || hyd_p_flag == 2) {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(pfac, dc, fctb)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t pond = weight[face_id];

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

          if (c_weight != NULL) {
            if (w_stride == 6) {
              /* (P_j - P_i)*/
              cs_real_t p_diff = (rhsv[jj][3] - rhsv[ii][3]);

              _compute_ani_weighting(&c_weight[ii*6],
                                     &c_weight[jj*6],
                                     p_diff,
                                     dc,
                                     pond,
                                     &rhsv[ii][0],
                                     &rhsv[jj][0]);
            }
            else {
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

    if (halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for private(dc, fctb, pfac)
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
        for (cs_lnum_t cidx = cell_cells_idx[ii];
             cidx < cell_cells_idx[ii+1];
             cidx++) {

          cs_lnum_t jj = cell_cells_lst[cidx];

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

    /* Contribution from coupled faces */

    if (cpl != NULL)
      cs_internal_coupling_lsq_scalar_gradient
        (cpl, c_weight, w_stride, rhsv);

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(extrab, \
                                      unddij, udbfs, umcbdd, pfac, dsij)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          if (cpl == NULL || !coupled_faces[face_id]) {

            cs_lnum_t ii = b_face_cells[face_id];

            extrab = pow((1. - isympa[face_id]*extrap*coefbp[face_id]), 2.0);
            unddij = 1. / b_dist[face_id];
            udbfs = 1. / b_face_surf[face_id];
            umcbdd = (1. - coefbp[face_id]) * unddij;

            for (cs_lnum_t ll = 0; ll < 3; ll++)
              dsij[ll] =   udbfs * b_face_normal[face_id][ll]
                       + umcbdd*diipb[face_id][ll];

            pfac =   (coefap[face_id]*inc + (coefbp[face_id] -1.)*rhsv[ii][3])
                 * unddij * extrab;

            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhsv[ii][ll] += dsij[ll] * pfac;

          } /* face without internal coupling */

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  }

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  else {  /* if hyd_p_flag == 1 */

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(dc, pfac, fctb)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_2_t poro = {
            i_poro_duq_0[is_porous*face_id],
            i_poro_duq_1[is_porous*face_id]
          };

          cs_real_t pond = weight[face_id];

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

          pfac =   (  rhsv[jj][3] - rhsv[ii][3]
                    + (cell_cen[ii][0] - i_face_cog[face_id][0]) * f_ext[ii][0]
                    + (cell_cen[ii][1] - i_face_cog[face_id][1]) * f_ext[ii][1]
                    + (cell_cen[ii][2] - i_face_cog[face_id][2]) * f_ext[ii][2]
                    + poro[0]
                    - (cell_cen[jj][0] - i_face_cog[face_id][0]) * f_ext[jj][0]
                    - (cell_cen[jj][1] - i_face_cog[face_id][1]) * f_ext[jj][1]
                    - (cell_cen[jj][2] - i_face_cog[face_id][2]) * f_ext[jj][2]
                    - poro[1])
                  / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            fctb[ll] = dc[ll] * pfac;

          if (c_weight != NULL) {
              cs_real_t denom = 1. / (  pond       *c_weight[ii]
                                      + (1. - pond)*c_weight[jj]);

              for (cs_lnum_t ll = 0; ll < 3; ll++)
                rhsv[ii][ll] += c_weight[jj] * denom * fctb[ll];

              for (cs_lnum_t ll = 0; ll < 3; ll++)
                rhsv[jj][ll] += c_weight[ii] * denom * fctb[ll];
          }
          else { // no cell weightening
            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhsv[ii][ll] += fctb[ll];

            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhsv[jj][ll] += fctb[ll];
          }
        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from extended neighborhood;
       We assume that the middle of the segment joining cell centers
       may replace the center of gravity of a fictitious face. */

    if (halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for private(dc, fctb, pfac)
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
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

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];

          pfac =   (  rhsv[jj][3] - rhsv[ii][3]
                    - 0.5 * dc[0] * f_ext[ii][0]
                    - 0.5 * dc[1] * f_ext[ii][1]
                    - 0.5 * dc[2] * f_ext[ii][2]
                    - 0.5 * dc[0] * f_ext[jj][0]
                    - 0.5 * dc[1] * f_ext[jj][1]
                    - 0.5 * dc[2] * f_ext[jj][2])
                  / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            fctb[ll] = dc[ll] * pfac;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ii][ll] += fctb[ll];

        }
      }

    } /* End for extended neighborhood */

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for private(extrab, \
                                      unddij, udbfs, umcbdd, pfac, dsij)
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          cs_real_t poro = b_poro_duq[is_porous*face_id];

          extrab = pow((1. - isympa[face_id]*extrap*coefbp[face_id]), 2.0);
          unddij = 1. / b_dist[face_id];
          udbfs = 1. / b_face_surf[face_id];
          umcbdd = (1. - coefbp[face_id]) * unddij;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dsij[ll] =   udbfs * b_face_normal[face_id][ll]
                       + umcbdd*diipb[face_id][ll];

          pfac
            =   (coefap[face_id]*inc
              + (  (coefbp[face_id] -1.)
                 * (  rhsv[ii][3]
                    + (b_face_cog[face_id][0] - cell_cen[ii][0]) * f_ext[ii][0]
                    + (b_face_cog[face_id][1] - cell_cen[ii][1]) * f_ext[ii][1]
                    + (b_face_cog[face_id][2] - cell_cen[ii][2]) * f_ext[ii][2]
                    + poro)))
              * unddij * extrab;

          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhsv[ii][ll] += dsij[ll] * pfac;

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  } /* End of test on hydrostatic pressure */

  /* Compute gradient */
  /*------------------*/

  if (hyd_p_flag == 1) {

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      grad[cell_id][0] =   cocg[cell_id][0][0] *rhsv[cell_id][0]
                         + cocg[cell_id][0][1] *rhsv[cell_id][1]
                         + cocg[cell_id][0][2] *rhsv[cell_id][2]
                         + f_ext[cell_id][0];
      grad[cell_id][1] =   cocg[cell_id][1][0] *rhsv[cell_id][0]
                         + cocg[cell_id][1][1] *rhsv[cell_id][1]
                         + cocg[cell_id][1][2] *rhsv[cell_id][2]
                         + f_ext[cell_id][1];
      grad[cell_id][2] =   cocg[cell_id][2][0] *rhsv[cell_id][0]
                         + cocg[cell_id][2][1] *rhsv[cell_id][1]
                         + cocg[cell_id][2][2] *rhsv[cell_id][2]
                         + f_ext[cell_id][2];
    }

  }
  else {

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      grad[cell_id][0] =   cocg[cell_id][0][0] *rhsv[cell_id][0]
                         + cocg[cell_id][0][1] *rhsv[cell_id][1]
                         + cocg[cell_id][0][2] *rhsv[cell_id][2];
      grad[cell_id][1] =   cocg[cell_id][1][0] *rhsv[cell_id][0]
                         + cocg[cell_id][1][1] *rhsv[cell_id][1]
                         + cocg[cell_id][1][2] *rhsv[cell_id][2];
      grad[cell_id][2] =   cocg[cell_id][2][0] *rhsv[cell_id][0]
                         + cocg[cell_id][2][1] *rhsv[cell_id][1]
                         + cocg[cell_id][2][2] *rhsv[cell_id][2];
    }

  }

  /* Synchronize halos */

  _sync_scalar_gradient_halo(m, CS_HALO_STANDARD, idimtr, grad);

   if (c_weight != NULL) {
     if (w_stride == 6) {
       BFT_FREE(_cocgb);
       BFT_FREE(_cocg);
     }
   }
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a vector if necessary. This function deals with the
 * standard or extended neighborhood.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   clip_mode      <-- type of clipping for the computation of the gradient
 *   verbosity      <-- output level
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   pvar           <-- variable
 *   gradv          <-> gradient of pvar (du_i/dx_j : gradv[][i][j])
 *   pvar           <-- variable
 *----------------------------------------------------------------------------*/

static void
_vector_gradient_clipping(const cs_mesh_t              *m,
                          const cs_mesh_quantities_t   *fvq,
                          cs_halo_type_t                halo_type,
                          int                           clip_mode,
                          int                           verbosity,
                          cs_real_t                     climgp,
                          const cs_real_3_t   *restrict pvar,
                          cs_real_33_t        *restrict gradv)
{
  int        g_id, t_id;
  cs_gnum_t  t_n_clip;
  cs_lnum_t  cell_id, cell_id1, cell_id2, face_id, i, j;
  cs_real_3_t dist, grad_dist1, grad_dist2;
  cs_real_t  dvar_sq, dist_sq1, dist_sq2;
  cs_real_t  global_min_factor, global_max_factor, factor1, factor2;
  cs_real_t  t_max_factor, t_min_factor;

  cs_gnum_t  n_clip = 0, n_g_clip =0;
  cs_real_t  min_factor = 1;
  cs_real_t  max_factor = 0;
  cs_real_t  clipp_coef_sq = climgp*climgp;
  cs_real_t  *restrict buf = NULL, *restrict clip_factor = NULL;
  cs_real_t  *restrict denom = NULL, *restrict denum = NULL;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  const cs_halo_t *halo = m->halo;

  if (clip_mode < 0)
    return;

  /* The gradient and the variable must be already synchronized */

  /* Allocate and initialize working buffers */

  if (clip_mode == 1)
    BFT_MALLOC(buf, 3*n_cells_ext, cs_real_t);
  else
    BFT_MALLOC(buf, 2*n_cells_ext, cs_real_t);

  denum = buf;
  denom = buf + n_cells_ext;

  if (clip_mode == 1)
    clip_factor = buf + 2*n_cells_ext;

  /* Initialization */

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    denum[cell_id] = 0;
    denom[cell_id] = 0;
    if (clip_mode == 1)
      clip_factor[cell_id] = (cs_real_t)DBL_MAX;
  }

  /* Remark:
     denum: holds the maximum l2 norm of the variation of the gradient squared
     denom: holds the maximum l2 norm of the variation of the variable squared */

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clip_mode == 0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(face_id, cell_id1, cell_id2, i, \
                                      dist, grad_dist1, grad_dist2, \
                                      dist_sq1, dist_sq2, dvar_sq)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cell_id1 = i_face_cells[face_id][0];
          cell_id2 = i_face_cells[face_id][1];

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++) {

            grad_dist1[i] =   gradv[cell_id1][i][0] * dist[0]
                            + gradv[cell_id1][i][1] * dist[1]
                            + gradv[cell_id1][i][2] * dist[2];

            grad_dist2[i] =   gradv[cell_id2][i][0] * dist[0]
                            + gradv[cell_id2][i][1] * dist[1]
                            + gradv[cell_id2][i][2] * dist[2];

          }

          dist_sq1 =   grad_dist1[0]*grad_dist1[0]
                     + grad_dist1[1]*grad_dist1[1]
                     + grad_dist1[2]*grad_dist1[2];

          dist_sq2 =   grad_dist2[0]*grad_dist2[0]
                     + grad_dist2[1]*grad_dist2[1]
                     + grad_dist2[2]*grad_dist2[2];

          dvar_sq =     (pvar[cell_id1][0]-pvar[cell_id2][0])
                      * (pvar[cell_id1][0]-pvar[cell_id2][0])
                    +   (pvar[cell_id1][1]-pvar[cell_id2][1])
                      * (pvar[cell_id1][1]-pvar[cell_id2][1])
                    +   (pvar[cell_id1][2]-pvar[cell_id2][2])
                      * (pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denum[cell_id2] = CS_MAX(denum[cell_id2], dist_sq2);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);
          denom[cell_id2] = CS_MAX(denom[cell_id2], dvar_sq);

        } /* End of loop on faces */

      } /* End of loop on threads */

    } /* End of loop on thread groups */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for private(cell_id2, i, dist, \
                                      grad_dist1, dist_sq1, dvar_sq)
      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
        for (cs_lnum_t cidx = cell_cells_idx[cell_id1];
             cidx < cell_cells_idx[cell_id1+1];
             cidx++) {

          cell_id2 = cell_cells_lst[cidx];

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++)
            grad_dist1[i] =   gradv[cell_id1][i][0] * dist[0]
                            + gradv[cell_id1][i][1] * dist[1]
                            + gradv[cell_id1][i][2] * dist[2];


          dist_sq1 =   grad_dist1[0]*grad_dist1[0]
                     + grad_dist1[1]*grad_dist1[1]
                     + grad_dist1[2]*grad_dist1[2];

          dvar_sq =     (pvar[cell_id1][0]-pvar[cell_id2][0])
                      * (pvar[cell_id1][0]-pvar[cell_id2][0])
                    +   (pvar[cell_id1][1]-pvar[cell_id2][1])
                      * (pvar[cell_id1][1]-pvar[cell_id2][1])
                    +   (pvar[cell_id1][2]-pvar[cell_id2][2])
                      * (pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);

        }
      }

    } /* End for extended halo */

  }

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clip_mode == 1) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(face_id, cell_id1, cell_id2, i, \
                                      dist, grad_dist1, dist_sq1, dvar_sq)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cell_id1 = i_face_cells[face_id][0];
          cell_id2 = i_face_cells[face_id][1];

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++)
            grad_dist1[i]
              = 0.5 * (  (gradv[cell_id1][i][0]+gradv[cell_id2][i][0])*dist[0]
                       + (gradv[cell_id1][i][1]+gradv[cell_id2][i][1])*dist[1]
                       + (gradv[cell_id1][i][2]+gradv[cell_id2][i][2])*dist[2]);

          dist_sq1 =   grad_dist1[0]*grad_dist1[0]
                     + grad_dist1[1]*grad_dist1[1]
                     + grad_dist1[2]*grad_dist1[2];

          dvar_sq =     (pvar[cell_id1][0]-pvar[cell_id2][0])
                      * (pvar[cell_id1][0]-pvar[cell_id2][0])
                    +   (pvar[cell_id1][1]-pvar[cell_id2][1])
                      * (pvar[cell_id1][1]-pvar[cell_id2][1])
                    +   (pvar[cell_id1][2]-pvar[cell_id2][2])
                      * (pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denum[cell_id2] = CS_MAX(denum[cell_id2], dist_sq1);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);
          denom[cell_id2] = CS_MAX(denom[cell_id2], dvar_sq);

        } /* End of loop on threads */

      } /* End of loop on thread groups */

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for private(cell_id2, i, dist, \
                                      grad_dist1, dist_sq1, dvar_sq)
      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
        for (cs_lnum_t cidx = cell_cells_idx[cell_id1];
             cidx < cell_cells_idx[cell_id1+1];
             cidx++) {

          cell_id2 = cell_cells_lst[cidx];

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++)
            grad_dist1[i]
              = 0.5 * (  (gradv[cell_id1][i][0]+gradv[cell_id2][i][0])*dist[0]
                       + (gradv[cell_id1][i][1]+gradv[cell_id2][i][1])*dist[1]
                       + (gradv[cell_id1][i][2]+gradv[cell_id2][i][2])*dist[2]);

          dist_sq1 =   grad_dist1[0]*grad_dist1[0]
                     + grad_dist1[1]*grad_dist1[1]
                     + grad_dist1[2]*grad_dist1[2];

          dvar_sq =     (pvar[cell_id1][0]-pvar[cell_id2][0])
                      * (pvar[cell_id1][0]-pvar[cell_id2][0])
                    +   (pvar[cell_id1][1]-pvar[cell_id2][1])
                      * (pvar[cell_id1][1]-pvar[cell_id2][1])
                    +   (pvar[cell_id1][2]-pvar[cell_id2][2])
                      * (pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);

        }
      }

    } /* End for extended neighborhood */

    /* Synchronize variable */

    if (halo != NULL) {
      cs_halo_sync_var(m->halo, halo_type, denom);
      cs_halo_sync_var(m->halo, halo_type, denum);
    }

  } /* End if clip_mode == 1 */

  /* Clipping of the gradient if denum/denom > climgp**2 */

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clip_mode == 0) {

#   pragma omp parallel private(t_min_factor, t_max_factor, t_n_clip, \
                                factor1, i, j)
    {
      t_n_clip = 0;
      t_min_factor = min_factor; t_max_factor = max_factor;

#     pragma omp for
      for (cell_id = 0; cell_id < n_cells; cell_id++) {

        if (denum[cell_id] > clipp_coef_sq * denom[cell_id]) {

          factor1 = sqrt(clipp_coef_sq * denom[cell_id]/denum[cell_id]);

          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++)
              gradv[cell_id][i][j] *= factor1;
          }

          t_min_factor = CS_MIN(factor1, t_min_factor);
          t_max_factor = CS_MAX(factor1, t_max_factor);
          t_n_clip++;

        } /* If clipping */

      } /* End of loop on cells */

#     pragma omp critical
      {
        min_factor = CS_MIN(min_factor, t_min_factor);
        max_factor = CS_MAX(max_factor, t_max_factor);
        n_clip += t_n_clip;
      }
    } /* End of omp parallel construct */

  }

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clip_mode == 1) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(face_id, cell_id1, cell_id2, \
                                      factor1, factor2, min_factor)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cell_id1 = i_face_cells[face_id][0];
          cell_id2 = i_face_cells[face_id][1];

          factor1 = 1.0;
          if (denum[cell_id1] > clipp_coef_sq * denom[cell_id1])
            factor1 = sqrt(clipp_coef_sq * denom[cell_id1]/denum[cell_id1]);

          factor2 = 1.0;
          if (denum[cell_id2] > clipp_coef_sq * denom[cell_id2])
            factor2 = sqrt(clipp_coef_sq * denom[cell_id2]/denum[cell_id2]);

          min_factor = CS_MIN(factor1, factor2);

          clip_factor[cell_id1] = CS_MIN(clip_factor[cell_id1], min_factor);
          clip_factor[cell_id2] = CS_MIN(clip_factor[cell_id2], min_factor);

        } /* End of loop on faces */

      } /* End of loop on threads */

    } /* End of loop on thread groups */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

#     pragma omp parallel for private(cell_id2, min_factor, factor2)
      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {

        min_factor = 1.0;

        for (cs_lnum_t cidx = cell_cells_idx[cell_id1];
             cidx < cell_cells_idx[cell_id1+1];
             cidx++) {

          cell_id2 = cell_cells_lst[cidx];
          factor2 = 1.0;

          if (denum[cell_id2] > clipp_coef_sq * denom[cell_id2])
            factor2 = sqrt(clipp_coef_sq * denom[cell_id2]/denum[cell_id2]);

          min_factor = CS_MIN(min_factor, factor2);

        }

        clip_factor[cell_id1] = CS_MIN(clip_factor[cell_id1], min_factor);

      } /* End of loop on cells */

    } /* End for extended neighborhood */

#   pragma omp parallel private(t_min_factor, t_max_factor, factor1, \
                                t_n_clip, i, j)
    {
      t_n_clip = 0;
      t_min_factor = min_factor; t_max_factor = max_factor;

#     pragma omp for
      for (cell_id = 0; cell_id < n_cells; cell_id++) {

        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++)
            gradv[cell_id][i][j] *= clip_factor[cell_id];
        }

        if (clip_factor[cell_id] < 0.99) {
          t_max_factor = CS_MAX(t_max_factor, clip_factor[cell_id]);
          t_min_factor = CS_MIN(t_min_factor, clip_factor[cell_id]);
          t_n_clip++;
        }

      } /* End of loop on cells */

#     pragma omp critical
      {
        min_factor = CS_MIN(min_factor, t_min_factor);
        max_factor = CS_MAX(max_factor, t_max_factor);
        n_clip += t_n_clip;
      }
    } /* End of omp parallel construct */

  } /* End if clip_mode == 1 */

  /* Update min/max and n_clip in case of parallelism */
  /*--------------------------------------------------*/

#if defined(HAVE_MPI)

  if (m->n_domains > 1) {

    assert(sizeof(cs_real_t) == sizeof(double));

    /* Global Max */

    MPI_Allreduce(&max_factor, &global_max_factor, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    max_factor = global_max_factor;

    /* Global min */

    MPI_Allreduce(&min_factor, &global_min_factor, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    min_factor = global_min_factor;

    /* Sum number of clippings */

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_GNUM,
                  MPI_SUM, cs_glob_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(HAVE_MPI) */

  /* Output warning if necessary */

  if (verbosity > 1)
    bft_printf(_(" Gradient of a vector limitation in %llu cells\n"
                 "   minimum factor = %14.5e; maximum factor = %14.5e\n"),
               (unsigned long long)n_clip, min_factor, max_factor);

  /* Synchronize the updated Gradient */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)gradv, 9);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)gradv);
  }

  BFT_FREE(buf);
}

/*----------------------------------------------------------------------------
 * Initialize the gradient of a vector for gradient reconstruction.
 *
 * A non-reconstructed gradient is computed at this stage.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_initialize_vector_gradient(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            const cs_internal_coupling_t *cpl,
                            cs_halo_type_t                halo_type,
                            int                           inc,
                            const cs_real_3_t   *restrict coefav,
                            const cs_real_33_t  *restrict coefbv,
                            const cs_real_3_t   *restrict pvar,
                            const cs_real_t     *restrict c_weight,
                            cs_real_33_t        *restrict grad)
{
  int g_id, t_id;
  cs_lnum_t face_id;
  cs_real_t pond;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Computation without reconstruction */
  /*------------------------------------*/

  /* Initialization */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        grad[cell_id][i][j] = 0.0;
    }
  }

  /* Interior faces contribution */

  for (g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for private(face_id, pond)
    for (t_id = 0; t_id < n_i_threads; t_id++) {

      for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        pond = weight[face_id];

        cs_real_t ktpond = (c_weight == NULL) ?
          pond :                    // no cell weightening
          pond * c_weight[cell_id1] // cell weightening active
            / (      pond * c_weight[cell_id1]
              + (1.0-pond)* c_weight[cell_id2]);

        /*
           Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                    + (1-\alpha_\ij) \varia_\cellj\f$
                   but for the cell \f$ \celli \f$ we remove
                   \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                   and for the cell \f$ \cellj \f$ we remove
                   \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
        */
        for (int i = 0; i < 3; i++) {
          cs_real_t pfaci = (1.0-ktpond) * (pvar[cell_id2][i] - pvar[cell_id1][i]);
          cs_real_t pfacj = - ktpond * (pvar[cell_id2][i] - pvar[cell_id1][i]);

          for (int j = 0; j < 3; j++) {
            grad[cell_id1][i][j] += pfaci * i_f_face_normal[face_id][j];
            grad[cell_id2][i][j] -= pfacj * i_f_face_normal[face_id][j];
          }
        }

      } /* End of loop on faces */

    } /* End of loop on threads */

  } /* End of loop on thread groups */

  /* Contribution from coupled faces */
  if (cpl != NULL)
    cs_internal_coupling_initialize_vector_gradient
      (cpl, c_weight, pvar, grad);

  /* Boundary face treatment */

  for (g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for private(face_id)
    for (t_id = 0; t_id < n_b_threads; t_id++) {

      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        if (cpl == NULL || !coupled_faces[face_id]) {

          cs_lnum_t cell_id = b_face_cells[face_id];


          /*
             Remark: for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
           */
          for (int i = 0; i < 3; i++) {
            cs_real_t pfac = inc*coefav[face_id][i];

            for (int k = 0; k < 3; k++) {
              if (i == k)
                pfac += (coefbv[face_id][i][k] - 1.0) * pvar[cell_id][k];
              else
                pfac += coefbv[face_id][i][k] * pvar[cell_id][k];
            }

            for (int j = 0; j < 3; j++)
              grad[cell_id][i][j] += pfac * b_f_face_normal[face_id][j];
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        grad[cell_id][i][j] *= dvol;
    }
  }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 9);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
  }
}

/*----------------------------------------------------------------------------
 * Reconstruct the gradient of a vector using a given gradient of
 * this vector (typically lsq).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   r_grad         --> gradient used for reconstruction
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_reconstruct_vector_gradient(const cs_mesh_t              *m,
                             const cs_mesh_quantities_t   *fvq,
                             const cs_internal_coupling_t *cpl,
                             cs_halo_type_t                halo_type,
                             const cs_real_33_t  *restrict coefbv,
                             cs_real_33_t        *restrict r_grad,
                             cs_real_33_t        *restrict grad)
{
  int g_id, t_id;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;
  const cs_real_33_t *restrict corr_grad_lin
    = (const cs_real_33_t *restrict)fvq->corr_grad_lin;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Initialize gradient */
  /*---------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        grad[cell_id][i][j] = grad[cell_id][i][j] * cell_f_vol[cell_id];
    }
  }

  /* Interior faces contribution */

  for (g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        for (int i = 0; i < 3; i++) {
          /* Reconstruction part */
          cs_real_t
          rfac = 0.5 *
              ( dofij[face_id][0]*( r_grad[cell_id1][i][0]
                                   +r_grad[cell_id2][i][0])
               +dofij[face_id][1]*( r_grad[cell_id1][i][1]
                                   +r_grad[cell_id2][i][1])
               +dofij[face_id][2]*( r_grad[cell_id1][i][2]
                                   +r_grad[cell_id2][i][2]));

          for (int j = 0; j < 3; j++) {
            grad[cell_id1][i][j] += rfac * i_f_face_normal[face_id][j];
            grad[cell_id2][i][j] -= rfac * i_f_face_normal[face_id][j];
          }
        }

      } /* End of loop on faces */

    } /* End of loop on threads */

  } /* End of loop on thread groups */

  /* Contribution from coupled faces */
  if (cpl != NULL)
    cs_internal_coupling_reconstruct_vector_gradient
      (cpl, r_grad, grad);

  /* Boundary face treatment */

  for (g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        if (cpl == NULL || !coupled_faces[face_id]) {
          cs_lnum_t cell_id = b_face_cells[face_id];

          /* Reconstruction part */
          for (int i = 0; i < 3; i++) {
            cs_real_t rfac = 0.;
            for (int k = 0; k < 3; k++) {
              cs_real_t vecfac = grad[cell_id][k][0] * diipb[face_id][0]
                               + grad[cell_id][k][1] * diipb[face_id][1]
                               + grad[cell_id][k][2] * diipb[face_id][2];
              rfac += coefbv[face_id][i][k] * vecfac;
            }

            for (int j = 0; j < 3; j++)
              grad[cell_id][i][j] += rfac * b_f_face_normal[face_id][j];
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        grad[cell_id][i][j] *= dvol;
    }

    if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_WARPED_CORRECTION) {
      cs_real_3_t gradpa;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          gradpa[j] = grad[cell_id][i][j];
          grad[cell_id][i][j] = 0.;
        }

        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            grad[cell_id][i][j] += corr_grad_lin[cell_id][j][k] * gradpa[k];
      }
    }
  }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 9);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
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
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   var_name       <-- variable's name
 *   gradient_info  <-- pointer to performance logging structure, or NULL
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_r_sweeps     --> >1: with reconstruction
 *   verbosity      --> verbosity level
 *   epsrgp         --> precision for iterative gradient calculation
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   grad           <-> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_iterative_vector_gradient(const cs_mesh_t              *m,
                           const cs_mesh_quantities_t   *fvq,
                           const cs_internal_coupling_t *cpl,
                           const char                   *var_name,
                           cs_gradient_info_t           *gradient_info,
                           cs_halo_type_t                halo_type,
                           int                           inc,
                           int                           n_r_sweeps,
                           int                           verbosity,
                           cs_real_t                     epsrgp,
                           const cs_real_3_t   *restrict coefav,
                           const cs_real_33_t  *restrict coefbv,
                           const cs_real_3_t   *restrict pvar,
                           const cs_real_t              *c_weight,
                           cs_real_33_t        *restrict grad)
{
  int isweep = 0;
  int g_id, t_id;
  cs_lnum_t face_id;
  cs_real_t l2_norm, l2_residual, vecfac, pond;

  cs_real_33_t *rhs;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  cs_real_33_t *restrict cocg = (cpl == NULL) ?
    fvq->cocg_it : cpl->cocg_it;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_33_t);

  /* Gradient reconstruction to handle non-orthogonal meshes */
  /*---------------------------------------------------------*/

  /* L2 norm */

  l2_norm = _l2_norm_1(9*n_cells, (cs_real_t *)grad);
  l2_residual = l2_norm;

  if (l2_norm > cs_math_epzero) {

    /* Iterative process */
    /*-------------------*/

    for (isweep = 1; isweep < n_r_sweeps && l2_residual > epsrgp*l2_norm; isweep++) {

      /* Computation of the Right Hand Side*/

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++)
            rhs[cell_id][i][j] = -grad[cell_id][i][j] * cell_f_vol[cell_id];
        }
      }

      /* Interior face treatment */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, pond)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id1 = i_face_cells[face_id][0];
            cs_lnum_t cell_id2 = i_face_cells[face_id][1];
            pond = weight[face_id];

            cs_real_t ktpond = (c_weight == NULL) ?
              pond :                     // no cell weightening
              pond  * c_weight[cell_id1] // cell weightening active
                / (      pond  * c_weight[cell_id1]
                  + (1.0-pond) * c_weight[cell_id2]);

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            for (int i = 0; i < 3; i++) {

              /* Reconstruction part */
              cs_real_t
              pfaci = 0.5 * ( ( grad[cell_id1][i][0] + grad[cell_id2][i][0])
                              * dofij[face_id][0]
                            + ( grad[cell_id1][i][1] + grad[cell_id2][i][1])
                              * dofij[face_id][1]
                            + ( grad[cell_id1][i][2] + grad[cell_id2][i][2])
                              * dofij[face_id][2]
                            );
              cs_real_t pfacj = pfaci;

              pfaci += (1.0-ktpond) * (pvar[cell_id2][i] - pvar[cell_id1][i]);
              pfacj -=      ktpond  * (pvar[cell_id2][i] - pvar[cell_id1][i]);

              for (int j = 0; j < 3; j++) {
                rhs[cell_id1][i][j] += pfaci * i_f_face_normal[face_id][j];
                rhs[cell_id2][i][j] -= pfacj * i_f_face_normal[face_id][j];
              }
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from coupled faces */
      if (cpl != NULL)
        cs_internal_coupling_iterative_vector_gradient
          (cpl, c_weight, grad, pvar, rhs);


      /* Boundary face treatment */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id, vecfac)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            if (cpl == NULL || !coupled_faces[face_id]) {

              cs_lnum_t cell_id = b_face_cells[face_id];

              for (int i = 0; i < 3; i++) {

                /*
                   Remark: for the cell \f$ \celli \f$ we remove
                           \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                 */

                cs_real_t pfac = inc*coefav[face_id][i];

                for (int k = 0; k < 3; k++) {
                  /* Reconstruction part */
                  vecfac = grad[cell_id][k][0] * diipb[face_id][0]
                         + grad[cell_id][k][1] * diipb[face_id][1]
                         + grad[cell_id][k][2] * diipb[face_id][2];
                  pfac += coefbv[face_id][i][k] * vecfac;

                  if (i == k)
                    pfac += (coefbv[face_id][i][k] - 1.0) * pvar[cell_id][k];
                  else
                    pfac += coefbv[face_id][i][k] * pvar[cell_id][k];
                }

                for (int j = 0; j < 3; j++)
                  rhs[cell_id][i][j] += pfac * b_f_face_normal[face_id][j];

              }
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Increment of the gradient */

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_real_t dvol;
        /* Is the cell disabled (for solid or porous)? Not the case if coupled */
        if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
          dvol = 1. / cell_f_vol[cell_id];
        else
          dvol = 0.;

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++)
            rhs[cell_id][i][j] *= dvol;
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++)
              grad[cell_id][i][j] += rhs[cell_id][i][k] * cocg[cell_id][k][j];
          }
        }
      }

      /* Periodicity and parallelism treatment */

      if (m->halo != NULL) {
        cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 9);
        if (cs_glob_mesh->n_init_perio > 0)
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

  if (gradient_info != NULL)
    _gradient_info_update_iter(gradient_info, isweep);

  BFT_FREE(rhs);
}

/*----------------------------------------------------------------------------
 * Initialize cocg for lsq vector gradient.
 *
 * parameters:
 *   cell_id          <-- cell id
 *   madj             <-- pointer to mesh adjacencies structure
 *   fvq              <-- pointer to associated finite volume quantities
 *   cocg             --> cocg
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 2 /* Bug with O3 or above with icc 18.0.1 20171018
                                at least on Xeon(R) Gold 6140 */
#endif

static void
_init_cocg_lsq(cs_lnum_t                     cell_id,
               const cs_mesh_adjacencies_t  *madj,
               const cs_mesh_quantities_t   *fvq,
               cs_real_t                     cocg[restrict 3][3])
{
  /* Short variable accesses */

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;

  cs_lnum_t s_id, e_id;
  cs_real_t dc[3];

  /* initialize cocg */

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++)
      cocg[ii][jj] = 0.;
  }

  /* Contribution from cell neighbors and extended cell neighbors */

  for (int adj_id = 0; adj_id < 2; adj_id++) {

    const cs_lnum_t   *restrict cell_cells;

    if (adj_id == 0) {
      s_id = madj->cell_cells_idx[cell_id];
      e_id = madj->cell_cells_idx[cell_id+1];
      cell_cells = (const cs_lnum_t *restrict)(madj->cell_cells);
    }
    else if (madj->cell_cells_e_idx != NULL) {
      s_id = madj->cell_cells_e_idx[cell_id];
      e_id = madj->cell_cells_e_idx[cell_id+1];
      cell_cells = (const cs_lnum_t *restrict)(madj->cell_cells_e);
    }
    else
      break;

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t cell_id1 = cell_cells[i];

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        dc[ii] = cell_cen[cell_id1][ii] - cell_cen[cell_id][ii];

      cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++)
          cocg[ii][jj] += dc[ii]*dc[jj]*ddc;
      }

    }

  }

  /* Contribution from boundary faces */

  const cs_lnum_t  *restrict cell_b_faces
    = (const cs_lnum_t  *restrict)(madj->cell_b_faces);

  s_id = madj->cell_b_faces_idx[cell_id];
  e_id = madj->cell_b_faces_idx[cell_id+1];

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    cs_lnum_t face_id = cell_b_faces[i];

    cs_real_t udbfs = 1. / fvq->b_face_surf[face_id];

    for (int ii = 0; ii < 3; ii++)
      dc[ii] = udbfs * b_face_normal[face_id][ii];

    for (int ii = 0; ii < 3; ii++) {
      for (int jj = 0; jj < 3; jj++)
        cocg[ii][jj] += dc[ii]*dc[jj];
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute cocg and RHS at boundaries for lsq vector gradient.
 *
 * parameters:
 *   cell_id          <-- cell id
 *   inc              <-- if 0, solve on increment; 1 otherwise
 *   madj             <-- pointer to mesh adjacencies structure
 *   fvq              <-- pointer to associated finite volume quantities
 *   _33_9_idx        <-- symmetric indexes mapping
 *   pvar             <-- variable
 *   coefav           <-- B.C. coefficients for boundary face normals
 *   coefbv           <-- B.C. coefficients for boundary face normals
 *   rhs              <-- right hand side
 *   cocgb_v          --> boundary cocg vector values
 *   rhsb_v           --> boundary RHS values
 *----------------------------------------------------------------------------*/

static void
_compute_cocgb_rhsb_lsq_v(cs_lnum_t                     cell_id,
                          const int                     inc,
                          const cs_mesh_adjacencies_t  *madj,
                          const cs_mesh_quantities_t   *fvq,
                          cs_lnum_t              _33_9_idx[const restrict 9][2],
                          const cs_real_3_t            *restrict pvar,
                          const cs_real_3_t            *restrict coefav,
                          const cs_real_33_t           *restrict coefbv,
                          const cs_real_33_t           *rhs,
                          cs_real_t                     cocgb_v[restrict 45],
                          cs_real_t                     rhsb_v[restrict 9])
{
  /* Short variable accesses */

  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)fvq->b_dist;

  cs_lnum_t s_id, e_id;
  cs_real_t cocg[3][3];

  /* initialize cocg and rhsfor lsq vector gradient */

  _init_cocg_lsq(cell_id, madj, fvq, cocg);

  for (int ll = 0; ll < 9; ll++) {

    /* index of row first coefficient */
    int ll_9 = ll*(ll+1)/2;

    for (int mm = 0; mm <= ll; mm++) {
      /* initialize */
      cocgb_v[ll_9+mm] = 0.;

      /* contribution of t[kk][qq] */
      int pp = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      /* derivative with respect to t[rr][ss] */
      int rr = _33_9_idx[mm][0];
      int ss = _33_9_idx[mm][1];

      /* part from cocg_s (BCs independant) */
      if (pp == rr)
        cocgb_v[ll_9+mm] = cocg[qq][ss];

      /* part already computed from rhs */
      rhsb_v[ll] = rhs[cell_id][pp][qq];
    }
  }

  s_id = madj->cell_b_faces_idx[cell_id];
  e_id = madj->cell_b_faces_idx[cell_id+1];

  const cs_lnum_t   *restrict cell_b_faces
    = (const cs_lnum_t *restrict)(madj->cell_b_faces);

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    cs_lnum_t face_id = cell_b_faces[i];

    /* build cocgb_v matrix */

    cs_real_t udbfs = 1. / b_face_surf[face_id];
    const cs_real_t *restrict iipbf = diipb[face_id];

    /* db = I'F / ||I'F|| */
    cs_real_t  nb[3];
    for (int ii = 0; ii < 3; ii++)
      nb[ii] = udbfs * b_face_normal[face_id][ii];

    cs_real_t db = 1./b_dist[face_id];
    cs_real_t db2 = db*db;

    /* A and (B - I) */
    cs_real_t a[3];
    cs_real_t bt[3][3];
    for (int ll = 0; ll < 3; ll++) {
      for (int pp = 0; pp < 3; pp++)
        bt[ll][pp] = coefbv[face_id][ll][pp];
    }
    for (int ll = 0; ll < 3; ll++) {
      a[ll] = inc*coefav[face_id][ll];
      bt[ll][ll] -= 1;
    }

    /* cocgb */

    for (int ll = 0; ll < 9; ll++) {

      /* contribution of t[kk][qq] */
      int kk = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      int ll_9 = ll*(ll+1)/2;
      for (int pp = 0; pp <= ll; pp++) {

        /* derivative with respect to t[rr][ss] */
        int rr = _33_9_idx[pp][0];
        int ss = _33_9_idx[pp][1];

        /* part from derivative of 1/2*|| B*t*IIp/db ||^2 */
        cs_real_t cocgv = 0.;
        for (int mm = 0; mm < 3; mm++)
          cocgv += bt[mm][kk]*bt[mm][rr];
        cocgb_v[ll_9+pp] += cocgv*(iipbf[qq]*iipbf[ss])*db2;

        /* part from derivative of -< t*nb , B*t*IIp/db > */
        cocgb_v[ll_9+pp] -= (  nb[ss]*bt[rr][kk]*iipbf[qq]
                             + nb[qq]*bt[kk][rr]*iipbf[ss])
                             *db;
      }
    }

    /* rhsb */

    for (int ll = 0; ll < 9; ll++) {
      int pp = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      /* part from derivative of < (B-1)*t*IIp/db , (A+(B-1)*v)/db > */
      cs_real_t rhsv = 0.;
      for (int rr = 0; rr < 3; rr++) {
        rhsv +=   bt[rr][pp]*diipb[face_id][qq]
                            *(a[rr]+ bt[rr][0]*pvar[cell_id][0]
                                   + bt[rr][1]*pvar[cell_id][1]
                                   + bt[rr][2]*pvar[cell_id][2]);
      }

      rhsb_v[ll] -= rhsv*db2;
    }

  }

  /* Crout factorization of 9x9 symmetric cocg at boundaries */

  _fact_crout_pp(9, cocgb_v);
}

/*----------------------------------------------------------------------------
 * Compute cocg and RHS at boundaries for lsq tensor gradient.
 *
 * parameters:
 *   cell_id          <-- cell id
 *   inc              <-- if 0, solve on increment; 1 otherwise
 *   madj             <-- pointer to mesh adjacencies structure
 *   fvq              <-- pointer to associated finite volume quantities
 *   _63_18_idx       <-- symmetric indexes mapping
 *   pvar             <-- variable
 *   coefat           <-- B.C. coefficients for boundary face normals
 *   coefbt           <-- B.C. coefficients for boundary face normals
 *   rhs              <-- right hand side
 *   cocgb_t          --> boundary cocg vector values
 *   rhsb_t           --> boundary RHS values
 *----------------------------------------------------------------------------*/

static void
_compute_cocgb_rhsb_lsq_t(cs_lnum_t                     cell_id,
                          const int                     inc,
                          const cs_mesh_adjacencies_t  *madj,
                          const cs_mesh_quantities_t   *fvq,
                          cs_lnum_t            _63_18_idx[const restrict 18][2],
                          const cs_real_6_t            *restrict pvar,
                          const cs_real_6_t            *restrict coefat,
                          const cs_real_66_t           *restrict coefbt,
                          const cs_real_63_t           *rhs,
                          cs_real_t                     cocgb_t[restrict 171],
                          cs_real_t                     rhsb_t[restrict 18])
{
  /* Short variable accesses */

  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)fvq->b_dist;

  cs_lnum_t s_id, e_id;
  cs_real_t cocg[3][3];

  /* initialize cocg and rhs for lsq tensor gradient */

  _init_cocg_lsq(cell_id, madj, fvq, cocg);

  for (int ll = 0; ll < 18; ll++) {

    /* index of row first coefficient */
    int ll_18 = ll*(ll+1)/2;

    for (int mm = 0; mm <= ll; mm++) {
      /* initialize */
      cocgb_t[ll_18+mm] = 0.;

      /* contribution of t[kk][qq] */
      int pp = _63_18_idx[ll][0];
      int qq = _63_18_idx[ll][1];

      /* derivative with respect to t[rr][ss] */
      int rr = _63_18_idx[mm][0];
      int ss = _63_18_idx[mm][1];

      /* part from cocg (BCs independant) */
      if (pp == rr)
        cocgb_t[ll_18+mm] = cocg[qq][ss];

      /* part already computed from rhs */
      rhsb_t[ll] = rhs[cell_id][pp][qq];
    }
  }

  s_id = madj->cell_b_faces_idx[cell_id];
  e_id = madj->cell_b_faces_idx[cell_id+1];

  const cs_lnum_t   *restrict cell_b_faces
    = (const cs_lnum_t *restrict)(madj->cell_b_faces);

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    cs_lnum_t face_id = cell_b_faces[i];

    /* build cocgb_v matrix */

    cs_real_t udbfs = 1. / b_face_surf[face_id];
    const cs_real_t *restrict iipbf = diipb[face_id];

    /* db = I'F / ||I'F|| */
    cs_real_t  nb[3];
    for (int ii = 0; ii < 3; ii++)
      nb[ii] = udbfs * b_face_normal[face_id][ii];

    cs_real_t db = 1./b_dist[face_id];
    cs_real_t db2 = db*db;

    /* A and (B - I) */
    cs_real_t a[6];
    cs_real_t bt[6][6];
    for (int ll = 0; ll < 6; ll++) {
      for (int pp = 0; pp < 6; pp++)
        bt[ll][pp] = coefbt[face_id][ll][pp];
    }
    for (int ll = 0; ll < 6; ll++) {
      a[ll] = inc*coefat[face_id][ll];
      bt[ll][ll] -= 1;
    }

    /* cocgb */

    for (int ll = 0; ll < 18; ll++) {

      /* contribution of t[kk][qq] */
      int kk = _63_18_idx[ll][0];
      int qq = _63_18_idx[ll][1];

      int ll_18 = ll*(ll+1)/2;
      for (int pp = 0; pp <= ll; pp++) {

        /* derivative with respect to t[rr][ss] */
        int rr = _63_18_idx[pp][0];
        int ss = _63_18_idx[pp][1];

        /* part from derivative of 1/2*|| B*t*IIp/db ||^2 */
        cs_real_t cocgt = 0.;
        for (int mm = 0; mm < 6; mm++)
          cocgt += bt[mm][kk]*bt[mm][rr];
        cocgb_t[ll_18+pp] += cocgt*(iipbf[qq]*iipbf[ss])*db2;

        /* part from derivative of -< t*nb , B*t*IIp/db > */
        cocgb_t[ll_18+pp] -= (  nb[ss]*bt[rr][kk]*iipbf[qq]
                              + nb[qq]*bt[kk][rr]*iipbf[ss])
                              *db;
      }
    }

    /* rhsb */

    for (int ll = 0; ll < 18; ll++) {
      int pp = _63_18_idx[ll][0];
      int qq = _63_18_idx[ll][1];

      /* part from derivative of < (B-1)*t*IIp/db , (A+(B-1)*v)/db > */
      cs_real_t rhst = 0.;
      for (int rr = 0; rr < 6; rr++) {
        cs_real_t tfac = a[rr];
        for (int kk = 0; kk < 6; kk++) {
          tfac += bt[rr][kk]*pvar[cell_id][kk];
        }
        rhst += bt[rr][pp]*diipb[face_id][qq]*tfac;
      }

      rhsb_t[ll] -= rhst*db2;
    }

  }

  /* Crout factorization of 18x18 symmetric cocg at boundaries */

  _fact_crout_pp(18, cocgb_t);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of a vector using least-squares reconstruction for
 * non-orthogonal meshes (n_r_sweeps > 1).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   gradv          --> gradient of pvar (du_i/dx_j : gradv[][i][j])
 *----------------------------------------------------------------------------*/

static void
_lsq_vector_gradient(const cs_mesh_t              *m,
                     const cs_mesh_adjacencies_t  *madj,
                     const cs_mesh_quantities_t   *fvq,
                     const cs_internal_coupling_t *cpl,
                     const cs_halo_type_t          halo_type,
                     const cs_int_t                inc,
                     const cs_real_3_t   *restrict coefav,
                     const cs_real_33_t  *restrict coefbv,
                     const cs_real_3_t   *restrict pvar,
                     const cs_real_t     *restrict c_weight,
                     cs_real_33_t        *restrict gradv)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;

  cs_real_33_t *restrict cocg = fvq->cocg_lsq;//FIXME for internal coupling, has to be recomputed

  cs_lnum_t  cell_id1, cell_id2, i, j, k;
  cs_real_t  pfac, ddc;
  cs_real_3_t  dc;
  cs_real_3_t  fctb;

  cs_real_33_t *rhs;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_33_t);

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)pvar, 3);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)pvar, 3);
  }

  /* Compute Right-Hand Side */
  /*-------------------------*/

# pragma omp parallel for private(i, j)
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        rhs[cell_id][i][j] = 0.0;
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for private(cell_id1, cell_id2,\
                                    i, j, pfac, dc, fctb, ddc)
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cell_id1 = i_face_cells[face_id][0];
        cell_id2 = i_face_cells[face_id][1];

        for (i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        if (c_weight != NULL) {
          cs_real_t pond = weight[face_id];
          cs_real_t denom = 1. / (  pond       *c_weight[cell_id1]
                                  + (1. - pond)*c_weight[cell_id2]);

          for (i = 0; i < 3; i++) {
            pfac =  (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

            for (j = 0; j < 3; j++) {
              fctb[j] = dc[j] * pfac;
              rhs[cell_id1][i][j] += c_weight[cell_id2] * denom * fctb[j];
              rhs[cell_id2][i][j] += c_weight[cell_id1] * denom * fctb[j];
            }
          }
        }
        else {
          for (i = 0; i < 3; i++) {
            pfac =  (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

            for (j = 0; j < 3; j++) {
              fctb[j] = dc[j] * pfac;
              rhs[cell_id1][i][j] += fctb[j];
              rhs[cell_id2][i][j] += fctb[j];
            }
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from extended neighborhood */

  if (halo_type == CS_HALO_EXTENDED) {

#   pragma omp parallel for private(cell_id2, dc, pfac, ddc, i, j)
    for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
      for (cs_lnum_t cidx = cell_cells_idx[cell_id1];
           cidx < cell_cells_idx[cell_id1+1];
           cidx++) {

        cell_id2 = cell_cells_lst[cidx];

        for (i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (i = 0; i < 3; i++) {

          pfac = (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

          for (j = 0; j < 3; j++) {
            rhs[cell_id1][i][j] += dc[j] * pfac;
          }
        }
      }
    }

  } /* End for extended neighborhood */

  /* Contribution from coupled faces */

  if (cpl != NULL)
    cs_internal_coupling_lsq_vector_gradient
      (cpl,
       c_weight,
       1, /* w_stride */
       pvar,
       rhs);

  /* Contribution from boundary faces */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for private(cell_id1, i, j, pfac)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        if (cpl == NULL || !coupled_faces[face_id]) {

          cell_id1 = b_face_cells[face_id];

          cs_real_3_t n_d_dist;
          /* Normal is vector 0 if the b_face_normal norm is too small */
          cs_math_3_normalise(b_face_normal[face_id], n_d_dist);

          cs_real_t d_b_dist = 1. / b_dist[face_id];

          /* Normal divided by b_dist */
          for (i = 0; i < 3; i++)
            n_d_dist[i] *= d_b_dist;

          for (i = 0; i < 3; i++) {
            pfac = (coefav[face_id][i]*inc
                 + ( coefbv[face_id][0][i] * pvar[cell_id1][0]
                   + coefbv[face_id][1][i] * pvar[cell_id1][1]
                   + coefbv[face_id][2][i] * pvar[cell_id1][2]
                   -                         pvar[cell_id1][i]));

            for (j = 0; j < 3; j++)
              rhs[cell_id1][i][j] += n_d_dist[j] * pfac;
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Compute gradient */
  /*------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) {

        gradv[cell_id][i][j] = 0.0;

        for (k = 0; k < 3; k++)
          gradv[cell_id][i][j] += rhs[cell_id][i][k] * cocg[cell_id][k][j];

      }
    }
  }

  /* Compute gradient on boundary cells */
  /*------------------------------------*/

  #pragma omp parallel
  {
    cs_lnum_t t_s_id, t_e_id;
    _thread_range(m->n_b_cells, &t_s_id, &t_e_id);

    /* Build indices bijection between [1-9] and [1-3]*[1-3] */

    cs_lnum_t _33_9_idx[9][2];
    int nn = 0;
    for (int ll = 0; ll < 3; ll++) {
      for (int mm = 0; mm < 3; mm++) {
        _33_9_idx[nn][0] = ll;
        _33_9_idx[nn][1] = mm;
        nn++;
      }
    }

    /* Loop on boundary cells */

    for (cs_lnum_t b_cell_id = t_s_id; b_cell_id < t_e_id; b_cell_id++) {

      cs_lnum_t cell_id = m->b_cells[b_cell_id];

      cs_real_t cocgb_v[45], rhsb_v[9], x[9];

      _compute_cocgb_rhsb_lsq_v
        (cell_id,
         inc,
         madj,
         fvq,
         _33_9_idx,
         (const cs_real_3_t *)pvar,
         (const cs_real_3_t *)coefav,
         (const cs_real_33_t *)coefbv,
         (const cs_real_33_t *)rhs,
         cocgb_v,
         rhsb_v);

      _fw_and_bw_ldtl_pp(cocgb_v,
                         9,
                         x,
                         rhsb_v);

      for (int kk = 0; kk < 9; kk++) {
        int ii = _33_9_idx[kk][0];
        int jj = _33_9_idx[kk][1];
        gradv[cell_id][ii][jj] = x[kk];
      }

    }

  }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)gradv, 9);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)gradv);
  }

  BFT_FREE(rhs);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of a tensor using least-squares reconstruction for
 * non-orthogonal meshes (n_r_sweeps > 1).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefat         <-- B.C. coefficients for boundary face normals
 *   coefbt         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   gradt          --> gradient of pvar (du_i/dx_j : gradv[][i][j])
 *----------------------------------------------------------------------------*/

static void
_lsq_tensor_gradient(const cs_mesh_t              *m,
                     const cs_mesh_adjacencies_t  *madj,
                     const cs_mesh_quantities_t   *fvq,
                     const cs_halo_type_t          halo_type,
                     const cs_int_t                inc,
                     const cs_real_6_t   *restrict coefat,
                     const cs_real_66_t  *restrict coefbt,
                     const cs_real_6_t   *restrict pvar,
                     const cs_real_t     *restrict c_weight,
                     cs_real_63_t        *restrict gradt)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  cs_real_33_t *restrict cocg = fvq->cocg_lsq;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;

  cs_real_63_t *rhs;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_63_t);

  /* Compute Right-Hand Side */
  /*-------------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 3; j++)
        rhs[cell_id][i][j] = 0.0;
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        cs_real_3_t dc, fctb;
        for (int i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        if (c_weight != NULL) {
          cs_real_t pond = weight[face_id];
          cs_real_t denom = 1. / (  pond       *c_weight[cell_id1]
                                  + (1. - pond)*c_weight[cell_id2]);

          for (int i = 0; i < 6; i++) {
            cs_real_t pfac =  (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

            for (int j = 0; j < 3; j++) {
              fctb[j] = dc[j] * pfac;
              rhs[cell_id1][i][j] += c_weight[cell_id2] * denom * fctb[j];
              rhs[cell_id2][i][j] += c_weight[cell_id1] * denom * fctb[j];
            }
          }
        }
        else {
          for (int i = 0; i < 6; i++) {
            cs_real_t pfac =  (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

            for (int j = 0; j < 3; j++) {
              fctb[j] = dc[j] * pfac;
              rhs[cell_id1][i][j] += fctb[j];
              rhs[cell_id2][i][j] += fctb[j];
            }
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from extended neighborhood */

  if (halo_type == CS_HALO_EXTENDED) {

#   pragma omp parallel for
    for (cs_lnum_t cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
      for (cs_lnum_t cidx = cell_cells_idx[cell_id1];
           cidx < cell_cells_idx[cell_id1+1];
           cidx++) {

        cs_lnum_t cell_id2 = cell_cells_lst[cidx];

        cs_real_3_t dc;
        for (int i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        cs_real_t ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (int i = 0; i < 6; i++) {

          cs_real_t pfac = (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

          for (int j = 0; j < 3; j++) {
            rhs[cell_id1][i][j] += dc[j] * pfac;
          }
        }
      }
    }

  } /* End for extended neighborhood */

  /* Contribution from boundary faces */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id1 = b_face_cells[face_id];

        cs_real_3_t n_d_dist;
        /* Normal is vector 0 if the b_face_normal norm is too small */
        cs_math_3_normalise(b_face_normal[face_id], n_d_dist);

        cs_real_t d_b_dist = 1. / b_dist[face_id];

        /* Normal divided by b_dist */
        for (int i = 0; i < 3; i++)
          n_d_dist[i] *= d_b_dist;

        for (int i = 0; i < 6; i++) {
          cs_real_t pfac = coefat[face_id][i]*inc - pvar[cell_id1][i];
          for (int j = 0; j < 6; j++)
            pfac += coefbt[face_id][j][i] * pvar[cell_id1][j];

          for (int j = 0; j < 3; j++)
            rhs[cell_id1][i][j] += pfac * n_d_dist[j];
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Compute gradient */
  /*------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 6; i++) {

        gradt[cell_id][i][j] = 0.0;

        for (int k = 0; k < 3; k++)
          gradt[cell_id][i][j] += rhs[cell_id][i][k] * cocg[cell_id][k][j];

      }
    }
  }

  /* Compute gradient on boundary cells */
  /*------------------------------------*/

  #pragma omp parallel
  {
    cs_lnum_t t_s_id, t_e_id;
    _thread_range(m->n_b_cells, &t_s_id, &t_e_id);

    /* Build indices bijection between [1-18] and [1-6]*[1-3] */

    cs_lnum_t _63_18_idx[18][2];
    int nn = 0;
    for (int ll = 0; ll < 6; ll++) {
      for (int mm = 0; mm < 3; mm++) {
        _63_18_idx[nn][0] = ll;
        _63_18_idx[nn][1] = mm;
        nn++;
      }
    }

    /* Loop on boundary cells */

    for (cs_lnum_t b_cell_id = t_s_id; b_cell_id < t_e_id; b_cell_id++) {

      cs_lnum_t cell_id = m->b_cells[b_cell_id];

      cs_real_t cocgb_t[171], rhsb_t[18], x[18];

      _compute_cocgb_rhsb_lsq_t
        (cell_id,
         inc,
         madj,
         fvq,
         _63_18_idx,
         (const cs_real_6_t *)pvar,
         (const cs_real_6_t *)coefat,
         (const cs_real_66_t *)coefbt,
         (const cs_real_63_t *)rhs,
         cocgb_t,
         rhsb_t);

      _fw_and_bw_ldtl_pp(cocgb_t,
                         18,
                         x,
                         rhsb_t);

      for (int kk = 0; kk < 18; kk++) {
        int ii = _63_18_idx[kk][0];
        int jj = _63_18_idx[kk][1];
        gradt[cell_id][ii][jj] = x[kk];
      }

    }

  }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)gradt, 18);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)gradt);
  }

  BFT_FREE(rhs);
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
 *   coefat         <-- B.C. coefficients for boundary face normals
 *   coefbt         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   grad          --> gradient of pvar (dts_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

static void
_initialize_tensor_gradient(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            int                           inc,
                            const cs_real_6_t   *restrict coefat,
                            const cs_real_66_t  *restrict coefbt,
                            const cs_real_6_t   *restrict pvar,
                            cs_real_63_t        *restrict grad)
{
  int g_id, t_id;
  cs_lnum_t face_id;
  cs_real_t pond;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  /* Computation without reconstruction */
  /*------------------------------------*/

  /* Initialization */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++)
        grad[cell_id][i][j] = 0.0;
    }
  }

  /* Interior faces contribution */

  for (g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for private(face_id, pond)
    for (t_id = 0; t_id < n_i_threads; t_id++) {

      for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        pond = weight[face_id];

        /*
           Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                    + (1-\alpha_\ij) \varia_\cellj\f$
                   but for the cell \f$ \celli \f$ we remove
                   \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                   and for the cell \f$ \cellj \f$ we remove
                   \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
        */
        for (int i = 0; i < 6; i++) {
          cs_real_t pfaci = (1.0-pond) * (pvar[cell_id2][i] - pvar[cell_id1][i]);
          cs_real_t pfacj = - pond * (pvar[cell_id2][i] - pvar[cell_id1][i]);
          for (int j = 0; j < 3; j++) {
            grad[cell_id1][i][j] += pfaci * i_f_face_normal[face_id][j];
            grad[cell_id2][i][j] -= pfacj * i_f_face_normal[face_id][j];
          }
        }

      } /* End of loop on faces */

    } /* End of loop on threads */

  } /* End of loop on thread groups */

  /* Boundary face treatment */

  for (g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for private(face_id)
    for (t_id = 0; t_id < n_b_threads; t_id++) {

      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t cell_id = b_face_cells[face_id];

        /*
           Remark: for the cell \f$ \celli \f$ we remove
                   \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
         */
        for (int i = 0; i < 6; i++) {
          cs_real_t pfac = inc*coefat[face_id][i];

          for (int k = 0; k < 6; k++) {
            if (i == k)
              pfac += (coefbt[face_id][i][k] - 1.0) * pvar[cell_id][k];
            else
              pfac += coefbt[face_id][i][k] * pvar[cell_id][k] ;

          }

          for (int j = 0; j < 3; j++)
            grad[cell_id][i][j] += pfac * b_f_face_normal[face_id][j];
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++)
        grad[cell_id][i][j] *= dvol;
    }
  }


  /* Periodicity and parallelism treatment */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 18);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                           halo_type,
                                           (cs_real_t *)grad);
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
 *   fvq            <-- pointer to associated finite volume quantities
 *   var_name       <-- variable's name
 *   gradient_info  <-- pointer to performance logging structure, or NULL
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_r_sweeps     --> >1: with reconstruction
 *   verbosity      --> verbosity level
 *   epsrgp         --> precision for iterative gradient calculation
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
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
                           const cs_real_6_t   *restrict coefat,
                           const cs_real_66_t  *restrict coefbt,
                           const cs_real_6_t   *restrict pvar,
                           cs_real_63_t        *restrict grad)
{
  int isweep = 0;
  int g_id, t_id;
  cs_lnum_t face_id;
  cs_real_t l2_norm, l2_residual, vecfac, pond;

  cs_real_63_t *rhs;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;
  cs_real_33_t *restrict cocg = fvq->cocg_it;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_63_t);

  /* Gradient reconstruction to handle non-orthogonal meshes */
  /*---------------------------------------------------------*/

  /* L2 norm */

  l2_norm = _l2_norm_1(18*n_cells, (cs_real_t *)grad);
  l2_residual = l2_norm ;

  if (l2_norm > cs_math_epzero) {

    /* Iterative process */
    /*-------------------*/

    for (isweep = 1; isweep < n_r_sweeps && l2_residual > epsrgp*l2_norm; isweep++) {

      /* Computation of the Right Hand Side*/

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        for (int i = 0; i < 6; i++) {
          for (int j = 0; j < 3; j++)
            rhs[cell_id][i][j] = - cell_f_vol[cell_id] * grad[cell_id][i][j];
        }
      }

      /* Interior face treatment */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, pond)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id1 = i_face_cells[face_id][0];
            cs_lnum_t cell_id2 = i_face_cells[face_id][1];
            pond = weight[face_id];

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            for (int i = 0; i < 6; i++) {

              /* Reconstruction part */
              cs_real_t
              pfaci = 0.5 * ( ( grad[cell_id1][i][0] + grad[cell_id2][i][0])
                              * dofij[face_id][0]
                            + ( grad[cell_id1][i][1] + grad[cell_id2][i][1])
                              * dofij[face_id][1]
                            + ( grad[cell_id1][i][2] + grad[cell_id2][i][2])
                              * dofij[face_id][2]
                            );
              cs_real_t pfacj = pfaci;

              pfaci += (1.0-pond) * (pvar[cell_id2][i] - pvar[cell_id1][i]);
              pfacj -=       pond * (pvar[cell_id2][i] - pvar[cell_id1][i]);
              for (int j = 0; j < 3; j++) {
                rhs[cell_id1][i][j] += pfaci * i_f_face_normal[face_id][j];
                rhs[cell_id2][i][j] -= pfacj * i_f_face_normal[face_id][j];
              }
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Boundary face treatment */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id, vecfac)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id = b_face_cells[face_id];

            for (int i = 0; i < 6; i++) {

              /*
                 Remark: for the cell \f$ \celli \f$ we remove
                         \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
               */

              cs_real_t pfac = inc*coefat[face_id][i];

              for (int k = 0; k < 6; k++) {
                /* Reconstruction part */
                vecfac = grad[cell_id][k][0] * diipb[face_id][0]
                       + grad[cell_id][k][1] * diipb[face_id][1]
                       + grad[cell_id][k][2] * diipb[face_id][2];
                pfac += coefbt[face_id][i][k] * vecfac;

                if (i == k)
                  pfac += (coefbt[face_id][i][k] - 1.0) * pvar[cell_id][k];
                else
                  pfac += coefbt[face_id][i][k] * pvar[cell_id][k];
              }

              for (int j = 0; j < 3; j++)
                rhs[cell_id][i][j] += pfac * b_f_face_normal[face_id][j];

            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Increment of the gradient */

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_real_t dvol;
        /* Is the cell disabled (for solid or porous)? Not the case if coupled */
        if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
          dvol = 1. / cell_f_vol[cell_id];
        else
          dvol = 0.;

        for (int i = 0; i < 6; i++) {
          for (int j = 0; j < 3; j++)
            rhs[cell_id][i][j] *= dvol;
        }

        for (int i = 0; i < 6; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++)
              grad[cell_id][i][j] += rhs[cell_id][i][k] * cocg[cell_id][k][j];
          }
        }
      }

      /* Periodicity and parallelism treatment */

      if (m->halo != NULL) {
        cs_halo_sync_var_strided(m->halo, halo_type, (cs_real_t *)grad, 18);
        if (cs_glob_mesh->n_init_perio > 0)
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

  if (gradient_info != NULL)
    _gradient_info_update_iter(gradient_info, isweep);

  BFT_FREE(rhs);
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
 *   cpl            <-> structure associated with internal coupling, or NULL
 *   idimtr         <-- 0 if ivar does not match a vector or tensor
 *                        or there is no periodicity of rotation
 *                      1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   f_ext          <-- exterior force generating pressure
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   r_grad         --> gradient used for reconstruction
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

static void
_reconstruct_scalar_gradient(const cs_mesh_t                 *m,
                             const cs_mesh_quantities_t      *fvq,
                             const cs_internal_coupling_t    *cpl,
                             int                              idimtr,
                             int                              hyd_p_flag,
                             const cs_real_3_t                f_ext[],
                             const cs_real_t                  coefbp[],
                             cs_real_3_t            *restrict r_grad,
                             cs_real_3_t            *restrict grad)
{
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;

  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const cs_real_33_t *restrict corr_grad_lin
    = (const cs_real_33_t *restrict)fvq->corr_grad_lin;

  int        g_id, t_id;

  cs_real_3_t  fexd;

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Initialize gradient */
  /*---------------------*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    grad[cell_id][0] = grad[cell_id][0] * cell_f_vol[cell_id];
    grad[cell_id][1] = grad[cell_id][1] * cell_f_vol[cell_id];
    grad[cell_id][2] = grad[cell_id][2] * cell_f_vol[cell_id];
  }

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  if (hyd_p_flag == 1) {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for private(fexd)
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t cell_id1 = i_face_cells[face_id][0];
          cs_lnum_t cell_id2 = i_face_cells[face_id][1];

          fexd[0] = 0.5 * (f_ext[cell_id1][0] + f_ext[cell_id2][0]);//FIXME not correct ...
          fexd[1] = 0.5 * (f_ext[cell_id1][1] + f_ext[cell_id2][1]);
          fexd[2] = 0.5 * (f_ext[cell_id1][2] + f_ext[cell_id2][2]);

          cs_real_t rfac =
                 weight[face_id]
                 * ( (cell_cen[cell_id1][0]-i_face_cog[face_id][0])*fexd[0]
                   + (cell_cen[cell_id1][1]-i_face_cog[face_id][1])*fexd[1]
                   + (cell_cen[cell_id1][2]-i_face_cog[face_id][2])*fexd[2])
              +  (1.0 - weight[face_id])
                 * ( (cell_cen[cell_id2][0]-i_face_cog[face_id][0])*fexd[0]
                   + (cell_cen[cell_id2][1]-i_face_cog[face_id][1])*fexd[1]
                   + (cell_cen[cell_id2][2]-i_face_cog[face_id][2])*fexd[2])
              + ( dofij[face_id][0] * (r_grad[cell_id1][0]+r_grad[cell_id2][0])
                + dofij[face_id][1] * (r_grad[cell_id1][1]+r_grad[cell_id2][1])
                + dofij[face_id][2] * (r_grad[cell_id1][2]+r_grad[cell_id2][2]))
                *0.5;

          for (int j = 0; j < 3; j++) {
            grad[cell_id1][j] += rfac * i_f_face_normal[face_id][j];
            grad[cell_id2][j] -= rfac * i_f_face_normal[face_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t cell_id = b_face_cells[face_id];

          /*
             Remark: for the cell \f$ \celli \f$ we remove
                     \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
           */

          /* Reconstruction part */
          cs_real_t
          rfac = coefbp[face_id]
                 * ( diipb[face_id][0] * (r_grad[cell_id][0] - f_ext[cell_id][0])
                   + diipb[face_id][1] * (r_grad[cell_id][1] - f_ext[cell_id][1])
                   + diipb[face_id][2] * (r_grad[cell_id][2] - f_ext[cell_id][2]));

          for (int j = 0; j < 3; j++) {
            grad[cell_id][j] += rfac * b_f_face_normal[face_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  } /* End of test on hydrostatic pressure */


  /* Standard case, without hydrostatic pressure */
  /*---------------------------------------------*/

  else {

    /* Contribution from interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t cell_id1 = i_face_cells[face_id][0];
          cs_lnum_t cell_id2 = i_face_cells[face_id][1];

          /* Reconstruction part */
          cs_real_t rfac = 0.5 *
                    (dofij[face_id][0]*(r_grad[cell_id1][0]+r_grad[cell_id2][0])
                    +dofij[face_id][1]*(r_grad[cell_id1][1]+r_grad[cell_id2][1])
                    +dofij[face_id][2]*(r_grad[cell_id1][2]+r_grad[cell_id2][2]));

          for (int j = 0; j < 3; j++) {
            grad[cell_id1][j] += rfac * i_f_face_normal[face_id][j];
            grad[cell_id2][j] -= rfac * i_f_face_normal[face_id][j];
          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

    /* Contribution from coupled faces */
    if (cpl != NULL)
      cs_internal_coupling_reconstruct_scalar_gradient
        (cpl, r_grad, grad);

    /* Contribution from boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          if (cpl == NULL || !coupled_faces[face_id]) {

            cs_lnum_t cell_id = b_face_cells[face_id];

            /* Reconstruction part */
            cs_real_t
            rfac = coefbp[face_id]
                   * ( diipb[face_id][0] * r_grad[cell_id][0]
                     + diipb[face_id][1] * r_grad[cell_id][1]
                     + diipb[face_id][2] * r_grad[cell_id][2] );

            for (int j = 0; j < 3; j++) {
              grad[cell_id][j] += rfac * b_f_face_normal[face_id][j];
            }

          }

        } /* loop on faces */

      } /* loop on threads */

    } /* loop on thread groups */

  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_f_vol[cell_id];
    else
      dvol = 0.;

    grad[cell_id][0] *= dvol;
    grad[cell_id][1] *= dvol;
    grad[cell_id][2] *= dvol;

    if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_WARPED_CORRECTION) {
      cs_real_3_t gradpa;
      for (int i = 0; i < 3; i++) {
        gradpa[i] = grad[cell_id][i];
        grad[cell_id][i] = 0.;
      }

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          grad[cell_id][i] += corr_grad_lin[cell_id][i][j] * gradpa[j];
    }
  }

  /* Synchronize halos */

  _sync_scalar_gradient_halo(m, CS_HALO_EXTENDED, idimtr, grad);
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
 *   cpl             <-- structure associated with internal coupling, or NULL
 *   var_name        <-- variable name
 *   gradient_info   <-- pointer to performance logging structure, or NULL
 *   nswrgp          <-- number of sweeps for gradient reconstruction
 *   idimtr          <-- 0 if ivar does not match a vector or tensor
 *                         or there is no periodicity of rotation
 *                       1 for velocity, 2 for Reynolds stress
 *   hyd_p_flag      <-- flag for hydrostatic pressure
 *   verbosity       <-- verbosity level
 *   inc             <-- if 0, solve on increment; 1 otherwise
 *   epsrgp          <-- relative precision for gradient reconstruction
 *   f_ext           <-- exterior force generating pressure
 *   coefap          <-- B.C. coefficients for boundary face normals
 *   coefbp          <-- B.C. coefficients for boundary face normals
 *   pvar            <-- variable
 *   c_weight        <-- weighted gradient coefficient variable
 *   grad            <-> gradient of pvar (halo prepared for periodicity
 *                       of rotation)
 *----------------------------------------------------------------------------*/

static void
_iterative_scalar_gradient(const cs_mesh_t                *m,
                           cs_mesh_quantities_t           *fvq,
                           const cs_internal_coupling_t   *cpl,
                           const char                     *var_name,
                           cs_gradient_info_t             *gradient_info,
                           int                             nswrgp,
                           int                             idimtr,
                           int                             hyd_p_flag,
                           int                             verbosity,
                           cs_real_t                       inc,
                           cs_real_t                       epsrgp,
                           const cs_real_3_t               f_ext[],
                           const cs_real_t                 coefap[],
                           const cs_real_t                 coefbp[],
                           const cs_real_t                 pvar[],
                           const cs_real_t                 c_weight[],
                           cs_real_3_t           *restrict grad)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  int has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  cs_real_33_t *restrict cocg = (cpl == NULL) ?
    fvq->cocg_it : cpl->cocg_it;

  cs_lnum_t  face_id;
  int        g_id, t_id;
  cs_real_t  rnorm;
  cs_real_3_t  fexd;
  cs_real_3_t *rhs;

  int n_sweeps = 0;
  cs_real_t l2_residual = 0.;

  if (nswrgp < 1) {
    if (gradient_info != NULL)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  bool  *coupled_faces = (cpl == NULL) ?
    NULL : (bool *)cpl->coupled_faces;

  /* Reconstruct gradients for non-orthogonal meshes */
  /*-------------------------------------------------*/

  /* Semi-implicit resolution on the whole mesh  */

  /* Compute normalization residual */

  rnorm = _l2_norm_1(3*n_cells, (cs_real_t *)grad);

  if (rnorm <= cs_math_epzero) {
    if (gradient_info != NULL)
      _gradient_info_update_iter(gradient_info, 0);
    return;
  }

  BFT_MALLOC(rhs, n_cells_ext, cs_real_3_t);

  /* Vector OijFij is computed in CLDijP */

  /* Start iterations */
  /*------------------*/

  for (n_sweeps = 1; n_sweeps < nswrgp; n_sweeps++) {

    /* Compute right hand side */

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      rhs[cell_id][0] = -grad[cell_id][0] * cell_f_vol[cell_id];
      rhs[cell_id][1] = -grad[cell_id][1] * cell_f_vol[cell_id];
      rhs[cell_id][2] = -grad[cell_id][2] * cell_f_vol[cell_id];
    }

    /* Case with hydrostatic pressure */
    /*--------------------------------*/

    if (hyd_p_flag == 1) {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id, fexd)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id1 = i_face_cells[face_id][0];
            cs_lnum_t cell_id2 = i_face_cells[face_id][1];

            cs_real_t ktpond = (c_weight == NULL) ?
              weight[face_id] :                     // no cell weighting
              weight[face_id]  * c_weight[cell_id1] // cell weighting active
                / (      weight[face_id]  * c_weight[cell_id1]
                  + (1.0-weight[face_id]) * c_weight[cell_id2]);

            // TODO add porous contribution
            fexd[0] = 0.5 * (f_ext[cell_id1][0] + f_ext[cell_id2][0]);
            fexd[1] = 0.5 * (f_ext[cell_id1][1] + f_ext[cell_id2][1]);
            fexd[2] = 0.5 * (f_ext[cell_id1][2] + f_ext[cell_id2][2]);

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            /* Reconstruction part */
            cs_real_t pfaci =
                     (i_face_cog[face_id][0]-cell_cen[cell_id1][0])
                    *(ktpond*f_ext[cell_id1][0]-weight[face_id]*fexd[0])
                   + (i_face_cog[face_id][1]-cell_cen[cell_id1][1])
                    *(ktpond*f_ext[cell_id1][1]-weight[face_id]*fexd[1])
                   + (i_face_cog[face_id][2]-cell_cen[cell_id1][2])
                    *(ktpond*f_ext[cell_id1][2]-weight[face_id]*fexd[2])
               +     (i_face_cog[face_id][0]-cell_cen[cell_id2][0])
                    *((1.0 - ktpond)*f_ext[cell_id2][0]-(1.-weight[face_id])*fexd[0])
                   + (i_face_cog[face_id][1]-cell_cen[cell_id2][1])
                    *((1.0 - ktpond)*f_ext[cell_id2][1]-(1.-weight[face_id])*fexd[1])
                   + (i_face_cog[face_id][2]-cell_cen[cell_id2][2])
                    *((1.0 - ktpond)*f_ext[cell_id2][2]-(1.-weight[face_id])*fexd[2])
               + ( dofij[face_id][0] * (grad[cell_id1][0]+grad[cell_id2][0])
                 + dofij[face_id][1] * (grad[cell_id1][1]+grad[cell_id2][1])
                 + dofij[face_id][2] * (grad[cell_id1][2]+grad[cell_id2][2]))*0.5;

            cs_real_t pfacj = pfaci;

            pfaci += (1.0-ktpond) * (pvar[cell_id2] - pvar[cell_id1]);
            pfacj -= ktpond * (pvar[cell_id2] - pvar[cell_id1]);

            for (int j = 0; j < 3; j++) {
              rhs[cell_id1][j] += pfaci * i_f_face_normal[face_id][j];
              rhs[cell_id2][j] -= pfacj * i_f_face_normal[face_id][j];
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from boundary faces */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id = b_face_cells[face_id];

            /*
               Remark: for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            /* Reconstruction part */
            cs_real_t
            pfac = coefap[face_id] * inc
                 + coefbp[face_id]
                   * ( diipb[face_id][0] * (grad[cell_id][0] - f_ext[cell_id][0])
                     + diipb[face_id][1] * (grad[cell_id][1] - f_ext[cell_id][1])
                     + diipb[face_id][2] * (grad[cell_id][2] - f_ext[cell_id][2])
                     + (b_face_cog[face_id][0]-cell_cen[cell_id][0])*f_ext[cell_id][0]
                     + (b_face_cog[face_id][1]-cell_cen[cell_id][1]) * f_ext[cell_id][1]
                     + (b_face_cog[face_id][2]-cell_cen[cell_id][2]) * f_ext[cell_id][2]);

            pfac += (coefbp[face_id] -1.0) * pvar[cell_id];

            rhs[cell_id][0] += pfac * b_f_face_normal[face_id][0];
            rhs[cell_id][1] += pfac * b_f_face_normal[face_id][1];
            rhs[cell_id][2] += pfac * b_f_face_normal[face_id][2];

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    } /* End of test on hydrostatic pressure */

    /* Standard case, without hydrostatic pressure */
    /*---------------------------------------------*/

    else {

      /* Contribution from interior faces */

      for (g_id = 0; g_id < n_i_groups; g_id++) {

#       pragma omp parallel for private(face_id)
        for (t_id = 0; t_id < n_i_threads; t_id++) {

          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            cs_lnum_t cell_id1 = i_face_cells[face_id][0];
            cs_lnum_t cell_id2 = i_face_cells[face_id][1];

            /*
               Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                        + (1-\alpha_\ij) \varia_\cellj\f$
                       but for the cell \f$ \celli \f$ we remove
                       \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
                       and for the cell \f$ \cellj \f$ we remove
                       \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
            */

            /* Reconstruction part */
            cs_real_t pfaci = 0.5 *
                     (dofij[face_id][0]*(grad[cell_id1][0]+grad[cell_id2][0])
                     +dofij[face_id][1]*(grad[cell_id1][1]+grad[cell_id2][1])
                     +dofij[face_id][2]*(grad[cell_id1][2]+grad[cell_id2][2]));
            cs_real_t pfacj = pfaci;

            cs_real_t ktpond = (c_weight == NULL) ?
              weight[face_id] :                     // no cell weighting
              weight[face_id]  * c_weight[cell_id1] // cell weighting active
                / (      weight[face_id]  * c_weight[cell_id1]
                  + (1.0-weight[face_id]) * c_weight[cell_id2]);

            pfaci += (1.0-ktpond) * (pvar[cell_id2] - pvar[cell_id1]);
            pfacj -=      ktpond  * (pvar[cell_id2] - pvar[cell_id1]);

            for (int j = 0; j < 3; j++) {
              rhs[cell_id1][j] += pfaci * i_f_face_normal[face_id][j];
              rhs[cell_id2][j] -= pfacj * i_f_face_normal[face_id][j];
            }

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

      /* Contribution from coupled faces */
      if (cpl != NULL)
        cs_internal_coupling_iterative_scalar_gradient
          (cpl,
           c_weight,
           grad,
           pvar,
           rhs);

      /* Contribution from boundary faces */

      for (g_id = 0; g_id < n_b_groups; g_id++) {

#       pragma omp parallel for private(face_id)
        for (t_id = 0; t_id < n_b_threads; t_id++) {

          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            if (cpl == NULL || !coupled_faces[face_id]) {

              cs_lnum_t cell_id = b_face_cells[face_id];

              /*
                 Remark: for the cell \f$ \celli \f$ we remove
                         \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
               */

              /* Reconstruction part */
              cs_real_t
              pfac = coefap[face_id] * inc
                   + coefbp[face_id]
                     * ( diipb[face_id][0] * grad[cell_id][0]
                       + diipb[face_id][1] * grad[cell_id][1]
                       + diipb[face_id][2] * grad[cell_id][2]
                       );

              pfac += (coefbp[face_id] -1.0) * pvar[cell_id];

              rhs[cell_id][0] += pfac * b_f_face_normal[face_id][0];
              rhs[cell_id][1] += pfac * b_f_face_normal[face_id][1];
              rhs[cell_id][2] += pfac * b_f_face_normal[face_id][2];

            } /* face without internal coupling */

          } /* loop on faces */

        } /* loop on threads */

      } /* loop on thread groups */

    }

    /* Increment gradient */
    /*--------------------*/

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cs_real_t dvol;
      /* Is the cell disabled (for solid or porous)? Not the case if coupled */
      if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
        dvol = 1. / cell_f_vol[cell_id];
      else
        dvol = 0.;

      rhs[cell_id][0] *= dvol;
      rhs[cell_id][1] *= dvol;
      rhs[cell_id][2] *= dvol;

      grad[cell_id][0] +=   rhs[cell_id][0] * cocg[cell_id][0][0]
                          + rhs[cell_id][1] * cocg[cell_id][1][0]
                          + rhs[cell_id][2] * cocg[cell_id][2][0];
      grad[cell_id][1] +=   rhs[cell_id][0] * cocg[cell_id][0][1]
                          + rhs[cell_id][1] * cocg[cell_id][1][1]
                          + rhs[cell_id][2] * cocg[cell_id][2][1];
      grad[cell_id][2] +=   rhs[cell_id][0] * cocg[cell_id][0][2]
                          + rhs[cell_id][1] * cocg[cell_id][1][2]
                          + rhs[cell_id][2] * cocg[cell_id][2][2];
    }

    /* Synchronize halos */

    _sync_scalar_gradient_halo(m, CS_HALO_STANDARD, idimtr, grad);

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

  if (gradient_info != NULL)
    _gradient_info_update_iter(gradient_info, n_sweeps);

  BFT_FREE(rhs);
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
 * \param[in]       var_name        variable name
 * \param[in]       gradient_info   performance logging structure, or NULL
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       tr_dim          2 for tensor with periodicity of rotation,
 *                                  0 otherwise
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
 * \param[in]       w_stride        stride for weighting coefficient
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in]       var             gradient's base variable
 * \param[in]       c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in]       cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_scalar(const char                    *var_name,
                 cs_gradient_info_t            *gradient_info,
                 cs_gradient_type_t             gradient_type,
                 cs_halo_type_t                 halo_type,
                 int                            inc,
                 bool                           recompute_cocg,
                 int                            n_r_sweeps,
                 int                            tr_dim,
                 int                            hyd_p_flag,
                 int                            w_stride,
                 int                            verbosity,
                 int                            clip_mode,
                 double                         epsilon,
                 double                         extrap,
                 double                         clip_coeff,
                 cs_real_t                      f_ext[][3],
                 const cs_real_t                bc_coeff_a[],
                 const cs_real_t                bc_coeff_b[],
                 const cs_real_t                var[restrict],
                 const cs_real_t                c_weight[restrict],
                 const cs_internal_coupling_t  *cpl,
                 cs_real_t                      grad[restrict][3])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_real_4_t  *restrict rhsv;

  cs_lnum_t n_b_faces = mesh->n_b_faces;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  static int last_fvm_count = 0;

  if (n_r_sweeps > 0) {
    int prev_fvq_count = last_fvm_count;
    last_fvm_count = cs_mesh_quantities_compute_count();
    if (last_fvm_count != prev_fvq_count)
      recompute_cocg = true;
  }

  /* Use Neumann BC's as default if not provided */

  cs_real_t *_bc_coeff_a = NULL;
  cs_real_t *_bc_coeff_b = NULL;

  if (bc_coeff_a == NULL) {
    BFT_MALLOC(_bc_coeff_a, n_b_faces, cs_real_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      _bc_coeff_a[i] = 0;
    bc_coeff_a = _bc_coeff_a;
  }
  if (bc_coeff_b == NULL) {
    BFT_MALLOC(_bc_coeff_b, n_b_faces, cs_real_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      _bc_coeff_b[i] = 1;
    bc_coeff_b = _bc_coeff_b;
  }

  /* Allocate work arrays */

  BFT_MALLOC(rhsv, n_cells_ext, cs_real_4_t);

  /* Compute gradient */

  if (gradient_type == CS_GRADIENT_ITER) {

    _initialize_scalar_gradient(mesh,
                                fvq,
                                cpl,
                                tr_dim,
                                hyd_p_flag,
                                inc,
                                (const cs_real_3_t *)f_ext,
                                bc_coeff_a,
                                bc_coeff_b,
                                var,
                                c_weight,
                                grad);

    _iterative_scalar_gradient(mesh,
                               fvq,
                               cpl,
                               var_name,
                               gradient_info,
                               n_r_sweeps,
                               tr_dim,
                               hyd_p_flag,
                               verbosity,
                               inc,
                               epsilon,
                               (const cs_real_3_t *)f_ext,
                               bc_coeff_a,
                               bc_coeff_b,
                               var,
                               c_weight,
                               grad);

  } else if (gradient_type == CS_GRADIENT_ITER_OLD) {

    _initialize_scalar_gradient_old(mesh,
                                    fvq,
                                    tr_dim,
                                    hyd_p_flag,
                                    inc,
                                    (const cs_real_3_t *)f_ext,
                                    bc_coeff_a,
                                    bc_coeff_b,
                                    var,
                                    c_weight,
                                    grad,
                                    rhsv);

    _iterative_scalar_gradient_old(mesh,
                                   fvq,
                                   var_name,
                                   gradient_info,
                                   recompute_cocg,
                                   n_r_sweeps,
                                   tr_dim,
                                   hyd_p_flag,
                                   verbosity,
                                   inc,
                                   epsilon,
                                   extrap,
                                   (const cs_real_3_t *)f_ext,
                                   bc_coeff_a,
                                   bc_coeff_b,
                                   grad,
                                   rhsv);

  } else if (gradient_type == CS_GRADIENT_LSQ) {

    _lsq_scalar_gradient(mesh,
                         fvq,
                         cpl,
                         halo_type,
                         recompute_cocg,
                         n_r_sweeps,
                         tr_dim,
                         hyd_p_flag,
                         w_stride,
                         inc,
                         extrap,
                         (const cs_real_3_t *)f_ext,
                         bc_coeff_a,
                         bc_coeff_b,
                         var,
                         c_weight,
                         grad,
                         rhsv);

  } else if (gradient_type == CS_GRADIENT_LSQ_ITER) {

    const cs_int_t  _imlini = 1;
    const cs_real_t _climin = 1.5;

    cs_real_3_t  *restrict r_grad;
    BFT_MALLOC(r_grad, n_cells_ext, cs_real_3_t);

    _lsq_scalar_gradient(mesh,
                         fvq,
                         cpl,
                         halo_type,
                         recompute_cocg,
                         n_r_sweeps,
                         tr_dim,
                         hyd_p_flag,
                         w_stride,
                         inc,
                         extrap,
                         (const cs_real_3_t *)f_ext,
                         bc_coeff_a,
                         bc_coeff_b,
                         var,
                         c_weight,
                         r_grad,
                         rhsv);

    _scalar_gradient_clipping(halo_type, _imlini, verbosity, tr_dim, _climin,
                              var, r_grad);

    _initialize_scalar_gradient(mesh,
                                fvq,
                                cpl,
                                tr_dim,
                                hyd_p_flag,
                                inc,
                                (const cs_real_3_t *)f_ext,
                                bc_coeff_a,
                                bc_coeff_b,
                                var,
                                c_weight,
                                grad);

    _reconstruct_scalar_gradient(mesh,
                                 fvq,
                                 cpl,
                                 tr_dim,
                                 hyd_p_flag,
                                 (const cs_real_3_t *)f_ext,
                                 bc_coeff_b,
                                 r_grad,
                                 grad);

    BFT_FREE(r_grad);

  }

  _scalar_gradient_clipping(halo_type, clip_mode, verbosity, tr_dim, clip_coeff,
                            var, grad);

  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION)
    cs_bad_cells_regularisation_vector(grad, 0);

  BFT_FREE(_bc_coeff_a);
  BFT_FREE(_bc_coeff_b);

  BFT_FREE(rhsv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_info   performance logging structure, or NULL
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in]       var             gradient's base variable
 * \param[in]       c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in]       cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
                                    (\f$ \der{u_i}{x_j} \f$ is grad[][i][j])
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_vector(const char                    *var_name,
                 cs_gradient_info_t            *gradient_info,
                 cs_gradient_type_t             gradient_type,
                 cs_halo_type_t                 halo_type,
                 int                            inc,
                 int                            n_r_sweeps,
                 int                            verbosity,
                 int                            clip_mode,
                 double                         epsilon,
                 double                         clip_coeff,
                 const cs_real_3_t              bc_coeff_a[],
                 const cs_real_33_t             bc_coeff_b[],
                 const cs_real_3_t    *restrict var,
                 const cs_real_t      *restrict c_weight,
                 const cs_internal_coupling_t  *cpl,
                 cs_real_33_t         *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  /* Use Neumann BC's as default if not provided */

  cs_real_3_t *_bc_coeff_a = NULL;
  cs_real_33_t *_bc_coeff_b = NULL;

  if (bc_coeff_a == NULL) {
    BFT_MALLOC(_bc_coeff_a, n_b_faces, cs_real_3_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _bc_coeff_a[i][j] = 0;
    }
    bc_coeff_a = (const cs_real_3_t *)_bc_coeff_a;
  }
  if (bc_coeff_b == NULL) {
    BFT_MALLOC(_bc_coeff_b, n_b_faces, cs_real_33_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          _bc_coeff_b[i][j][j] = 1;
      }
    }
    bc_coeff_b = (const cs_real_33_t *)_bc_coeff_b;
  }

  /* Compute gradient */

  if (  gradient_type == CS_GRADIENT_ITER
     || gradient_type == CS_GRADIENT_ITER_OLD) {

    _initialize_vector_gradient(mesh,
                                fvq,
                                cpl,
                                halo_type,
                                inc,
                                bc_coeff_a,
                                bc_coeff_b,
                                var,
                                c_weight,
                                grad);

    /* If reconstructions are required */

    if (n_r_sweeps > 1)
      _iterative_vector_gradient(mesh,
                                 fvq,
                                 cpl,
                                 var_name,
                                 gradient_info,
                                 halo_type,
                                 inc,
                                 n_r_sweeps,
                                 verbosity,
                                 epsilon,
                                 bc_coeff_a,
                                 bc_coeff_b,
                                 (const cs_real_3_t *)var,
                                 c_weight,
                                 grad);

  } else if (gradient_type == CS_GRADIENT_LSQ) {

    _lsq_vector_gradient(mesh,
                         cs_glob_mesh_adjacencies,
                         fvq,
                         cpl,
                         halo_type,
                         inc,
                         bc_coeff_a,
                         bc_coeff_b,
                         (const cs_real_3_t *)var,
                         c_weight,
                         grad);

  } else if (gradient_type == CS_GRADIENT_LSQ_ITER) {

    /* Clipping algorithm and clipping factor */

    const cs_int_t  _imlini = 1;
    const cs_real_t _climin = 1.5;

    cs_real_33_t  *restrict r_gradv;
    BFT_MALLOC(r_gradv, n_cells_ext, cs_real_33_t);

    _lsq_vector_gradient(mesh,
                         cs_glob_mesh_adjacencies,
                         fvq,
                         cpl,
                         halo_type,
                         inc,
                         bc_coeff_a,
                         bc_coeff_b,
                         (const cs_real_3_t *)var,
                         c_weight,
                         r_gradv);

    _vector_gradient_clipping(mesh,
                              fvq,
                              halo_type,
                              _imlini,
                              verbosity,
                              _climin,
                              (const cs_real_3_t *)var,
                              r_gradv);

    _initialize_vector_gradient(mesh,
                                fvq,
                                cpl,
                                halo_type,
                                inc,
                                bc_coeff_a,
                                bc_coeff_b,
                                var,
                                c_weight,
                                grad);

    _reconstruct_vector_gradient(mesh,
                                 fvq,
                                 cpl,
                                 halo_type,
                                 bc_coeff_b,
                                 r_gradv,
                                 grad);

    BFT_FREE(r_gradv);
  }

  _vector_gradient_clipping(mesh,
                            fvq,
                            halo_type,
                            clip_mode,
                            verbosity,
                            clip_coeff,
                            (const cs_real_3_t *)var,
                            grad);

  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION)
    cs_bad_cells_regularisation_tensor((cs_real_9_t *)grad, 0);

  BFT_FREE(_bc_coeff_a);
  BFT_FREE(_bc_coeff_b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_info   performance logging structure, or NULL
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in]       var             gradient's base variable
 * \param[out]      grad            gradient
                                      (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_tensor(const char                *var_name,
                 cs_gradient_info_t        *gradient_info,
                 cs_gradient_type_t         gradient_type,
                 cs_halo_type_t             halo_type,
                 int                        inc,
                 int                        n_r_sweeps,
                 int                        verbosity,
                 int                        clip_mode,
                 double                     epsilon,
                 double                     clip_coeff,
                 const cs_real_6_t          bc_coeff_a[],
                 const cs_real_66_t         bc_coeff_b[],
                 const cs_real_6_t      *restrict var,
                 cs_real_63_t           *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  /* Use Neumann BC's as default if not provided */

  cs_real_6_t *_bc_coeff_a = NULL;
  cs_real_66_t *_bc_coeff_b = NULL;

  if (bc_coeff_a == NULL) {
    BFT_MALLOC(_bc_coeff_a, n_b_faces, cs_real_6_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 6; j++)
        _bc_coeff_a[i][j] = 0;
    }
    bc_coeff_a = (const cs_real_6_t *)_bc_coeff_a;
  }
  if (bc_coeff_b == NULL) {
    BFT_MALLOC(_bc_coeff_b, n_b_faces, cs_real_66_t);
    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      for (cs_lnum_t j = 0; j < 6; j++) {
        for (cs_lnum_t k = 0; k < 6; k++)
          _bc_coeff_b[i][j][j] = 1;
      }
    }
    bc_coeff_b = (const cs_real_66_t *)_bc_coeff_b;
  }

  /* Compute gradient */

  if (   gradient_type == CS_GRADIENT_ITER
      || gradient_type == CS_GRADIENT_ITER_OLD) {

    _initialize_tensor_gradient(mesh,
                                fvq,
                                halo_type,
                                inc,
                                bc_coeff_a,
                                bc_coeff_b,
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
                                 bc_coeff_a,
                                 bc_coeff_b,
                                 (const cs_real_6_t *)var,
                                 grad);

  } else if (gradient_type == CS_GRADIENT_LSQ) {

    /* If NO reconstruction are required */

    if (n_r_sweeps <= 1)
      _initialize_tensor_gradient(mesh,
                                  fvq,
                                  halo_type,
                                  inc,
                                  bc_coeff_a,
                                  bc_coeff_b,
                                  var,
                                  grad);

    /* Reconstruction with least squares method */

    else
      _lsq_tensor_gradient(mesh,
                           cs_glob_mesh_adjacencies,
                           fvq,
                           halo_type,
                           inc,
                           bc_coeff_a,
                           bc_coeff_b,
                           (const cs_real_6_t *)var,
                           NULL, /* c_weight */
                           grad);

  }

  BFT_FREE(_bc_coeff_a);
  BFT_FREE(_bc_coeff_b);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the steady balance due to porous modelling for the pressure
 * gradient
 *----------------------------------------------------------------------------*/

void CS_PROCF (grdpor, GRDPOR)
(
 const cs_int_t   *const inc          /* <-- 0 or 1: increment or not         */
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  const cs_halo_t  *halo = m->halo;

  const cs_real_t *restrict cell_f_vol = mq->cell_f_vol;
  cs_real_2_t *i_f_face_factor = mq->i_f_face_factor;
  cs_real_t *b_f_face_factor = mq->b_f_face_factor;
  cs_real_t *i_massflux = cs_field_by_name("inner_mass_flux")->val;
  cs_real_t *b_massflux = cs_field_by_name("boundary_mass_flux")->val;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)mq->i_f_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)mq->b_f_face_normal;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict i_f_face_surf = mq->i_f_face_surf;
  const cs_real_t *restrict i_face_surf = mq->i_face_surf;
  const cs_real_t *restrict b_f_face_surf = mq->b_f_face_surf;
  const cs_real_t *restrict b_face_surf = mq->b_face_surf;

  const int *restrict c_disable_flag = mq->c_disable_flag;
  int has_dc = mq->has_disable_flag; /* Has cells disabled? */

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  /*Additional terms due to porosity */
  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  if (f_i_poro_duq_0 == NULL)
    return;

  cs_real_t *i_poro_duq_0 = f_i_poro_duq_0->val;
  cs_real_t *i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
  cs_real_t *b_poro_duq = cs_field_by_name("b_poro_duq")->val;
  cs_real_3_t *c_poro_div_duq
    = (cs_real_3_t *restrict)cs_field_by_name("poro_div_duq")->val;

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      c_poro_div_duq[cell_id][i] = 0.;
  }

  if (*inc == 1) {

    /* Inner faces corrections */
    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_3_t normal;

          cs_math_3_normalise(i_face_normal[face_id], normal);

          cs_real_t *vel_i = &(CS_F_(vel)->val_pre[3*ii]);
          cs_real_t *vel_j = &(CS_F_(vel)->val_pre[3*jj]);

          cs_real_t veli_dot_n =    (1. - i_f_face_factor[face_id][0])
                                  * cs_math_3_dot_product(vel_i, normal);
          cs_real_t velj_dot_n =    (1. - i_f_face_factor[face_id][1])
                                  * cs_math_3_dot_product(vel_j, normal);

          cs_real_t d_f_surf = 0.;
          /* Is the cell disabled (for solid or porous)? Not the case if coupled */
          if (has_dc *  c_disable_flag[has_dc * ii] == 0
              && has_dc * c_disable_flag[has_dc * jj] == 0)
            d_f_surf = 1. / CS_MAX(i_f_face_surf[face_id],
                                   cs_math_epzero * i_face_surf[face_id]);

          i_poro_duq_0[face_id] = veli_dot_n * i_massflux[face_id] * d_f_surf;
          i_poro_duq_1[face_id] = velj_dot_n * i_massflux[face_id] * d_f_surf;

          for (cs_lnum_t i = 0; i < 3; i++) {
            c_poro_div_duq[ii][i] +=   i_poro_duq_0[face_id]
                                     * i_f_face_normal[face_id][i];
            c_poro_div_duq[jj][i] -=   i_poro_duq_1[face_id]
                                    * i_f_face_normal[face_id][i];
          }
        }
      }

    }

    /* Boundary faces corrections */
    for (int g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          cs_real_3_t normal;

          cs_math_3_normalise(b_face_normal[face_id], normal);

          cs_real_t *vel_i = &(CS_F_(vel)->val_pre[3*ii]);

          cs_real_t veli_dot_n =   (1. - b_f_face_factor[face_id])
                                 * cs_math_3_dot_product(vel_i, normal);

          cs_real_t d_f_surf = 0.;
          /* Is the cell disabled (for solid or porous)? Not the case if coupled */
          if (has_dc * c_disable_flag[has_dc * ii] == 0)
            d_f_surf = 1. / CS_MAX(b_f_face_surf[face_id],
                                   cs_math_epzero * b_face_surf[face_id]);

          b_poro_duq[face_id] = veli_dot_n * b_massflux[face_id] * d_f_surf;

          for (cs_lnum_t i = 0; i < 3; i++)
            c_poro_div_duq[ii][i] +=   b_poro_duq[face_id]
                                     * b_f_face_normal[face_id][i];
        }

        /* Finalisation of cell terms */
        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
          /* Is the cell disabled (for solid or porous)? Not the case if coupled */
          cs_real_t dvol = 0.;
          if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
            dvol = 1. / cell_f_vol[cell_id];

          for (cs_lnum_t i = 0; i < 3; i++)
            c_poro_div_duq[cell_id][i] *= dvol;
        }
      }
    }

    /* Handle parallelism and periodicity */
    if (halo != NULL)
      cs_halo_sync_var_strided(halo,
                               CS_HALO_STANDARD,
                               (cs_real_t *)c_poro_div_duq,
                               3);

  }
  else {
#   pragma omp parallel for
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_poro_duq_0[face_id] = 0.;
      i_poro_duq_1[face_id] = 0.;
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdcel, CGDCEL)
(
 const cs_int_t   *const f_id,        /* <-- field id                         */
 const cs_int_t   *const imrgra,      /* <-- gradient computation mode        */
 const cs_int_t   *const inc,         /* <-- 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* <-- 1 or 0: recompute COCG or not    */
 const cs_int_t   *const n_r_sweeps,  /* <-- >1: with reconstruction          */
 const cs_int_t   *const idimtr,      /* <-- 0, 1, 2: scalar, vector, tensor
                                             in case of rotation              */
 const cs_int_t   *const iphydp,      /* <-- use hydrosatatic pressure        */
 const cs_int_t   *const ipond,       /* <-- >0: weighted gradient computation*/
 const cs_int_t   *const iwarnp,      /* <-- verbosity level                  */
 const cs_int_t   *const imligp,      /* <-- type of clipping                 */
 const cs_real_t  *const epsrgp,      /* <-- precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* <-- extrapolate gradient at boundary */
 const cs_real_t  *const climgp,      /* <-- clipping coefficient             */
       cs_real_3_t       f_ext[],     /* <-- exterior force generating the
                                             hydrostatic pressure             */
 const cs_real_t         coefap[],    /* <-- boundary condition term          */
 const cs_real_t         coefbp[],    /* <-- boundary condition term          */
       cs_real_t         pvar[],      /* <-- gradient's base variable         */
       cs_real_t         ktvar[],     /* <-- gradient coefficient variable    */
       cs_real_3_t       grad[]       /* <-> gradient                         */
)
{
  cs_real_t *c_weight = (*ipond > 0) ? ktvar : NULL;

  bool recompute_cocg = (*iccocg) ? true : false;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  char var_name[32];
  if (*f_id > -1) {
    cs_field_t *f = cs_field_by_id(*f_id);
    snprintf(var_name, 31, "%s", f->name);
  }
  else
    strcpy(var_name, "Work array");
  var_name[31] = '\0';

  /* Choose gradient type */

  cs_gradient_type_by_imrgra(*imrgra,
                             &gradient_type,
                             &halo_type);

  /* Check if given field has internal coupling  */
  cs_internal_coupling_t  *cpl = NULL;
  if (*f_id > -1) {
    const int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      const cs_field_t *f = cs_field_by_id(*f_id);
      int coupl_id = cs_field_get_key_int(f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }
  }

  /* Compute gradient */

  cs_gradient_scalar(var_name,
                     gradient_type,
                     halo_type,
                     *inc,
                     recompute_cocg,
                     *n_r_sweeps,
                     *idimtr,
                     *iphydp,
                     1,             /* w_stride */
                     *iwarnp,
                     *imligp,
                     *epsrgp,
                     *extrap,
                     *climgp,
                     f_ext,
                     coefap,
                     coefbp,
                     pvar,
                     c_weight,
                     cpl,
                     grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of vector field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdvec, CGDVEC)
(
 const cs_int_t         *const f_id,
 const cs_int_t         *const imrgra,    /* <-- gradient computation mode    */
 const cs_int_t         *const inc,       /* <-- 0 or 1: increment or not     */
 const cs_int_t         *const n_r_sweeps,    /* <-- >1: with reconstruction      */
 const cs_int_t         *const iwarnp,    /* <-- verbosity level              */
 const cs_int_t         *const imligp,    /* <-- type of clipping             */
 const cs_real_t        *const epsrgp,    /* <-- precision for iterative
                                                 gradient calculation         */
 const cs_real_t        *const climgp,    /* <-- clipping coefficient         */
 const cs_real_3_t             coefav[],  /* <-- boundary condition term      */
 const cs_real_33_t            coefbv[],  /* <-- boundary condition term      */
       cs_real_3_t             pvar[],    /* <-- gradient's base variable     */
       cs_real_33_t            grad[]     /* <-> gradient of the variable
                                                 (du_i/dx_j : grad[][i][j])  */
)
{
  char var_name[32];

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(*imrgra,
                             &gradient_type,
                             &halo_type);

  if (*f_id > -1)
    snprintf(var_name, 31, "Field %2d", *f_id);
  else
    strcpy(var_name, "Work array");
  var_name[31] = '\0';

  /* Check if given field has internal coupling  */
  cs_internal_coupling_t  *cpl = NULL;
  if (*f_id > -1) {
    const int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      const cs_field_t *f = cs_field_by_id(*f_id);
      int coupl_id = cs_field_get_key_int(f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }
  }

  cs_gradient_vector(var_name,
                     gradient_type,
                     halo_type,
                     *inc,
                     *n_r_sweeps,
                     *iwarnp,
                     *imligp,
                     *epsrgp,
                     *climgp,
                     coefav,
                     coefbv,
                     pvar,
                     NULL, /* weighted gradient */
                     cpl,
                     grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of tensor field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdts, CGDTS)
(
 const cs_int_t         *const f_id,
 const cs_int_t         *const imrgra,    /* <-- gradient computation mode    */
 const cs_int_t         *const inc,       /* <-- 0 or 1: increment or not     */
 const cs_int_t         *const n_r_sweeps,    /* <-- >1: with reconstruction  */
 const cs_int_t         *const iwarnp,    /* <-- verbosity level              */
 const cs_int_t         *const imligp,    /* <-- type of clipping             */
 const cs_real_t        *const epsrgp,    /* <-- precision for iterative
                                                 gradient calculation         */
 const cs_real_t        *const climgp,    /* <-- clipping coefficient         */
 const cs_real_6_t             coefav[],  /* <-- boundary condition term      */
 const cs_real_66_t            coefbv[],  /* <-- boundary condition term      */

       cs_real_6_t             pvar[],    /* <-- gradient's base variable     */
       cs_real_63_t            grad[]    /* <-> gradient of the variable
                                                 (du_i/dx_j : gradv[][i][j])  */
)
{
  char var_name[32];

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(*imrgra,
                             &gradient_type,
                             &halo_type);

  if (*f_id > -1)
    snprintf(var_name, 31, "Field %2d", *f_id);
  else
    strcpy(var_name, "Work array");
  var_name[31] = '\0';

  cs_gradient_tensor(var_name,
                     gradient_type,
                     halo_type,
                     *inc,
                     *n_r_sweeps,
                     *iwarnp,
                     *imligp,
                     *epsrgp,
                     *climgp,
                     coefav,
                     coefbv,
                     pvar,
                     grad);
}

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
  assert(cs_glob_mesh != NULL);

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
  int ii;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Total elapsed time for all gradient computations:  %.3f s\n"),
                _gradient_t_tot.wall_nsec*1e-9);

  /* Free system info */

  for (ii = 0; ii < cs_glob_gradient_n_systems; ii++) {
    _gradient_info_dump(cs_glob_gradient_systems[ii]);
    _gradient_info_destroy(&(cs_glob_gradient_systems[ii]));
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  BFT_FREE(cs_glob_gradient_systems);

  cs_glob_gradient_n_systems = 0;
  cs_glob_gradient_n_max_systems = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       tr_dim          2 for tensor with periodicity of rotation,
 *                                  0 otherwise
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
 * \param[in]       w_stride        stride for weighting coefficient
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in, out]  cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar(const char                *var_name,
                   cs_gradient_type_t         gradient_type,
                   cs_halo_type_t             halo_type,
                   int                        inc,
                   bool                       recompute_cocg,
                   int                        n_r_sweeps,
                   int                        tr_dim,
                   int                        hyd_p_flag,
                   int                        w_stride,
                   int                        verbosity,
                   int                        clip_mode,
                   double                     epsilon,
                   double                     extrap,
                   double                     clip_coeff,
                   cs_real_3_t                f_ext[],
                   const cs_real_t            bc_coeff_a[],
                   const cs_real_t            bc_coeff_b[],
                   cs_real_t        *restrict var,
                   cs_real_t        *restrict c_weight,
                   cs_internal_coupling_t    *cpl,
                   cs_real_3_t      *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  /* Synchronize variable */

  if (mesh->halo != NULL) {

    if (tr_dim > 0)
      cs_halo_sync_component(mesh->halo, halo_type,
                             CS_HALO_ROTATION_IGNORE, var);
    else
      cs_halo_sync_var(mesh->halo, halo_type, var);

    if (c_weight != NULL) {
      if (w_stride == 6) {
        cs_halo_sync_var_strided(mesh->halo, halo_type, c_weight, 6);
        cs_halo_perio_sync_var_sym_tens(mesh->halo, halo_type, c_weight);
      }
      else
        cs_halo_sync_var(mesh->halo, halo_type, c_weight);
    }

    if (hyd_p_flag == 1) {
      cs_halo_sync_var_strided(mesh->halo, halo_type, (cs_real_t *)f_ext, 3);
      cs_halo_perio_sync_var_vect(mesh->halo, halo_type, (cs_real_t *)f_ext, 3);
    }

  }

  _gradient_scalar(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   recompute_cocg,
                   n_r_sweeps,
                   tr_dim,
                   hyd_p_flag,
                   w_stride,
                   verbosity,
                   clip_mode,
                   epsilon,
                   extrap,
                   clip_coeff,
                   f_ext,
                   bc_coeff_a,
                   bc_coeff_b,
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
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in, out]  cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
                                    (\f$ \der{u_i}{x_j} \f$ is grad[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector(const char                *var_name,
                   cs_gradient_type_t         gradient_type,
                   cs_halo_type_t             halo_type,
                   int                        inc,
                   int                        n_r_sweeps,
                   int                        verbosity,
                   int                        clip_mode,
                   double                     epsilon,
                   double                     clip_coeff,
                   const cs_real_3_t          bc_coeff_a[],
                   const cs_real_33_t         bc_coeff_b[],
                   cs_real_3_t      *restrict var,
                   cs_real_t        *restrict c_weight,
                   cs_internal_coupling_t    *cpl,
                   cs_real_33_t     *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true) {
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (mesh->halo != NULL) {

    cs_halo_sync_var_strided(mesh->halo, halo_type, (cs_real_t *)var, 3);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_vect(mesh->halo, halo_type, (cs_real_t *)var, 3);

    if (c_weight != NULL)
      cs_halo_sync_var(mesh->halo, halo_type, c_weight);

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
                   bc_coeff_a,
                   bc_coeff_b,
                   (const cs_real_3_t *)var,
                   (const cs_real_t *)c_weight,
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
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad           gradient
                                    (\f$ \der{u_i}{x_j} \f$ is grad[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor(const char                *var_name,
                   cs_gradient_type_t         gradient_type,
                   cs_halo_type_t             halo_type,
                   int                        inc,
                   int                        n_r_sweeps,
                   int                        verbosity,
                   int                        clip_mode,
                   double                     epsilon,
                   double                     clip_coeff,
                   const cs_real_6_t          bc_coeff_a[],
                   const cs_real_66_t         bc_coeff_b[],
                   cs_real_6_t      *restrict var,
                   cs_real_63_t     *restrict grad)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  bool update_stats = true;

  t0 = cs_timer_time();

  if (update_stats == true) {
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (mesh->halo != NULL) {
    cs_halo_sync_var_strided(mesh->halo, halo_type, (cs_real_t *)var, 6);
    if (mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens(mesh->halo, halo_type, (cs_real_t *)var);
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
                   bc_coeff_a,
                   bc_coeff_b,
                   (const cs_real_6_t *)var,
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
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * This variant of the \ref cs_gradient_scalar function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       tr_dim          2 for tensor with periodicity of rotation,
 *                                  0 otherwise
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
 * \param[in]       w_stride        stride for weighting coefficient
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in]       var             gradient's base variable
 * \param[in]       c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in, out]  cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_synced_input(const char                 *var_name,
                                cs_gradient_type_t          gradient_type,
                                cs_halo_type_t              halo_type,
                                int                         inc,
                                bool                        recompute_cocg,
                                int                         n_r_sweeps,
                                int                         tr_dim,
                                int                         hyd_p_flag,
                                int                         w_stride,
                                int                         verbosity,
                                int                         clip_mode,
                                double                      epsilon,
                                double                      extrap,
                                double                      clip_coeff,
                                cs_real_t                   f_ext[][3],
                                const cs_real_t             bc_coeff_a[],
                                const cs_real_t             bc_coeff_b[],
                                const cs_real_t             var[restrict],
                                const cs_real_t             c_weight[restrict],
                                const cs_internal_coupling_t  *cpl,
                                cs_real_t                   grad[restrict][3])
{
  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  bool update_stats = true;

  if (hyd_p_flag == 1) {
    if (cs_glob_mesh->halo != NULL) {
      cs_halo_sync_var_strided(cs_glob_mesh->halo, halo_type,
                               (cs_real_t *)f_ext, 3);
      if (cs_glob_mesh->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(cs_glob_mesh->halo, halo_type,
                                    (cs_real_t *)f_ext, 3);
    }
  }

  t0 = cs_timer_time();

  if (update_stats == true)
    gradient_info = _find_or_add_system(var_name, gradient_type);

  _gradient_scalar(var_name,
                   gradient_info,
                   gradient_type,
                   halo_type,
                   inc,
                   recompute_cocg,
                   n_r_sweeps,
                   tr_dim,
                   hyd_p_flag,
                   w_stride,
                   verbosity,
                   clip_mode,
                   epsilon,
                   extrap,
                   clip_coeff,
                   f_ext,
                   bc_coeff_a,
                   bc_coeff_b,
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
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        weighted gradient coefficient variable,
 *                                  or NULL
 * \param[in, out]  cpl             structure associated with internal coupling,
 *                                  or NULL
 * \param[out]      grad            gradient
                                    (\f$ \der{u_i}{x_j} \f$ is grad[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_synced_input(const char                *var_name,
                                cs_gradient_type_t         gradient_type,
                                cs_halo_type_t             halo_type,
                                int                        inc,
                                int                        n_r_sweeps,
                                int                        verbosity,
                                int                        clip_mode,
                                double                     epsilon,
                                double                     clip_coeff,
                                const cs_real_t            bc_coeff_a[][3],
                                const cs_real_t            bc_coeff_b[][3][3],
                                const cs_real_t            var[restrict][3],
                                const cs_real_t            c_weight[restrict],
                                const cs_internal_coupling_t    *cpl,
                                cs_real_33_t     *restrict grad)
{
  cs_gradient_info_t *gradient_info = NULL;
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
                   bc_coeff_a,
                   bc_coeff_b,
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
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_synced_input(const char                *var_name,
                                cs_gradient_type_t         gradient_type,
                                cs_halo_type_t             halo_type,
                                int                        inc,
                                int                        n_r_sweeps,
                                int                        verbosity,
                                int                        clip_mode,
                                double                     epsilon,
                                double                     clip_coeff,
                                const cs_real_t            bc_coeff_a[][6],
                                const cs_real_t            bc_coeff_b[][6][6],
                                const cs_real_t            var[restrict][6],
                                cs_real_63_t     *restrict grad)
{
  cs_gradient_info_t *gradient_info = NULL;
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
                   bc_coeff_a,
                   bc_coeff_b,
                   var,
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
  *gradient_type = CS_GRADIENT_ITER;

  switch (CS_ABS(imrgra)) {
  case 0:
    *gradient_type = CS_GRADIENT_ITER;
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
    *gradient_type = CS_GRADIENT_LSQ_ITER;
    break;
  case 5:
  case 6:
    *gradient_type = CS_GRADIENT_LSQ_ITER;
    *halo_type = CS_HALO_EXTENDED;
    break;
  case 10:
    *gradient_type = CS_GRADIENT_ITER_OLD;
    break;
  default: break;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
