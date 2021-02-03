/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  BFT headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"
#include "cs_equation.h"
#include "cs_evaluate.h"
#include "cs_fp_exception.h"
#include "cs_iter_algo.h"
#include "cs_matrix_default.h"
#include "cs_navsto_coupling.h"
#include "cs_parall.h"
#include "cs_saddle_itsol.h"
#include "cs_sles.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_monolithic_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic_sles.c
 *
 * \brief Functions dedicated to to the linear algebra settings and operations
 *        in case of CDO face-based schemes with a monolithic velocity-pressure
 *        coupling
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_MONOLITHIC_SLES_DBG      0

/* GKB advanced settings */

#define CS_GKB_TRUNCATION_THRESHOLD       5

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* This structure follow notations given in the article entitled
 * "An iterative generalized Golub-Kahan algorithm for problems in structural
 *  mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 *
 * M space is isomorphic to the velocity space (size = 3.n_faces)
 * N space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */
  cs_real_t                gamma;

  /* Size of spaces */
  cs_lnum_t                n_u_dofs; /* Size of the space M */
  cs_lnum_t                n_p_dofs; /* Size of the space N */

  /* Vector transformation */
  cs_real_t               *b_tilda;  /* Modified RHS */
  cs_real_t               *u_tilda;  /* Modified velocity unknown */

  /* Auxiliary vectors */
  cs_real_t               *q;        /* vector iterates in space N */
  cs_real_t               *d;        /* vector iterates in space N */
  cs_real_t               *d__v;     /* buffer in space N */
  cs_real_t               *dt_q;     /* buffer in space M */
  cs_real_t               *m__v;     /* vector iterates in space M */
  cs_real_t               *v;        /* vector iterates in space M */

  /* Orthogonalization coefficients */
  cs_real_t                alpha;
  cs_real_t                beta;
  cs_real_t                zeta;

  /* Store z_size zeta coefficients */
  int                      z_size;
  cs_real_t               *zeta_array;
  cs_real_t                zeta_square_sum;

  cs_iter_algo_info_t     *info;     /* Information related to the convergence
                                        of the algorithm */

} cs_gkb_builder_t;

/* This structure is used to manage the Uzawa algorithm and its variants
 *
 * U space is isomorphic to the velocity space (size = 3.n_faces)
 * P space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */
  cs_real_t               gamma;

  /* Size of spaces */
  cs_lnum_t               n_u_dofs; /* Size of the space U */
  cs_lnum_t               n_p_dofs; /* Size of the space P */

  /* Vector transformation */
  cs_real_t              *b_tilda;  /* Modified RHS (size U) */

  /* Auxiliary scaling coefficient */
  cs_real_t               alpha;

  /* Auxiliary vectors */
  cs_real_t              *inv_mp;   /* reciprocal of the pressure mass matrix */
  cs_real_t              *res_p;    /* buffer in space P */
  cs_real_t              *d__v;     /* buffer in space P */
  cs_real_t              *gk;       /* buffer in space P */
  cs_real_t              *dzk;      /* buffer in space U */
  cs_real_t              *rhs;      /* buffer in space U */

  cs_iter_algo_info_t    *info;     /* Information related to the convergence
                                       of the algorithm */

} cs_uza_builder_t;

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_range_set_t         *cs_shared_range_set;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a norm for a scalar-valued cell-based array "a"
 *        The parallel synchronization is performed inside this function
 *
 * \param[in]    a    array of size n_cells
 *
 * \return the computed norm
 */
/*----------------------------------------------------------------------------*/

static inline double
_get_cbscal_norm(cs_real_t  *a)
{
  double norm2 = cs_dot_xx(cs_shared_quant->n_cells, a);

  cs_parall_sum(1, CS_DOUBLE, &norm2);
  assert(norm2 > -DBL_MIN);

  return sqrt(norm2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a norm for a face-based "a" v with 3*n_faces elements
 *        The parallel synchronization is performed inside this function
 *
 * \param[in]    a    array of size 3*n_faces
 *
 * \return the computed norm
 */
/*----------------------------------------------------------------------------*/

static inline double
_get_fbvect_norm(cs_real_t  *a)
{
  double norm2 = cs_evaluate_3_square_wc2x_norm(a,
                                                cs_shared_connect->c2f,
                                                cs_shared_quant->pvol_fc);

  assert(norm2 > -DBL_MIN);
  return sqrt(norm2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dot product between two arrays on face unknowns.
 *         One assumes that input arrays are in a "scattered" distribution
 *         So the size should be 3*n_faces.
 *
 * \param[in]       size   size of arrays
 * \param[in, out]  x      first array
 * \param[in, out]  y      second array
 *
 * \return the computed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_face_gdot(cs_lnum_t    size,
           cs_real_t    x[],
           cs_real_t    y[])
{
  CS_UNUSED(size);       /* Avoid a compilation warning in during compilation */
  const cs_range_set_t  *rset = cs_shared_range_set;

  assert(size == rset->n_elts[1]);
  assert(size == 3*cs_shared_quant->n_faces);

  /* x and y are scattered arrays. One assumes that values are synchronized
     across ranks (for instance by using a cs_interface_set_sum()) */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE,/* type */
                      1,           /* stride (treated as scalar up to now) */
                      x,           /* in: size = n_sles_scatter_elts */
                      x);          /* out: size = n_sles_gather_elts */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE,/* type */
                      1,           /* stride (treated as scalar up to now) */
                      y,           /* in: size = n_sles_scatter_elts */
                      y);          /* out: size = n_sles_gather_elts */

  cs_real_t  result = cs_gdot(rset->n_elts[0], x, y);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       x,
                       x);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       y,
                       y);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform a matrix-vector multiplication (matvec is not allocated
 *         when given as input parameter)
 *
 * \param[in]      mat        matrix
 * \param[in]      rset       pointer to a cs_range_set_t structure
 * \param[in]      vec_len    size of the vector length
 * \param[in, out] vec        vector
 * \param[out]     p_matvec   resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_vector_multiply(const cs_matrix_t         *mat,
                        const cs_range_set_t      *rset,
                        cs_lnum_t                  vec_len,
                        cs_real_t                 *vec,
                        cs_real_t                **p_matvec)
{
  if (mat == NULL || vec == NULL)
    return;

  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(mat);

  /* Handle the input array
   * n_rows = n_gather_elts <= n_scatter_elts = n_dofs (mesh view) <= n_cols
   */
  cs_real_t  *vecx = NULL;
  if (n_cols > vec_len) {
    BFT_MALLOC(vecx, n_cols, cs_real_t);
    memcpy(vecx, vec, sizeof(cs_real_t)*vec_len);
  }
  else
    vecx = vec;

  /* scatter view to gather view */
  if (rset != NULL)
    cs_range_set_gather(rset,
                        CS_REAL_TYPE,  /* type */
                        1,             /* stride */
                        vecx,          /* in:  size=n_sles_scatter_elts */
                        vecx);         /* out: size=n_sles_gather_elts */

  /* Handle the output array */
  cs_real_t  *matvec = NULL;
  BFT_MALLOC(matvec, n_cols, cs_real_t);
  memset(matvec, 0, sizeof(cs_real_t)*n_cols);

  cs_matrix_vector_multiply(CS_HALO_ROTATION_IGNORE, mat, vecx, matvec);

  /* gather to scatter view (i.e. algebraic to mesh view) */
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       vecx, vec);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       matvec, matvec);

  /* Free allocated memory if needed */
  if (vecx != vec) {
    assert(n_cols > vec_len);
    BFT_FREE(vecx);
  }

  /* return the resulting array */
  *p_matvec = matvec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the matrix for an approximation of the Schur complement based
 *         on the inverse of the diagonal of the velocity block
 *
 * \param[in]   nsp     pointer to a cs_navsto_param_t structure
 * \param[in]   a       (MSR) matrix for the velocity block
 * \param[out]  diagK   double pointer to the diagonal coefficients
 * \param[out]  xtraK   double pointer to the extra-diagonal coefficients
 *
 * \return a pointer to a the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_diag_schur_approximation(const cs_navsto_param_t   *nsp,
                          const cs_matrix_t         *a,
                          cs_uza_builder_t          *uza,
                          cs_real_t                **p_diagK,
                          cs_real_t                **p_xtraK)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Compute scaling coefficients */
  const cs_real_t  rho0 = nsp->mass_density->ref_value;
  cs_real_t  *visc_val = NULL;
  int  visc_stride = 0;
  if (nsp->turbulence->model->iturb == CS_TURB_NONE) {
    BFT_MALLOC(visc_val, 1, cs_real_t);
    visc_val[0] = nsp->lam_viscosity->ref_value;
  }
  else {
    visc_stride = 1;
    BFT_MALLOC(visc_val, n_cells, cs_real_t);
    cs_property_eval_at_cells(ts->t_cur, nsp->tot_viscosity, visc_val);
  }

  cs_real_t  alpha = 1/ts->dt[0];
  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    alpha = 0.01*nsp->lam_viscosity->ref_value;
  uza->alpha = rho0*alpha;

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    uza->inv_mp[ip] = visc_val[visc_stride*ip]/quant->cell_vol[ip];

  BFT_FREE(visc_val);

  /* Synchronize the diagonal values for A */
  const cs_real_t  *diagA = NULL;
  cs_real_t  *_diagA = NULL;

  if (cs_glob_n_ranks > 1) {

    size_t  size = 3*quant->n_faces;
    BFT_MALLOC(_diagA, size, cs_real_t);
    cs_range_set_scatter(cs_shared_range_set,
                         CS_REAL_TYPE,
                         1,     /* treated as scalar-valued up to now */
                         cs_matrix_get_diagonal(a), /* gathered view */
                         _diagA);

    diagA = _diagA; /* scatter view (synchronized)*/
  }
  else {
    diagA = cs_matrix_get_diagonal(a);
    assert(m->periodicity == NULL); /* TODO */
  }

  /* Native format for the Schur approximation matrix */
  cs_real_t   *diagK = NULL;
  cs_real_t   *xtraK = NULL;

  BFT_MALLOC(diagK, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtraK, 2*n_i_faces, cs_real_t);

  memset(diagK, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtraK, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = diagA + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1/a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */
    cs_real_t  *_xtraK = xtraK + 2*f_id;
    _xtraK[0] = _xtraK[1] = contrib;

    /* Diagonal contributions */
    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diagK[cell_i] -= contrib;
    diagK[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/
  const cs_real_t  *diagA_shift = diagA + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diagA_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1/a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */
    diagK[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */
  cs_lnum_t  db_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */
  cs_lnum_t  eb_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */
  cs_matrix_t  *K = cs_matrix_native(false, /* symmetry */
                                     db_size,
                                     eb_size);

  cs_matrix_set_coefficients(K, false, /* symmetry */
                             db_size, eb_size,
                             n_i_faces, i_face_cells,
                             diagK, xtraK);

  /* Return arrays (to be freed when the algorithm is converged) */
  *p_diagK = diagK;
  *p_xtraK = xtraK;

  BFT_FREE(_diagA);

  return K;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the matrix for an approximation of the Schur complement based
 *         on the inverse of the sum of the absolute values of the velocity
 *         block
 *
 * \param[in]   nsp     pointer to a cs_navsto_param_t structure
 * \param[in]   a       (MSR) matrix for the velocity block
 * \param[out]  diagK   double pointer to the diagonal coefficients
 * \param[out]  xtraK   double pointer to the extra-diagonal coefficients
 *
 * \return a pointer to a the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_invlumped_schur_approximation(const cs_navsto_param_t     *nsp,
                               const cs_equation_param_t   *eqp,
                               cs_cdofb_monolithic_sles_t  *msles,
                               const cs_matrix_t           *a,
                               cs_uza_builder_t            *uza,
                               cs_real_t                  **p_diagK,
                               cs_real_t                  **p_xtraK)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Compute scaling coefficients */
  const cs_real_t  rho0 = nsp->mass_density->ref_value;
  cs_real_t  *visc_val = NULL;
  int  visc_stride = 0;
  if (nsp->turbulence->model->iturb == CS_TURB_NONE) {
    BFT_MALLOC(visc_val, 1, cs_real_t);
    visc_val[0] = nsp->lam_viscosity->ref_value;
  }
  else {
    visc_stride = 1;
    BFT_MALLOC(visc_val, n_cells, cs_real_t);
    cs_property_eval_at_cells(ts->t_cur, nsp->tot_viscosity, visc_val);
  }

  cs_real_t  alpha = 1/ts->dt[0];
  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    alpha = 0.01*nsp->lam_viscosity->ref_value;
  uza->alpha = rho0*alpha;

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    uza->inv_mp[ip] = visc_val[visc_stride*ip]/quant->cell_vol[ip];

  BFT_FREE(visc_val);

  /* Compute A^-1 lumped */

  /* Modify the tolerance in order to be more accurate on this step */
  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":inv_lumped") + 1, char);
  sprintf(system_name, "%s:inv_lumped", eqp->name);

  cs_param_sles_t  *slesp0 = cs_param_sles_create(-1, system_name);

  cs_param_sles_copy_from(eqp->sles_param, slesp0);
  slesp0->eps = 1e-1; /* Only a coarse approximation is needed */
  slesp0->n_max_iter = 50;

  for (cs_lnum_t i = 0; i < uza->n_u_dofs; i++)
    uza->rhs[i] = 1;

  cs_real_t  *invA_lumped = NULL;
  BFT_MALLOC(invA_lumped, uza->n_u_dofs, cs_real_t);
  memset(invA_lumped, 0, sizeof(cs_real_t)*uza->n_u_dofs);

  uza->info->n_inner_iter
    += (uza->info->last_inner_iter =
        cs_equation_solve_scalar_system(uza->n_u_dofs,
                                        slesp0,
                                        a,
                                        cs_shared_range_set,
                                        1,     /* no normalization */
                                        false, /* rhs_redux --> already done */
                                        msles->sles,
                                        invA_lumped,
                                        uza->rhs));
  /* Partial memory free */
  BFT_FREE(system_name);
  cs_param_sles_free(&slesp0);

  /* Native format for the Schur approximation matrix */
  cs_real_t   *diagK = NULL;
  cs_real_t   *xtraK = NULL;

  BFT_MALLOC(diagK, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtraK, 2*n_i_faces, cs_real_t);

  memset(diagK, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtraK, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = invA_lumped + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */
    cs_real_t  *_xtraK = xtraK + 2*f_id;
    _xtraK[0] = _xtraK[1] = contrib;

    /* Diagonal contributions */
    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diagK[cell_i] -= contrib;
    diagK[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/
  cs_real_t  *diagA_shift = invA_lumped + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diagA_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */
    diagK[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */
  cs_lnum_t  db_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */
  cs_lnum_t  eb_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */
  cs_matrix_t  *K = cs_matrix_native(false, /* symmetry */
                                     db_size,
                                     eb_size);

  cs_matrix_set_coefficients(K, false, /* symmetry */
                             db_size, eb_size,
                             n_i_faces, i_face_cells,
                             diagK, xtraK);

  /* Return arrays (to be freed when the algorithm is converged) */
  *p_diagK = diagK;
  *p_xtraK = xtraK;

  BFT_FREE(invA_lumped);

  return K;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_diag_schur_sbp(const cs_navsto_param_t       *nsp,
                const cs_saddle_system_t      *ssys,
                cs_saddle_block_precond_t     *sbp)
{
  CS_UNUSED(nsp);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t  b11_size = ssys->x1_size;

  /* Synchronize the diagonal values for the block m11 */
  cs_matrix_t  *m11 = ssys->m11_matrices[0];
  const cs_real_t  *diag_m11 = NULL;
  cs_real_t  *_diag_m11 = NULL;

  if (cs_shared_range_set != NULL) {

    BFT_MALLOC(_diag_m11, b11_size, cs_real_t);
    cs_range_set_scatter(cs_shared_range_set,
                         CS_REAL_TYPE,
                         1,     /* treated as scalar-valued up to now */
                         cs_matrix_get_diagonal(m11), /* gathered view */
                         _diag_m11);

    diag_m11 = _diag_m11; /* scatter view (synchronized)*/

  }
  else {

    diag_m11 = cs_matrix_get_diagonal(m11);
    assert(m->periodicity == NULL && cs_glob_n_ranks == 1);

  }

  /* Native format for the Schur approximation matrix */
  cs_real_t   *diagK = NULL;
  cs_real_t   *xtraK = NULL;

  BFT_MALLOC(diagK, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtraK, 2*n_i_faces, cs_real_t);

  memset(diagK, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtraK, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = diag_m11 + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1./a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */
    cs_real_t  *_xtraK = xtraK + 2*f_id;
    _xtraK[0] = _xtraK[1] = contrib;

    /* Diagonal contributions */
    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diagK[cell_i] -= contrib;
    diagK[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/
  const cs_real_t  *diag_shift = diag_m11 + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diag_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1./a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */
    diagK[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */
  cs_lnum_t  db_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */
  cs_lnum_t  eb_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */
  sbp->schur_matrix = cs_matrix_native(false, /* symmetry */
                                       db_size,
                                       eb_size);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             db_size, eb_size,
                             n_i_faces, i_face_cells,
                             diagK, xtraK);

  /* Return arrays (to be freed when the algorithm is converged) */
  sbp->schur_diag = diagK;
  sbp->schur_xtra = xtraK;

  BFT_FREE(_diag_m11);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a scaled mass matrix (on the pressure space) and a scaling
 *        coefficient for the compatible Laplacian
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_scaled_mass_sbp(const cs_navsto_param_t       *nsp,
                 const cs_saddle_system_t      *ssys,
                 cs_saddle_block_precond_t     *sbp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;

  assert(ssys->x2_size == n_cells);
  BFT_MALLOC(sbp->mass22_diag, n_cells, cs_real_t);

  /* Compute scaling coefficients */
  if (nsp->turbulence->model->iturb == CS_TURB_NONE) {

    const cs_real_t  visc_val = nsp->lam_viscosity->ref_value;
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n_cells; i2++)
      sbp->mass22_diag[i2] = visc_val/quant->cell_vol[i2];

  }
  else {

    cs_property_eval_at_cells(ts->t_cur, nsp->tot_viscosity, sbp->mass22_diag);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n_cells; i2++)
      sbp->mass22_diag[i2] /= quant->cell_vol[i2];

  }

  const cs_real_t  rho0 = nsp->mass_density->ref_value;
  cs_real_t  alpha = 1/ts->dt[0];
  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    alpha = 0.01*nsp->lam_viscosity->ref_value;
  sbp->schur_scaling = rho0*alpha;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_elman_schur_sbp(const cs_navsto_param_t       *nsp,
                 const cs_saddle_system_t      *ssys,
                 cs_saddle_block_precond_t     *sbp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t  b11_size = ssys->x1_size;
  const cs_lnum_t  schur_size = ssys->x2_size;

  /* Native format for the Schur approximation matrix */
  cs_real_t   *diagK = NULL;
  cs_real_t   *xtraK = NULL;

  BFT_MALLOC(diagK, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtraK, 2*n_i_faces, cs_real_t);

  memset(diagK, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtraK, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */
    cs_real_t  *_xtraK = xtraK + 2*f_id;
    _xtraK[0] = _xtraK[1] = contrib;

    /* Diagonal contributions */
    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diagK[cell_i] -= contrib;
    diagK[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */
    diagK[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */
  cs_lnum_t  db_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */
  cs_lnum_t  eb_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */
  sbp->schur_matrix = cs_matrix_native(false, /* symmetry */
                                       db_size,
                                       eb_size);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             db_size, eb_size,
                             n_i_faces, i_face_cells,
                             diagK, xtraK);

  /* Return arrays (to be freed when the algorithm is converged) */
  sbp->schur_diag = diagK;
  sbp->schur_xtra = xtraK;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_invlumped_schur_sbp(const cs_navsto_param_t       *nsp,
                     const cs_saddle_system_t      *ssys,
                     cs_saddle_block_precond_t     *sbp)
{
  CS_UNUSED(nsp);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t  b11_size = ssys->x1_size;

  /* Compute m11^-1 lumped */

  /* Modify the tolerance in order to be less accurate on this step */
  cs_param_sles_t  *slesp0 = cs_param_sles_create(-1, "schur:inv_lumped");

  cs_param_sles_copy_from(sbp->m11_slesp, slesp0);
  slesp0->eps = 1e-3;
  slesp0->n_max_iter = 50;

  cs_real_t  *rhs = NULL;
  BFT_MALLOC(rhs, b11_size, cs_real_t);
  for (cs_lnum_t i = 0; i < b11_size; i++) rhs[i] = 1;

  cs_real_t  *inv_lumped = NULL;
  BFT_MALLOC(inv_lumped, b11_size, cs_real_t);
  memset(inv_lumped, 0, sizeof(cs_real_t)*b11_size);

  cs_equation_solve_scalar_system(b11_size,
                                  slesp0,
                                  ssys->m11_matrices[0],
                                  ssys->rset,
                                  1,     /* no normalization */
                                  false, /* rhs_redux --> already done */
                                  sbp->m11_sles,
                                  inv_lumped,
                                  rhs);
  /* Partial memory free */
  BFT_FREE(rhs);
  cs_param_sles_free(&slesp0);

  /* Native format for the Schur approximation matrix */
  cs_real_t   *diagK = NULL;
  cs_real_t   *xtraK = NULL;

  BFT_MALLOC(diagK, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtraK, 2*n_i_faces, cs_real_t);

  memset(diagK, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtraK, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = inv_lumped + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */
    cs_real_t  *_xtraK = xtraK + 2*f_id;
    _xtraK[0] = _xtraK[1] = contrib;

    /* Diagonal contributions */
    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diagK[cell_i] -= contrib;
    diagK[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/
  cs_real_t  *diag_shift = inv_lumped + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diag_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */
    diagK[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */
  cs_lnum_t  db_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */
  cs_lnum_t  eb_size[4] = {1, 1, 1, 1}; /* 1, 1, 1, 1*1 */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */
  sbp->schur_matrix = cs_matrix_native(false, /* symmetry */
                                       db_size,
                                       eb_size);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             db_size, eb_size,
                             n_i_faces, i_face_cells,
                             diagK, xtraK);

  /* Return arrays (to be freed when the algorithm is converged) */
  sbp->schur_diag = diagK;
  sbp->schur_xtra = xtraK;

  sbp->m11_inv_diag = inv_lumped;
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the main iterative solver for the velocity block
 *
 * \param[in]      model    type of model related to the Navsto system
 * \param[in]      nslesp   set of parameter for the monolithic SLES
 * \param[in]      slesp    set of parameters for the velocity SLES
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_set_petsc_main_solver(const cs_navsto_param_model_t   model,
                       const cs_navsto_param_sles_t   *nslesp,
                       const cs_param_sles_t          *slesp,
                       KSP                             ksp)
{
  if (model == CS_NAVSTO_MODEL_STOKES)
    KSPSetType(ksp, KSPFCG);

  else { /* Advection is present, so one needs a more genric iterative solver */

    const int  n_max_restart = 30;

    KSPSetType(ksp, KSPFGMRES);
    KSPGMRESSetRestart(ksp, n_max_restart);

  }

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  max_it;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
  KSPSetTolerances(ksp,
                   nslesp->il_algo_rtol,   /* relative convergence tolerance */
                   abstol,                  /* absolute convergence tolerance */
                   dtol,                    /* divergence tolerance */
                   nslesp->n_max_il_algo_iter); /* max number of iterations */

  switch (slesp->resnorm_type) {

  case CS_PARAM_RESNORM_NORM2_RHS: /* Try to have "true" norm */
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;
  case CS_PARAM_RESNORM_NONE:
    KSPSetNormType(ksp, KSP_NORM_NONE);
    break;

  default:
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;
  }

}

#if defined(PETSC_HAVE_HYPRE)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when BoomerAMG from the HYPRE library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_boomeramg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_no_CF", "");
#else
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_no_CF","");
#endif
}
#endif  /* PETSC_HAVE_HYPRE */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when GAMG from the PETSc library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_gamg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_square_graph", "4");
#else
  PetscOptionsSetValue("-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue("-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_velocity_gamg_square_graph", "4");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generate IndexSet for the PETSc FieldSplit preconditioner
 *
 * \param[in, out]  isp     IndexSet for the pressure DoFs
 * \param[in, out]  isv     IndexSet for the velocity DoFs
 */
/*----------------------------------------------------------------------------*/

static void
_build_is_for_fieldsplit(IS   *isp,
                         IS   *isv)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_range_set_t  *rset = cs_shared_range_set;

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  *indices = NULL;

  PetscMalloc1(3*n_faces, &indices);

  /* IndexSet for the velocity DoFs */
  if (rset->n_elts[0] == rset->n_elts[1]) {

    for (PetscInt i = 0; i < 3*n_faces; i++)
      indices[i] = rset->g_id[i];
    ISCreateGeneral(PETSC_COMM_WORLD, 3*n_faces, indices, PETSC_COPY_VALUES,
                    isv);

  }
  else {

    PetscInt  n_velocity_elts = 0;
    for (PetscInt i = 0; i < 3*n_faces; i++) {
      cs_gnum_t  g_id = rset->g_id[i];
      if (g_id >= rset->l_range[0] && g_id < rset->l_range[1])
        indices[n_velocity_elts++] = g_id;
    }
    ISCreateGeneral(PETSC_COMM_WORLD, n_velocity_elts, indices,
                    PETSC_COPY_VALUES, isv);

  }

  /* IndexSet for the velocity DoFs
   * Pressure unknowns are located at cell centers so the treatment should be
   * the same in sequential and parallel computation */
  for (PetscInt i = 0; i < n_cells; i++)
    indices[i] = rset->g_id[i + 3*n_faces];
  ISCreateGeneral(PETSC_COMM_WORLD, n_cells, indices, PETSC_COPY_VALUES,
                  isp);

  PetscFree(indices);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in]   slesp      set of parameters for the linear algebra
 */
/*----------------------------------------------------------------------------*/

static PCType
_petsc_get_pc_type(const cs_param_sles_t    *slesp)
{
  PCType  pc_type = PCNONE;

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_NONE:
    return PCNONE;

  case CS_PARAM_PRECOND_DIAG:
    return PCJACOBI;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
  case CS_PARAM_PRECOND_BJACOB_SGS:
    return PCBJACOBI;

  case CS_PARAM_PRECOND_SSOR:
    return PCSOR;

  case CS_PARAM_PRECOND_ICC0:
    return PCICC;

  case CS_PARAM_PRECOND_ILU0:
    return PCILU;

  case CS_PARAM_PRECOND_AS:
    return PCASM;

  case CS_PARAM_PRECOND_AMG:
    {
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG:
        return PCGAMG;
        break;

      case CS_PARAM_AMG_PETSC_PCMG:
        return PCMG;
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
        return PCHYPRE;
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Switch to MG since BoomerAMG is not available.\n",
                      __func__);
        return PCMG;
#endif

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid AMG type for the PETSc library.", __func__);
        break;

      } /* End of switch on the AMG type */

    } /* AMG as preconditioner */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Preconditioner not interfaced with PETSc.", __func__);
    break;
  }

  return pc_type;
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a GMRES
 *
 * \param[in]      slesp     pointer to a set of SLES settings
 * \param[in]      rtol      relative tolerance to set
 * \param[in]      max_it    max number of iterations
 * \param[in, out] u_ksp     pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_set_velocity_ksp(const cs_param_sles_t   *slesp,
                  PetscReal                rtol,
                  PetscInt                 max_it,
                  KSP                      u_ksp)
{
  PC u_pc;
  KSPGetPC(u_ksp, &u_pc);
  PCType  pc_type = _petsc_get_pc_type(slesp);

  switch (slesp->resnorm_type) {

  case CS_PARAM_RESNORM_NORM2_RHS: /* Try to have "true" norm */
    KSPSetNormType(u_ksp, KSP_NORM_UNPRECONDITIONED);
    break;
  case CS_PARAM_RESNORM_NONE:
    KSPSetNormType(u_ksp, KSP_NORM_NONE);
    break;
  default:
    KSPSetNormType(u_ksp, KSP_NORM_UNPRECONDITIONED);
    break;

  }

  /* Set the solver */
  switch (slesp->solver) {

  case CS_PARAM_ITSOL_NONE:
    KSPSetType(u_ksp, KSPPREONLY);
    break;
  case CS_PARAM_ITSOL_FCG:
  case CS_PARAM_ITSOL_CG:
    KSPSetType(u_ksp, KSPCG);
    break;
  case CS_PARAM_ITSOL_BICG:      /* Improved Bi-CG stab */
    KSPSetType(u_ksp, KSPIBCGS);
    break;
  case CS_PARAM_ITSOL_BICGSTAB2: /* BiCGstab2 */
    KSPSetType(u_ksp, KSPBCGSL);
    break;
  case CS_PARAM_ITSOL_MUMPS:     /* Direct solver (factorization) */
#if defined(PETSC_HAVE_MUMPS)
    KSPSetType(u_ksp, KSPPREONLY);
    PCSetType(u_pc, PCLU);
    PCFactorSetMatSolverType(u_pc, MATSOLVERMUMPS);
#else
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS not interfaced with this installation of PETSc.",
              __func__);
#endif
    break;
  case CS_PARAM_ITSOL_GMRES:
    KSPSetType(u_ksp, KSPGMRES);
    break;
  case CS_PARAM_ITSOL_FGMRES:
    KSPSetType(u_ksp, KSPFGMRES);
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:     /* Direct solver (factorization) */
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver. Try mumps.",
              __func__);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver.", __func__);
    break;

  } /* Switch on solver */

  if (slesp->solver != CS_PARAM_ITSOL_MUMPS)
    PCSetType(u_pc, pc_type);

  /* Additional settings for the preconditioner */
  switch (slesp->precond) {

  case CS_PARAM_PRECOND_AMG:
    switch (slesp->amg_type) {

    case CS_PARAM_AMG_PETSC_GAMG:
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      _setup_velocity_gamg();
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(u_pc, "boomeramg");
      _setup_velocity_boomeramg();

#if 0 /* JB: TEST TO PERFORM IN 3D*/
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#else
  PetscOptionsSetValue(
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#endif
#endif  /* if 0 */

#else
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      _setup_velocity_gamg();
#endif
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
      break;

    } /* AMG type */
    break;

  default:
    break; /* Nothing else to do */

  } /* Switch on preconditioner */

  /* Set tolerance and number of iterations */
  PetscReal _rtol, abstol, dtol;
  PetscInt  _max_it;
  KSPGetTolerances(u_ksp, &_rtol, &abstol, &dtol, &_max_it);
  KSPSetTolerances(u_ksp,
                   rtol,        /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   max_it);     /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_additive_amg_hook(void     *context,
                   KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the velocity block */
  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */
  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_ADDITIVE);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
   * Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the KSP used as preconditioner for the velocity block */
  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative block preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_multiplicative_hook(void     *context,
                     KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  /* Set the main iterative solver for the velocity block */
  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */
  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_MULTIPLICATIVE);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */
  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of diagonal Schur preconditioner by block
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_diag_schur_hook(void     *context,
                 KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the velocity block */
  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */
  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */
  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of upper Schur preconditioner by block
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_upper_schur_hook(void     *context,
                  KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the velocity block */
  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */
  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */
  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

#if PETSC_VERSION_GE(3,11,0)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB as a solver.
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_hook(void     *context,
          KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPPREONLY);

  /* Apply modifications to the KSP structure */
  PC up_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc, 10*nslesp->il_algo_rtol);
  PCFieldSplitSetGKBMaxit(up_pc, nslesp->n_max_il_algo_iter);
  PCFieldSplitSetGKBNu(up_pc, 0);
  PCFieldSplitSetGKBDelay(up_pc, CS_GKB_TRUNCATION_THRESHOLD);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB as a preconditioner.
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_precond_hook(void     *context,
                  KSP       ksp)
{
  IS  isv, isp;

  cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the velocity block */
  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */
  PC up_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc,  slesp->eps);
  PCFieldSplitSetGKBMaxit(up_pc, slesp->n_max_iter);
  PCFieldSplitSetGKBNu(up_pc, 0);
  PCFieldSplitSetGKBDelay(up_pc, CS_GKB_TRUNCATION_THRESHOLD);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  /* Set KSP tolerances */
  PetscInt  max_it = 50;
  PetscReal  rtol = 1e-2;
  _set_velocity_ksp(slesp, rtol, max_it, up_subksp[0]);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* GKB available only if version >= 3.11 */

#if defined(PETSC_HAVE_MUMPS)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of MUMPS via PETSc
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_mumps_hook(void     *context,
            KSP       ksp)
{
  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  PC  pc;
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCLU);
  PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

  PetscReal rtol, abstol, dtol;
  PetscInt  max_it;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
  KSPSetTolerances(ksp,
                   slesp->eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   slesp->n_max_iter); /* max number of iterations */

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* PETSC_HAVE_MUMPS */

#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a GKB builder structure
 *
 * \param[in]  nsp        pointer to a cs_navsto_param_t structure
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 *
 * \return a pointer to a new allocated GKB builder
 */
/*----------------------------------------------------------------------------*/

static cs_gkb_builder_t *
_init_gkb_builder(const cs_navsto_param_t    *nsp,
                  cs_real_t                   gamma,
                  cs_lnum_t                   n_u_dofs,
                  cs_lnum_t                   n_p_dofs)
{
  cs_gkb_builder_t  *gkb = NULL;

  BFT_MALLOC(gkb, 1, cs_gkb_builder_t);

  gkb->gamma = gamma;
  gkb->n_u_dofs = n_u_dofs;
  gkb->n_p_dofs = n_p_dofs;

  /* Vector transformation */
  BFT_MALLOC(gkb->u_tilda, n_u_dofs, cs_real_t);
  /* Rk: b_tilda stores quantities in space M and N alternatively */
  assert(n_u_dofs >= n_p_dofs);
  BFT_MALLOC(gkb->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */
  BFT_MALLOC(gkb->v, n_u_dofs, cs_real_t);
  memset(gkb->v, 0, n_u_dofs*sizeof(cs_real_t));

  BFT_MALLOC(gkb->q, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->dt_q, n_u_dofs, cs_real_t);
  BFT_MALLOC(gkb->m__v, n_u_dofs, cs_real_t);

  /* Orthogonalization coefficients */
  gkb->alpha = gkb->beta = gkb->zeta = 0.;

  /* Convergence members */
  if (gamma < 1)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD + 1;
  else if (gamma < 10)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD;
  else if (gamma < 100)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 1);
  else if (gamma < 1e3)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 2);
  else if (gamma < 1e4)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 3);
  else
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 4);

  BFT_MALLOC(gkb->zeta_array, gkb->z_size, cs_real_t);
  memset(gkb->zeta_array, 0, gkb->z_size*sizeof(cs_real_t));

  gkb->zeta_square_sum = 0.;

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  gkb->info = cs_iter_algo_define(nslesp->il_algo_verbosity,
                                  nslesp->n_max_il_algo_iter,
                                  nslesp->il_algo_atol,
                                  nslesp->il_algo_rtol,
                                  nslesp->il_algo_dtol);

  return gkb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a GKB builder structure
 *
 * \param[in, out]  p_gkb   double pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_gkb_builder(cs_gkb_builder_t   **p_gkb)
{
  cs_gkb_builder_t  *gkb = *p_gkb;

  if (gkb == NULL)
    return;

  BFT_FREE(gkb->b_tilda);
  BFT_FREE(gkb->u_tilda);

  BFT_FREE(gkb->q);
  BFT_FREE(gkb->d);
  BFT_FREE(gkb->d__v);
  BFT_FREE(gkb->dt_q);
  BFT_FREE(gkb->m__v);
  BFT_FREE(gkb->v);

  BFT_FREE(gkb->zeta_array);

  BFT_FREE(gkb->info);

  BFT_FREE(gkb);
  *p_gkb = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a Uzawa builder structure
 *
 * \param[in]  nsp        pointer to a cs_navsto_param_t structure
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 * \param[in]  quant      pointer to additional mesh quantities
 *
 * \return a pointer to a new allocated Uzawa builder
 */
/*----------------------------------------------------------------------------*/

static cs_uza_builder_t *
_init_uzawa_builder(const cs_navsto_param_t      *nsp,
                    cs_real_t                     gamma,
                    cs_lnum_t                     n_u_dofs,
                    cs_lnum_t                     n_p_dofs,
                    const cs_cdo_quantities_t    *quant)
{
  cs_uza_builder_t  *uza = NULL;

  BFT_MALLOC(uza, 1, cs_uza_builder_t);

  uza->alpha = 0;
  uza->gamma = gamma;
  uza->n_u_dofs = n_u_dofs;
  uza->n_p_dofs = n_p_dofs;

  BFT_MALLOC(uza->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */
  BFT_MALLOC(uza->inv_mp, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->res_p, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->rhs, n_u_dofs, cs_real_t);

  uza->gk = NULL;
  uza->dzk = NULL;
  if (nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_CG) {

    /* Since gk is used as a variable in a cell system, one has to take into
       account the space for synchronization */
    cs_lnum_t  size = n_p_dofs;
    if (cs_glob_n_ranks > 1)
      size = CS_MAX(n_p_dofs, cs_glob_mesh->n_cells_with_ghosts);
    BFT_MALLOC(uza->gk, size, cs_real_t);

    BFT_MALLOC(uza->dzk, n_u_dofs, cs_real_t);

  }
  else  {

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      uza->inv_mp[ip] = 1./quant->cell_vol[ip];

  }

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  uza->info = cs_iter_algo_define(nslesp->il_algo_verbosity,
                                  nslesp->n_max_il_algo_iter,
                                  nslesp->il_algo_atol,
                                  nslesp->il_algo_rtol,
                                  nslesp->il_algo_dtol);

  return uza;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a Uzawa builder structure
 *
 * \param[in, out]  p_uza   double pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_uza_builder(cs_uza_builder_t   **p_uza)
{
  cs_uza_builder_t  *uza = *p_uza;

  if (uza == NULL)
    return;

  BFT_FREE(uza->b_tilda);

  BFT_FREE(uza->inv_mp);
  BFT_FREE(uza->res_p);
  BFT_FREE(uza->d__v);
  BFT_FREE(uza->rhs);
  BFT_FREE(uza->gk);
  BFT_FREE(uza->dzk);

  BFT_FREE(uza->info);

  BFT_FREE(uza);
  *p_uza = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the divergence operator and store the result in div_v
 *
 * \param[in]      div_op  pointer to the values of divergence operator
 * \param[in]      v       vector to apply in velocity space
 * \param[in, out] div_v   resulting vector in pressure space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op(const cs_real_t   *div_op,
              const cs_real_t   *v,
              cs_real_t         *div_v)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t _div_v = 0;
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
      _div_v += cs_math_3_dot_product(div_op + 3*j, v + 3*c2f->ids[j]);
    div_v[c_id] = _div_v;

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the gradient operator (which is the transpose of the
 *         divergence operator) and store the result in dt_q
 *
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in]      q        vector to apply in pressure space
 * \param[in, out] dt_q     resulting vector in velocity space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op_transpose(const cs_real_t   *div_op,
                        const cs_real_t   *q,
                        cs_real_t         *dt_q)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

  memset(dt_q, 0, 3*quant->n_faces*sizeof(cs_real_t));

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  qc = q[c_id];
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_real_t  *_div_f = div_op + 3*j;

      cs_real_t  *_dt_q = dt_q + 3*c2f->ids[j];
#     pragma omp critical
      {
        _dt_q[0] += qc * _div_f[0];
        _dt_q[1] += qc * _div_f[1];
        _dt_q[2] += qc * _div_f[2];
      }

    } /* Loop on cell faces */

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Transform the initial saddle-point problem. The velocity unknown
 *         is modified and is stored in u_tilda as well as the RHS related to
 *         the mass equation and stored in b_tilda
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      nslesp   pointer to SLES settings for NavSto
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] u_f      initial velocity on faces
 * \param[in, out] b_f      right-hand side (scatter/gather if needed) on faces
 * \param[in, out] b_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_transform_gkb_system(const cs_matrix_t              *matrix,
                      const cs_equation_param_t      *eqp,
                      const cs_navsto_param_sles_t   *nslesp,
                      const cs_real_t                *div_op,
                      cs_gkb_builder_t               *gkb,
                      cs_sles_t                      *sles,
                      cs_real_t                      *u_f,
                      cs_real_t                      *b_f,
                      cs_real_t                      *b_c)
{
  assert(gkb != NULL);

  cs_real_t  normalization = 1.0; /* TODO */

  /* Modifiy the tolerance in order to be more accurate on this step */
  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":gkb_transfo") + 1, char);
  sprintf(system_name, "%s:gkb_transfo", eqp->name);

  cs_param_sles_t  *slesp = cs_param_sles_create(-1, system_name);

  cs_param_sles_copy_from(eqp->sles_param, slesp);

  slesp->eps = nslesp->il_algo_rtol;

  if (gkb->gamma > 0) {

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->b_tilda[ip] = gkb->gamma*b_c[ip]/cs_shared_quant->cell_vol[ip];

    /* Solve Dt.b_tilda */
    _apply_div_op_transpose(div_op, gkb->b_tilda, gkb->dt_q);

#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
      gkb->b_tilda[iu] = b_f[iu] + gkb->dt_q[iu];

  }
  else
    memcpy(gkb->b_tilda, b_f, gkb->n_u_dofs*sizeof(cs_real_t));

  /* Compute M^-1.(b_f + gamma. Bt.N^-1.b_c) */
  gkb->info->n_inner_iter
    += (gkb->info->last_inner_iter
        = cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                          slesp,
                                          matrix,
                                          cs_shared_range_set,
                                          normalization,
                                          true, /* rhs_redux, */
                                          sles,
                                          gkb->v,
                                          gkb->b_tilda));

  /* Compute the initial u_tilda := u_f - M^-1.(b_f + gamma. Bt.N^-1.b_c) */
# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    gkb->u_tilda[iu] = u_f[iu] - gkb->v[iu];

  /* Compute b_tilda := b_c - div(M^-1.b_f) */
  _apply_div_op(div_op, gkb->v, gkb->d__v);

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
    gkb->b_tilda[ip] = b_c[ip] - gkb->d__v[ip];

  /* Free memory */
  BFT_FREE(system_name);
  cs_param_sles_free(&slesp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the GKB algorithm
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] p_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_init_gkb_algo(const cs_matrix_t             *matrix,
               const cs_equation_param_t     *eqp,
               const cs_real_t               *div_op,
               cs_gkb_builder_t              *gkb,
               cs_sles_t                     *sles,
               cs_real_t                     *p_c)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  size = quant->n_cells;

  double beta2 = 0.0;

  /* Compute beta := ||b_tilta||_N^-1 and q := N^-1(b_tilda)/beta */
# pragma omp parallel reduction(+:beta2) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_real_t  *_w = quant->cell_vol + s_id;
    const cs_real_t  *_b = gkb->b_tilda + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_real_t  *_q = gkb->q + s_id;
    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_beta2 = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _beta2 = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const  cs_real_t  b_ov_w = _b[j]/_w[j];
          _beta2 += b_ov_w*_b[j];
          _q[j] = b_ov_w;

        } /* Loop on block_size */

        s_beta2 += _beta2;

      } /* Loop on blocks */

      beta2 += s_beta2;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel synchronization */
  cs_parall_sum(1, CS_DOUBLE, &beta2);

  /* Keep the value of beta = ||b||_{N^-1} */
  assert(beta2 > -DBL_MIN);
  gkb->beta = sqrt(beta2);

  /* Store M^-1.(b_f + gamma. Bt.N^-1.b_c) in b_tilda which is not useful
   * anymore */
  memcpy(gkb->b_tilda, gkb->v, gkb->n_u_dofs*sizeof(cs_real_t));

  if (fabs(gkb->beta) > FLT_MIN) {
    const  cs_real_t  inv_beta = 1./gkb->beta;
# pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      gkb->q[i] *= inv_beta;
  }
  else {
    gkb->info->cvg = CS_SLES_CONVERGED;
    return;
  }

  /* Solve M.w = Dt.q */
  _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

  if (cs_shared_range_set->ifs != NULL)
    cs_interface_set_sum(cs_shared_range_set->ifs,
                         gkb->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         gkb->dt_q);

  cs_real_t  normalization = 1.0; /* TODO */

  gkb->info->n_inner_iter
    += (gkb->info->last_inner_iter =
        cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                        eqp->sles_param,
                                        matrix,
                                        cs_shared_range_set,
                                        normalization,
                                        false, /* rhs_redux */
                                        sles,
                                        gkb->v,
                                        gkb->dt_q));

  gkb->alpha = _face_gdot(gkb->n_u_dofs, gkb->v, gkb->dt_q);
  assert(gkb->alpha > -DBL_MIN);
  gkb->alpha = sqrt(gkb->alpha);

  const double ov_alpha = 1./gkb->alpha;

  gkb->zeta = gkb->beta * ov_alpha;

  /* Initialize auxiliary vectors and first update of the solution vectors */

# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
    gkb->v[iu] *= ov_alpha;
    gkb->u_tilda[iu] = gkb->zeta * gkb->v[iu];
    gkb->m__v[iu] = ov_alpha * gkb->dt_q[iu];
  }

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
    gkb->d[ip] = gkb->q[ip] * ov_alpha;
    p_c[ip] = -gkb->zeta * gkb->d[ip];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more GKB iteration
 *
 * \param[in, out] gkb     pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_cvg_test(cs_gkb_builder_t           *gkb)
{
  /* Update the sum of square of zeta values (used for renormalization) */
  cs_real_t  z2 = gkb->zeta*gkb->zeta;

  gkb->zeta_square_sum += z2;
  gkb->zeta_array[gkb->info->n_algo_iter % gkb->z_size] = z2;

  /* Increment the number of Picard iterations */
  gkb->info->n_algo_iter += 1;

  /* Compute the relative energy norm. The normalization arises from an
     iterative estimation of the initial error in the energy norm */
  const cs_real_t  prev_res = gkb->info->res;

  int  n = gkb->z_size;
  if (gkb->info->n_algo_iter < gkb->z_size)
    n = gkb->info->n_algo_iter;

  cs_real_t  err2_energy = 0.;
  for (int i = 0; i < n; i++)
    err2_energy += gkb->zeta_array[i];

  double  tau = (gkb->gamma > 0) ?
    gkb->gamma*gkb->info->rtol : gkb->info->rtol;

  gkb->info->res = sqrt(err2_energy);

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nGKB.It%02d-- err2 = %6.4e ?<? tau * square_sum %6.4e\n",
                gkb->info->n_algo_iter, err2_energy, tau*gkb->zeta_square_sum);
#endif

  if (err2_energy < tau * gkb->zeta_square_sum)
    gkb->info->cvg = CS_SLES_CONVERGED;

  else if (gkb->info->n_algo_iter >= gkb->info->n_max_algo_iter)
    gkb->info->cvg = CS_SLES_MAX_ITERATION;

  else if (gkb->info->res > gkb->info->dtol * prev_res)
    gkb->info->cvg = CS_SLES_DIVERGED;

  else
    gkb->info->cvg = CS_SLES_ITERATING;

  if (gkb->info->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "### GKB.It%d %5.3e %5d %6d z2:%6.4e renorm:%6.4e cvg:%d\n",
                  gkb->info->n_algo_iter, gkb->info->res,
                  gkb->info->last_inner_iter, gkb->info->n_inner_iter,
                  z2, sqrt(gkb->zeta_square_sum), gkb->info->cvg);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration in case of an Uzawa
 *         CG (conjugate gradient variant). The residual criterion has to be
 *         computed before calling this function.
 *
 * \param[in, out] uza     pointer to a Uzawa builder structure
 *
 * \return true (one moe iteration) otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_uza_cg_cvg_test(cs_uza_builder_t           *uza)
{
  /* Increment the number of algo. iterations */
  uza->info->n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */
  const cs_real_t  prev_res = uza->info->res;
  const double  tau = fmax(uza->info->rtol*uza->info->res0, uza->info->atol);

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nUZA-CG.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->info->n_algo_iter, uza->info->res, uza->info->rtol);
#endif


  if (uza->info->res < tau)
    uza->info->cvg = CS_SLES_CONVERGED;

  else if (uza->info->n_algo_iter >= uza->info->n_max_algo_iter)
    uza->info->cvg = CS_SLES_MAX_ITERATION;

  else if (uza->info->res > uza->info->dtol * prev_res)
    uza->info->cvg = CS_SLES_DIVERGED;

  else
    uza->info->cvg = CS_SLES_ITERATING;

  if (uza->info->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "<UZACG.It%02d> res %5.3e | %4d %6d cvg%d | fit.eps %5.3e\n",
                  uza->info->n_algo_iter, uza->info->res,
                  uza->info->last_inner_iter, uza->info->n_inner_iter,
                  uza->info->cvg, tau);

  if (uza->info->cvg == CS_SLES_ITERATING)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration
 *
 * \param[in, out] uza     pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_uza_cvg_test(cs_uza_builder_t           *uza)
{
  /* Increment the number of algo. iterations */
  uza->info->n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */
  const cs_real_t  prev_res = uza->info->res;

  cs_real_t  res_square = cs_dot_wxx(uza->n_p_dofs, uza->inv_mp, uza->res_p);
  cs_parall_sum(1, CS_DOUBLE, &res_square);
  assert(res_square > -DBL_MIN);
  uza->info->res = sqrt(res_square);

  double  tau = (uza->gamma > 0) ?
    uza->info->rtol/sqrt(uza->gamma) : uza->info->rtol;

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nUZA.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->info->n_algo_iter, uza->info->res, uza->info->rtol);
#endif

  if (uza->info->res < tau)
    uza->info->cvg = CS_SLES_CONVERGED;

  else if (uza->info->n_algo_iter >= uza->info->n_max_algo_iter)
    uza->info->cvg = CS_SLES_MAX_ITERATION;

  else if (uza->info->res > uza->info->dtol * prev_res)
    uza->info->cvg = CS_SLES_DIVERGED;

  else
    uza->info->cvg = CS_SLES_ITERATING;

  if (uza->info->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "### UZA.It%02d-- %5.3e %5d %6d cvg:%d\n",
                  uza->info->n_algo_iter, uza->info->res,
                  uza->info->last_inner_iter, uza->info->n_inner_iter,
                  uza->info->cvg);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration in case of an incremental
 *         formulation
 *
 * \param[in]      delta_u_l2   value of the weighted L2 norm of delta_u
 * \param[in, out] uza          pointer to a Uzawa builder structure
 *
 * \return true if one more iteration is needed otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_uza_incr_cvg_test(cs_real_t                   delta_u_l2,
                   cs_uza_builder_t           *uza)
{
  /* Increment the number of algo. iterations */
  uza->info->n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */
  const cs_real_t  prev_res = uza->info->res;

  cs_real_t  res_square = cs_dot_wxx(uza->n_p_dofs, uza->inv_mp, uza->d__v);
  cs_parall_sum(1, CS_DOUBLE, &res_square);
  assert(res_square > -DBL_MIN);
  cs_real_t  divu_l2 = sqrt(res_square);
  uza->info->res = fmax(delta_u_l2, divu_l2);

  if (uza->info->n_algo_iter == 1) { /* First call */
    uza->info->res0 = uza->info->res;
    uza->info->tol = fmax(uza->info->atol, uza->info->rtol*uza->info->res0);

    if (uza->info->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### UZAi.res0: %5.3e tol: %5.3e\n",
                    uza->info->res0, uza->info->tol);
  }

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nUZAi.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->info->n_algo_iter, uza->info->res, uza->info->rtol);
#endif

  if (uza->info->res < uza->info->tol)
    uza->info->cvg = CS_SLES_CONVERGED;

  else if (uza->info->n_algo_iter >= uza->info->n_max_algo_iter)
    uza->info->cvg = CS_SLES_MAX_ITERATION;

  else if (uza->info->res > uza->info->dtol * prev_res)
    uza->info->cvg = CS_SLES_DIVERGED;

  else
    uza->info->cvg = CS_SLES_ITERATING;

  if (uza->info->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "### UZAi.It%02d %5.3e %5d %6d cvg:%d div:%5.3e, du:%5.3e\n",
                  uza->info->n_algo_iter, uza->info->res,
                  uza->info->last_inner_iter, uza->info->n_inner_iter,
                  uza->info->cvg, divu_l2, delta_u_l2);

  if (uza->info->cvg == CS_SLES_ITERATING)
    return true;
  else
    return false;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdofb_monolithic_sles_t structure
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdofb_monolithic_sles_t *
cs_cdofb_monolithic_sles_create(void)
{
  cs_cdofb_monolithic_sles_t  *msles = NULL;

  BFT_MALLOC(msles, 1, cs_cdofb_monolithic_sles_t);

  msles->block_matrices = NULL;
  msles->compatible_laplacian = NULL;
  msles->div_op = NULL;

  msles->graddiv_coef = 0.;

  msles->sles = NULL;
  msles->schur_sles = NULL;

  msles->n_faces = 0;
  msles->n_cells = 0;

  msles->u_f = NULL;
  msles->p_c = NULL;
  msles->b_f = NULL;
  msles->b_c = NULL;

  return msles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the rhs
 *
 * \param[in]       n_cells    local number of cells
 * \param[in]       n_faces    local number of faces
 * \param[in, out]  msles      pointer to the structure to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init(cs_lnum_t                     n_cells,
                              cs_lnum_t                     n_faces,
                              cs_cdofb_monolithic_sles_t   *msles)
{
  if (msles == NULL)
    return;

  msles->n_cells = n_cells;
  msles->n_faces = n_faces;

  cs_lnum_t  full_size = 3*n_faces + n_cells;
  BFT_MALLOC(msles->b_f, full_size, cs_real_t);
  msles->b_c = msles->b_f + 3*n_faces;

  /* Set rhs to zero */
#if defined(HAVE_OPENMP)
# pragma omp parallel if (full_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < full_size; i++)
    msles->b_f[i] = 0.;
#else
  memset(msles->b_f, 0, full_size*sizeof(cs_real_t));
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset to zero rhs and clean the cs_sles_t structure
 *
 * \param[in, out]  msles   pointer to the structure to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_reset(cs_cdofb_monolithic_sles_t   *msles)
{
  if (msles == NULL)
    return;

  for (int i = 0; i < msles->n_row_blocks*msles->n_row_blocks; i++)
    cs_matrix_destroy(&(msles->block_matrices[i]));

  cs_matrix_destroy(&msles->compatible_laplacian);

  cs_sles_free(msles->sles);
  cs_sles_free(msles->schur_sles);

  cs_lnum_t  full_size = 3*msles->n_faces + msles->n_cells;

  /* Set rhs to zero (b_f and b_c are stored consecutively) */
#if defined(HAVE_OPENMP)
# pragma omp parallel if (full_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < full_size; i++)
    msles->b_f[i] = 0.;
#else
  memset(msles->b_f, 0, full_size*sizeof(cs_real_t));
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a part of the structure
 *
 * \param[in, out]  msles   pointer to the structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_clean(cs_cdofb_monolithic_sles_t   *msles)
{
  if (msles == NULL)
    return;

  for (int i = 0; i < msles->n_row_blocks*msles->n_row_blocks; i++)
    cs_matrix_destroy(&(msles->block_matrices[i]));

  cs_sles_free(msles->sles);
  cs_sles_free(msles->schur_sles);

  /* b_f and b_c are stored consecutively */
  BFT_FREE(msles->b_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdofb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_free(cs_cdofb_monolithic_sles_t   **p_msles)
{
  cs_cdofb_monolithic_sles_t  *msles = *p_msles;

  if (msles == NULL)
    return;

  BFT_FREE(msles->block_matrices);
  BFT_FREE(msles->div_op);
  /* other pointer are shared, thus no free at this stage */

  BFT_FREE(msles);
  *p_msles = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  connect  pointer to cdo connectivities
 * \param[in]  quant    pointer to additional mesh quantities
 * \param[in]  rset     pointer to a \ref cs_range_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_set_shared(const cs_cdo_connect_t        *connect,
                                    const cs_cdo_quantities_t     *quant,
                                    const cs_range_set_t          *rset)
{
  assert(rset != NULL);

  /* Assign static const pointers */
  cs_shared_connect = connect;
  cs_shared_quant = quant;
  cs_shared_range_set = rset;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage.
 *         nsp is not declared as cnost to avoid compilation warnings but
 *         it should be modified at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(cs_navsto_param_t    *nsp,
                             void                 *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  assert(nsp != NULL && nsc != NULL);

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  cs_param_sles_t  *mom_slesp = mom_eqp->sles_param;
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  mom_slesp->field_id = field_id;
  if (mom_slesp->amg_type == CS_PARAM_AMG_NONE)
    mom_slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER;

  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  switch (nslesp->strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
    /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
     * Notice that zeta can be equal to 0 */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_MINRES:
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
    {
      cs_equation_param_set_sles(mom_eqp);

      /* Set the solver for the compatible Laplacian (the related SLES is
         defined using the system name instead of the field id since this is an
         auxiliary system) */
      int ier = cs_param_sles_set(false, nslesp->schur_sles_param);

      if (ier == -1)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: The requested class of solvers is not available"
                  " for the system %s\n Please modify your settings.",
                  __func__, nslesp->schur_sles_param->name);
    }
    break;

  case CS_NAVSTO_SLES_UZAWA_CG:
    {
      /* Set solver and preconditioner for solving A */
      cs_equation_param_set_sles(mom_eqp);

      /* Set the solver for the compatible Laplacian (the related SLES is
         defined using the system name instead of the field id since this is an
         auxiliary system) */
      int ier = cs_param_sles_set(false, nslesp->schur_sles_param);

      if (ier == -1)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: The requested class of solvers is not available"
                  " for the system %s\n Please modify your settings.",
                  __func__, nslesp->schur_sles_param->name);
    }
    break;

  case CS_NAVSTO_SLES_UZAWA_AL:
     /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
      * Notice that zeta can be equal to 0 */
    cs_equation_param_set_sles(mom_eqp);
    break;

#if defined(HAVE_PETSC)
#if PETSC_VERSION_GE(3,11,0)    /* Golub-Kahan Bi-diagonalization */
  case CS_NAVSTO_SLES_GKB_PETSC:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_hook,
                         (void *)nsp);
    break;

  case CS_NAVSTO_SLES_GKB_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_precond_hook,
                         (void *)nsp);
    break;
#else  /* PETSC version < 3.11 */
  case CS_NAVSTO_SLES_GKB_PETSC:
  case CS_NAVSTO_SLES_GKB_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc 3.11.x or greater is required with this option.\n",
              __func__, mom_eqp->name);
    break;
#endif
#else  /* no HAVE_PETSC */
  case CS_NAVSTO_SLES_GKB_PETSC:
  case CS_NAVSTO_SLES_GKB_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please use a version of Code_Saturne built with PETSc.",
              __func__, mom_eqp->name);
    break;

#endif

  case CS_NAVSTO_SLES_MUMPS:
#if defined(HAVE_MUMPS)
    if (mom_slesp->solver != CS_PARAM_ITSOL_MUMPS &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_LDLT &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_FLOAT &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT)
      mom_slesp->solver = CS_PARAM_ITSOL_MUMPS;

    cs_sles_mumps_define(field_id,
                         NULL,
                         mom_slesp,
                         cs_user_sles_mumps_hook,
                         NULL);
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _mumps_hook,
                         (void *)mom_eqp);
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc with MUMPS is required with this option.\n",
              __func__, mom_eqp->name);
#endif  /* PETSC_HAVE_MUMPS */
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " Neither PETSc nor MUMPS is available.\n",
              __func__, mom_eqp->name);
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
    break;

#if defined(HAVE_PETSC)
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _additive_amg_hook,
                         (void *)nsp);
    break;

  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _multiplicative_hook,
                         (void *)nsp);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _diag_schur_hook,
                         (void *)nsp);
    break;

  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _upper_schur_hook,
                         (void *)nsp);
    break;

#else
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please use a version of Code_Saturne built with PETSc.",
              __func__, mom_eqp->name);
    break;
#endif /* HAVE_PETSC */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

  /* Define the level of verbosity for SLES structure */
  if (mom_slesp->verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, mom_slesp->verbosity);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         Navier-Stokes equation with a CDO face-based approach.
 *         The full system is treated as one block and then sent to PETSc
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_solve(const cs_navsto_param_t       *nsp,
                          const cs_equation_param_t     *eqp,
                          cs_cdofb_monolithic_sles_t    *msles)
{
  const cs_matrix_t  *matrix = msles->block_matrices[0];
  const cs_lnum_t  n_faces = msles->n_faces;
  const cs_lnum_t  n_cells = msles->n_cells;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_scatter_elts = 3*n_faces + n_cells;

  /* De-interlace the velocity array and the rhs for the face DoFs */
  cs_real_t  *xsol = NULL;
  BFT_MALLOC(xsol, n_cols, cs_real_t);

  cs_real_t  *b = NULL;
  BFT_MALLOC(b, n_scatter_elts, cs_real_t);

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, xsol, b)                                \
  firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    xsol[f            ] = msles->u_f[3*f];
    xsol[f +   n_faces] = msles->u_f[3*f+1];
    xsol[f + 2*n_faces] = msles->u_f[3*f+2];

    b[f            ] = msles->b_f[3*f];
    b[f +   n_faces] = msles->b_f[3*f+1];
    b[f + 2*n_faces] = msles->b_f[3*f+2];

  }

  /* Add the pressure related elements */
  memcpy(xsol + 3*n_faces, msles->p_c, n_cells*sizeof(cs_real_t));
  memcpy(b + 3*n_faces, msles->b_c, n_cells*sizeof(cs_real_t));

  const cs_range_set_t  *rset = cs_shared_range_set;
  int  n_iters = 0;
  double  residual = DBL_MAX;

  /* Prepare solving (handle parallelism) */
  cs_equation_prepare_system(1,     /* stride */
                             n_scatter_elts,
                             matrix,
                             rset,
                             true,  /* rhs_redux */
                             xsol, b);

  /* Solve the linear solver */
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  const cs_param_sles_t  *sles_param = eqp->sles_param;
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */

  cs_real_t  rtol = sles_param->eps;

  if (nslesp->strategy == CS_NAVSTO_SLES_UPPER_SCHUR_GMRES              ||
      nslesp->strategy == CS_NAVSTO_SLES_DIAG_SCHUR_GMRES               ||
      nslesp->strategy == CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK  ||
      nslesp->strategy == CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK)
    rtol = nslesp->il_algo_rtol;

  cs_sles_convergence_state_t  code = cs_sles_solve(msles->sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    rtol,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %d |"
                  " residual % -8.4e\n",
                  eqp->name, code, n_iters, residual);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       xsol, xsol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 1
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_MONOLITHIC_DBG,
                        xsol, b, 3*n_faces);
#endif

  /* Interlace xsol --> u_f and p_c */
# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, xsol) firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    msles->u_f[3*f]   = xsol[f];
    msles->u_f[3*f+1] = xsol[f +   n_faces];
    msles->u_f[3*f+2] = xsol[f + 2*n_faces];

  }

  /* Copy the part of the solution array related to the pressure in cells */
  memcpy(msles->p_c, xsol + 3*n_faces, n_cells*sizeof(cs_real_t));

  /* Free what can be freed at this stage */
  BFT_FREE(xsol);
  BFT_FREE(b);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach. The system is
 *        split into blocks to enable more efficient preconditioning
 *        techniques. The main iterative solver is a Krylov solver such as GCR,
 *        GMRES or MINRES
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_krylov_block_precond(const cs_navsto_param_t       *nsp,
                                         const cs_equation_param_t     *eqp,
                                         cs_cdofb_monolithic_sles_t    *msles)
{
  /*  Sanity checks */
  if (msles == NULL)
    return 0;
  if (msles->n_row_blocks != 1) /* Only this case is handled up to now */
    bft_error(__FILE__, __LINE__, 0,
              "%s: Only one block for the velocity is possible up to now.",
              __func__);

  /* Set the structure to manage the iterative solver */
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_iter_algo_info_t
    *saddle_info = cs_iter_algo_define(nslesp->il_algo_verbosity,
                                       nslesp->n_max_il_algo_iter,
                                       nslesp->il_algo_atol,
                                       nslesp->il_algo_rtol,
                                       nslesp->il_algo_dtol);

  /* Set the saddle-point system */
  /* --------------------------- */

  cs_saddle_system_t  *ssys = NULL;
  BFT_MALLOC(ssys, 1, cs_saddle_system_t);

  ssys->n_m11_matrices = 1;
  ssys->m11_matrices = msles->block_matrices;
  ssys->x1_size = 3*msles->n_faces;
  ssys->max_x1_size =
    CS_MAX(ssys->x1_size, cs_matrix_get_n_columns(msles->block_matrices[0]));
  ssys->rhs1 = msles->b_f;

  ssys->x2_size = msles->n_cells;
  ssys->rhs2 = msles->b_c;

  ssys->m21_stride = 3;
  ssys->m21_unassembled = msles->div_op;
  ssys->m21_adjacency = cs_shared_connect->c2f;

  ssys->rset = cs_shared_range_set;

  /* u_f is allocated to 3*n_faces (the size of the scatter view but during the
     resolution process one need a vector at least of size n_cols of the matrix
     m11. */
  cs_real_t  *xu = NULL;
  BFT_MALLOC(xu, ssys->max_x1_size, cs_real_t);
  memcpy(xu, msles->u_f, ssys->x1_size*sizeof(cs_real_t));

  switch (nsp->sles_param->strategy) {

  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
    {
      /* Define block preconditionning */
      cs_saddle_block_precond_t  *sbp =
        cs_saddle_block_precond_create(CS_PARAM_PRECOND_BLOCK_DIAG,
                                       nslesp->schur_approximation,
                                       eqp->sles_param,
                                       msles->sles);

      /* Schur preconditionning */
      cs_param_sles_t  *schur_slesp = nslesp->schur_sles_param;

      sbp->schur_slesp = schur_slesp;
      if (msles->schur_sles == NULL)
        /* This sles structure should have been defined by name */
        sbp->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

      /* Compute the schur approximation matrix */
      switch (nslesp->schur_approximation) {

      case CS_PARAM_SCHUR_DIAG_INVERSE:
        _diag_schur_sbp(nsp, ssys, sbp);
        break;
      case CS_PARAM_SCHUR_ELMAN:
        _elman_schur_sbp(nsp, ssys, sbp);
        break;
      case CS_PARAM_SCHUR_IDENTITY:
        break; /* Nothing to do */
      case CS_PARAM_SCHUR_LUMPED_INVERSE:
        _invlumped_schur_sbp(nsp, ssys, sbp);
        break;
      case CS_PARAM_SCHUR_MASS_SCALED:
        _scaled_mass_sbp(nsp, ssys, sbp);
        break; /* Nothing to do */
      case CS_PARAM_SCHUR_MASS_SCALED_DIAG_INVERSE:
        _scaled_mass_sbp(nsp, ssys, sbp);
        _diag_schur_sbp(nsp, ssys, sbp);
        break;
      case CS_PARAM_SCHUR_MASS_SCALED_LUMPED_INVERSE:
        _scaled_mass_sbp(nsp, ssys, sbp);
        _invlumped_schur_sbp(nsp, ssys, sbp);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid Schur approximation.", __func__);
      }

      /* Call the inner linear algorithm */
      cs_saddle_minres(ssys, sbp, xu, msles->p_c, saddle_info);

      cs_saddle_block_precond_free(&sbp);
    }
    break;

  case CS_NAVSTO_SLES_MINRES:   /* No block preconditioning */
    cs_saddle_minres(ssys, NULL, xu, msles->p_c, saddle_info);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy to solve the system.\n"
              "Please used a Krylov-based iterative solver.", __func__);
    break;
  }

  memcpy(msles->u_f, xu, ssys->x1_size*sizeof(cs_real_t));

  /* Free the saddle-point system */
  BFT_FREE(xu);
  BFT_FREE(ssys); /* only shared pointers inside */

  int n_algo_iter = saddle_info->n_algo_iter;

  if (nslesp->il_algo_verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  " -cvg- inner_algo: cumulated_iters: %d\n",
                  saddle_info->n_inner_iter);


  BFT_FREE(saddle_info);

  return n_algo_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_gkb_solve(const cs_navsto_param_t       *nsp,
                              const cs_equation_param_t     *eqp,
                              cs_cdofb_monolithic_sles_t    *msles)
{
  assert(nsp != NULL);
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  /* Sanity checks */
  assert(nslesp->strategy == CS_NAVSTO_SLES_GKB_SATURNE);
  assert(cs_shared_range_set != NULL);

  const cs_real_t  *vol = cs_shared_quant->cell_vol;
  const cs_real_t  *div_op = msles->div_op;
  const cs_matrix_t  *matrix = msles->block_matrices[0];
  const cs_real_t  gamma = msles->graddiv_coef;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the GKB builder structure */
  cs_gkb_builder_t  *gkb = _init_gkb_builder(nsp,
                                             gamma,
                                             3*msles->n_faces,
                                             msles->n_cells);

  /* Transformation of the initial saddle-point system */
  _transform_gkb_system(matrix, eqp, nslesp, div_op, gkb, msles->sles,
                        u_f, b_f, b_c);

  /* Initialization */
  _init_gkb_algo(matrix, eqp, div_op, gkb, msles->sles, p_c);

  /* Main loop */
  /* ========= */

  while (gkb->info->cvg == CS_SLES_ITERATING) {

    /* Compute g (store as an update of d__v), q */
    _apply_div_op(div_op, gkb->v, gkb->d__v);

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d__v[ip] /= vol[ip];
      gkb->d__v[ip] -= gkb->alpha * gkb->q[ip];
    }

    /* Compute beta */
    gkb->beta = cs_dot_wxx(gkb->n_p_dofs, vol, gkb->d__v);
    cs_parall_sum(1, CS_DOUBLE, &(gkb->beta));
    assert(gkb->beta > -DBL_MIN);
    gkb->beta = sqrt(gkb->beta);

    const double  ov_beta = 1./gkb->beta;

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->q[ip] = ov_beta*gkb->d__v[ip];

    /* Solve M.w_tilda = Dt.q */
    _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

    if (cs_shared_range_set->ifs != NULL)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           gkb->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           gkb->dt_q);

    /* Prepare update of m__v:
     *  m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */
#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->m__v[iu] *= -gkb->beta;
      gkb->m__v[iu] +=  gkb->dt_q[iu];
    }

    cs_real_t  normalization = gkb->alpha; /* TODO */
    gkb->info->n_inner_iter
      += (gkb->info->last_inner_iter =
          cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                          eqp->sles_param,
                                          matrix,
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux */
                                          msles->sles,
                                          gkb->v,
                                          gkb->m__v));

    /* Compute alpha */
    gkb->alpha = _face_gdot(gkb->n_u_dofs, gkb->v, gkb->m__v);
    assert(gkb->alpha > -DBL_MIN);
    gkb->alpha = sqrt(gkb->alpha);

    const double ov_alpha = 1./gkb->alpha;

    /* zeta(k+1) = -beta/alpha * zeta(k) */
    gkb->zeta *= -gkb->beta * ov_alpha;

    /* Update vectors and solutions */
#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->v[iu] *= ov_alpha;
      gkb->u_tilda[iu] += gkb->zeta * gkb->v[iu];
      /* Last step: m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */
      gkb->m__v[iu] *= ov_alpha;
    }

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d[ip] = ov_alpha * (gkb->q[ip] - gkb->beta*gkb->d[ip]);
      p_c[ip] += -gkb->zeta * gkb->d[ip];
    }

    /* Update error norm and test if one needs one more iteration */
    _gkb_cvg_test(gkb);

  }

  /* Return to the initial velocity formulation
   * u: = u_tilda + M^-1.(b_f + gamma.N^-1.b_c)
   * where M^-1.(b_f + gamma.N^-1.b_c) is stored in b_tilda */
# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    u_f[iu] = gkb->u_tilda[iu] + gkb->b_tilda[iu];

  int n_inner_iter = gkb->info->n_inner_iter;

  /* Last step: Free temporary memory */
  _free_gkb_builder(&gkb);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the preconditioned Uzawa-CG algorithm to solve the saddle-point
 *         problem arising from CDO-Fb schemes for Stokes, Oseen and
 *         Navier-Stokes with a monolithic coupling
 *         This algorithm is based on Koko's paper "Uzawa conjugate gradient
 *         method for the Stokes problem: Matlab implementation with P1-iso-P2/
 *         P1 finite element"
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_cg_solve(const cs_navsto_param_t       *nsp,
                                   const cs_equation_param_t     *eqp,
                                   cs_cdofb_monolithic_sles_t    *msles)
{
  int  _n_iter;

  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_CG);
  assert(cs_shared_range_set != NULL);

  const cs_real_t  *B_op = msles->div_op;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the Uzawa-CG builder structure */
  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               0, /* grad-div scaling */
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               cs_shared_quant);

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  /* The Schur complement approximation (B.A^-1.Bt) is build and stored in the
     native format */
  cs_matrix_t  *A = msles->block_matrices[0];

  cs_matrix_t  *K = NULL;
  cs_real_t  *diagK = NULL, *xtraK = NULL;

  switch (nsp->sles_param->schur_approximation) {

  case CS_PARAM_SCHUR_DIAG_INVERSE:
    K = _diag_schur_approximation(nsp, A, uza, &diagK, &xtraK);
    break;
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
    K = _invlumped_schur_approximation(nsp, eqp, msles, A,
                                       uza, &diagK, &xtraK);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid Schur approximation.",
              __func__);
  }

  cs_param_sles_t  *schur_slesp = nslesp->schur_sles_param;

  if (msles->schur_sles == NULL) /* has been defined by name */
    msles->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

  /* Compute the first RHS: A.u0 = rhs = b_f - B^t.p_0 to solve */
  _apply_div_op_transpose(B_op, p_c, uza->rhs);

  if (cs_shared_range_set->ifs != NULL) {

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->rhs);

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         b_f);

  }

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
    uza->rhs[iu] = b_f[iu] - uza->rhs[iu];

  /* Initial residual for rhs */
  double normalization = _get_fbvect_norm(uza->rhs);

  /* Compute the first velocity guess */

  /* Modify the tolerance in order to be more accurate on this step */
  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":init_guess") + 1, char);
  sprintf(system_name, "%s:init_guess", eqp->name);

  cs_param_sles_t  *slesp0 = cs_param_sles_create(-1, system_name);

  cs_param_sles_copy_from(eqp->sles_param, slesp0);
  slesp0->eps = fmin(1e-6, 0.05*nslesp->il_algo_rtol);

  _n_iter =
    cs_equation_solve_scalar_system(uza->n_u_dofs,
                                    slesp0,
                                    A,
                                    cs_shared_range_set,
                                    normalization,
                                    false, /* rhs_redux --> already done */
                                    msles->sles,
                                    u_f,
                                    uza->rhs);
  uza->info->n_inner_iter += _n_iter;
  uza->info->last_inner_iter += _n_iter;

  /* Partial memory free */
  BFT_FREE(system_name);
  cs_param_sles_free(&slesp0);

  /* Set pointers used in this algorithm */
  cs_real_t  *gk = uza->gk;
  cs_real_t  *dk = uza->res_p;
  cs_real_t  *rk = uza->d__v;    /* P space */
  cs_real_t  *wk = uza->b_tilda; /* U space */
  cs_real_t  *dwk = uza->dzk;    /* P space */
  cs_real_t  *zk = uza->rhs;     /* P or U space */

  /* Compute the first residual rk0 (in fact the velocity divergence) */
  _apply_div_op(B_op, u_f, rk);

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    rk[ip] = b_c[ip] - rk[ip];

  double div_l2_norm = _get_cbscal_norm(rk);

  /* Compute g0 as
   *   Solve K.zk = r0
   *   g0 = alpha zk + nu Mp^-1 r0 */
  memset(zk, 0, sizeof(cs_real_t)*uza->n_p_dofs);

  _n_iter = cs_equation_solve_scalar_cell_system(uza->n_p_dofs,
                                                 schur_slesp,
                                                 K,
                                                 div_l2_norm,
                                                 msles->schur_sles,
                                                 zk,
                                                 rk);
  uza->info->n_inner_iter += _n_iter;
  uza->info->last_inner_iter += _n_iter;

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    gk[ip] = uza->alpha*zk[ip] + uza->inv_mp[ip]*rk[ip];

  /* dk0 <-- gk0 */
  memcpy(dk, gk, uza->n_p_dofs*sizeof(cs_real_t));

  uza->info->res0 = cs_gdot(uza->n_p_dofs, rk, gk);
  uza->info->res = uza->info->res0;

  /* Main loop knowing g0, r0, d0, u0, p0 */
  /* ------------------------------------ */

  while (_uza_cg_cvg_test(uza)) {

    /* Sensitivity step: Compute wk as the solution of A.wk = B^t.dk */

    /* Define the rhs for this system */
    _apply_div_op_transpose(B_op, dk, uza->rhs);
    if (cs_shared_range_set->ifs != NULL)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

    normalization = _get_fbvect_norm(uza->rhs);

    /* Solve A.wk = B^t.dk (should be -B^t this implies a sign modification
       during the update step) */
    memset(wk, 0, sizeof(cs_real_t)*uza->n_u_dofs);

    uza->info->n_inner_iter
      += (uza->info->last_inner_iter =
          cs_equation_solve_scalar_system(uza->n_u_dofs,
                                          eqp->sles_param,
                                          A,
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux -->already done */
                                          msles->sles,
                                          wk,
                                          uza->rhs));

    _apply_div_op(B_op, wk, dwk); /* -B -w --> dwk has the right sign */

    normalization = _get_cbscal_norm(dwk);

    /* Solve K.zk = dwk */
    memset(zk, 0, sizeof(cs_real_t)*uza->n_p_dofs);

    _n_iter = cs_equation_solve_scalar_cell_system(uza->n_p_dofs,
                                                   schur_slesp,
                                                   K,
                                                   normalization,
                                                   msles->schur_sles,
                                                   zk,
                                                   dwk);
    uza->info->n_inner_iter += _n_iter;
    uza->info->last_inner_iter += _n_iter;

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      zk[ip] = uza->alpha*zk[ip] + uza->inv_mp[ip]*dwk[ip];

    /* Updates
     *  - Compute the rho_factor = <rk,dk> / <dk, dwk>
     *  - u(k+1) = u(k) + rho_factor * wk  --> --wk
     *  - p(k+1) = p(k) - rho_factor * dk
     *  - gk     = gk   - rho_factor * zk
     */
    double rho_factor_denum = cs_gdot(uza->n_p_dofs, dk, dwk);
    assert(fabs(rho_factor_denum) > 0);
    double  rho_factor = cs_gdot(uza->n_p_dofs, rk, dk) / rho_factor_denum;

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      u_f[iu] = u_f[iu] + rho_factor*wk[iu]; /* --wk */

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      p_c[ip] = p_c[ip] - rho_factor*dk[ip];
      gk[ip]  = gk[ip]  - rho_factor*zk[ip];
      rk[ip]  = rk[ip]  - rho_factor*dwk[ip];
    }

    /* Conjugate gradient direction: update dk */

    double  beta_num = cs_gdot(uza->n_p_dofs, rk, gk);
    double  beta_factor = beta_num/uza->info->res;

    uza->info->res = beta_num;

    /* dk <-- gk + beta_factor * dk */
#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      dk[ip] = gk[ip] + beta_factor*dk[ip];

  } /* End of main loop */

  /* Cumulated sum of iterations to solve the Schur complement and the velocity
     block (save before freeing the uzawa structure). */
  int n_inner_iter = uza->info->n_inner_iter;

  /* Last step: Free temporary memory */
  BFT_FREE(diagK);
  BFT_FREE(xtraK);
  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the preconditioned Uzawa-CG algorithm to solve the saddle-point
 *         problem arising from CDO-Fb schemes for Stokes, Oseen and
 *         Navier-Stokes with a monolithic coupling
 *         This algorithm is based on the EDF report H-E40-1991-03299-FR
 *         devoted the numerical algorithms used in the code N3S.
 *         Specifically a Cahout-Chabard preconditioning is used to approximate
 *         the Schur complement.
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_n3s_solve(const cs_navsto_param_t       *nsp,
                                    const cs_equation_param_t     *eqp,
                                    cs_cdofb_monolithic_sles_t    *msles)
{
  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_CG);
  assert(cs_shared_range_set != NULL);

  const cs_real_t  *B_op = msles->div_op;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the Uzawa-CG builder structure */
  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               0, /* grad-div scaling */
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               cs_shared_quant);

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  /* The Schur complement approximation (B.A^-1.Bt) is build and stored in the
     native format */
  cs_matrix_t  *A = msles->block_matrices[0];

  cs_real_t  *diagK = NULL, *xtraK = NULL;
  cs_matrix_t  *K = NULL;

  switch (nsp->sles_param->schur_approximation) {

  case CS_PARAM_SCHUR_DIAG_INVERSE:
    K = _diag_schur_approximation(nsp, A, uza, &diagK, &xtraK);
    break;
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
    K = _invlumped_schur_approximation(nsp, eqp, msles, A,
                                       uza, &diagK, &xtraK);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid Schur approximation.",
              __func__);
  }

  cs_param_sles_t  *schur_slesp = nslesp->schur_sles_param;

  if (msles->schur_sles == NULL) /* has been defined by name */
    msles->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

  /* Compute the first RHS: A.u0 = rhs = b_f - B^t.p_0 to solve */
  _apply_div_op_transpose(B_op, p_c, uza->rhs);

  if (cs_shared_range_set->ifs != NULL) {

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->rhs);

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         b_f);

  }

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
    uza->rhs[iu] = b_f[iu] - uza->rhs[iu];

  /* Initial residual for rhs */
  double normalization = _get_fbvect_norm(uza->rhs);

  /* Compute the first velocity guess */

  /* Modify the tolerance in order to be more accurate on this step */
  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":init_guess") + 1, char);
  sprintf(system_name, "%s:init_guess", eqp->name);

  cs_param_sles_t  *slesp0 = cs_param_sles_create(-1, system_name);

  cs_param_sles_copy_from(eqp->sles_param, slesp0);
  slesp0->eps = 0.05*nslesp->il_algo_rtol;

  uza->info->n_inner_iter
    += (uza->info->last_inner_iter =
        cs_equation_solve_scalar_system(uza->n_u_dofs,
                                        slesp0,
                                        A,
                                        cs_shared_range_set,
                                        normalization,
                                        false, /* rhs_redux --> already done */
                                        msles->sles,
                                        u_f,
                                        uza->rhs));
  /* Partial memory free */
  BFT_FREE(system_name);
  cs_param_sles_free(&slesp0);

  /* Compute the first residual rk0 (in fact the velocity divergence) */
  _apply_div_op(B_op, u_f, uza->d__v);

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    uza->d__v[ip] = uza->d__v[ip] - b_c[ip];

  double div_l2_norm = _get_cbscal_norm(uza->d__v);

  uza->info->res0 = div_l2_norm;
  uza->info->res = div_l2_norm;

  /* Set pointers used in this algorithm */
  cs_real_t  *gk = uza->gk;
  cs_real_t  *dk = uza->res_p;
  cs_real_t  *rk = uza->d__v;
  cs_real_t  *zk = uza->b_tilda;
  cs_real_t  *dzk = uza->dzk;

  /* Compute gk0 as
   *   Solve K.gk00 = rk0
   *   gk0 = alpha gk00 + nu Mp^-1 rk0 */
  memset(gk, 0, sizeof(cs_real_t)*msles->n_cells);

  uza->info->n_inner_iter
    += (uza->info->last_inner_iter =
        cs_equation_solve_scalar_cell_system(uza->n_p_dofs,
                                             schur_slesp,
                                             K,
                                             div_l2_norm,
                                             msles->schur_sles,
                                             gk,
                                             rk));

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    uza->gk[ip] = uza->alpha*uza->gk[ip] + uza->inv_mp[ip]*rk[ip];

  double  rho_factor_num = cs_gdot(uza->n_p_dofs, rk, gk);

  /* dk0 <-- gk0 */
  memcpy(dk, gk, uza->n_p_dofs*sizeof(cs_real_t));

  /* Initialize z0 and compute the rhs for A.z0 = B^t.dk0 */
  memset(zk, 0, sizeof(cs_real_t)*3*msles->n_faces);

  _apply_div_op_transpose(B_op, dk, uza->rhs);

  if (cs_shared_range_set->ifs != NULL)
    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->rhs);

  normalization =_get_fbvect_norm(uza->rhs);

  uza->info->n_inner_iter
    += (uza->info->last_inner_iter =
        cs_equation_solve_scalar_system(uza->n_u_dofs,
                                        eqp->sles_param,
                                        A,
                                        cs_shared_range_set,
                                        normalization,
                                        false, /* rhs_redux --> already done */
                                        msles->sles,
                                        zk,
                                        uza->rhs));

  /* First update (different from the next one)
   * k = 0
   *  - Compute the rho_factor = <r(k),c(k)> / <c(k), div(z(k))>
   *  - u(k+1) = u(k) - rho_factor * z(k)
   *  - p(k+1) = p(k) + rho_factor * d(k)
   *  - r(k+1) = r(k) - rho_factor * div(z(k))
   */

  _apply_div_op(B_op, zk, dzk);

  double  rho_factor_denum = cs_gdot(uza->n_p_dofs, dzk, gk);
  if (rho_factor_denum < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Orthogonality <gk, div(zk)> \approx %5.3e\n"
              " This induces a division by zero during initialization!",
              __func__, rho_factor_denum);
  double rho_factor = rho_factor_num/rho_factor_denum;

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
    u_f[iu] = u_f[iu] - rho_factor*zk[iu];

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
    p_c[ip] = p_c[ip] + rho_factor*dk[ip];
    rk[ip]  = rk[ip]  - rho_factor*dzk[ip];
  }

  /* Main loop */
  /* --------- */

  while (uza->info->cvg == CS_SLES_ITERATING) {

    /* Solve K.gk = rk */
    uza->info->n_inner_iter
      += (uza->info->last_inner_iter =
          cs_equation_solve_scalar_cell_system(uza->n_p_dofs,
                                               schur_slesp,
                                               K,
                                               uza->info->res, /* ||r(k-1)|| */
                                               msles->schur_sles,
                                               gk,
                                               rk));

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      uza->gk[ip] = uza->alpha*uza->gk[ip] + uza->inv_mp[ip]*rk[ip];

    double  prev_rho_factor_num = rho_factor_num;
    rho_factor_num = cs_gdot(uza->n_p_dofs, rk, gk);
    assert(fabs(prev_rho_factor_num) > 0);
    double  lambda_factor = rho_factor_num / prev_rho_factor_num;

    /* dk <-- gk + lambda_factor * dk */
#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      dk[ip] = gk[ip] + lambda_factor*dk[ip];

    /* Compute the rhs for A.zk = B^t.dk */
    _apply_div_op_transpose(B_op, dk, uza->rhs);

    if (cs_shared_range_set->ifs != NULL)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

    normalization = _get_fbvect_norm(uza->rhs);

    /* Solve A.zk = B^t.dk */
    uza->info->n_inner_iter
      += (uza->info->last_inner_iter =
          cs_equation_solve_scalar_system(uza->n_u_dofs,
                                          eqp->sles_param,
                                          A,
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux -->already done */
                                          msles->sles,
                                          zk,
                                          uza->rhs));

    /* Updates
     *  - Compute the rho_factor = <r(k),c(k)> / <c(k), div(z(k))>
     *  - u(k+1) = u(k) - rho_factor * z(k)
     *  - p(k+1) = p(k) + rho_factor * chi(k)
     *  - r(k+1) = r(k) - rho_factor * div(z(k))
     */

    _apply_div_op(B_op, zk, dzk);

    rho_factor_denum = cs_gdot(uza->n_p_dofs, dzk, gk);
    if (fabs(rho_factor_denum) < FLT_MIN)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Orthogonality <gk, div(zk)> \approx %5.3e\n"
                " This induces a division by zero!",
                __func__, rho_factor_denum);
    rho_factor = rho_factor_num/rho_factor_denum;

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      u_f[iu] = u_f[iu] - rho_factor*zk[iu];

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      p_c[ip] = p_c[ip] + rho_factor*dk[ip];
      rk[ip]  = rk[ip]  - rho_factor*dzk[ip];
    }

    /* Update error norm and test if one needs one more iteration */
    _uza_cg_cvg_test(uza);

  } /* End of main loop */


  int n_inner_iter = uza->info->n_inner_iter;

  /* Last step: Free temporary memory */
  BFT_FREE(diagK);
  BFT_FREE(xtraK);
  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the Uzawa algorithm with an Augmented Lagrangian technique to
 *         solve the saddle-point problem arising from CDO-Fb schemes for
 *         Stokes, Oseen and Navier-Stokes with a monolithic coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_al_solve(const cs_navsto_param_t       *nsp,
                                   const cs_equation_param_t     *eqp,
                                   cs_cdofb_monolithic_sles_t    *msles)
{
  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_AL);
  assert(cs_shared_range_set != NULL);

  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_real_t  *div_op = msles->div_op;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the ALU builder structure */
  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               gamma,
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               cs_shared_quant);

  /* Transformation of the initial right-hand side */
  cs_real_t  *btilda_c = uza->d__v;
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    btilda_c[ip] = uza->inv_mp[ip]*b_c[ip];

  _apply_div_op_transpose(div_op, btilda_c, uza->b_tilda);

  if (cs_shared_range_set->ifs != NULL) {

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->b_tilda);

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         b_f);

  }

  /* Update the modify right-hand side: b_tilda = b_f + gamma*Dt.W^-1.b_c */
# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->b_tilda[iu] *= gamma;
    uza->b_tilda[iu] += b_f[iu];
  }

  /* Main loop */
  /* ========= */

  cs_param_sles_t  *slesp = cs_param_sles_create(-1, eqp->sles_param->name);
  cs_param_sles_copy_from(eqp->sles_param, slesp);

  while (uza->info->cvg == CS_SLES_ITERATING) {

    /* Compute the RHS for the Uzawa system: rhs = b_tilda - Dt.p_c */
    _apply_div_op_transpose(div_op, p_c, uza->rhs);

    if (cs_shared_range_set->ifs != NULL)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
      uza->rhs[iu] *= -1;
      uza->rhs[iu] += uza->b_tilda[iu];
    }

    /* Solve AL.u_f = rhs */
    cs_real_t  normalization = 1.; /* No need to change this */
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_CS) {

      /* Inexact algorithm: adaptive tolerance criterion */
      cs_real_t  eps = fmin(1e-2, 0.1*uza->info->res);
      slesp->eps = fmax(eps, eqp->sles_param->eps);
      if (uza->info->verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT, "### UZA.It%02d-- eps=%5.3e\n",
                      uza->info->n_algo_iter, slesp->eps);

    }

    uza->info->n_inner_iter
      += (uza->info->last_inner_iter =
          cs_equation_solve_scalar_system(uza->n_u_dofs,
                                          slesp,
                                          msles->block_matrices[0],
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux */
                                          msles->sles,
                                          u_f,
                                          uza->rhs));

    /* Update p_c = p_c - gamma * (D.u_f - b_c) */
    _apply_div_op(div_op, u_f, uza->d__v);

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      uza->res_p[ip] = uza->d__v[ip] - b_c[ip];
      p_c[ip] += gamma * uza->inv_mp[ip] * uza->res_p[ip];
    }

    /* Update error norm and test if one needs one more iteration */
    _uza_cvg_test(uza);

  }

  int n_inner_iter = uza->info->n_inner_iter;

  /* Last step: Free temporary memory */
  _free_uza_builder(&uza);
  cs_param_sles_free(&slesp);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the Uzawa algorithm with an Augmented Lagrangian technique in an
 *         incremental way to solve the saddle-point problem arising from
 *         CDO-Fb schemes for Stokes, Oseen and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_al_incr_solve(const cs_navsto_param_t       *nsp,
                                        const cs_equation_param_t     *eqp,
                                        cs_cdofb_monolithic_sles_t    *msles)
{
  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_AL);
  assert(cs_shared_range_set != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_real_t  *div_op = msles->div_op;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the ALU builder structure */
  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               gamma,
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               quant);

  /* Transformation of the initial right-hand side */
  cs_real_t  *btilda_c = uza->d__v;
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    btilda_c[ip] = uza->inv_mp[ip]*b_c[ip];

  _apply_div_op_transpose(div_op, btilda_c, uza->b_tilda);

  if (cs_shared_range_set->ifs != NULL) {

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->b_tilda);

    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         b_f);

  }

  /* Update the modify right-hand side: b_tilda = b_f + gamma*Dt.W^-1.b_c */
# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->b_tilda[iu] *= gamma;
    uza->b_tilda[iu] += b_f[iu];
  }

  /* Initialization */
  /* ============== */

  /* Compute the RHS for the Uzawa system: rhs = b_tilda - Dt.p_c */
  _apply_div_op_transpose(div_op, p_c, uza->rhs);

  if (cs_shared_range_set->ifs != NULL)
    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->rhs);

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->rhs[iu] *= -1;
    uza->rhs[iu] += uza->b_tilda[iu];
  }

  /* Solve AL.u_f = rhs */
  /* Modifiy the tolerance in order to be more accurate on this step */
  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":alu0") + 1, char);
  sprintf(system_name, "%s:alu0", eqp->name);

  cs_param_sles_t  *slesp = cs_param_sles_create(-1, system_name);
  cs_param_sles_copy_from(eqp->sles_param, slesp);
  slesp->eps = nsp->sles_param->il_algo_rtol;

  cs_real_t  normalization = cs_evaluate_3_square_wc2x_norm(uza->rhs,
                                                            c2f,
                                                            quant->pvol_fc);

  if (fabs(normalization) > 0)
    normalization = sqrt(normalization);
  else
    normalization = 1.0;

  uza->info->n_inner_iter
    += (uza->info->last_inner_iter =
        cs_equation_solve_scalar_system(uza->n_u_dofs,
                                        slesp,
                                        msles->block_matrices[0],
                                        cs_shared_range_set,
                                        normalization,
                                        false, /* rhs_redux */
                                        msles->sles,
                                        u_f,
                                        uza->rhs));

  /* Partia free */
  BFT_FREE(system_name);
  cs_param_sles_free(&slesp);

  /* Main loop */
  /* ========= */

  cs_real_t  *delta_u = uza->b_tilda;
  cs_real_t  delta_u_l2 = normalization;

  /* Compute the divergence of u since this is a stopping criteria */
  _apply_div_op(div_op, u_f, uza->d__v);

  /* Update p_c = p_c - gamma * (D.u_f - b_c). Recall that B = -div
   * Compute the RHS for the Uzawa system: rhs = -gamma*B^T.W^-1.(B.u - g) */
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
    uza->d__v[ip] -= b_c[ip];
    uza->res_p[ip] = uza->inv_mp[ip] * uza->d__v[ip];
    p_c[ip] += gamma * uza->res_p[ip];
  }

  while (_uza_incr_cvg_test(delta_u_l2, uza)) {

    /* Continue building the RHS */
    _apply_div_op_transpose(div_op, uza->res_p, uza->rhs);

    if (cs_shared_range_set->ifs != NULL)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      uza->rhs[iu] *= -gamma;

    /* Solve AL.u_f = rhs */
    memset(delta_u, 0, sizeof(cs_real_t)*uza->n_u_dofs);

    uza->info->n_inner_iter
      += (uza->info->last_inner_iter =
          cs_equation_solve_scalar_system(uza->n_u_dofs,
                                          eqp->sles_param,
                                          msles->block_matrices[0],
                                          cs_shared_range_set,
                                          delta_u_l2, /* normalization */
                                          false, /* rhs_redux */
                                          msles->sles,
                                          delta_u,
                                          uza->rhs));

    delta_u_l2 = cs_evaluate_3_square_wc2x_norm(delta_u,
                                                c2f, quant->pvol_fc);
    if (fabs(delta_u_l2) > 0)
      delta_u_l2 = sqrt(delta_u_l2);
    else
      delta_u_l2 = 1.0;

    /* Update the velocity */
#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      u_f[iu] += delta_u[iu];

    /* Update the divergence */
    _apply_div_op(div_op, u_f, uza->d__v);

    /* Update p_c = p_c - gamma * (D.u_f - b_c). Recall that B = -div
     * Prepare the computation of the RHS for the Uzawa system:
     * rhs = -gamma*B^T.W^-1.(B.u - g) */
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      uza->d__v[ip] -= b_c[ip];
      uza->res_p[ip] = uza->inv_mp[ip] * uza->d__v[ip];
      p_c[ip] += gamma * uza->res_p[ip];
    }

  } /* End of Uzawa iterations */

  int n_inner_iter = uza->info->n_inner_iter;

  /* Last step: Free temporary memory */
  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
