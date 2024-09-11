/*============================================================================
 * Matrix building
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

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
#include "cs_dispatch.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_gradient.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_building.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_matrix_building.c

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a scalar field.
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 *                               (implicit part)
 * \param[in]     rovsdt        working array
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

static void
_sym_matrix_scalar(const cs_mesh_t            *m,
                   int                         idiffp,
                   double                      thetap,
                   const cs_field_bc_coeffs_t *bc_coeffs,
                   const cs_real_t             rovsdt[],
                   const cs_real_t             i_visc[],
                   const cs_real_t             b_visc[],
                   cs_real_t         *restrict da,
                   cs_real_t         *restrict xa)
{
  const cs_real_t *cofbfp = (const cs_real_t *)bc_coeffs->bf;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    da[c_id] = rovsdt[c_id];
  });

  /* Computation of extradiagonal terms and contribution to the diagonal */

  if (idiffp) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      xa[f_id] = -thetap*i_visc[f_id];

      if (ii < n_cells)
        cs_dispatch_sum(&da[ii], -xa[f_id], i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&da[jj], -xa[f_id], i_sum_type);

    });

  }

  else {

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      xa[f_id] = 0.;
    });

  }

  /* Contribution of border faces to the diagonal */

  if (idiffp) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = b_face_cells[f_id];

      cs_dispatch_sum(&da[ii],
                      thetap*b_visc[f_id]*cofbfp[f_id],
                      b_sum_type);

    });
  }

  ctx.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a scalar field
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension 2 time the number of internal
 * faces).
 *
 * template parameters:
 *   is_thermal        true for the temperature, otherwise false
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convective term by Cp
 *                               - 1 do multiply the convective term by Cp
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     rovsdt        working array
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <bool is_thermal>
static void
_matrix_scalar(const cs_mesh_t            *m,
               int                         iconvp,
               int                         idiffp,
               double                      thetap,
               const cs_field_bc_coeffs_t *bc_coeffs,
               const cs_real_t             rovsdt[],
               const cs_real_t             i_massflux[],
               const cs_real_t             b_massflux[],
               const cs_real_t             i_visc[],
               const cs_real_t             b_visc[],
               const cs_real_t             xcpp[],
               cs_real_t         *restrict da,
               cs_real_2_t       *restrict xa)
{
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    da[c_id] = rovsdt[c_id];
  });

  /* Contribution of the extra-diagonal terms to the diagonal */

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    cs_real_t _i_massflux = i_massflux[f_id];
    cs_real_t flui =  0.5*(_i_massflux - cs_math_fabs(_i_massflux));
    cs_real_t fluj = -0.5*(_i_massflux + cs_math_fabs(_i_massflux));

    cs_real_t cpi = 1.0, cpj = 1.0;
    /* When solving the temperature, the convective part is multiplied by Cp */
    if (is_thermal) {
      cpi = xcpp[ii];
      cpj = xcpp[jj];
    }

    /* Computation of extradiagonal terms */

    xa[f_id][0] = thetap*(iconvp*cpi*flui -idiffp*i_visc[f_id]);
    xa[f_id][1] = thetap*(iconvp*cpj*fluj -idiffp*i_visc[f_id]);

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */

    cs_real_t ifac = xa[f_id][0] + iconvp * (1.-thetap) * cpi * _i_massflux;
    cs_real_t jfac = xa[f_id][1] - iconvp * (1.-thetap) * cpj * _i_massflux;

    if (ii < n_cells)
      cs_dispatch_sum(&da[ii], -ifac, i_sum_type);
    if (jj < n_cells)
      cs_dispatch_sum(&da[jj], -jfac, i_sum_type);

   });

  /* Contribution of border faces to the diagonal */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = b_face_cells[f_id];

    cs_real_t _b_massflux = b_massflux[f_id];
    cs_real_t flui = 0.5*(_b_massflux - cs_math_fabs(_b_massflux));

    cs_real_t cpi = 1.0;
    /* When solving the temperature, the convective part is multiplied by Cp */
    if (is_thermal)
      cpi = xcpp[ii];

    /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
     *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
     *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
     */

    cs_real_t bfac = iconvp * cpi *(  flui * thetap * (coefbp[f_id] - 1.)
                                    - (1. - thetap) * _b_massflux)
                   + idiffp * thetap * b_visc[f_id] * cofbfp[f_id];

    cs_dispatch_sum(&da[ii], bfac, b_sum_type);

  });

  ctx.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a vector field
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * template parameters:
 *   stride        3 for vectors, 6 for tensors
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sym_matrix_strided(const cs_mesh_t            *m,
                    int                         idiffp,
                    cs_real_t                   thetap,
                    const cs_field_bc_coeffs_t *bc_coeffs_v,
                    const cs_real_t             fimp[][stride][stride],
                    const cs_real_t             i_visc[],
                    const cs_real_t             b_visc[],
                    cs_real_t        (*restrict da)[stride][stride],
                    cs_real_t        (*restrict xa))
{
  using b_t = cs_real_t[stride][stride];
  const b_t *cofbfp = (const b_t *)bc_coeffs_v->bf;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++)
      for (cs_lnum_t j = 0; j < stride; j++)
        da[c_id][i][j] = fimp[c_id][i][j];
  });

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    /* Computation of extradiagonal terms */

    xa[f_id] = -thetap*idiffp*i_visc[f_id];

    /* Contribution of the extra-diagonal terms to the diagonal */

    if (ii < n_cells)
      for (cs_lnum_t i = 0; i < stride; i++)
        cs_dispatch_sum(&da[ii][i][i], -xa[f_id], i_sum_type);

    if (jj < n_cells)
      for (cs_lnum_t i = 0; i < stride; i++)
        cs_dispatch_sum(&da[jj][i][i], -xa[f_id], i_sum_type);

  });

  /* Contribution of border faces to the diagonal */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = b_face_cells[f_id];
    cs_real_t bfac[stride*stride];

    for (cs_lnum_t i = 0; i < stride; i++)
      for (cs_lnum_t j = 0; j < stride; j++)
        bfac[stride*i+j] = thetap * idiffp * b_visc[f_id] * cofbfp[f_id][i][j];

    cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                   bfac,
                                   b_sum_type);

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a vector field
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 2 time the number of internal
 * faces).
 *
 * template parameters:
 *   eb_size       extra-diagonal block size
 *                 - 1 when blocks are of the form aij.I
 *                 - 3 for full blocks
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     mq            pointer to mesh quantities structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     fimp          part of the diagonal
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride, cs_lnum_t eb_size>
static void
_matrix_strided(const cs_mesh_t            *m,
                const cs_mesh_quantities_t *mq,
                int                         iconvp,
                int                         idiffp,
                double                      thetap,
                const cs_field_bc_coeffs_t *bc_coeffs_v,
                const cs_real_t             fimp[][stride][stride],
                const cs_real_t             i_massflux[],
                const cs_real_t             b_massflux[],
                const cs_real_t             i_visc[],
                const cs_real_t             b_visc[],
                cs_real_t        (*restrict da)[stride][stride],
                cs_real_t        (*restrict xa)[2])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;
  const cs_real_3_t *i_face_u_normal
    = (const cs_real_3_t *)mq->i_face_u_normal;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)mq->b_face_u_normal;
  cs_real_2_t *i_f_face_factor;
  cs_real_t *b_f_face_factor;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Discontinuous porous treatment */
  cs_real_2_t *_i_f_face_factor = nullptr;
  cs_real_t *_b_f_face_factor = nullptr;
  int is_p = 0;

  if (cs_glob_porous_model == 3 && stride == 3) {
    i_f_face_factor = mq->i_f_face_factor;
    b_f_face_factor = mq->b_f_face_factor;
    is_p = 1;
  }
  else {
    CS_MALLOC_HD(_i_f_face_factor, 1, cs_real_2_t, cs_alloc_mode);
    CS_MALLOC_HD(_b_f_face_factor, 1, cs_real_t, cs_alloc_mode);
    i_f_face_factor = _i_f_face_factor;
    b_f_face_factor = _b_f_face_factor;
    i_f_face_factor[0][0] = i_f_face_factor[0][1] = 1.0;
    b_f_face_factor[0] = 1.0;
  }

  cs_real_332_t *_xa = (cs_real_332_t *) xa; //FIXME why 332 and use 233...

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < stride; j++) {
        da[c_id][i][j] = fimp[c_id][i][j];
      }
    }
  });

  if ((stride == 3 && eb_size == 1) || stride == 6) {
    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      const cs_lnum_t _p = is_p*f_id;
      cs_real_t _i_massflux = i_massflux[f_id];
      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t flu[2] = {0.5*iconvp*(_i_massflux - cs_math_fabs(_i_massflux)),
                         -0.5*iconvp*(_i_massflux + cs_math_fabs(_i_massflux))};

      /* Computation of extradiagonal terms */
      /*
       * X_ij = - theta f_j (m_ij)^-
       * X_ji = - theta f_i (m_ij)^+
       */

      xa[f_id][0] = thetap * (flu[0] - idiffp * i_visc[f_id])
        * i_f_face_factor[_p][1]; //FIXME also diffusion? MF thinks so

      xa[f_id][1] = thetap * (flu[1] - idiffp * i_visc[f_id])
        * i_f_face_factor[_p][0];

      /* Contribution of the extra-diagonal terms to the diagonal */

      /* D_ii =  theta f_i (m_ij)^+ - m_ij
       *      = -X_ij - (1-theta)*m_ij
       *      = -X_ji - m_ij
       * D_jj = -theta f_j (m_ij)^- + m_ij
       *      = -X_ji + (1-theta)*m_ij
       *      = -X_ij + m_ij
       */

      cs_real_t ifac = xa[f_id][1] + iconvp * i_massflux[f_id];
      cs_real_t jfac = xa[f_id][0] - iconvp * i_massflux[f_id];

      if (ii < n_cells)
        for (cs_lnum_t i = 0; i < stride; i++)
          cs_dispatch_sum(&da[ii][i][i], -ifac, i_sum_type);

      if (jj < n_cells)
        for (cs_lnum_t i = 0; i < stride; i++)
          cs_dispatch_sum(&da[jj][i][i], -jfac, i_sum_type);

    });
  }
  else {
    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      const cs_lnum_t _p = is_p*f_id;
      const cs_real_t *n = i_face_u_normal[f_id];
      const cs_real_t _i_massflux = i_massflux[f_id];
      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t flu[2] = {
         0.5 * iconvp * (_i_massflux - cs_math_fabs(_i_massflux))
          - idiffp*i_visc[f_id],
        -0.5 * iconvp * (_i_massflux + cs_math_fabs(_i_massflux))
          - idiffp*i_visc[f_id]
      };

      /* Computation of extradiagonal terms */
      /*
       * X_ij = - theta f_j (m_ij)^-
       * X_ji = - theta f_i (m_ij)^+
       */

      /* Diagonal part:
       * the n(x)n term is multiplied by i_f_face_factor and (1 - n(x)n) by 1
       * XA_ij <= XA_ik n_k n_j (factor - 1) + XA_ij
       * XA_ij used to be diagonal: XA_ik n_k n_j = XA_ii n_i n_j*/

      for (cs_lnum_t i = 0; i < eb_size; i++) {
        for (cs_lnum_t j = 0; j < eb_size; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);

          _xa[f_id][0][i][j] = thetap * flu[0]
            * (d_ij + (i_f_face_factor[_p][1] - 1.) * n[i] * n[j]);

          //FIXME also diffusion? MF thinks so
          _xa[f_id][1][i][j] = thetap * flu[1]
            * (d_ij + (i_f_face_factor[_p][0] - 1.) * n[i] * n[j]);
        }
      }

      /* D_ii =  theta (m_ij)^+ - m_ij
       *      = -X_ij - (1-theta)*m_ij
       *      = -X_ji - m_ij
       * D_jj = -theta (m_ij)^- + m_ij
       *      = -X_ji + (1-theta)*m_ij
       *      = -X_ij + m_ij
       */

      cs_real_t vfaci[eb_size*eb_size], vfacj[eb_size*eb_size];

      for (cs_lnum_t i = 0; i < eb_size; i++) {
        for (cs_lnum_t j = 0; j < eb_size; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);
          cs_real_t diag = d_ij * iconvp * _i_massflux;

          vfaci[stride*i+j] = - diag -_xa[f_id][1][i][j];

          vfacj[stride*i+j] =   diag -_xa[f_id][0][i][j];

        }
      }

      if (ii < n_cells)
        cs_dispatch_sum<eb_size*eb_size>(reinterpret_cast<cs_real_t*>(da[ii]),
                                         vfaci,
                                         i_sum_type);

      if (jj < n_cells)
        cs_dispatch_sum<eb_size*eb_size>(reinterpret_cast<cs_real_t*>(da[jj]),
                                         vfacj,
                                         i_sum_type);

    });
  }

  /* Contribution of border faces to the diagonal */

  if (stride == 3) {

    const cs_real_33_t *coefbp = (const cs_real_33_t *)bc_coeffs_v->b;
    const cs_real_33_t *cofbfp = (const cs_real_33_t *)bc_coeffs_v->bf;

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      const cs_real_t *n = b_face_u_normal[f_id];
      const cs_lnum_t _p = is_p*f_id;
      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t flu[2] = {
        /* (m_ij)^+ */
        iconvp * 0.5 * (b_massflux[f_id] + cs_math_fabs(b_massflux[f_id])),
        /* (m_ij)^- */
        iconvp * 0.5 * (b_massflux[f_id] - cs_math_fabs(b_massflux[f_id]))
      };

      cs_real_t n_b_n
        = cs_math_3_33_3_dot_product(n, coefbp[f_id], n);
      cs_real_t n_bf_n
        = cs_math_3_33_3_dot_product(n, cofbfp[f_id], n);

      cs_real_t bfac[stride*stride];

      for (cs_lnum_t i = 0; i < stride; i++) {

        for (cs_lnum_t j = 0; j < stride; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);
          /* D = theta (m_f)^+.1 + theta B (m_f)^- - m_f.1
           * NB: stop here because the first two terms maybe scaled
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */

          bfac[stride*i+j]
            =  thetap * (
                   d_ij * flu[0]
                   + flu[1] * coefbp[f_id][i][j]
                   + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]
                   + (flu[0] + flu[1] * n_b_n + idiffp * b_visc[f_id] * n_bf_n)
                   * (b_f_face_factor[_p] - 1.) * n[i] * n[j]
                   )
             - iconvp * d_ij * b_massflux[f_id];
        }
      }

      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     bfac,
                                     b_sum_type);

    });
  }
  else if (stride == 6) {

    const cs_real_66_t *coefbp = (const cs_real_66_t *)bc_coeffs_v->b;
    const cs_real_66_t *cofbfp = (const cs_real_66_t *)bc_coeffs_v->bf;

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      cs_real_t _b_massflux = b_massflux[f_id];
      cs_lnum_t ii = b_face_cells[f_id];
      cs_real_t flui = 0.5*(_b_massflux - cs_math_fabs(_b_massflux));

      cs_real_t bfac[stride*stride];

      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);
          /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */

          cs_real_t diag = d_ij * iconvp
            * (thetap * flui + (1. - thetap) * iconvp * b_massflux[f_id]);

          bfac[stride*i+j]
            = - diag + thetap * (  iconvp * flui * coefbp[f_id][i][j]
                                 + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]);
        }
      }

      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     bfac,
                                     b_sum_type);
    });

  }

  ctx.wait();

  CS_FREE_HD(_i_f_face_factor);
  CS_FREE_HD(_b_f_face_factor);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a vector field with a
 * tensorial diffusivity (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 3x3 the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sym_matrix_anisotropic_diffusion_strided
  (const cs_mesh_t            *m,
   int                         idiffp,
   double                      thetap,
   const cs_field_bc_coeffs_t *bc_coeffs_v,
   const cs_real_t             fimp[][stride][stride],
   const cs_real_t             i_visc[][stride][stride],
   const cs_real_t             b_visc[],
   cs_real_t        (*restrict da)[stride][stride],
   cs_real_t        (*restrict xa)[stride][stride])
{
  using b_t = cs_real_t[stride][stride];
  const b_t *cofbfp = (const b_t *)bc_coeffs_v->bf;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++)
      for (cs_lnum_t j = 0; j < stride; j++)
        da[c_id][i][j] = fimp[c_id][i][j];
  });

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    /* Computation of extradiagonal terms */

    cs_real_t vfac[stride*stride];
    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < stride; j++) {
        vfac[stride*i+j] = thetap*idiffp*i_visc[f_id][i][j];
        xa[f_id][i][j] = -vfac[stride*i+j];
      }
    }

    /* Contribution of the extra-diagonal terms to the diagonal */

    if (ii < n_cells)
      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     vfac,
                                     i_sum_type);

    if (jj < n_cells)
      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[jj]),
                                     vfac,
                                     i_sum_type);
  });

  /* Contribution of border faces to the diagonal */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = b_face_cells[f_id];
    cs_real_t bfac[stride*stride];

    for (cs_lnum_t i = 0; i < stride; i++)
      for (cs_lnum_t j = 0; j < stride; j++)
        bfac[stride*i+j] = thetap * idiffp * b_visc[f_id] * cofbfp[f_id][i][j];

    cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                   bfac,
                                   b_sum_type);

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a vector field with a
 * tensorial diffusivity.
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 2 times 3x3 the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     mq            pointer to mesh quantities structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     fimp          part of the diagonal
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_matrix_anisotropic_diffusion_strided
(
 const cs_mesh_t            *m,
 const cs_mesh_quantities_t *mq,
 int                         iconvp,
 int                         idiffp,
 double                      thetap,
 const cs_field_bc_coeffs_t *bc_coeffs_v,
 const cs_real_t             fimp[][stride][stride],
 const cs_real_t             i_massflux[],
 const cs_real_t             b_massflux[],
 const cs_real_t             i_visc[][stride][stride],
 const cs_real_t             b_visc[],
 cs_real_t        (*restrict da)[stride][stride],
 cs_real_t        (*restrict xa)[2][stride][stride]
)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;
  const cs_real_3_t *i_face_u_normal
    = (const cs_real_3_t *)mq->i_face_u_normal;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)mq->b_face_u_normal;

  cs_real_2_t *i_f_face_factor = nullptr;
  cs_real_t *b_f_face_factor = nullptr;
  cs_real_2_t *_i_f_face_factor = nullptr;
  cs_real_t *_b_f_face_factor = nullptr;
  int is_p = 0;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Discontinuous porous treatment */

  /* Is it porous? */
  if (cs_glob_porous_model == 3 && stride == 3) {
    i_f_face_factor = mq->i_f_face_factor;
    b_f_face_factor = mq->b_f_face_factor;
    is_p = 1;
  }
  else {
    CS_MALLOC_HD(_i_f_face_factor, 1, cs_real_2_t, cs_alloc_mode);
    CS_MALLOC_HD(_b_f_face_factor, 1, cs_real_t, cs_alloc_mode);
    i_f_face_factor = _i_f_face_factor;
    b_f_face_factor = _b_f_face_factor;
    i_f_face_factor[0][0] = i_f_face_factor[0][1] = 1.0;
    b_f_face_factor[0] = 1.0;
  }

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++)
      for (cs_lnum_t j = 0; j < stride; j++)
        da[c_id][i][j] = fimp[c_id][i][j];
  });

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    const cs_lnum_t _p = is_p*f_id;
    const cs_real_t *n = i_face_u_normal[f_id];
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    cs_real_t flu[2] = {
       iconvp * 0.5*(i_massflux[f_id] - cs_math_fabs(i_massflux[f_id])),
      -iconvp * 0.5*(i_massflux[f_id] + cs_math_fabs(i_massflux[f_id]))
    };

    /* Computation of extradiagonal terms */

    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < stride; j++) {
        cs_real_t d_ij = ((i == j) ? 1. : 0.);

        xa[f_id][0][i][j] = thetap*(  d_ij * flu[0]
                                    + (i_f_face_factor[_p][0] - 1.)
                                      * n[i] * n[j] * flu[0]
                                    //FIXME also diffusion? MF thinks so
                                    - idiffp*i_visc[f_id][i][j]);

        xa[f_id][1][i][j] = thetap*(  d_ij * flu[1]
                                    + (i_f_face_factor[_p][1] - 1.)
                                      * n[i] * n[j] * flu[1]
                                    - idiffp*i_visc[f_id][i][j]);
      }
    }

    /* Contribution of the extra-diagonal terms to the diagonal */

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */

    cs_real_t vfaci[stride*stride], vfacj[stride*stride];

    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < stride; j++) {
        cs_real_t d_ij = ((i == j) ? 1. : 0.);

        cs_real_t diag = d_ij * iconvp * (1. - thetap) * i_massflux[f_id];

        vfaci[stride*i+j] = - diag - xa[f_id][0][i][j];

        vfacj[stride*i+j] =   diag - xa[f_id][1][i][j];
      }
    }

    if (ii < n_cells)
      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     vfaci,
                                     i_sum_type);

    if (jj < n_cells)
      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[jj]),
                                     vfacj,
                                     i_sum_type);

  });

  /* Contribution of border faces to the diagonal */

  if (stride == 3) {

    const cs_real_33_t *coefbp = (const cs_real_33_t *)bc_coeffs_v->b;
    const cs_real_33_t *cofbfp = (const cs_real_33_t *)bc_coeffs_v->bf;

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      const cs_lnum_t _p = is_p*f_id;
      const cs_real_t *n = b_face_u_normal[f_id];
      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t flu[2] = {
         /* (m_ij)^+ */
         iconvp * 0.5 * (b_massflux[f_id] + cs_math_fabs(b_massflux[f_id])),
         /* (m_ij)^- */
         iconvp * 0.5 * (b_massflux[f_id] - cs_math_fabs(b_massflux[f_id]))
      };

      cs_real_t n_b_n
        = cs_math_3_33_3_dot_product(n, coefbp[f_id], n);
      cs_real_t n_bf_n
        = cs_math_3_33_3_dot_product(n, cofbfp[f_id], n);

      cs_real_t bfac[stride*stride];

      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);
          /* D = theta (m_f)^+.1 + theta B (m_f)^- - m_f.1
           * NB: stop here because the first two terms maybe scaled
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */
          bfac[stride*i+j] =
            thetap * (
                d_ij * flu[0]
                + flu[1] * coefbp[f_id][i][j]
                + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]
                + (flu[0] + flu[1] * n_b_n + idiffp * b_visc[f_id] * n_bf_n)
                * (b_f_face_factor[_p] - 1.) * n[i] * n[j]
                )
            - iconvp * d_ij * b_massflux[f_id];
        }
      }

      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     bfac,
                                     b_sum_type);

    });
  }
  else if (stride == 6) {

    const cs_real_66_t *coefbp = (const cs_real_66_t *)bc_coeffs_v->b;
    const cs_real_66_t *cofbfp = (const cs_real_66_t *)bc_coeffs_v->bf;

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = b_face_cells[f_id];
      cs_real_t flui = 0.5*(b_massflux[f_id] - cs_math_fabs(b_massflux[f_id]));
      cs_real_t bfac[stride*stride];

      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++) {
          cs_real_t d_ij = ((i == j) ? 1. : 0.);
          /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */

          bfac[stride*i+j] = - d_ij * iconvp
                  * (thetap * flui + (1. - thetap) * iconvp * b_massflux[f_id])
                    + thetap * (  iconvp * flui * coefbp[f_id][i][j]
                                + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]);
        }
      }

      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     bfac,
                                     b_sum_type);


    });

  }

  ctx.wait();

  CS_FREE_HD(_i_f_face_factor);
  CS_FREE_HD(_b_f_face_factor);
}

#endif /* cplusplus */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_scalar(int                         iconvp,
                         int                         idiffp,
                         int                         ndircp,
                         int                         isym,
                         double                      thetap,
                         int                         imucpp,
                         const cs_field_bc_coeffs_t *bc_coeffs,
                         const cs_real_t             rovsdt[],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         const cs_real_t             xcpp[],
                         cs_real_t                   da[],
                         cs_real_t                   xa[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    _sym_matrix_scalar(m,
                       idiffp,
                       thetap,
                       bc_coeffs,
                       rovsdt,
                       i_visc,
                       b_visc,
                       da,
                       xa);
  }

  /* Non-symmetric matrix */
  else {
    if (imucpp == 0)
      _matrix_scalar<false>(m,
                            iconvp,
                            idiffp,
                            thetap,
                            bc_coeffs,
                            rovsdt,
                            i_massflux,
                            b_massflux,
                            i_visc,
                            b_visc,
                            xcpp,
                            da,
                            (cs_real_2_t*) xa);
    else
      _matrix_scalar<true>(m,
                           iconvp,
                           idiffp,
                           thetap,
                           bc_coeffs,
                           rovsdt,
                           i_massflux,
                           b_massflux,
                           i_visc,
                           b_visc,
                           xcpp,
                           da,
                           (cs_real_2_t*) xa);
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const cs_real_t epsi = 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      da[c_id] *= (1. + epsi);
    });
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      da[c_id] += (cs_real_t)c_disable_flag[c_id];
    });
  }
  ctx.wait();

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_vector (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_vector(int                         iconvp,
                         int                         idiffp,
                         int                         tensorial_diffusion,
                         int                         ndircp,
                         int                         isym,
                         cs_lnum_t                   eb_size,
                         double                      thetap,
                         const cs_field_bc_coeffs_t *bc_coeffs_v,
                         const cs_real_t             fimp[][3][3],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         cs_real_t                   da[][3][3],
                         cs_real_t                   xa[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* scalar diffusion or right anisotropic diffusion */
  if (tensorial_diffusion == 1) {
    /* Symmetric matrix */
    if (isym == 1) {
      assert(eb_size == 1);
      _sym_matrix_strided<3>(m,
                             idiffp,
                             thetap,
                             bc_coeffs_v,
                             fimp,
                             i_visc,
                             b_visc,
                             da,
                             xa);

    /* Non-symmetric matrix */
    }
    else {
      if (eb_size == 1)
        _matrix_strided<3,1>(m,
                             mq,
                             iconvp,
                             idiffp,
                             thetap,
                             bc_coeffs_v,
                             fimp,
                             i_massflux,
                             b_massflux,
                             i_visc,
                             b_visc,
                             da,
                             (cs_real_2_t*) xa);
      else
        _matrix_strided<3,3>(m,
                             mq,
                             iconvp,
                             idiffp,
                             thetap,
                             bc_coeffs_v,
                             fimp,
                             i_massflux,
                             b_massflux,
                             i_visc,
                             b_visc,
                             da,
                             (cs_real_2_t*) xa);
    }
  }
  /* left tensor diffusion */
  else {

    /* Symmetric matrix */
    if (isym == 1) {
      _sym_matrix_anisotropic_diffusion_strided<3>(m,
                                                   idiffp,
                                                   thetap,
                                                   bc_coeffs_v,
                                                   fimp,
                                                   (const cs_real_33_t *)i_visc,
                                                   b_visc,
                                                   da,
                                                   (cs_real_33_t *) xa);

    /* Non-symmetric matrix */
    } else {
      _matrix_anisotropic_diffusion_strided<3>(m,
                                               mq,
                                               iconvp,
                                               idiffp,
                                               thetap,
                                               bc_coeffs_v,
                                               fimp,
                                               i_massflux,
                                               b_massflux,
                                               (const cs_real_33_t *)i_visc,
                                               b_visc,
                                               da,
                                               (cs_real_332_t *) xa);
    }

  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if (ndircp <= 0) {
    const cs_real_t epsi = 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        da[c_id][i][i] = (1. + epsi) * da[c_id][i][i];
    });
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        da[c_id][i][i] += (cs_real_t)c_disable_flag[c_id];
    });
  }
  ctx.wait();

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_tensor (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_tensor(int                         iconvp,
                         int                         idiffp,
                         int                         tensorial_diffusion,
                         int                         ndircp,
                         int                         isym,
                         double                      thetap,
                         const cs_field_bc_coeffs_t *bc_coeffs_ts,
                         const cs_real_t             fimp[][6][6],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         cs_real_t                   da[][6][6],
                         cs_real_t                   xa[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* scalar diffusion or right anisotropic diffusion */
  if (tensorial_diffusion == 1) {
    /* Symmetric matrix */
    if (isym == 1) {
      _sym_matrix_strided<6>(m,
                             idiffp,
                             thetap,
                             bc_coeffs_ts,
                             fimp,
                             i_visc,
                             b_visc,
                             da,
                             xa);

    /* Non-symmetric matrix */
    } else {
      _matrix_strided<6,1>(m,
                           mq,
                           iconvp,
                           idiffp,
                           thetap,
                           bc_coeffs_ts,
                           fimp,
                           i_massflux,
                           b_massflux,
                           i_visc,
                           b_visc,
                           da,
                           (cs_real_2_t*) xa);
    }
  }
  /* left tensor diffusion */
  else {
    /* Symmetric matrix */
    if (isym == 1) {
      _sym_matrix_anisotropic_diffusion_strided<6>(m,
                                                   idiffp,
                                                   thetap,
                                                   bc_coeffs_ts,
                                                   fimp,
                                                   (const cs_real_66_t *)i_visc,
                                                   b_visc,
                                                   da,
                                                   (cs_real_66_t *)xa);

    /* Non-symmetric matrix */
    }
    else {
      _matrix_anisotropic_diffusion_strided<6>(m,
                                               mq,
                                               iconvp,
                                               idiffp,
                                               thetap,
                                               bc_coeffs_ts,
                                               fimp,
                                               i_massflux,
                                               b_massflux,
                                               (const cs_real_66_t *)i_visc,
                                               b_visc,
                                               da,
                                               (cs_real_662_t *)xa);
    }
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const cs_real_t epsi = 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 6; i++)
        da[c_id][i][i] = (1. + epsi)*da[c_id][i][i];
    });
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 6; i++)
        da[c_id][i][i] += (cs_real_t)c_disable_flag[c_id];
    });
  }
  ctx.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_time_step(const cs_mesh_t            *m,
                    int                         iconvp,
                    int                         idiffp,
                    int                         isym,
                    const cs_field_bc_coeffs_t *bc_coeffs,
                    const cs_real_t             i_massflux[],
                    const cs_real_t             b_massflux[],
                    const cs_real_t             i_visc[],
                    const cs_real_t             b_visc[],
                    cs_real_t         *restrict da)
{
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  /* 1. Initialization */

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    da[c_id] = 0.;
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = n_cells; c_id < n_cells_ext; c_id++) {
      da[c_id] = 0.;
    }
  }

  /* 2. Computation of extradiagonal terms unnecessary */

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  if (isym == 2) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t fluj =-0.5*(i_massflux[f_id] + fabs(i_massflux[f_id]));
          cs_real_t flui = 0.5*(i_massflux[f_id] - fabs(i_massflux[f_id]));

          cs_real_t xaifa2 = iconvp*fluj -idiffp*i_visc[f_id];
          cs_real_t xaifa1 = iconvp*flui -idiffp*i_visc[f_id];
          da[ii] -= xaifa2;
          da[jj] -= xaifa1;

        }
      }
    }

  } else {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t f_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             f_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             f_id++) {

          cs_lnum_t ii = i_face_cells[f_id][0];
          cs_lnum_t jj = i_face_cells[f_id][1];

          cs_real_t flui = 0.5*(i_massflux[f_id] - fabs(i_massflux[f_id]));

          cs_real_t xaifa1 = iconvp*flui -idiffp*i_visc[f_id];
          da[ii] -= xaifa1;
          da[jj] -= xaifa1;

        }
      }
    }

  }

  /* 4. Contribution of border faces to the diagonal */

# pragma omp parallel for if (m->n_b_faces > CS_THR_MIN)
  for (int t_id = 0; t_id < n_b_threads; t_id++) {
    for (cs_lnum_t f_id = b_group_index[t_id*2];
         f_id < b_group_index[t_id*2 + 1];
         f_id++) {

      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t flui =  0.5*(b_massflux[f_id] - fabs(b_massflux[f_id]));
      cs_real_t fluj = -0.5*(b_massflux[f_id] + fabs(b_massflux[f_id]));

      da[ii] +=   iconvp*(-fluj + flui*coefbp[f_id])
                + idiffp*b_visc[f_id]*cofbfp[f_id];

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
