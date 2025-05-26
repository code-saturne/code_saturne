/*============================================================================
 * Matrix building
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

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

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "alge/cs_blas.h"
#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh.h"
#include "base/cs_field.h"
#include "alge/cs_gradient.h"
#include "base/cs_ext_neighborhood.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_timer.h"

#include "base/cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_matrix_building.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_matrix_building.cpp

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
 * \param[in]     ctx           dispatch context
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
 * \param[out]    da             diagonal part of the matrix
 * \param[out]    ea             extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

static void
_sym_coeffs_msr(const cs_mesh_t            *m,
                cs_dispatch_context        &ctx,
                int                         idiffp,
                double                      thetap,
                const cs_field_bc_coeffs_t *bc_coeffs,
                const cs_real_t             rovsdt[],
                const cs_real_t             i_visc[],
                const cs_real_t             b_visc[],
                cs_real_t         *restrict da,
                cs_real_t         *restrict ea)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  if (idiffp) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t _da = rovsdt[c_id];

      /* Loop on interior faces */
      const cs_lnum_t s_id_i = c2c_idx[c_id];
      const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

      for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
        const cs_lnum_t f_id = cell_i_faces[cidx];
        cs_real_t _xa = -thetap*i_visc[f_id];
        ea[cidx] = _xa;
        _da -= _xa;
      }

      da[c_id] = _da;

    });

  }

  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      da[c_id] = rovsdt[c_id];

      /* Loop on interior faces */
      const cs_lnum_t s_id_i = c2c_idx[c_id];
      const cs_lnum_t e_id_i = c2c_idx[c_id + 1];
      for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
        ea[cidx] = 0.;
      }

    });

  }

  /* Contribution of border faces to the diagonal */

  if (idiffp) {
    const cs_real_t *cofbfp = (const cs_real_t *)bc_coeffs->bf;
    cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

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
 * \brief Build the diffusion matrix for a scalar field.
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     ctx           dispatch context
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$, or null
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da             diagonal part of the matrix
 * \param[out]    ea             extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sym_coeffs_msr(const cs_mesh_t            *m,
                cs_dispatch_context        &ctx,
                int                         idiffp,
                double                      thetap,
                const cs_field_bc_coeffs_t *bc_coeffs,
                const cs_real_t             fimp[][stride][stride],
                const cs_real_t             i_visc[],
                const cs_real_t             b_visc[],
                cs_real_t                   da[][stride][stride],
                cs_real_t         *restrict ea)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  if (idiffp) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      auto &_da = da[c_id];

      if (fimp != nullptr) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < stride; j++)
            _da[i][j] = fimp[c_id][i][j];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < stride; j++)
            _da[i][j] = 0.;
        }
      }

      /* Loop on interior faces */
      const cs_lnum_t s_id_i = c2c_idx[c_id];
      const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

      for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
        const cs_lnum_t f_id = cell_i_faces[cidx];
        cs_real_t _xa = -thetap*i_visc[f_id];
        ea[cidx] = _xa;
        for (cs_lnum_t i = 0; i < stride; i++) {
          _da[i][i] -= _xa;
        }
      }

    });

  }

  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      auto &_da = da[c_id];

      if (fimp != nullptr) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < stride; j++)
            _da[i][j] = fimp[c_id][i][j];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < stride; i++) {
          for (cs_lnum_t j = 0; j < stride; j++)
            _da[i][j] = 0.;
        }
      }

      /* Loop on interior faces */
      const cs_lnum_t s_id_i = c2c_idx[c_id];
      const cs_lnum_t e_id_i = c2c_idx[c_id + 1];
      for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
        ea[cidx] = 0.;
      }

    });

  }

  /* Contribution of boundary faces to the diagonal */

  if (idiffp) {
    using b_t = cs_real_t[stride][stride];
    const b_t *cofbfp = (const b_t *)bc_coeffs->bf;
    cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = b_face_cells[f_id];
      cs_real_t bfac[stride*stride];

      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++)
          bfac[stride*i+j] = thetap * b_visc[f_id] * cofbfp[f_id][i][j];
      }

      cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                     bfac,
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
 * and an extra diagonal part).
 *
 * template parameters:
 *   is_thermal        true for the temperature, otherwise false
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     ctx           dispatch context
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
 * \param[out]    ea            extra-diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <bool is_thermal>
static void
_coeffs_msr(const cs_mesh_t            *m,
            cs_dispatch_context        &ctx,
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
            cs_real_t         *restrict ea)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2c = ma->cell_cells;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* iconvp should always be 1 here; when it is 0, _sym_coeffs_scalar_msr
     should be called instead */

  cs_assert(iconvp);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t _da = rovsdt[c_id];

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];

      cs_real_t _i_massflux = i_massflux[f_id];

      if (c2f_sgn[cidx] > 0) {

        cs_real_t flui =  0.5*(_i_massflux - cs::abs(_i_massflux));

        // When solving the temperature, multiply convective part by Cp
        cs_real_t cpi = (is_thermal) ? xcpp[c_id] : 1.0;

        // Computation of extradiagonal terms
        cs_real_t xai = thetap*(cpi*flui -idiffp*i_visc[f_id]);

        // D_ii =  theta (m_ij)^+ - m_ij
        //      = -X_ij - (1-theta)*m_ij

        ea[cidx] = xai;
        _da -= xai + (1.-thetap) * cpi * _i_massflux;

      }
      else {

        cs_real_t fluj = -0.5*(_i_massflux + cs::abs(_i_massflux));

        // When solving the temperature, multiply convective part by Cp
        cs_real_t cpj = (is_thermal) ? xcpp[c2c[cidx]] : 1.0;

        // Computation of extradiagonal terms
        cs_real_t xaj = thetap*(cpj*fluj -idiffp*i_visc[f_id]);

        // D_jj = -theta (m_ij)^- + m_ij
        //      = -X_ji + (1-theta)*m_ij

        ea[cidx] = xaj;
        _da -= xaj - (1.-thetap) * cpj * _i_massflux;

      }
    }

    da[c_id] = _da;

  });

  /* Contribution of border faces to the diagonal */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_lnum_t ii = b_face_cells[f_id];

    cs_real_t _b_massflux = b_massflux[f_id];
    cs_real_t flui = 0.5*(_b_massflux - cs::abs(_b_massflux));

    // When solving the temperature, multiply convective part by Cp
    cs_real_t cpi = (is_thermal) ? xcpp[ii] : 1.0;

    // D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
    //      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
    //      = theta*(B -1)*(m_f)^- - (1-theta)*m_f

    cs_real_t bfac =          cpi *(  flui * thetap * (coefbp[f_id] - 1.)
                                    - (1. - thetap) * _b_massflux)
                   + idiffp * thetap * b_visc[f_id] * cofbfp[f_id];

    cs_dispatch_sum(&da[ii], bfac, b_sum_type);

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a strided field.
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part.
 *
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     ctx           dispatch context
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$, or null
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da             diagonal part of the matrix
 * \param[out]    ea             extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_coeffs_msr(const cs_mesh_t            *m,
            cs_dispatch_context        &ctx,
            int                         idiffp,
            double                      thetap,
            const cs_field_bc_coeffs_t *bc_coeffs,
            const cs_real_t             fimp[][stride][stride],
            const cs_real_t             i_massflux[],
            const cs_real_t             b_massflux[],
            const cs_real_t             i_visc[],
            const cs_real_t             b_visc[],
            cs_real_t                   da[][stride][stride],
            cs_real_t         *restrict ea)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  using b_t = cs_real_t[stride][stride];
  const b_t *coefbp = (const b_t *)bc_coeffs->b;
  const b_t *cofbfp = (const b_t *)bc_coeffs->bf;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    auto &_da = da[c_id];

    if (fimp != nullptr) {
      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++)
          _da[i][j] = fimp[c_id][i][j];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++)
          _da[i][j] = 0.;
      }
    }

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];

      cs_real_t _i_massflux = i_massflux[f_id];

      if (c2f_sgn[cidx] > 0) {

        cs_real_t flui =  0.5*(_i_massflux - cs::abs(_i_massflux));

        // Computation of extradiagonal terms
        cs_real_t xai = thetap*(flui -idiffp*i_visc[f_id]);

        // D_ii =  theta (m_ij)^+ - m_ij
        //      = -X_ij - (1-theta)*m_ij

        ea[cidx] = xai;

        // alternative: cs_real_t ifac = xaj + i_massflux;
        cs_real_t ifac = xai + (1.-thetap) * _i_massflux;
        for (cs_lnum_t i = 0; i < stride; i++) {
          _da[i][i] -= ifac;
        }

      }
      else {

        cs_real_t fluj = -0.5*(_i_massflux + cs::abs(_i_massflux));

        // When solving the temperature, multiply convective part by Cp

        // Computation of extradiagonal terms
        cs_real_t xaj = thetap*(fluj -idiffp*i_visc[f_id]);

        // D_jj = -theta (m_ij)^- + m_ij
        //      = -X_ji + (1-theta)*m_ij

        ea[cidx] = xaj;

        // alternative: cs_real_t jfac = xai - i_massflux[f_id];
        cs_real_t jfac = xaj - (1.-thetap) * _i_massflux;
        for (cs_lnum_t i = 0; i < stride; i++) {
          _da[i][i] -= jfac;
        }

      }
    }

  });

  /* Contribution of border faces to the diagonal */

  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    cs_real_t _b_massflux = b_massflux[f_id];
    cs_lnum_t ii = b_face_cells[f_id];
    cs_real_t flui = 0.5*(_b_massflux - cs::abs(_b_massflux));

    cs_real_t bfac[stride*stride];

    for (cs_lnum_t i = 0; i < stride; i++) {
      for (cs_lnum_t j = 0; j < stride; j++) {
        cs_real_t d_ij = ((i == j) ? 1. : 0.);
        /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */

        cs_real_t diag = d_ij
          * (thetap * flui + (1. - thetap) * b_massflux[f_id]);

        bfac[stride*i+j]
          = - diag + thetap * (  flui * coefbp[f_id][i][j]
                               + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]);
      }
    }

    cs_dispatch_sum<stride*stride>(reinterpret_cast<cs_real_t*>(da[ii]),
                                   bfac,
                                   b_sum_type);

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the convection-diffusion matrix for a strided field.
 * (non-symmetric matrix), where the diagonal terms are isotropic.
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     mq            pointer to mesh quantities structure
 * \param[in]     ctx           dispatch context
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$, or null
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da             diagonal part of the matrix
 * \param[out]    ea             extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

static void
_coeffs_msr_porous(const cs_mesh_t            *m,
                   const cs_mesh_quantities_t *mq,
                   cs_dispatch_context        &ctx,
                   int                         idiffp,
                   double                      thetap,
                   const cs_field_bc_coeffs_t *bc_coeffs,
                   const cs_real_t             fimp[][3][3],
                   const cs_real_t             i_massflux[],
                   const cs_real_t             b_massflux[],
                   const cs_real_t             i_visc[],
                   const cs_real_t             b_visc[],
                   cs_real_t                   da[][3][3],
                   cs_real_t                   ea[][3][3])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  using b_t = cs_real_t[3][3];
  const b_t *coefbp = (const b_t *)bc_coeffs->b;
  const b_t *cofbfp = (const b_t *)bc_coeffs->bf;

  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;

  // Without the porous model's face factors,
  // The extra-diagonal terms are isotropic, do there is no
  // need to use this version.

  assert(cs_glob_porous_model == 3);

  const cs_real_2_t *restrict i_f_face_factor = mq->i_f_face_factor;
  const cs_real_t *restrict b_f_face_factor = mq->b_f_face_factor;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    auto &_da = da[c_id];

    if (fimp != nullptr) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          _da[i][j] = fimp[c_id][i][j];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          _da[i][j] = 0.;
      }
    }

    /* Computation of extradiagonal terms */
    /*
     * X_ij = - theta f_j (m_ij)^-
     * X_ji = - theta f_i (m_ij)^+
     */

    /* Diagonal part:
     * the n(x)n term is multiplied by i_f_face_factor and (1 - n(x)n) by 1
     * XA_ij <= XA_ik n_k n_j (factor - 1) + XA_ij
     * XA_ij used to be diagonal: XA_ik n_k n_j = XA_ii n_i n_j*/

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];

      const cs_nreal_t *restrict n = i_face_u_normal[f_id];
      cs_real_t _i_massflux = i_massflux[f_id];

      if (c2f_sgn[cidx] > 0) {

        cs_real_t flui = 0.5*(_i_massflux - cs::abs(_i_massflux))
                         - idiffp*i_visc[f_id];
        cs_real_t ffm1 = (i_f_face_factor[f_id][1] - 1.);

        // Computation of extradiagonal terms
        for (cs_lnum_t k = 0; k < 3; k++) {
          for (cs_lnum_t l = 0; l < 3; l++) {
            cs_real_t d_kl = ((k == l) ? 1. : 0.);
            ea[cidx][k][l] = thetap * flui * (d_kl + ffm1*n[k]*n[l]);
            _da[k][l] += -d_kl*_i_massflux - (1.-thetap)*ea[cidx][k][l];
          }
        }

      }
      else {

        cs_real_t fluj = -0.5*(_i_massflux + cs::abs(_i_massflux))
                         - idiffp*i_visc[f_id];
        cs_real_t ffm1 = (i_f_face_factor[f_id][0] - 1.);

        // Computation of extradiagonal terms
        for (cs_lnum_t k = 0; k < 3; k++) {
          for (cs_lnum_t l = 0; l < 3; l++) {
            cs_real_t d_kl = ((k == l) ? 1. : 0.);
            ea[cidx][k][l] = thetap * fluj * (d_kl + ffm1*n[k]*n[l]);
            _da[k][l] += d_kl*_i_massflux - (1.-thetap)*ea[cidx][k][l];
          }
        }

      }
    }

  });

  /* Contribution of border faces to the diagonal */

  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

    const cs_nreal_t *n = b_face_u_normal[f_id];
    cs_lnum_t ii = b_face_cells[f_id];

    cs_real_t ffm1 = (b_f_face_factor[f_id] - 1.);

    cs_real_t flu[2] = {
      /* (m_ij)^+ */
      0.5 * (b_massflux[f_id] + cs::abs(b_massflux[f_id])),
      /* (m_ij)^- */
      0.5 * (b_massflux[f_id] - cs::abs(b_massflux[f_id]))
    };

    cs_real_t n_b_n
      = cs_math_3_33_3_dot_product(n, coefbp[f_id], n);
    cs_real_t n_bf_n
      = cs_math_3_33_3_dot_product(n, cofbfp[f_id], n);

    cs_real_t bfac[3*3];

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t d_ij = ((i == j) ? 1. : 0.);

        /* D = theta (m_f)^+.1 + theta B (m_f)^- - m_f.1
         * NB: stop here because the first two terms maybe scaled
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */

        bfac[3*i+j]
          =    thetap
             * (  d_ij * flu[0]
                + flu[1] * coefbp[f_id][i][j]
                + idiffp * b_visc[f_id] * cofbfp[f_id][i][j]
                + (flu[0] + flu[1] * n_b_n + idiffp * b_visc[f_id] * n_bf_n)
                  * ffm1*n[i]*n[j])
             - d_ij * b_massflux[f_id];

      }
    }

    cs_dispatch_sum<3*3>(reinterpret_cast<cs_real_t*>(da[ii]),
                         bfac,
                         b_sum_type);
  });

  ctx.wait();
}

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

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

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
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

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
    cs_real_t flui =  0.5*(_i_massflux - cs::abs(_i_massflux));
    cs_real_t fluj = -0.5*(_i_massflux + cs::abs(_i_massflux));

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
    cs_real_t flui = 0.5*(_b_massflux - cs::abs(_b_massflux));

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
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     ctx           Reference to dispatch context
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
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 *                              (initialized to \f$ \tens{f_s}^{imp} \f$)
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sym_matrix_strided(const cs_mesh_t            *m,
                    cs_dispatch_context        &ctx,
                    int                         idiffp,
                    cs_real_t                   thetap,
                    const cs_field_bc_coeffs_t *bc_coeffs_v,
                    const cs_real_t             i_visc[],
                    const cs_real_t             b_visc[],
                    cs_real_t        (*restrict da)[stride][stride],
                    cs_real_t        (*restrict xa))
{
  using b_t = cs_real_t[stride][stride];
  const b_t *cofbfp = (const b_t *)bc_coeffs_v->bf;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

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
 * \param[in]     ctx           Reference to dispatch context
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
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 *                              (initialized to \f$ \tens{f_s}^{imp} \f$)
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride, cs_lnum_t eb_size>
static void
_matrix_strided(const cs_mesh_t            *m,
                const cs_mesh_quantities_t *mq,
                cs_dispatch_context        &ctx,
                int                         iconvp,
                int                         idiffp,
                double                      thetap,
                const cs_field_bc_coeffs_t *bc_coeffs_v,
                const cs_real_t             i_massflux[],
                const cs_real_t             b_massflux[],
                const cs_real_t             i_visc[],
                const cs_real_t             b_visc[],
                cs_real_t        (*restrict da)[stride][stride],
                cs_real_t        (*restrict xa)[2])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  cs_real_2_t *i_f_face_factor;
  cs_real_t *b_f_face_factor;

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

  if ((stride == 3 && eb_size == 1) || stride == 6) {
    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      const cs_lnum_t _p = is_p*f_id;
      cs_real_t _i_massflux = i_massflux[f_id];
      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t flu[2] = {0.5*iconvp*(_i_massflux - cs::abs(_i_massflux)),
                         -0.5*iconvp*(_i_massflux + cs::abs(_i_massflux))};

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
      const cs_nreal_t *n = i_face_u_normal[f_id];
      const cs_real_t _i_massflux = i_massflux[f_id];
      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t flu[2] = {
         0.5 * iconvp * (_i_massflux - cs::abs(_i_massflux))
          - idiffp*i_visc[f_id],
        -0.5 * iconvp * (_i_massflux + cs::abs(_i_massflux))
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

      const cs_nreal_t *n = b_face_u_normal[f_id];
      const cs_lnum_t _p = is_p*f_id;
      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t flu[2] = {
        /* (m_ij)^+ */
        iconvp * 0.5 * (b_massflux[f_id] + cs::abs(b_massflux[f_id])),
        /* (m_ij)^- */
        iconvp * 0.5 * (b_massflux[f_id] - cs::abs(b_massflux[f_id]))
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
      cs_real_t flui = 0.5*(_b_massflux - cs::abs(_b_massflux));

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
 * \param[in]     ctx           Reference to dispatch context
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
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 *                              (initialized to \f$ \tens{f_s}^{imp} \f$)
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sym_matrix_anisotropic_diffusion_strided
  (const cs_mesh_t            *m,
   cs_dispatch_context        &ctx,
   int                         idiffp,
   double                      thetap,
   const cs_field_bc_coeffs_t *bc_coeffs_v,
   const cs_real_t             i_visc[][stride][stride],
   const cs_real_t             b_visc[],
   cs_real_t        (*restrict da)[stride][stride],
   cs_real_t        (*restrict xa)[stride][stride])
{
  using b_t = cs_real_t[stride][stride];
  const b_t *cofbfp = (const b_t *)bc_coeffs_v->bf;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Initialization */

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
 * \param[in]     ctx           Reference to dispatch context
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
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 *                              (initialized to \f$ \tens{f_s}^{imp} \f$)
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_matrix_anisotropic_diffusion_strided
(
 const cs_mesh_t            *m,
 const cs_mesh_quantities_t *mq,
 cs_dispatch_context        &ctx,
 int                         iconvp,
 int                         idiffp,
 double                      thetap,
 const cs_field_bc_coeffs_t *bc_coeffs_v,
 const cs_real_t             i_massflux[],
 const cs_real_t             b_massflux[],
 const cs_real_t             i_visc[][stride][stride],
 const cs_real_t             b_visc[],
 cs_real_t        (*restrict da)[stride][stride],
 cs_real_t        (*restrict xa)[2][stride][stride]
)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;

  cs_real_2_t *i_f_face_factor = nullptr;
  cs_real_t *b_f_face_factor = nullptr;
  cs_real_2_t *_i_f_face_factor = nullptr;
  cs_real_t *_b_f_face_factor = nullptr;
  int is_p = 0;

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

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    const cs_lnum_t _p = is_p*f_id;
    const cs_nreal_t *n = i_face_u_normal[f_id];
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    cs_real_t flu[2] = {
       iconvp * 0.5*(i_massflux[f_id] - cs::abs(i_massflux[f_id])),
      -iconvp * 0.5*(i_massflux[f_id] + cs::abs(i_massflux[f_id]))
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
      const cs_nreal_t *n = b_face_u_normal[f_id];
      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t flu[2] = {
         /* (m_ij)^+ */
         iconvp * 0.5 * (b_massflux[f_id] + cs::abs(b_massflux[f_id])),
         /* (m_ij)^- */
         iconvp * 0.5 * (b_massflux[f_id] - cs::abs(b_massflux[f_id]))
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
      cs_real_t flui = 0.5*(b_massflux[f_id] - cs::abs(b_massflux[f_id]));
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \param[in, out]  a             pointer to matrix structure
 * \param[in]       f             pointer to field, or null
 * \param[in]       iconvp        indicator
 *                                 - 1 advection
 *                                 - 0 otherwise
 * \param[in]       idiffp        indicator
 *                                 - 1 diffusion
 *                                 - 0 otherwise
 * \param[in]       ndircp        number of Dirichlet BCs
 * \param[in]       thetap        time scheme parameter
 * \param[in]       relaxp        relaxation coefficient (if < 1)
 * \param[in]       imucp         1 for temperature (with Cp), 0 otherwise
 * \param[in]       bc_coeffs     boundary condition structure
 * \param[in]       rovsdt        implicit terms (rho / dt)
 * \param[in]       i_massflux    mass flux at interior faces
 * \param[in]       b_massflux    mass flux at border faces
 * \param[in]       i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                 at interior faces for the matrix
 * \param[in]       b_visc        \f$ S_\fib \f$
 *                                 at border faces for the matrix
 * \param[in]       xcpp          Cp per cell, or null
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_compute_coeffs(cs_matrix_t                 *a,
                         const cs_field_t            *f,
                         int                          iconvp,
                         int                          idiffp,
                         int                          ndircp,
                         double                       thetap,
                         double                       relaxp,
                         int                          imucpp,
                         const cs_field_bc_coeffs_t  *bc_coeffs,
                         const cs_real_t              rovsdt[],
                         const cs_real_t              i_massflux[],
                         const cs_real_t              b_massflux[],
                         const cs_real_t              i_visc[],
                         const cs_real_t              b_visc[],
                         const cs_real_t              xcpp[])
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_matrix_assembler_t  *ma = cs_matrix_get_assembler(a);

  cs_matrix_type_t m_type = cs_matrix_get_type(a);
  cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);
  bool symmetric = (iconvp == 1) ? false : true;
  bool need_xa = cs_matrix_get_need_xa(a);

  /* Case with matrix assembler or non-msr format;
     we use the legacy assembly here, but should move to a row-based
     (i.e. gather) algorithm so as to be able to call assembler in a
     multithreaded manner (and later even a device-accelerated manner). */

  if (   ma != nullptr
      || m_type != CS_MATRIX_MSR
      || need_xa) {

    const cs_lnum_t n_edges = m->n_i_faces;
    const cs_lnum_2_t *edges = m->i_face_cells;

    const int isym = (iconvp == 1) ? 2 : 1;

    cs_real_t *da, *xa;
    CS_MALLOC_HD(da, m->n_cells_with_ghosts, cs_real_t, amode);
    CS_MALLOC_HD(xa, m->n_i_faces*(cs_lnum_t)isym, cs_real_t, amode);

    cs_matrix_wrapper(iconvp,
                      idiffp,
                      ndircp,
                      isym,
                      thetap,
                      imucpp,
                      bc_coeffs,
                      rovsdt,
                      i_massflux,  b_massflux,
                      i_visc, b_visc,
                      xcpp,
                      da, xa);

    /* For steady computations, the diagonal is relaxed */
    if (relaxp < 1) {
      cs_dispatch_context ctx;
      ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);
      cs_real_t rf = 1. / relaxp;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        da[i] *= rf;
      });
    }

    if (ma != nullptr) {

      cs_lnum_t s0 = 2;
      cs_lnum_t s1 = 1;

      if (isym == 1) {
        s0 = 1;
        s1 = 0;
      }

      cs_matrix_assembler_values_t *mav
        = cs_matrix_assembler_values_init(a, 1, 1);
      assert(n_cells == cs_matrix_get_n_rows(a));

      const cs_gnum_t *r_g_id = cs_matrix_get_block_row_g_id(a, m->halo);

      const cs_lnum_t block_size = 800;
      cs_gnum_t g_row_id[800];
      cs_gnum_t g_col_id[800];
      cs_real_t val[1600];

      /* Diagonal values */

      cs_matrix_assembler_values_add_g(mav, n_cells, r_g_id, r_g_id, da);

      /* Extradiagonal values based on internal faces */

      cs_lnum_t jj = 0;

      for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
        cs_lnum_t i0 = edges[ii][0];
        cs_lnum_t i1 = edges[ii][1];
        if (i0 < n_cells) {
          g_row_id[jj] = r_g_id[i0];
          g_col_id[jj] = r_g_id[i1];
          val[jj] = xa[ii*s0];
          jj++;
        }
        if (i1 < n_cells) {
          g_row_id[jj] = r_g_id[i1];
          g_col_id[jj] = r_g_id[i0];
          val[jj] = xa[ii*s0+s1];
          jj++;
        }
        if (jj >= block_size - 1) {
          cs_matrix_assembler_values_add_g(mav, jj,
                                           g_row_id, g_col_id, val);
          jj = 0;
        }
      }
      cs_matrix_assembler_values_add_g(mav, jj,
                                       g_row_id, g_col_id, val);
      jj = 0;

      /* Set extended contribution for domain coupling */

      if (f != nullptr) {
        int k_cpl = cs_field_key_id("coupling_entity");
        int coupling_id = cs_field_get_key_int(f, k_cpl);

        if (coupling_id > -1)
          cs_internal_coupling_matrix_add_values(f, 1, 1, r_g_id, mav);
      }

      /* Finalize assembly */

      cs_matrix_assembler_values_finalize(&mav);

    }

    else {

      /* As arrays are transferred, we assume the array type is
         neither "native", nor "dist", as the asscoaited formats
         do not currently have a "transfer coefficients" method.*/

      cs_assert(   m_type != CS_MATRIX_NATIVE
                && m_type != CS_MATRIX_DIST);

      cs_matrix_transfer_coefficients(a,
                                      symmetric,
                                      1,
                                      1,
                                      n_edges,
                                      edges,
                                      &da,
                                      &xa);

    }

    /* Free remaining local (non-transferred) arrays */

    CS_FREE(xa);
    CS_FREE(da);

    if (cs_glob_timer_kernels_flag > 0) {
      std::chrono::high_resolution_clock::time_point
        t_stop = std::chrono::high_resolution_clock::now();
      std::chrono::microseconds elapsed
        = std::chrono::duration_cast
        <std::chrono::microseconds>(t_stop - t_start);
      printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
             elapsed.count());
    }

    cs_matrix_default_set_tuned(a);

    return;
  }

  /* Common case: direct assigment of matrix coefficients
     ---------------------------------------------------- */

  cs_mesh_adjacencies_update_cell_i_faces();

  cs_dispatch_context ctx;
  ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);

  cs_real_t *da, *ea;
  cs_matrix_get_coefficients_msr_w(a,
                                   symmetric,
                                   1, 1,
                                   &da,
                                   &ea);

  /* Symmetric matrix */
  if (symmetric) {
    _sym_coeffs_msr(m,
                    ctx,
                    idiffp,
                    thetap,
                    bc_coeffs,
                    rovsdt,
                    i_visc,
                    b_visc,
                    da,
                    ea);
  }

  /* Non-symmetric matrix */
  else {
    if (imucpp == 0)
      _coeffs_msr<false>(m,
                         ctx,
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
                         ea);
    else
      _coeffs_msr<true>(m,
                        ctx,
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
                        ea);
  }

  /* For steady computations, the diagonal is relaxed */
  if (relaxp < 1) {
    cs_real_t rf = 1. / relaxp;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      da[i] *= rf;
    });
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    constexpr cs_real_t epsi_p1 = 1. + 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      da[c_id] *= epsi_p1;
    });
  }

  /* If a whole row of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      da[c_id] += (cs_real_t)c_disable_flag[c_id];

      /* FIXME: additional precaution if diagonal is still 0 due to
         all surrounding cells being disabled used here as it was
         present in cs_equation_iterative_solve, seems redundant
         with check above, as diagonal cannot be 0 if rovsdt > 0,
         (unless extradiagonal contributions bring it to 0 ?) */
      if (fabs(da[c_id]) < DBL_MIN)
        da[c_id] += 1;
    });
  }
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
           elapsed.count());
  }

  cs_matrix_default_set_tuned(a);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in, out]  a                   pointer to matrix structure
 * \param[in]       f                    pointer to field, or null
 * \param[in]       iconvp               indicator
 *                                         - 1 advection
 *                                         - 0 otherwise
 * \param[in]       idiffp               indicator
 *                                         - 1 diffusion
 *                                         - 0 otherwise
 * \param[in]       tensorial_diffusion  indicator
 * \param[in]       ndircp               number of Dirichlet BCs
 * \param[in]       thetap               time scheme parameter
 * \param[in]       relaxp               relaxation coefficient (if < 1)
 * \param[in]       eb_size              extra-diagonal block size
 *                                       (1 or 3 for stride 3, 1 for stride 6)
 * \param[in]       bc_coeffs            boundary conditions structure
 * \param[in]       fimp                 implicit terms, or null
 * \param[in]       i_massflux           mass flux at interior faces
 * \param[in]       b_massflux           mass flux at border faces
 * \param[in]       i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                 at interior faces for the matrix
 * \param[in]       b_visc        \f$ S_\fib \f$
 *                                 at boundary faces for the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_matrix_compute_coeffs
(
  cs_matrix_t                 *a,
  const cs_field_t            *f,
  int                          iconvp,
  int                          idiffp,
  int                          tensorial_diffusion,
  int                          ndircp,
  cs_lnum_t                    eb_size,
  double                       thetap,
  double                       relaxp,
  const cs_field_bc_coeffs_t  *bc_coeffs,
  const cs_real_t              fimp[][stride][stride],
  const cs_real_t              i_massflux[],
  const cs_real_t              b_massflux[],
  const cs_real_t              i_visc[],
  const cs_real_t              b_visc[]
)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_matrix_assembler_t  *ma = cs_matrix_get_assembler(a);

  cs_matrix_type_t m_type = cs_matrix_get_type(a);
  cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);
  bool symmetric = (iconvp == 1) ? false : true;
  bool need_xa = cs_matrix_get_need_xa(a);

  using b_t = cs_real_t[stride][stride];

  /* Case with matrix assembler or non-msr format;
     we use the legacy assembly here, but should move to a row-based
     (i.e. gather) algorithm so as to be able to call assembler in a
     multithreaded manner (and later even a device-accelerated manner). */

  if (   ma != nullptr
      || m_type != CS_MATRIX_MSR
      || tensorial_diffusion != 1
      || need_xa
      || (cs_glob_porous_model == 3 && stride == 3)) { // should be test on vel ?

    const cs_lnum_t n_edges = m->n_i_faces;
    const cs_lnum_2_t *edges = m->i_face_cells;
    const int isym = (iconvp == 1) ? 2 : 1;

    b_t *da;
    cs_real_t *xa;
    cs_lnum_t nnd = m->n_i_faces * (cs_lnum_t)isym * eb_size * eb_size;
    CS_MALLOC_HD(da, m->n_cells_with_ghosts, b_t, amode);
    CS_MALLOC_HD(xa, nnd, cs_real_t, amode);

    cs_matrix_wrapper<stride>(iconvp,
                      idiffp,
                      tensorial_diffusion,
                      ndircp,
                      isym,
                      eb_size,
                      thetap,
                      bc_coeffs,
                      fimp,
                      i_massflux,  b_massflux,
                      i_visc, b_visc,
                      da, xa);

    cs_real_t *da_p = reinterpret_cast<cs_real_t *>(da);

    /* For steady computations, the diagonal is relaxed */
    if (relaxp < 1) {
      cs_dispatch_context ctx;
      ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);
      cs_real_t rf = 1. / relaxp;
      ctx.parallel_for(n_cells*stride*stride, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        da_p[i] *= rf;
      });
    }

    if (ma != nullptr) {

      cs_lnum_t s0 = 2;
      cs_lnum_t s1 = 1;

      if (isym == 1) {
        s0 = 1;
        s1 = 0;
      }

      cs_matrix_assembler_values_t *mav
        = cs_matrix_assembler_values_init(a, stride, eb_size);
      assert(n_cells == cs_matrix_get_n_rows(a));

      const cs_gnum_t *r_g_id = cs_matrix_get_block_row_g_id(a, m->halo);

      const cs_lnum_t block_size = 800;
      cs_gnum_t g_row_id[800];
      cs_gnum_t g_col_id[800];
      cs_real_t val[1600];

      /* Diagonal values */

      cs_matrix_assembler_values_add_g(mav, n_cells, r_g_id, r_g_id, da_p);

      /* Extradiagonal values based on internal faces */

      cs_lnum_t jj = 0;

      for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
        cs_lnum_t i0 = edges[ii][0];
        cs_lnum_t i1 = edges[ii][1];
        if (i0 < n_cells) {
          g_row_id[jj] = r_g_id[i0];
          g_col_id[jj] = r_g_id[i1];
          val[jj] = xa[ii*s0];
          jj++;
        }
        if (i1 < n_cells) {
          g_row_id[jj] = r_g_id[i1];
          g_col_id[jj] = r_g_id[i0];
          val[jj] = xa[ii*s0+s1];
          jj++;
        }
        if (jj >= block_size - 1) {
          cs_matrix_assembler_values_add_g(mav, jj,
                                           g_row_id, g_col_id, val);
          jj = 0;
        }
      }
      cs_matrix_assembler_values_add_g(mav, jj,
                                       g_row_id, g_col_id, val);
      jj = 0;

      /* Set extended contribution for domain coupling */

      if (f != nullptr) {
        int k_cpl = cs_field_key_id("coupling_entity");
        int coupling_id = cs_field_get_key_int(f, k_cpl);

        if (coupling_id > -1)
          cs_internal_coupling_matrix_add_values(f,
                                                 stride, eb_size,
                                                 r_g_id, mav);
      }

      /* Finalize assembly */

      cs_matrix_assembler_values_finalize(&mav);

    }

    else {

      /* As arrays are transferred, we assume the array type is
         neither "native", nor "dist", as the associated formats
         do not currently have a "transfer coefficients" method.*/

      cs_assert(   m_type != CS_MATRIX_NATIVE
                && m_type != CS_MATRIX_DIST);

      cs_matrix_transfer_coefficients(a,
                                      symmetric,
                                      stride,
                                      eb_size,
                                      n_edges,
                                      edges,
                                      &(da_p),
                                      &xa);

      da = reinterpret_cast<b_t *>(da_p);  // should be nullptr
    }

    /* Free remaining local (non-transferred) arrays */

    CS_FREE(xa);
    CS_FREE(da);

    if (cs_glob_timer_kernels_flag > 0) {
      std::chrono::high_resolution_clock::time_point
        t_stop = std::chrono::high_resolution_clock::now();
      std::chrono::microseconds elapsed
        = std::chrono::duration_cast
        <std::chrono::microseconds>(t_stop - t_start);
      printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
             elapsed.count());
    }

    cs_matrix_default_set_tuned(a);

    return;
  }

  /* Common case: direct assigment of matrix coefficients
     ---------------------------------------------------- */

  cs_mesh_adjacencies_update_cell_i_faces();

  cs_dispatch_context ctx;
  ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);

  cs_real_t *da_p, *ea_p;
  cs_matrix_get_coefficients_msr_w(a,
                                   symmetric,
                                   stride, eb_size,
                                   &da_p,
                                   &ea_p);

  b_t *da = reinterpret_cast<b_t *>(da_p);

  /* scalar diffusion or right anisotropic diffusion */
  if (tensorial_diffusion == 1) {
    assert(eb_size == 1);

    /* Symmetric matrix */
    if (symmetric) {
      _sym_coeffs_msr<stride>(m,
                              ctx,
                              idiffp,
                              thetap,
                              bc_coeffs,
                              fimp,
                              i_visc,
                              b_visc,
                              da,
                              ea_p);
    }

    /* Non-symmetric matrix */
    else {
      if (eb_size == 1)
        _coeffs_msr<stride>(m,
                            ctx,
                            idiffp,
                            thetap,
                            bc_coeffs,
                            fimp,
                            i_massflux,
                            b_massflux,
                            i_visc,
                            b_visc,
                            da,
                            ea_p);
      else {
        if (stride == 3) {
          assert(cs_glob_porous_model == 3 && stride == 3); // test on vel ?
          b_t *restrict ea = reinterpret_cast<b_t *>(ea_p);
          _coeffs_msr_porous(m,
                             mq,
                             ctx,
                             idiffp,
                             thetap,
                             bc_coeffs,
                             (const cs_real_33_t *)fimp,
                             i_massflux,
                             b_massflux,
                             i_visc,
                             b_visc,
                             (cs_real_33_t *)da,
                             (cs_real_33_t *)ea);
        }
      }
    }

  }
  else { // tensorial_diffusion > 1
    cs_assert(0);  // Not handled direcly yet; handled above
  }

  /* For steady computations, the diagonal is relaxed */
  if (relaxp < 1) {
    cs_real_t rf = 1. / relaxp;
    ctx.parallel_for(n_cells*stride*stride, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      da_p[i] *= rf;
    });
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    constexpr cs_real_t epsi_p1 = 1. + 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        da[c_id][i][i] = epsi_p1 * da[c_id][i][i];
    });
  }

  /* If a whole row of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++) {
        da[c_id][i][i] += (cs_real_t)c_disable_flag[c_id];

        /* FIXME: additional precaution if diagonal is still 0 due to
           all surrounding cells being disabled used here as it was
           present in cs_equation_iterative_solve, seems redundant
           with check above, as diagonal cannot be 0 if rovsdt > 0,
           (unless extradiagonal contributions bring it to 0 ?) */
        if (fabs(da[c_id][i][i]) < DBL_MIN)
          da[c_id][i][i] += 1;
      }
    });
  }
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
           elapsed.count());
  }

  cs_matrix_default_set_tuned(a);
}

// Force instanciation

template void
cs_matrix_compute_coeffs(cs_matrix_t                 *a,
                         const cs_field_t            *f,
                         int                          iconvp,
                         int                          idiffp,
                         int                          tensorial_diffusion,
                         int                          ndircp,
                         cs_lnum_t                    eb_size,
                         double                       thetap,
                         double                       relaxp,
                         const cs_field_bc_coeffs_t  *bc_coeffs,
                         const cs_real_t              fimp[][3][3],
                         const cs_real_t              i_massflux[],
                         const cs_real_t              b_massflux[],
                         const cs_real_t              i_visc[],
                         const cs_real_t              b_visc[]);

template void
cs_matrix_compute_coeffs(cs_matrix_t                 *a,
                         const cs_field_t            *f,
                         int                          iconvp,
                         int                          idiffp,
                         int                          tensorial_diffusion,
                         int                          ndircp,
                         cs_lnum_t                    eb_size,
                         double                       thetap,
                         double                       relaxp,
                         const cs_field_bc_coeffs_t  *bc_coeffs,
                         const cs_real_t              fimp[][6][6],
                         const cs_real_t              i_massflux[],
                         const cs_real_t              b_massflux[],
                         const cs_real_t              i_visc[],
                         const cs_real_t              b_visc[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper(int                         iconvp,
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
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

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

      /* FIXME: additional precaution if diagonal is still 0 due to
         all surrounding cells being disabled used here as it was
         present in cs_equation_iterative_solve, seems redundant
         with check above, as diagonal cannot be 0 if rovsdt > 0,
         (unless extradiagonal contributions bring it to 0 ?) */
      if (fabs(da[c_id]) < DBL_MIN)
        da[c_id] += 1;
    });
  }
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
           elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the matrix for a vector or tensor field
 *
 * The advection (if present) is upwind.
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal part (stride*stride blocks)
 * and an extra diagonal part.
 *
 * \tparam  stride 3 for vectors, 6 for tensors
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$, or null
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_matrix_wrapper(int                          iconvp,
                  int                          idiffp,
                  int                          tensorial_diffusion,
                  int                          ndircp,
                  int                          isym,
                  cs_lnum_t                    eb_size,
                  double                       thetap,
                  const cs_field_bc_coeffs_t  *bc_coeffs_v,
                  const cs_real_t              fimp[][stride][stride],
                  const cs_real_t              i_massflux[],
                  const cs_real_t              b_massflux[],
                  const cs_real_t              i_visc[],
                  const cs_real_t              b_visc[],
                  cs_real_t                    da[][stride][stride],
                  cs_real_t                    xa[])
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  if (fimp != nullptr) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        for (cs_lnum_t j = 0; j < stride; j++)
          da[c_id][i][j] = fimp[c_id][i][j];
    });
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        for (cs_lnum_t j = 0; j < stride; j++)
          da[c_id][i][j] = 0.;
    });
  }

  /* scalar diffusion or right anisotropic diffusion */
  if (tensorial_diffusion == 1) {
    /* Symmetric matrix */
    if (isym == 1) {
      assert(eb_size == 1);
      _sym_matrix_strided<stride>(m,
                                  ctx,
                                  idiffp,
                                  thetap,
                                  bc_coeffs_v,
                                  i_visc,
                                  b_visc,
                                  da,
                                  xa);

    /* Non-symmetric matrix */
    }
    else {
      if (eb_size == 1)
        _matrix_strided<stride,1>(m,
                                  mq,
                                  ctx,
                                  iconvp,
                                  idiffp,
                                  thetap,
                                  bc_coeffs_v,
                                  i_massflux,
                                  b_massflux,
                                  i_visc,
                                  b_visc,
                                  da,
                                  (cs_real_2_t*) xa);
      else {
        if (stride == 3)
          _matrix_strided<stride,stride>(m,
                                         mq,
                                         ctx,
                                         iconvp,
                                         idiffp,
                                         thetap,
                                         bc_coeffs_v,
                                         i_massflux,
                                         b_massflux,
                                         i_visc,
                                         b_visc,
                                         da,
                                         (cs_real_2_t*) xa);
        else
          cs_assert(0); // Only for vectors in this case
      }
    }
  }
  /* left tensor diffusion */
  else {

    using b_t = cs_real_t[stride][stride];

    /* Symmetric matrix */
    if (isym == 1) {
      _sym_matrix_anisotropic_diffusion_strided<stride>
        (m,
         ctx,
         idiffp,
         thetap,
         bc_coeffs_v,
         (const b_t *)i_visc,
         b_visc,
         da,
         (b_t *) xa);

    /* Non-symmetric matrix */
    }
    else {
      using xa_t = cs_real_t[2][stride][stride];

      _matrix_anisotropic_diffusion_strided<stride>
        (m,
         mq,
         ctx,
         iconvp,
         idiffp,
         thetap,
         bc_coeffs_v,
         i_massflux,
         b_massflux,
         (const b_t *)i_visc,
         b_visc,
         da,
         (xa_t *) xa);
    }

  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if (ndircp <= 0) {
    const cs_real_t epsi = 1.e-7;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        da[c_id][i][i] = (1. + epsi) * da[c_id][i][i];
    });
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++) {
        da[c_id][i][i] += (cs_real_t)c_disable_flag[c_id];

        /* FIXME: additional precaution if diagonal is still 0 due to
           all surrounding cells being disabled used here as it was
           present in cs_equation_iterative_solve, seems redundant
           with check above, as diagonal cannot be 0 if rovsdt > 0,
           (unless extradiagonal contributions bring it to 0 ?) */
        if (fabs(da[c_id][i][i]) < DBL_MIN)
          da[c_id][i][i] += 1;
      }
    });
  }
  ctx.wait();

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s = %ld\n", cs_glob_rank_id, __func__,
           elapsed.count());
  }
}

// Force instanciation

template void
cs_matrix_wrapper(int                          iconvp,
                  int                          idiffp,
                  int                          tensorial_diffusion,
                  int                          ndircp,
                  int                          isym,
                  cs_lnum_t                    eb_size,
                  double                       thetap,
                  const cs_field_bc_coeffs_t  *bc_coeffs_v,
                  const cs_real_t              fimp[][3][3],
                  const cs_real_t              i_massflux[],
                  const cs_real_t              b_massflux[],
                  const cs_real_t              i_visc[],
                  const cs_real_t              b_visc[],
                  cs_real_t                    da[][3][3],
                  cs_real_t                    xa[]);

template void
cs_matrix_wrapper(int                          iconvp,
                  int                          idiffp,
                  int                          tensorial_diffusion,
                  int                          ndircp,
                  int                          isym,
                  cs_lnum_t                    eb_size,
                  double                       thetap,
                  const cs_field_bc_coeffs_t  *bc_coeffs_v,
                  const cs_real_t              fimp[][6][6],
                  const cs_real_t              i_massflux[],
                  const cs_real_t              b_massflux[],
                  const cs_real_t              i_visc[],
                  const cs_real_t              b_visc[],
                  cs_real_t                    da[][6][6],
                  cs_real_t                    xa[]);

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  cs_dispatch_context ctx;

  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* With convection
     --------------- */

  if (iconvp == 1) {

    /* Contribution of the extra-diagonal terms to the diagonal */

    if (isym == 2) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t _da = 0.;

        /* Loop on interior faces */
        const cs_lnum_t s_id_i = c2c_idx[c_id];
        const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

        for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
          const cs_lnum_t f_id = cell_i_faces[cidx];

          cs_real_t _i_massflux = i_massflux[f_id];

          if (c2f_sgn[cidx] > 0) {
            cs_real_t fluj = -0.5*(_i_massflux + cs::abs(_i_massflux));
            _da -= fluj -idiffp*i_visc[f_id];
          }
          else {
            cs_real_t flui =  0.5*(_i_massflux - cs::abs(_i_massflux));
            _da -= flui -idiffp*i_visc[f_id];
          }
        }

        da[c_id] = _da;
      });

    }
    else { // if (isym == 1)

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t _da = 0.;

        /* Loop on interior faces */
        const cs_lnum_t s_id_i = c2c_idx[c_id];
        const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

        for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
          const cs_lnum_t f_id = cell_i_faces[cidx];
          cs_real_t _i_massflux = i_massflux[f_id];
          cs_real_t flui =  0.5*(_i_massflux - cs::abs(_i_massflux));
          _da -= flui -idiffp*i_visc[f_id];
        }

        da[c_id] = _da;
      });

    }

    /* Contribution of boundary faces to the diagonal */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {

      cs_lnum_t ii = b_face_cells[f_id];

      cs_real_t _b_massflux = b_massflux[f_id];
      cs_real_t flui =  0.5*(_b_massflux - cs::abs(_b_massflux));
      cs_real_t fluj = -0.5*(_b_massflux + cs::abs(_b_massflux));

      cs_real_t bfac =   (-fluj + flui*coefbp[f_id])
                       + idiffp*b_visc[f_id]*cofbfp[f_id];

      cs_dispatch_sum(&da[ii], bfac, b_sum_type);

    });

  }

  /* Without convection
     ------------------ */

  else {  // if iconvp == 0

    /* Contribution of the extra-diagonal terms to the diagonal */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t _da = 0.;

      if (idiffp > 0) {
        /* Loop on interior faces */
        const cs_lnum_t s_id_i = c2c_idx[c_id];
        const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

        for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
          const cs_lnum_t f_id = cell_i_faces[cidx];
          _da -= -idiffp*i_visc[f_id];
        }
      }

      da[c_id] = _da;
    });

    /* Contribution of boundary faces to the diagonal */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      cs_lnum_t ii = b_face_cells[f_id];
      cs_real_t bfac = idiffp*b_visc[f_id]*cofbfp[f_id];

      cs_dispatch_sum(&da[ii], bfac, b_sum_type);
    });

  }

  ctx.wait();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
