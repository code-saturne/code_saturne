/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_blas.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_param_sles_setup.h"
#include "alge/cs_saddle_solver_setup.h"
#include "base/cs_array.h"
#include "base/cs_fp_exception.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_timer.h"
#include "cdo/cs_cdo_blas.h"
#include "cdo/cs_cdo_solve.h"
#include "cdo/cs_equation.h"
#include "cdo/cs_saddle_system.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cdo/cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdofb_monolithic_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic_sles.cpp
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

/* Redefined the name of functions from cs_param_sles to get shorter names */

#define _petsc_cmd  cs_param_sles_setup_petsc_cmd

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_mesh_t  *cs_shared_mesh;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define several structures such as the cs_range_set_t,
 *        cs_interface_set_t, cs_matrix_assembler_t and cs_matrix_structure_t
 *        when the full saddle-point matrix is built into one block.
 *        A variant, activated with add_pressure_diag, is available in order to
 *        enforce the pressure.
 *
 * \param[in, out] block              pointer to a block structure
 * \param[in]      add_pressure_diag  true or false (pressure diagonal block)
 */
/*----------------------------------------------------------------------------*/

static void
_build_shared_structures_full_system(cs_cdo_system_block_t  *block,
                                     bool                    add_pressure_diag)
{
  /* Compute the range set for an array of size 3*n_faces + n_cells
   * velocity is attached to faces (one for each component) and pressure
   * to cells
   *
   * Storage for the global numbering: Vel_X | Vel_Y | Vel_Z | Pressure */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  size = 3*n_faces + m->n_cells;

  assert(block->type == CS_CDO_SYSTEM_BLOCK_EXTERN);
  cs_cdo_system_xblock_t *xb = (cs_cdo_system_xblock_t *)block->block_pointer;

  /* 1. Build the interface set and the range set structures */

  cs_interface_set_t *ifs = connect->face_ifs;

  if (ifs != nullptr)
    xb->interface_set = cs_interface_set_dup_blocks(ifs, n_faces, 3);

  xb->range_set = cs_range_set_create(xb->interface_set,
                                      nullptr, /* halo */
                                      size,
                                      false, /* TODO: add balance option */
                                      1,     /* tr_ignore */
                                      0);    /* g_id_base */

  /* 2. Build the matrix assembler structure */

  const cs_adjacency_t  *f2f = connect->f2f;
  const cs_adjacency_t  *f2c = connect->f2c;

  /* The second parameter is set to "true" meaning that the diagonal is stored
   * separately --> MSR storage
   * Create the matrix assembler structure */

  xb->matrix_assembler = cs_matrix_assembler_create(xb->range_set->l_range,
                                                    true);

  /* First loop to count max size of the buffer used to fill the matrix
   * structure. +1 to take into account the diagonal term. */

  int  max_sten = 0;
  for (cs_lnum_t f = 0; f < n_faces; f++) {
    int  sten
      = 9*(f2f->idx[f+1]-f2f->idx[f] + 1) + 6*(f2c->idx[f+1]-f2c->idx[f]);
    max_sten = cs::max(max_sten, sten);
  }

  cs_gnum_t *grows = nullptr, *gcols = nullptr;
  CS_MALLOC(grows, max_sten, cs_gnum_t);
  CS_MALLOC(gcols, max_sten, cs_gnum_t);

  /*
   *   | A_xx  |       |       | Bt_x  |
   *   |-------|-------|-------|-------|
   *   |       | A_yy  |       | Bt_y  |
   *   |-------|-------|-------|-------|
   *   |       |       | A_zz  | Bt_z  |
   *   |-------|-------|-------|-------|
   *   | B_x   | B_y   | B_z   |  0    |
   *
   *  Each block A_.. is n_faces * n_faces
   *  Each block B_.  is n_cells * n_faces
   */

  /* Only on faces (B_x is build in the same time as Bt_x for pressure DoFs) */

  for (cs_lnum_t frow_id = 0; frow_id < n_faces; frow_id++) {

    const cs_lnum_t  start = f2f->idx[frow_id];
    const cs_lnum_t  end = f2f->idx[frow_id+1];

    /* A face-face entry corresponds to 3x3 block + the diagonal which is not
       taken into account in the face --> face connectivity. The B and Bt
       operators have the same sparsity. 3x1 entries for the c2f
       connectivity. This is multiply by two since one considers B and Bt. */

    int  n_entries = (end-start + 1)*9
                   + 6*(f2c->idx[frow_id+1]-f2c->idx[frow_id]);

    const cs_gnum_t  grow_ids[3]
      = { xb->range_set->g_id[frow_id],               /* x-component */
          xb->range_set->g_id[frow_id +   n_faces],   /* y-component */
          xb->range_set->g_id[frow_id + 2*n_faces] }; /* z-component */

    int shift = 0;

    /* Diagonal term is excluded in this connectivity. Add it "manually" */

    for (int i = 0; i < 3; i++) {
      const cs_gnum_t  grow_id = grow_ids[i];
      for (int j = 0; j < 3; j++) {
        grows[shift] = grow_id;
        gcols[shift] = grow_ids[j];
        shift++;
      }
    }

    /* Extra diagonal couples */

    for (cs_lnum_t idx = start; idx < end; idx++) {

      const cs_lnum_t  fcol_id = f2f->ids[idx];
      const cs_gnum_t  gcol_ids[3]
        = { xb->range_set->g_id[fcol_id],               /* x-component */
            xb->range_set->g_id[fcol_id + n_faces],     /* y-component */
            xb->range_set->g_id[fcol_id + 2*n_faces] }; /* z-component */

      for (int i = 0; i < 3; i++) {
        const cs_gnum_t  grow_id = grow_ids[i];
        for (int j = 0; j < 3; j++) {
          grows[shift] = grow_id;
          gcols[shift] = gcol_ids[j];
          shift++;
        }
      }

    } /* Loop on extra-diag. entries */

    /* Loop on pressure-related  entries */

    for (cs_lnum_t idx = f2c->idx[frow_id]; idx < f2c->idx[frow_id+1]; idx++) {

      const cs_lnum_t  ccol_id = f2c->ids[idx];
      const cs_gnum_t  gcol_id = xb->range_set->g_id[3*n_faces + ccol_id];

      for (int i = 0; i < 3; i++) { /* x,y,z-component */

        grows[shift] = grow_ids[i];
        gcols[shift] = gcol_id;
        shift++;

        /* Its transposed B_x, B_y, B_z */

        grows[shift] = gcol_id;
        gcols[shift] = grow_ids[i];
        shift++;

      }

    } /* Loop on pressure related DoFs */

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  if (add_pressure_diag) {

    const cs_gnum_t  *cell_g_ids = xb->range_set->g_id + 3*n_faces;

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  m->n_cells, cell_g_ids, cell_g_ids);

  }

  /* 3. Build the matrix structure */

  cs_matrix_assembler_compute(xb->matrix_assembler);

  xb->matrix_structure
    = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR,
                                                xb->matrix_assembler);

  /* Free temporary buffers */

  CS_FREE(grows);
  CS_FREE(gcols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a diagonal scaled mass matrix
 *
 * \param[in] nsp    set of parameters for the Navier-Stokes system
 * \param[in] pty    property related to the (2,2) block
 * \param[in] gamma  scaling of the augmentation term
 *
 * \return a pointer to a newly allocated array
 */
/*----------------------------------------------------------------------------*/

static cs_real_t *
_get_scaled_diag_m22(const cs_navsto_param_t  *nsp,
                     const cs_property_t      *pty,
                     double                    gamma)
{
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_cdo_quantities_t  *cdoq = cs_shared_quant;
  const cs_lnum_t  n_cells = cdoq->n_cells;

  cs_real_t *m22_mass_diag = nullptr;
  CS_MALLOC(m22_mass_diag, n_cells, cs_real_t);

  /* Compute scaling coefficients */

  if (nsp->turbulence->model->model == CS_TURB_NONE) {

    const cs_real_t  scaling = gamma + pty->ref_value;
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      m22_mass_diag[i] = scaling/cdoq->cell_vol[i];

  }
  else {

    cs_property_eval_at_cells(ts->t_cur, pty, m22_mass_diag);

    if (fabs(gamma) > 0) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        m22_mass_diag[i] += gamma;
        m22_mass_diag[i] /= cdoq->cell_vol[i];
      }
    }
    else {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        m22_mass_diag[i] /= cdoq->cell_vol[i];
    }

  }

  return m22_mass_diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and define the matrix aiming at approximating the Schur
 *        complement. This approximation is based on a diagonal approximation
 *        of the inverse of the (1,1)-matrix
 *
 * \param[in]  schur_mat_class  type of Schur matrix to define
 * \param[in]  m11_inv_approx   diagonal approx. of the inverse of m11
 * \param[in]  m21_adj          adjacency related to the m21
 * \param[in]  m21_val          array associated to the m21 matrix (unassembled)
 * \param[out] p_diag_smat      diagonal coeffs for the Schur matrix
 * \param[out] p_xtra_smat      extra-diagonal coeffs for the Schur matrix
 *
 * \return a pointer to the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_schur_approx_diag_inv_m11(cs_param_solver_class_t  mat_class,
                           const cs_real_t         *m11_inv_approx,
                           const cs_adjacency_t    *m21_adj,
                           const cs_real_t         *m21_val,
                           cs_real_t              **p_diag_smat,
                           cs_real_t              **p_xtra_smat)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_2_t *restrict i_face_cells = mesh->i_face_cells;

  /* Native format for the Schur approximation matrix */

  cs_real_t *diag_smat = nullptr;
  cs_real_t *xtra_smat = nullptr;

  CS_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  CS_MALLOC(xtra_smat, 2 * n_i_faces, cs_real_t);

  cs_array_real_fill_zero(n_cells_ext, diag_smat);
  cs_array_real_fill_zero(2*n_i_faces, xtra_smat);

  const cs_lnum_t n_elts = m21_adj->n_elts;

  /* Add diagonal contributions from interior and boundary faces */

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) {
    for (cs_lnum_t j = m21_adj->idx[i]; j < m21_adj->idx[i+1]; j++) {

      const cs_lnum_t f_id = m21_adj->ids[j];
      const cs_real_t *m11_inv_ff = m11_inv_approx + 3*f_id;

      const cs_real_t *m21_vals = m21_val + 3*j;

      cs_real_t _m11_inv_m12[3];
      for (short int k = 0; k < 3; k++)
        _m11_inv_m12[k] = m11_inv_ff[k] * m21_vals[k];

      cs_real_t _m21_m11_inv_m21 =
        cs_math_3_dot_product(m21_vals, _m11_inv_m12);

      diag_smat[i] += _m21_m11_inv_m21;
    }
  }

  /* Add extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *m11_inv_ff = m11_inv_approx + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += m11_inv_ff[k] * nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

  } /* Loop on interior faces */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  cs_matrix_t *smat = nullptr;

  if (mat_class == CS_PARAM_SOLVER_CLASS_HYPRE)
    smat = cs_matrix_external("HYPRE_ParCSR",
                              false, /* symmetry */
                              1, 1);
  else
    smat = cs_matrix_msr();

  cs_matrix_set_coefficients(smat, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Associate some mesh quantities to the matrix (useful for grid
     coarsening) */

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

  cs_matrix_set_mesh_association(smat,
                                 ma->cell_cells_idx,
                                 ma->cell_i_faces,
                                 ma->cell_i_faces_sgn,
                                 (const cs_real_3_t *)quant->cell_centers,
                                 quant->cell_vol,
                                 quant->i_face_u_normal,
                                 quant->i_face_surf);

  /* Return the associated matrix and set the pointers to return
   * Return arrays (to be freed when the algorithm is converged) */

  *p_diag_smat = diag_smat;
  *p_xtra_smat = xtra_smat;

  return smat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context structure associated to a GKB algorithm
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] ctx     additional memebers specific to the GKB algo.
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_init_context(cs_saddle_solver_t              *solver,
                  cs_saddle_solver_context_gkb_t  *ctx)
{
  assert(ctx != nullptr && solver != nullptr);

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  assert(n1_dofs == 3*quant->n_faces);
  assert(n2_dofs == quant->n_cells);

  /* Orthogonalization coefficients */

  ctx->alpha = ctx->beta = ctx->zeta = 0.;

  if (solver->do_setup == false)
    return;

  /* Buffers of size n2_dofs */

  CS_MALLOC(ctx->q, n2_dofs, cs_real_t);
  CS_MALLOC(ctx->d, n2_dofs, cs_real_t);
  CS_MALLOC(ctx->m21v, n2_dofs, cs_real_t);
  CS_MALLOC(ctx->inv_m22, n2_dofs, cs_real_t);

  ctx->m22 = quant->cell_vol;   /* shared pointer */
  for (cs_lnum_t i = 0; i < n2_dofs; i++)
    ctx->inv_m22[i] = 1./quant->cell_vol[i];

  /* Buffers of size n1_dofs */

  CS_MALLOC(ctx->m12q, n1_dofs, cs_real_t);
  CS_MALLOC(ctx->x1_tilda, n1_dofs, cs_real_t);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  max_b11_size = cs::max(cs_matrix_get_n_columns(m11), n1_dofs);

  CS_MALLOC(ctx->w, max_b11_size, cs_real_t);
  CS_MALLOC(ctx->v, max_b11_size, cs_real_t);

  /* Rk: rhs_tilda stores quantities in space X1 and X2 alternatively */

  CS_MALLOC(ctx->rhs_tilda, cs::max(n1_dofs, n2_dofs), cs_real_t);

  /* Convergence members (energy norm estimation) */

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_gkb_t *ctxp
    = (const cs_param_saddle_context_gkb_t *)saddlep->context;
  const double  gamma = ctxp->augmentation_scaling;
  const int  tt = ctxp->truncation_threshold;

  if (gamma < 1)
    ctx->zeta_size = tt + 1;
  else if (gamma < 10)
    ctx->zeta_size = tt;
  else if (gamma < 100)
    ctx->zeta_size = cs::max(1, tt - 1);
  else if (gamma < 1e3)
    ctx->zeta_size = cs::max(1, tt - 2);
  else if (gamma < 1e4)
    ctx->zeta_size = cs::max(1, tt - 3);
  else
    ctx->zeta_size = cs::max(1, tt - 4);

  CS_MALLOC(ctx->zeta_array, ctx->zeta_size, cs_real_t);
  for (int i = 0; i < ctx->zeta_size; i++)
    ctx->zeta_array[i] = 0.;

  ctx->zeta_square_sum = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context structure associated to an ALU algorithm
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] ctx     additional memebers specific to ALU algo.
 */
/*----------------------------------------------------------------------------*/

static void
_alu_init_context(const cs_navsto_param_t         *nsp,
                  cs_saddle_solver_t              *solver,
                  cs_saddle_solver_context_alu_t  *ctx)
{
  CS_NO_WARN_IF_UNUSED(nsp);
  assert(ctx != nullptr);

  if (solver->do_setup == false)
    return;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  assert(solver->n2_scatter_dofs == quant->n_cells);
  assert(solver->n1_scatter_dofs == 3*quant->n_faces);

  /* Buffers of size n2_scatter_dofs */

  CS_MALLOC(ctx->inv_m22, solver->n2_scatter_dofs, cs_real_t);

# pragma omp parallel for if (solver->n2_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < solver->n2_scatter_dofs; i2++)
    ctx->inv_m22[i2] = 1./quant->cell_vol[i2];

  CS_MALLOC(ctx->res2, solver->n2_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->m21x1, solver->n2_scatter_dofs, cs_real_t);

  /* Buffers of size n1_scatter_dofs */

  CS_MALLOC(ctx->b1_tilda, solver->n1_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->rhs, solver->n1_scatter_dofs, cs_real_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context structure associated to a Uzawa-CG algorithm
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] ctx     additional members specific to the Uzawa-CG algo.
 */
/*----------------------------------------------------------------------------*/

static void
_uzawa_cg_init_context(const cs_navsto_param_t              *nsp,
                       cs_saddle_solver_t                   *solver,
                       cs_saddle_solver_context_uzawa_cg_t  *ctx)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  assert(ctx != nullptr);
  assert(solver->n2_scatter_dofs == cs_shared_quant->n_cells);
  assert(solver->n1_scatter_dofs == 3 * cs_shared_quant->n_faces);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

  ctx->m11 = cs_cdo_system_get_matrix(sh, 0);
  ctx->b11_max_size = cs::max(cs_matrix_get_n_columns(ctx->m11),
                              solver->n1_scatter_dofs);

  /* Buffers of size n1_scatter_dofs */

  CS_MALLOC(ctx->b1_tilda, solver->n1_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->rhs, solver->n1_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->dzk, solver->n1_scatter_dofs, cs_real_t);

  /* Buffers of size n2_scatter_dofs */

  CS_MALLOC(ctx->res2, solver->n2_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->m21x1, solver->n2_scatter_dofs, cs_real_t);

  /* Since gk is used as a variable in a cell system, one has to take into
     account extra-space for synchronization */

  cs_lnum_t  size = solver->n2_scatter_dofs;
  if (cs_glob_n_ranks > 1)
    size = cs::max(size, connect->n_cells_with_ghosts);
  CS_MALLOC(ctx->gk, size, cs_real_t);

  // No scaling for an augmentation term with Uzawa-CG

  ctx->inv_m22 = _get_scaled_diag_m22(nsp, ctx->pty_22, 0.);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context structure associated to a SIMPLE algorithm
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] ctx     additional members specific to the Uzawa-CG algo.
 */
/*----------------------------------------------------------------------------*/

static void
_simple_init_context(cs_saddle_solver_t                *solver,
                     cs_saddle_solver_context_simple_t *ctx)
{
  assert(ctx != nullptr);
  assert(solver->n2_scatter_dofs == cs_shared_quant->n_cells);
  assert(solver->n1_scatter_dofs == 3 * cs_shared_quant->n_faces);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

  ctx->m11 = cs_cdo_system_get_matrix(sh, 0);
  ctx->b11_max_size = cs::max(cs_matrix_get_n_columns(ctx->m11),
                              solver->n1_scatter_dofs);

  /* Buffers of size n1_scatter_dofs */

  CS_MALLOC(ctx->b1_tilda, solver->n1_scatter_dofs, cs_real_t);
  CS_MALLOC(ctx->rhs, solver->n1_scatter_dofs, cs_real_t);

  /* Buffers of size n2_scatter_dofs */

  CS_MALLOC(ctx->m21x1, solver->n2_scatter_dofs, cs_real_t);

  ctx->inv_m22 = nullptr;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set pointers to shared structures
 *
 * \param[in] mesh     pointer to the mesh structure
 * \param[in] connect  pointer to additional CDO connectivities
 * \param[in] quant    pointer to additional CDO mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_sharing(const cs_mesh_t            *mesh,
                                      const cs_cdo_connect_t     *connect,
                                      const cs_cdo_quantities_t  *quant)
{
  /* Assign static const pointers */

  cs_shared_mesh = mesh;
  cs_shared_connect = connect;
  cs_shared_quant = quant;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the system helper for a CDO-Fb scheme solving the
 *        Navier-Stokes equation using a monolithic approach for the
 *        velocity-pressure coupling
 *
 * \param[in]      nsp      Navier-Stokes paremeters
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] sc       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_system_helper(const cs_navsto_param_t  *nsp,
                                            const cs_param_saddle_t  *saddlep,
                                            cs_cdofb_monolithic_t    *sc)
{
  const cs_cdo_quantities_t  *cdoq = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_param_sles_t  *b11_slesp = saddlep->block11_sles_param;

  /* Define the system helper */

  cs_cdo_system_helper_t *sh = nullptr;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_GKB:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
      cs_lnum_t  block_sizes[2];
      block_sizes[0] = 3*cdoq->n_faces, block_sizes[1] = cdoq->n_cells;

      /* Create the system helper */

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       2,
                                       block_sizes,
                                       2);

      /* Add a first block for the (1,1)-block and then define the underpinning
         structures */

      /* Choose the right class of matrix to avoid copy.
       * The way to perform the assembly may change if an external librairy is
       * used for solving the linear system.
       */

      cs_cdo_system_matrix_class_t
        matclass = cs_cdo_system_which_matrix_class(b11_slesp->solver_class);

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_dblock(sh, 0,                /* block id */
                                 matclass,
                                 cs_flag_primal_face , /* location */
                                 cdoq->n_faces,        /* n_elements */
                                 3,                    /* stride */
                                 true,                 /* interlaced */
                                 true);                /* unrolled */

      cs_cdo_system_build_block(sh, 0); /* build structures */

      /* Add a second block for the (1,0) and (0,1) blocks and then define the
         underpinning structures. The (0,1) block needs to be transposed before
         using it */

      cs_cdo_system_block_t  *bdiv =
        cs_cdo_system_add_ublock(sh, 1,               /* block_id */
                                 connect->c2f,        /* adjacency */
                                 cs_flag_primal_face, /* column location */
                                 cdoq->n_faces,       /* n_elements */
                                 3,                   /* stride */
                                 true);               /* interlaced */

      cs_cdo_system_dblock_t *a_db = (cs_cdo_system_dblock_t *)a->block_pointer;
      cs_cdo_system_ublock_t *b_ub
        = (cs_cdo_system_ublock_t *)bdiv->block_pointer;

      /* Define the bdiv block by hand */

      b_ub->adjacency = connect->c2f;            /* shared pointer */
      b_ub->values = sc->block21_op;             /* shared pointer */
      assert(b_ub->values != nullptr);
      b_ub->shared_structures = true;
      b_ub->range_set = a_db->range_set;         /* shared pointer */
      b_ub->interface_set = a_db->interface_set; /* shared pointer */
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    {
      cs_lnum_t block_size = 3*cdoq->n_faces + cdoq->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       2); /* Only one block for the system
                                            * matrix. The 2nd block is for the
                                            * algebraic transformation */

      /* Add the only block and then define the underpinning structures */

      cs_cdo_system_block_t
        *a = cs_cdo_system_add_xblock(sh, 0,       /* block id */
                                      block_size); /* n_dofs */

      _build_shared_structures_full_system(a, false);

      /* Add a second block for the (1,0) and (0,1) blocks and then define the
         underpinning structures. The (0,1) block needs to be transposed before
         using it */

      cs_cdo_system_block_t  *bdiv =
        cs_cdo_system_add_ublock(sh, 1,               /* block_id */
                                 connect->c2f,        /* adjacency */
                                 cs_flag_primal_face, /* column location */
                                 cdoq->n_faces,       /* n_elements */
                                 3,                   /* stride */
                                 true);               /* interlaced */

      cs_cdo_system_ublock_t *b_ub
        = (cs_cdo_system_ublock_t *)bdiv->block_pointer;

      /* Define the bdiv block by hand */

      cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

      b_ub->adjacency = connect->c2f;  /* shared pointer */
      b_ub->values = sc->block21_op;   /* shared pointer */
      assert(b_ub->values != nullptr);
      b_ub->shared_structures = true;
      b_ub->range_set = rset;          /* shared pointer */
      b_ub->interface_set
        = (cs_interface_set_t *)rset->ifs; /* shared pointer */
    }
    break;

  default:
    /* CS_PARAM_SADDLE_SOLVER_FGMRES
     * CS_PARAM_SADDLE_SOLVER_MUMPS */
    {
      cs_lnum_t block_size = 3*cdoq->n_faces + cdoq->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       1);

      /* Add the only block and then define the underpinning structures */

      cs_cdo_system_block_t
        *a = cs_cdo_system_add_xblock(sh, 0,       /* block id */
                                      block_size); /* n_dofs */

      /* Fill the xblock (with diagonal pressure block) */

      if (nsp->model_flag & CS_NAVSTO_MODEL_WITH_SOLIDIFICATION)
        _build_shared_structures_full_system(a, true);
      else
        _build_shared_structures_full_system(a, false);
    }
    break;

  } /* Switch on saddle-point solver */

  assert(sh != nullptr);
  sc->system_helper = sh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the saddle solver and its context for a CDO-Fb scheme solving
 *        the Navier-Stokes equation using a monolithic approach for the
 *        velocity-pressure coupling
 *
 * \param[in]      nsp      set of parameters for the Navier-Stokes system
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] sc       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_solver(const cs_navsto_param_t  *nsp,
                                     const cs_param_saddle_t  *saddlep,
                                     cs_cdofb_monolithic_t    *sc)
{
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  n_cells = cs_shared_quant->n_cells;
  const cs_param_sles_t  *b11_slesp = saddlep->block11_sles_param;

  cs_sles_t *sles = cs_sles_find_or_add(b11_slesp->field_id, nullptr);

  cs_saddle_solver_t  *solver = cs_saddle_solver_add(n_faces, 3,
                                                     n_cells, 1,
                                                     saddlep,
                                                     sc->system_helper,
                                                     sles); /* main sles */

  sc->saddle_solver = solver;

  /* Set the solve function pointer */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
      sc->solve = cs_cdofb_monolithic_sles_alu;

      cs_saddle_solver_context_alu_create(solver);

      cs_saddle_solver_context_alu_t *ctx
        = (cs_saddle_solver_context_alu_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfvp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_vector;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    {
      sc->solve = cs_cdofb_monolithic_sles_notay;

      cs_saddle_solver_context_notay_create(solver);

      cs_saddle_solver_context_notay_t *ctx
        = (cs_saddle_solver_context_notay_t *)solver->context;

      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    if (saddlep->solver_class == CS_PARAM_SOLVER_CLASS_PETSC)
      sc->solve = cs_cdofb_monolithic_sles_full_system;

    else {

      sc->solve = cs_cdofb_monolithic_sles_gkb_inhouse;

      cs_saddle_solver_context_gkb_create(solver);

      cs_saddle_solver_context_gkb_t *ctx
        = (cs_saddle_solver_context_gkb_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfvp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_vector;

    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
    {
      sc->solve = cs_cdofb_monolithic_sles_block_krylov;

      cs_saddle_solver_context_block_pcd_create(m->n_cells_with_ghosts, solver);

      cs_saddle_solver_context_block_pcd_t *ctx
        = (cs_saddle_solver_context_block_pcd_t *)solver->context;

      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_vector;

      if (nsp->turbulence->model->model == CS_TURB_NONE)
        ctx->pty_22 = nsp->lam_viscosity;
      else
        ctx->pty_22 = nsp->tot_viscosity;

      const cs_real_t  rho0 = nsp->mass_density->ref_value;

      /* alpha coefficient is related to time */

      cs_real_t  alpha;
      if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
        alpha = 0.01*nsp->lam_viscosity->ref_value;
      else
        alpha = 1/ts->dt[0];

      ctx->schur_scaling = rho0 * alpha;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
      sc->solve = cs_cdofb_monolithic_sles_uzawa_cg;

      cs_saddle_solver_context_uzawa_cg_create(m->n_cells_with_ghosts, solver);

      cs_saddle_solver_context_uzawa_cg_t *ctx
        = (cs_saddle_solver_context_uzawa_cg_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfvp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_vector;

      if (nsp->turbulence->model->model == CS_TURB_NONE)
        ctx->pty_22 = nsp->lam_viscosity;
      else
        ctx->pty_22 = nsp->tot_viscosity;

      const cs_real_t  rho0 = nsp->mass_density->ref_value;

      /* zeta coefficient is related to time */

      cs_real_t  zeta;
      if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
        zeta = 0.01*nsp->lam_viscosity->ref_value;
      else
        zeta = 1/ts->dt[0];

      ctx->alpha = rho0 * zeta; /* similar to a schur scaling */
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
      sc->solve = cs_cdofb_monolithic_sles_simple;

      cs_saddle_solver_context_simple_create(m->n_cells_with_ghosts, solver);

      cs_saddle_solver_context_simple_t *ctx
        = (cs_saddle_solver_context_simple_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfvp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_vector;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_vector;

      if (nsp->turbulence->model->model == CS_TURB_NONE)
        ctx->pty_22 = nsp->lam_viscosity;
      else
        ctx->pty_22 = nsp->tot_viscosity;
    }
    break;

  default:
    /* CS_PARAM_SADDLE_SOLVER_FGMRES
     * CS_PARAM_SADDLE_SOLVER_MUMPS */
    sc->solve = cs_cdofb_monolithic_sles_full_system;
    break;

  } /* Switch on saddle-point solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Augmented Lagrangian-Uzawa algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_alu(const cs_navsto_param_t  *nsp,
                             cs_saddle_solver_t       *solver,
                             cs_real_t                *u_f,
                             cs_real_t                *p_c)
{
  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  /* Sanity checks */

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_ALU)
    bft_error(__FILE__, __LINE__, 0,
              "%s: ALU algorithm is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

#if defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_system_helper_t  *sh = solver->system_helper;

  assert(sh != nullptr);
  assert(sh->n_blocks == 2);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  cs_iter_algo_type_t  type = CS_ITER_ALGO_DEFAULT | CS_ITER_ALGO_TWO_LEVEL;
  solver->algo = cs_iter_algo_create_with_settings(type,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;

  /* 1. Build the uzawa context */
  /* -------------------------- */

  cs_saddle_solver_context_alu_t *ctx
    = (cs_saddle_solver_context_alu_t *)solver->context;
  assert(ctx != nullptr);

  _alu_init_context(nsp, solver, ctx);

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  cs_saddle_solver_alu_incr(solver, u_f, p_c);

  /* 3. Monitoring and output */
  /* ------------------------ */

  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d (inner:%4d) | residual:% -8.4e\n",
                  __func__, saddlep->name, algo_ctx->cvg_status,
                  n_iters, algo_ctx->n_inner_iter, algo_ctx->res);

  /* Memory cleaning */

  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach. The system is
 *        split into a velocity block and the (unassembled) divergence operator
 *        Block preconditioning using a Schur approximation on a Krylov solver
 *        such as the GCR or MINRES is available.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_block_krylov(const cs_navsto_param_t  *nsp,
                                      cs_saddle_solver_t       *solver,
                                      cs_real_t                *u_f,
                                      cs_real_t                *p_c)
{
  if (solver == nullptr)
    return 0;

  /* 0. Initialization and checkings */
  /* ------------------------------- */

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_GCR &&
      saddlep->solver != CS_PARAM_SADDLE_SOLVER_MINRES)
    bft_error(__FILE__, __LINE__, 0,
              "%s: GCR or MINRES is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

#if defined(DEBUG) && !defined(NDEBUG)
  assert(sh != nullptr);
  assert(sh->n_blocks == 2);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  cs_iter_algo_type_t  algo_type =
    CS_ITER_ALGO_DEFAULT | CS_ITER_ALGO_TWO_LEVEL;

  solver->algo = cs_iter_algo_create_with_settings(algo_type,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_iter_algo_default_t *algo_ctx
    = static_cast<cs_iter_algo_default_t *>(solver->algo->context);

  /* 1. Build the block preconditioner */
  /* --------------------------------- */

  cs_saddle_solver_context_block_pcd_t *ctx
    = static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  assert(ctx != nullptr);

  /* Update the context after the matrix building */

  ctx->m11 = cs_cdo_system_get_matrix(sh, 0);
  ctx->b11_max_size = cs::max(cs_matrix_get_n_columns(ctx->m11),
                              solver->n1_scatter_dofs);

  /* Prepare the solution array at faces. It has to be allocated to a greater
   * size in case of parallelism in order to allow for a correct matrix-vector
   * product */

  cs_real_t *x1 = nullptr;

  if (cs_glob_n_ranks > 1) {
    CS_MALLOC(x1, ctx->b11_max_size, cs_real_t);
    cs_array_real_copy(solver->n1_scatter_dofs, u_f, x1);
  }
  else
    x1 = u_f;

  /* Prepare the context structure according to the choice of block
     preconditioner. In particular, define the Schur complement approximation
     if needed */

  const double  gamma = cs_param_saddle_get_augmentation_coef(saddlep);
  const cs_param_sles_t  *schur_slesp = saddlep->schur_sles_param;
  int  n_xtra_iters = 0;

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
    ctx->m11_inv_diag
      = cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

    algo_ctx->n_inner_iter += n_xtra_iters;

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    ctx->m22_mass_diag = _get_scaled_diag_m22(nsp, ctx->pty_22, gamma);
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
    ctx->m22_mass_diag = _get_scaled_diag_m22(nsp, ctx->pty_22, gamma);

    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    ctx->m11_inv_diag
      = cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

    algo_ctx->n_inner_iter += n_xtra_iters;

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));

    ctx->m22_mass_diag = _get_scaled_diag_m22(nsp, ctx->pty_22, gamma);
    break;

  default:
    /* Do nothing else */
    break;

  }

  /* Augmentation of rhs */

  if (gamma > 0.) {

    const cs_lnum_t n1_dofs = solver->n1_scatter_dofs;
    const cs_lnum_t n2_dofs = solver->n2_scatter_dofs;
    const cs_real_t *cell_vol = cs_shared_quant->cell_vol;

    cs_real_t *rhs1 = sh->rhs_array[0];
    cs_real_t *rhs2 = sh->rhs_array[1];

    cs_real_t *btilda_f = nullptr, *btilda_c = nullptr;

    CS_MALLOC(btilda_f, ctx->b11_max_size, cs_real_t);
    CS_MALLOC(btilda_c, solver->n2_scatter_dofs, cs_real_t);

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < solver->n2_scatter_dofs; i2++)
      btilda_c[i2] = gamma/cell_vol[i2]*rhs2[i2];

    cs_saddle_system_b12_matvec(sh, btilda_c, btilda_f,
                                true); /* reset btilda_f */

  /* Build augmented b1 = b1 + gamma*m12.W^-1.b_c */

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      rhs1[i1] += btilda_f[i1];

    /* Free arrays */
    CS_FREE(btilda_c);
    CS_FREE(btilda_f);
  }

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_MINRES:
    cs_saddle_solver_minres(solver, x1, p_c);
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
    cs_saddle_solver_gcr(solver, x1, p_c);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid saddle solver", __func__);
    break;
  }

  /* Copy back to the original array the velocity values at faces */

  if (cs_glob_n_ranks > 1) {
    cs_array_real_copy(solver->n1_scatter_dofs, x1, u_f);
    CS_FREE(x1);
  }

  /* 3. Monitoring and output */
  /* ------------------------ */

  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d (inner:%4d) | residual:% -8.4e\n",
                  __func__, saddlep->name, algo_ctx->cvg_status,
                  n_iters, algo_ctx->n_inner_iter, algo_ctx->res);

  /* Memory cleaning */

  cs_saddle_solver_context_block_pcd_clean(ctx);
  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_full_system(const cs_navsto_param_t  *nsp,
                                     cs_saddle_solver_t       *solver,
                                     cs_real_t                *u_f,
                                     cs_real_t                *p_c)
{
  CS_NO_WARN_IF_UNUSED(nsp);

  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_MUMPS &&
      saddlep->solver != CS_PARAM_SADDLE_SOLVER_FGMRES)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Full system solver is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

  /* Prepare the solution and RHS arrays */

#if defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_system_helper_t  *sh = solver->system_helper;

  assert(sh != nullptr);
  assert(sh->n_blocks == 1);

  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  /* Solve the saddle-point problem */

  solver->algo = cs_iter_algo_create_with_settings(CS_ITER_ALGO_DEFAULT,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_saddle_solver_sles_full_system(solver, u_f, p_c);

  /* Monitoring */

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;
  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code=%-d | n_iter:%d | residual:% -8.4e\n",
                  __func__, saddlep->name,
                  algo_ctx->cvg_status, n_iters, algo_ctx->res);

  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Golub-Kahan Bidiagonalization algorithm.
 *        In-house implementation. The PETSc implementation is also available
 *        but appears less efficient in our tests.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_gkb_inhouse(const cs_navsto_param_t  *nsp,
                                     cs_saddle_solver_t       *solver,
                                     cs_real_t                *u_f,
                                     cs_real_t                *p_c)
{
  CS_NO_WARN_IF_UNUSED(nsp);

  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  /* Sanity checks */

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_GKB)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: GKB algorithm is expected.\n"
              "%s: Please check your settings.\n",
              __func__,
              __func__);

#if defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_system_helper_t  *sh = solver->system_helper;

  assert(sh != nullptr);
  assert(sh->n_blocks == 2);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  cs_iter_algo_type_t  type = CS_ITER_ALGO_DEFAULT | CS_ITER_ALGO_TWO_LEVEL;
  solver->algo = cs_iter_algo_create_with_settings(type,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;

  /* 1. Build the GKB context */
  /* ------------------------- */

  cs_saddle_solver_context_gkb_t *ctx
    = (cs_saddle_solver_context_gkb_t *)solver->context;
  assert(ctx != nullptr);

  _gkb_init_context(solver, ctx);

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  cs_saddle_solver_gkb_inhouse(solver, u_f, p_c);

  /* 3. Monitoring and output */
  /* ------------------------ */

  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d (inner:%4d) | residual:% -8.4e\n",
                  __func__, saddlep->name, algo_ctx->cvg_status,
                  n_iters, algo_ctx->n_inner_iter, algo_ctx->res);

  /* Memory cleaning */

  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Notay's algebraic transformation.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_notay(const cs_navsto_param_t  *nsp,
                               cs_saddle_solver_t       *solver,
                               cs_real_t                *u_f,
                               cs_real_t                *p_c)
{
  CS_NO_WARN_IF_UNUSED(nsp);

  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Notay's transformation is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

  solver->algo = cs_iter_algo_create_with_settings(CS_ITER_ALGO_DEFAULT,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  /* Prepare the solution and RHS arrays */

#if defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_system_helper_t  *sh = solver->system_helper;

  assert(sh != nullptr);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  /* Solve the saddle-point problem */

  cs_saddle_solver_notay(solver, u_f, p_c);

  /* Monitoring */

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;
  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code=%-d | n_iter:%d | residual:% -8.4e\n",
                  __func__, saddlep->name,
                  algo_ctx->cvg_status, n_iters, algo_ctx->res);

  /* Memory cleaning. The Notay's context is simple. No need to be cleaned. */

  cs_iter_algo_free(&solver->algo);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Uzawa-CG algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_uzawa_cg(const cs_navsto_param_t  *nsp,
                                  cs_saddle_solver_t       *solver,
                                  cs_real_t                *u_f,
                                  cs_real_t                *p_c)
{
  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_UZAWA_CG)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Uzawa-CG algorithm is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

#if defined(DEBUG) && !defined(NDEBUG)
  assert(sh != nullptr);
  assert(sh->n_blocks == 2);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  cs_iter_algo_type_t  type = CS_ITER_ALGO_DEFAULT | CS_ITER_ALGO_TWO_LEVEL;

  solver->algo = cs_iter_algo_create_with_settings(type,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;

  /* 0. Partial initialization of the context */
  /* ---------------------------------------- */

  cs_saddle_solver_context_uzawa_cg_t *ctx
    = (cs_saddle_solver_context_uzawa_cg_t *)solver->context;
  assert(ctx != nullptr);

  /* Update the context after the matrix building */

  _uzawa_cg_init_context(nsp, solver, ctx);

  /* Prepare the solution array at faces. It has to be allocated to a greater
   * size in case of parallelism in order to allow for a correct matrix-vector
   * product */

  cs_real_t *x1 = nullptr;

  if (cs_glob_n_ranks > 1) {
    CS_MALLOC(x1, ctx->b11_max_size, cs_real_t);
    cs_array_real_copy(solver->n1_scatter_dofs, u_f, x1);
  }
  else
    x1 = u_f;


  /* 1. Build the schur approximation */
  /* -------------------------------- */

  const cs_param_sles_t  *schur_slesp = saddlep->schur_sles_param;
  int  n_xtra_iters = 0;

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
    {
      ctx->m11_inv_diag =
        cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

      algo_ctx->n_inner_iter += n_xtra_iters;

      ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                     ctx->m11_inv_diag,
                                                     ctx->m21_adj,
                                                     ctx->m21_val,
                                                     &(ctx->schur_diag),
                                                     &(ctx->schur_xtra));

    }
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      ctx->m11_inv_diag =
        cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

      algo_ctx->n_inner_iter += n_xtra_iters;

      ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                     ctx->m11_inv_diag,
                                                     ctx->m21_adj,
                                                     ctx->m21_val,
                                                     &(ctx->schur_diag),
                                                     &(ctx->schur_xtra));

    }
    break;

  default:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    /* Do nothing else */
    break;

  } /* Switch on the type of Schur approximation */

  /* 2. Solve the saddle-point system */
  /* -------------------------------- */

  cs_saddle_solver_uzawa_cg(solver, x1, p_c);

  /* Copy back to the original array the velocity values at faces */

  if (cs_glob_n_ranks > 1) {
    cs_array_real_copy(solver->n1_scatter_dofs, x1, u_f);
    CS_FREE(x1);
  }

  /* 3. Monitoring and output */
  /* ------------------------ */

  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d (inner:%4d) | residual:% -8.4e\n",
                  __func__, saddlep->name, algo_ctx->cvg_status,
                  n_iters, algo_ctx->n_inner_iter, algo_ctx->res);

  /* Memory cleaning */

  cs_saddle_solver_context_uzawa_cg_clean(ctx);
  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the SIMPLE algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_simple(const cs_navsto_param_t  *nsp,
                                cs_saddle_solver_t       *solver,
                                cs_real_t                *u_f,
                                cs_real_t                *p_c)
{
  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_SIMPLE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: SINGLE algorithm is expected.\n"
              "%s: Please check your settings.\n", __func__, __func__);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

#if defined(DEBUG) && !defined(NDEBUG)
  assert(sh != nullptr);
  assert(sh->n_blocks == 2);
  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  cs_iter_algo_type_t  type = CS_ITER_ALGO_DEFAULT | CS_ITER_ALGO_TWO_LEVEL;

  solver->algo = cs_iter_algo_create_with_settings(type,
                                                   saddlep->verbosity,
                                                   saddlep->cvg_param);

  cs_iter_algo_default_t *algo_ctx
    = (cs_iter_algo_default_t *)solver->algo->context;

  /* 0. Partial initialization of the context */
  /* ---------------------------------------- */

  cs_saddle_solver_context_simple_t *ctx
    = (cs_saddle_solver_context_simple_t *)solver->context;
  assert(ctx != nullptr);

  /* Update the context after the matrix building */

  _simple_init_context(solver, ctx);

  /* Prepare the solution array at faces. It has to be allocated to a greater
   * size in case of parallelism in order to allow for a correct matrix-vector
   * product */

  cs_real_t *x1 = nullptr;

  if (cs_glob_n_ranks > 1) {
    CS_MALLOC(x1, ctx->b11_max_size, cs_real_t);
    cs_array_real_copy(solver->n1_scatter_dofs, u_f, x1);
  }
  else
    x1 = u_f;


  /* 1. Build the schur approximation */
  /* -------------------------------- */

  const cs_param_sles_t  *schur_slesp = saddlep->schur_sles_param;
  int  n_xtra_iters = 0;

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                   ctx->m11_inv_diag,
                                                   ctx->m21_adj,
                                                   ctx->m21_val,
                                                   &(ctx->schur_diag),
                                                   &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
    {
      ctx->m11_inv_diag =
        cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

      algo_ctx->n_inner_iter += n_xtra_iters;

      ctx->schur_matrix = _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                                     ctx->m11_inv_diag,
                                                     ctx->m21_adj,
                                                     ctx->m21_val,
                                                     &(ctx->schur_diag),
                                                     &(ctx->schur_xtra));
    }
    break;

  default:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    /* Do nothing else */
    break;

  } /* Switch on the type of Schur approximation */

  /* 2. Solve the saddle-point system */
  /* -------------------------------- */

  cs_saddle_solver_simple(solver, x1, p_c);

  /* Copy back to the original array the velocity values at faces */

  if (cs_glob_n_ranks > 1) {
    cs_array_real_copy(solver->n1_scatter_dofs, x1, u_f);
    CS_FREE(x1);
  }

  /* 3. Monitoring and output */
  /* ------------------------ */

  int  n_iters = algo_ctx->n_algo_iter;

  cs_saddle_solver_update_monitoring(solver, n_iters);

  /* Output information about the convergence */

  if (saddlep->verbosity > 0 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT, "\n  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d (inner:%4d) | residual:% -8.4e\n",
                  __func__, saddlep->name, algo_ctx->cvg_status,
                  n_iters, algo_ctx->n_inner_iter, algo_ctx->res);

  /* Memory cleaning */

  cs_saddle_solver_context_simple_clean(ctx);
  cs_iter_algo_free(&(solver->algo));

  return n_iters;
}

/*----------------------------------------------------------------------------*/

#undef _petsc_cmd
END_C_DECLS
