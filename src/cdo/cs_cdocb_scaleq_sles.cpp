/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO cell-based schemes with a scaleq equations
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

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  BFT headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_blas.h"
#include "cs_cdo_blas.h"
#include "cs_cdo_solve.h"
#include "cs_equation.h"
#include "cs_cdocb_scaleq.h"
#include "cs_fp_exception.h"
#include "cs_matrix_default.h"
#include "cs_parall.h"
#include "cs_saddle_solver.h"
#include "cs_saddle_system.h"
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

#include "cs_cdocb_scaleq_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdocb_scaleq_sles.c
 *
 * \brief Functions dedicated to to the linear algebra settings and operations
 *        in case of CDO cell-based schemes with a scalar-valued equations
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOCB_SCALEQ_SLES_DBG      0

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
 *
 * \param[in, out] block  pointer to a block structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_shared_structures_full_system(cs_cdo_system_block_t  *block)
{
  /* Compute the range set for an array of size n_faces + n_cells
   * flux is attached to faces and potential to cells
   *
   * Storage for the global numbering: flux | potential
   */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  size = n_faces + m->n_cells;

  assert(block->type == CS_CDO_SYSTEM_BLOCK_EXT);
  cs_cdo_system_xblock_t *xb = (cs_cdo_system_xblock_t *)block->block_pointer;

  /* 1. Build the interface set and the range set structures */

  cs_interface_set_t  *ifs = connect->face_ifs;

  xb->interface_set = ifs;
  xb->range_set     = cs_range_set_create(xb->interface_set,
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

  int  max_stencil = 0;
  for (cs_lnum_t f = 0; f < n_faces; f++) {
    int  sten =
      f2f->idx[f+1]-f2f->idx[f] + 1 + 2*(f2c->idx[f+1] - f2c->idx[f]);
    max_stencil = CS_MAX(max_stencil, sten);
  }

  cs_gnum_t *grows = nullptr, *gcols = nullptr;
  BFT_MALLOC(grows, max_stencil, cs_gnum_t);
  BFT_MALLOC(gcols, max_stencil, cs_gnum_t);

  /*
   *   | A | Bt |
   *   |---|----|
   *   | B | 0  |
   *
   *  The block A is n_faces * n_faces
   *  The block B is n_cells * n_faces
   */

  /* Only on faces (B is build in the same time as Bt for potential DoFs) */

  for (cs_lnum_t frow_id = 0; frow_id < n_faces; frow_id++) {

    const cs_lnum_t  start = f2f->idx[frow_id];
    const cs_lnum_t  end = f2f->idx[frow_id+1];

    /* A face-face entry corresponds to an entry + the diagonal which is not
       taken into account in the face --> face connectivity. The B and Bt
       operators have the same sparsity related to the c2f connectivity. This
       is multiply by two since one considers B and Bt. */

    int  n_extra = end - start;
    int  n_entries = n_extra + 1 + 2*(f2c->idx[frow_id+1]-f2c->idx[frow_id]);

    const cs_gnum_t  grow_id = xb->range_set->g_id[frow_id];

    int shift = 0;

    /* Diagonal term is excluded in this connectivity. Add it "manually" */

    grows[shift] = grow_id;
    gcols[shift] = grow_id;
    shift++;

    /* Extra diagonal couples */

    for (cs_lnum_t idx = start; idx < end; idx++) {

      grows[shift] = grow_id;
      gcols[shift] = xb->range_set->g_id[f2f->ids[idx]];
      shift++;

    } /* Loop on extra-diag. entries */

    /* Loop on grad/div entries */

    for (cs_lnum_t idx = f2c->idx[frow_id]; idx < f2c->idx[frow_id+1]; idx++) {

      const cs_lnum_t  ccol_id = f2c->ids[idx];
      const cs_gnum_t  gcol_id = xb->range_set->g_id[n_faces + ccol_id];

      grows[shift] = grow_id;
      gcols[shift] = gcol_id;
      shift++;

      /* Its transposed B_x, B_y, B_z */

      grows[shift] = gcol_id;
      gcols[shift] = grow_id;
      shift++;

    } /* Loop on grad/div entries */

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  /* 3. Build the matrix structure */

  cs_matrix_assembler_compute(xb->matrix_assembler);

  xb->matrix_structure =
    cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR,
                                              xb->matrix_assembler);

  /* Free temporary buffers */

  BFT_FREE(grows);
  BFT_FREE(gcols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a diagonal scaled mass matrix
 *
 * \param[in] pty  property related to the (2,2) block
 *
 * \return a pointer to a newly allocated array
 */
/*----------------------------------------------------------------------------*/

static cs_real_t *
_get_scaled_diag_m22(const cs_property_t  *pty)
{
  const cs_cdo_quantities_t  *cdoq = cs_shared_quant;
  const cs_lnum_t  n_cells = cdoq->n_cells;

  cs_real_t *m22_mass_diag = nullptr;
  BFT_MALLOC(m22_mass_diag, n_cells, cs_real_t);

  /* Compute scaling coefficients */

  const cs_real_t  pty_scaling = pty->ref_value;
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    m22_mass_diag[i] = pty_scaling/cdoq->cell_vol[i];

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
 * \param[out] p_diag_smat      diagonal coeffs for the Schur matrix
 * \param[out] p_xtra_smat      extra-diagonal coeffs for the Schur matrix
 *
 * \return a pointer to the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_schur_approx_diag_inv_m11(cs_param_solver_class_t  mat_class,
                                  const cs_real_t         *m11_inv_approx,
                                  cs_real_t              **p_diag_smat,
                                  cs_real_t              **p_xtra_smat)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)mesh->b_face_cells;

  /* Native format for the Schur approximation matrix */

  cs_real_t *diag_smat = nullptr;
  cs_real_t *xtra_smat = nullptr;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  cs_array_real_fill_zero(n_cells_ext, diag_smat);
  cs_array_real_fill_zero(2*n_i_faces, xtra_smat);

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  m11_inv_ff = m11_inv_approx[f_id];

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = -m11_inv_ff;

    /* Diagonal contributions */

    cs_lnum_t  cell_i = i_face_cells[f_id][0];
    cs_lnum_t  cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] += m11_inv_ff;
    diag_smat[cell_j] += m11_inv_ff;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/

  const cs_real_t  *_shift = m11_inv_approx + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *m11_inv_ff = _shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += m11_inv_ff[k] * nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  cs_matrix_t *smat = nullptr;

  if (mat_class == CS_PARAM_SOLVER_CLASS_HYPRE)
    smat = cs_matrix_external("HYPRE_ParCSR",
                              false, /* symmetry */
                              1, 1);
  else
    smat = cs_matrix_msr(false, /* symmetry */
                         1, 1);

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
                                 (const cs_real_t *)quant->cell_vol,
                                 (const cs_real_3_t *)quant->i_face_normal);

  /* Return the associated matrix and set the pointers to return
   * Return arrays (to be freed when the algorithm is converged) */

  *p_diag_smat = diag_smat;
  *p_xtra_smat = xtra_smat;

  return smat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context structure associated to an ALU algorithm
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] ctx     additional memebers specific to ALU algo.
 */
/*----------------------------------------------------------------------------*/

static void
_alu_init_context(cs_saddle_solver_t              *solver,
                  cs_saddle_solver_context_alu_t  *ctx)
{
  assert(ctx != nullptr);

  if (solver->do_setup == false)
    return;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  assert(solver->n2_scatter_dofs == quant->n_cells);
  assert(solver->n1_scatter_dofs == quant->n_faces);

  /* Buffers of size n2_scatter_dofs */

  BFT_MALLOC(ctx->inv_m22, solver->n2_scatter_dofs, cs_real_t);

# pragma omp parallel for if (solver->n2_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < solver->n2_scatter_dofs; i2++)
    ctx->inv_m22[i2] = 1./quant->cell_vol[i2];

  BFT_MALLOC(ctx->res2, solver->n2_scatter_dofs, cs_real_t);
  BFT_MALLOC(ctx->m21x1, solver->n2_scatter_dofs, cs_real_t);

  /* Buffers of size n1_scatter_dofs */

  BFT_MALLOC(ctx->b1_tilda, solver->n1_scatter_dofs, cs_real_t);
  BFT_MALLOC(ctx->rhs, solver->n1_scatter_dofs, cs_real_t);
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

  assert(n1_dofs == quant->n_faces);
  assert(n2_dofs == quant->n_cells);

  /* Orthogonalization coefficients */

  ctx->alpha = ctx->beta = ctx->zeta = 0.;

  if (solver->do_setup == false)
    return;

  /* Buffers of size n2_dofs */

  BFT_MALLOC(ctx->q, n2_dofs, cs_real_t);
  BFT_MALLOC(ctx->d, n2_dofs, cs_real_t);
  BFT_MALLOC(ctx->m21v, n2_dofs, cs_real_t);
  BFT_MALLOC(ctx->inv_m22, n2_dofs, cs_real_t);

  ctx->m22 = quant->cell_vol;   /* shared pointer */
  for (cs_lnum_t i = 0; i < n2_dofs; i++)
    ctx->inv_m22[i] = 1./quant->cell_vol[i];

  /* Buffers of size n1_dofs */

  BFT_MALLOC(ctx->m12q, n1_dofs, cs_real_t);
  BFT_MALLOC(ctx->x1_tilda, n1_dofs, cs_real_t);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  max_b11_size = CS_MAX(cs_matrix_get_n_columns(m11), n1_dofs);

  BFT_MALLOC(ctx->w, max_b11_size, cs_real_t);
  BFT_MALLOC(ctx->v, max_b11_size, cs_real_t);

  /* Rk: rhs_tilda stores quantities in space X1 and X2 alternatively */

  BFT_MALLOC(ctx->rhs_tilda, CS_MAX(n1_dofs, n2_dofs), cs_real_t);

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
    ctx->zeta_size = CS_MAX(1, tt - 1);
  else if (gamma < 1e3)
    ctx->zeta_size = CS_MAX(1, tt - 2);
  else if (gamma < 1e4)
    ctx->zeta_size = CS_MAX(1, tt - 3);
  else
    ctx->zeta_size = CS_MAX(1, tt - 4);

  BFT_MALLOC(ctx->zeta_array, ctx->zeta_size, cs_real_t);
  for (int i = 0; i < ctx->zeta_size; i++)
    ctx->zeta_array[i] = 0.;

  ctx->zeta_square_sum = 0.;
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
cs_cdocb_scaleq_sles_init_sharing(const cs_mesh_t            *mesh,
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
 * \brief Define the system helper for a CDO-Cb scheme solving a scalar-valued
 *        equation (saddle-point system)
 *
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] eqc      pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_sles_init_system_helper(const cs_param_saddle_t  *saddlep,
                                        cs_cdocb_scaleq_t        *eqc)
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
    {
      cs_lnum_t  block_sizes[2];
      block_sizes[0] = cdoq->n_faces, block_sizes[1] = cdoq->n_cells;

      /* Create the system helper */

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       2,
                                       block_sizes,
                                       2);

      /* Add a first block for the (1,1)-block and then define the underpinning
         structures */

      /* Choose the right class of matrix to avoid copy.
       * The way to perform the assembly may change if an external library is
       * used for solving the linear system.
       */

      cs_cdo_system_matrix_class_t
        matclass = cs_cdo_system_which_matrix_class(b11_slesp->solver_class);

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_dblock(sh, 0,                /* block id */
                                 matclass,
                                 cs_flag_primal_face , /* location */
                                 cdoq->n_faces,        /* n_elements */
                                 1,                    /* stride */
                                 false,                /* interlaced */
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
                                 1,                   /* stride */
                                 true);               /* interlaced */

      cs_cdo_system_dblock_t *a_db = (cs_cdo_system_dblock_t *)a->block_pointer;
      cs_cdo_system_ublock_t *b_ub
        = (cs_cdo_system_ublock_t *)bdiv->block_pointer;

      /* Define the bdiv block by hand */

      b_ub->adjacency = connect->c2f;            /* shared pointer */
      b_ub->values = eqc->block21_op;             /* shared pointer */
      assert(b_ub->values != nullptr);
      b_ub->shared_structures = true;
      b_ub->range_set = a_db->range_set;         /* shared pointer */
      b_ub->interface_set = a_db->interface_set; /* shared pointer */
    }
    break;

  default:
    /* CS_PARAM_SADDLE_SOLVER_FGMRES
     * CS_PARAM_SADDLE_SOLVER_NOTAY
     * CS_PARAM_SADDLE_SOLVER_MUMPS */
    {
      cs_lnum_t block_size = cdoq->n_faces + cdoq->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       1);

      /* Add the only block and then define the underpinning structures */

      cs_cdo_system_block_t
        *a = cs_cdo_system_add_xblock(sh, 0,       /* block id */
                                      block_size); /* n_dofs */

      /* Fill the xblock (with diagonal pressure block) */

      _build_shared_structures_full_system(a);
    }
    break;

  } /* Switch on saddle-point solver */

  assert(sh != nullptr);
  eqc->system_helper = sh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the saddle solver and its context for a CDO-Cb scheme solving
 *        a scalar-valued equation (saddle-point problem)
 *
 * \param[in]      eqp      set of equation parameters
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] eqc      pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_sles_init_solver(const cs_equation_param_t  *eqp,
                                 const cs_param_saddle_t    *saddlep,
                                 cs_cdocb_scaleq_t          *eqc)
{
  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  n_cells = cs_shared_quant->n_cells;
  const cs_param_sles_t  *b11_slesp = saddlep->block11_sles_param;

  cs_sles_t *sles = cs_sles_find_or_add(b11_slesp->field_id, nullptr);

  cs_saddle_solver_t  *solver = cs_saddle_solver_add(n_faces, 1,
                                                     n_cells, 1,
                                                     saddlep,
                                                     eqc->system_helper,
                                                     sles); /* main sles */

  eqc->saddle_solver = solver;

  /* Set the solve function pointer */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
      eqc->solve = cs_cdocb_scaleq_sles_alu;

      cs_saddle_solver_context_alu_create(solver);

      cs_saddle_solver_context_alu_t *ctx
        = (cs_saddle_solver_context_alu_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfsp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_scalar;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_scalar;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    {
      eqc->solve = cs_cdocb_scaleq_sles_notay;

      cs_saddle_solver_context_notay_create(solver);

      cs_saddle_solver_context_notay_t *ctx
        = (cs_saddle_solver_context_notay_t *)solver->context;

      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_scalar;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    if (saddlep->solver_class == CS_PARAM_SOLVER_CLASS_PETSC)
      eqc->solve = cs_cdocb_scaleq_sles_full_system;

    else {

      eqc->solve = cs_cdocb_scaleq_sles_gkb_inhouse;

      cs_saddle_solver_context_gkb_create(solver);

      cs_saddle_solver_context_gkb_t *ctx
        = (cs_saddle_solver_context_gkb_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfsp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_scalar;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_scalar;

    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
    {
      eqc->solve = cs_cdocb_scaleq_sles_block_krylov;

      cs_saddle_solver_context_block_pcd_create(m->n_cells_with_ghosts, solver);

      cs_saddle_solver_context_block_pcd_t *ctx
        = (cs_saddle_solver_context_block_pcd_t *)solver->context;

      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_scalar;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_scalar;

      /* Scaling coefficient for the Schur complement approximation */

      ctx->pty_22 = eqp->diffusion_property;
      ctx->schur_scaling = 1.;  /* better value ? */
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
      cs_saddle_solver_context_uzawa_cg_create(m->n_cells_with_ghosts, solver);

      cs_saddle_solver_context_uzawa_cg_t *ctx
        = (cs_saddle_solver_context_uzawa_cg_t *)solver->context;

      ctx->square_norm_b11 = cs_cdo_blas_square_norm_pfsp;
      ctx->m12_vector_multiply = cs_saddle_solver_m12_multiply_scalar;
      ctx->m21_vector_multiply = cs_saddle_solver_m21_multiply_scalar;

      eqc->solve = cs_cdocb_scaleq_sles_uzawa_cg;
    }
    break;

  default:
    /* CS_PARAM_SADDLE_SOLVER_FGMRES
     * CS_PARAM_SADDLE_SOLVER_MUMPS */
    eqc->solve = cs_cdocb_scaleq_sles_full_system;
    break;

  } /* Switch on saddle-point solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes. Solve
 *        the saddle-point system using the Augmented Lagrangian-Uzawa
 *        algorithm.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_alu(cs_saddle_solver_t  *solver,
                         cs_real_t           *flux,
                         cs_real_t           *pot)
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

  _alu_init_context(solver, ctx);

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  cs_saddle_solver_alu_incr(solver, flux, pot);

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
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes. The
 *        system is split into a flux block and the (unassembled) divergence
 *        operator.
 *        Block preconditioning using a Schur approximation on a Krylov solver
 *        such as the GCR or MINRES is available.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_block_krylov(cs_saddle_solver_t  *solver,
                                  cs_real_t           *flux,
                                  cs_real_t           *pot)
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
    = (cs_iter_algo_default_t *)solver->algo->context;

  /* 1. Build the block preconditioner */
  /* --------------------------------- */

  cs_saddle_solver_context_block_pcd_t *ctx
    = (cs_saddle_solver_context_block_pcd_t *)solver->context;
  assert(ctx != nullptr);

  /* Update the context after the matrix building */

  ctx->m11 = cs_cdo_system_get_matrix(sh, 0);
  ctx->b11_max_size = CS_MAX(cs_matrix_get_n_columns(ctx->m11),
                             solver->n1_scatter_dofs);

  /* Prepare the solution array at faces. It has to be allocated to a greater
   * size in case of parallelism in order to allow for a correct matrix-vector
   * product */

  cs_real_t *x1 = nullptr;

  if (cs_glob_n_ranks > 1) {
    BFT_MALLOC(x1, ctx->b11_max_size, cs_real_t);
    cs_array_real_copy(solver->n1_scatter_dofs, flux, x1);
  }
  else
    x1 = flux;

  /* Prepare the context structure according to the choice of block
     preconditioner. In particular, define the Schur complement approximation
     if needed */

  const cs_param_sles_t  *schur_slesp = saddlep->schur_sles_param;
  int  n_xtra_iters = 0;

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);

    ctx->schur_matrix =
      _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                 ctx->m11_inv_diag,
                                 &(ctx->schur_diag),
                                 &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
    {
      cs_real_t  *m11_inv_lumped =
        cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

      algo_ctx->n_inner_iter += n_xtra_iters;

      ctx->schur_matrix =
        _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                   m11_inv_lumped,
                                   &(ctx->schur_diag),
                                   &(ctx->schur_xtra));

      BFT_FREE(m11_inv_lumped);
    }
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    ctx->m22_mass_diag = _get_scaled_diag_m22(ctx->pty_22);
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
    ctx->m22_mass_diag = _get_scaled_diag_m22(ctx->pty_22);

    ctx->m11_inv_diag = cs_saddle_system_b11_inv_diag(ctx->b11_max_size, sh);
    ctx->schur_matrix =
      _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                 ctx->m11_inv_diag,
                                 &(ctx->schur_diag),
                                 &(ctx->schur_xtra));
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      cs_real_t  *m11_inv_lumped =
        cs_saddle_solver_m11_inv_lumped(solver,
                                        ctx->m11,
                                        ctx->b11_range_set,
                                        ctx->xtra_sles,
                                        &n_xtra_iters);

      algo_ctx->n_inner_iter += n_xtra_iters;

      ctx->schur_matrix =
        _schur_approx_diag_inv_m11(schur_slesp->solver_class,
                                   m11_inv_lumped,
                                   &(ctx->schur_diag),
                                   &(ctx->schur_xtra));

      ctx->m22_mass_diag = _get_scaled_diag_m22(ctx->pty_22);

      BFT_FREE(m11_inv_lumped);
    }
    break;

  default:
    /* Do nothing else */
    break;

  }

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_MINRES:
    cs_saddle_solver_minres(solver, x1, pot);
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
    cs_saddle_solver_gcr(solver, x1, pot);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid saddle solver", __func__);
    break;
  }

  /* Copy back to the original array the velocity values at faces */

  if (cs_glob_n_ranks > 1) {
    cs_array_real_copy(solver->n1_scatter_dofs, x1, flux);
    BFT_FREE(x1);
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
 * \brief Solve the saddle-point linear system arising from the discretization
 *        of the scalar-valued CDO cell-based scheme.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_full_system(cs_saddle_solver_t  *solver,
                                 cs_real_t           *flux,
                                 cs_real_t           *pot)
{
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

  cs_saddle_solver_sles_full_system(solver, flux, pot);

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
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Golub-Kahan Bidiagonalization algorithm.
 *        In-house implementation. The PETSc implementation is also available
 *        but appears less efficient in our tests.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_gkb_inhouse(cs_saddle_solver_t  *solver,
                                 cs_real_t           *flux,
                                 cs_real_t           *pot)
{
  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  /* Sanity checks */

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_GKB)
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

  /* 1. Build the GKB context */
  /* ------------------------- */

  cs_saddle_solver_context_gkb_t *ctx
    = (cs_saddle_solver_context_gkb_t *)solver->context;
  assert(ctx != nullptr);

  _gkb_init_context(solver, ctx);

  /* 2. Solve the saddle-point problem */
  /* --------------------------------- */

  cs_saddle_solver_gkb_inhouse(solver, flux, pot);

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
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Notay's algebraic transformation.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_notay(cs_saddle_solver_t  *solver,
                           cs_real_t           *flux,
                           cs_real_t           *pot)
{
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
  assert(sh->n_blocks == 1);

  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
#endif

  /* Solve the saddle-point problem */

  cs_saddle_solver_notay(solver, flux, pot);

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
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Uzawa-CG algorithm.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_uzawa_cg(cs_saddle_solver_t  *solver,
                              cs_real_t           *flux,
                              cs_real_t           *pot)
{
  if (solver == nullptr)
    return 0;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_UZAWA_CG)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Uzawa-CG algorithm is expected.\n"
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

  /* 0. Partial initialization of the context */
  /* ---------------------------------------- */

  cs_saddle_solver_context_uzawa_cg_t *ctx
    = (cs_saddle_solver_context_uzawa_cg_t *)solver->context;
  assert(ctx != nullptr);

  /* Prepare the solution array at faces. It has to be allocated to a greater
   * size in case of parallelism in order to allow for a correct matrix-vector
   * product */

  cs_real_t *x1 = nullptr;

  if (cs_glob_n_ranks > 1) {
    BFT_MALLOC(x1, ctx->b11_max_size, cs_real_t);
    cs_array_real_copy(solver->n1_scatter_dofs, flux, x1);
  }
  else
    x1 = flux;

  /* 1. Build the schur approximation */
  /* -------------------------------- */

  /* TODO */

  /* 2. Solve the saddle-point system */
  /* -------------------------------- */

  cs_saddle_solver_uzawa_cg(solver, x1, pot);

  /* Copy back to the original array the velocity values at faces */

  if (cs_glob_n_ranks > 1) {
    cs_array_real_copy(solver->n1_scatter_dofs, x1, flux);
    BFT_FREE(x1);
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
END_C_DECLS
