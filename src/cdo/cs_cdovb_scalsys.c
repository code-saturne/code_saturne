/*============================================================================
 * Build an algebraic CDO vertex-based system of equations. These equations
 * corresponds to scalar-valued unsteady convection diffusion reaction
 * equations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_assembly.h"
#include "cs_cdo_system.h"
#include "cs_cdovb_priv.h"
#include "cs_cdovb_scaleq.h"
#include "cs_matrix.h"
#include "cs_parameters.h"
#include "cs_sles.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_scalsys.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdovb_scalsys.c

  \brief Build an algebraic CDO vertex-based system for a set of coupled
         unsteady convection-diffusion-reaction of scalar-valued equations with
         source terms
*/

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_SCALSYS_DBG     0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme.
 *
 * \param[in]      cur2prev     do a "current to previous" operation ?
 * \param[in]      n_eqs        number of equations
 * \param[in, out] blocks       array of the core members for an equation
 * \param[in, out] scalsys      pointer to a cs_cdovb_scalsys_t structure
 * \param[in, out] fields       array of pointers to the associated fields
 * \param[in, out] sh           pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_scalsys_build_t)(bool                           cur2prev,
                           int                            n_equations,
                           cs_equation_core_t           **blocks,
                           cs_cdovb_scalsys_t            *scalsys,
                           cs_field_t                   **fields,
                           cs_cdo_system_helper_t        *sh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function pointer to perform the assembly step.
 *        Add the block attached to the block (row_id, col_id) in the full
 *        (coupled) system
 *
 * \param[in]      row_id     id of the row in the coupled system
 * \param[in]      csys       pointer to a cellwise view of the system
 * \param[in, out] sh         pointer to the system helper of the coupled sys.
 * \param[in, out] eqc        context for this kind of discretization
 * \param[in, out] asb        pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_scalsys_asb_t)(int                           row_id,
                         const cs_cell_sys_t          *csys,
                         cs_cdo_system_helper_t       *sh,
                         cs_equation_builder_t        *eqb,
                         cs_cdo_assembly_t            *asb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising for a system of equations relying on
 *         scalar-valued CDO vertex-based schemes
 *
 * \param[in]      n_eqs     number of equations constituting the system
 * \param[in]      n_dofs    local number of DoFs (may be != n_gather_elts)
 * \param[in]      sysp      parameter settings
 * \param[in, out] sh        pointer to the system helper structure
 * \param[in, out] fields    array of field pointers (one for each eq.)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

typedef int
(cs_cdovb_scalsys_solve_t)(int                                 n_eqs,
                           cs_lnum_t                           n_dofs,
                           const cs_equation_system_param_t   *sysp,
                           cs_cdo_system_helper_t             *sh,
                           cs_field_t                        **fields);

/*=============================================================================
 * Structure definitions
 *============================================================================*/

/* Structure related to system of equations for CDO vertex-based
   discretization */

struct _cs_cdovb_scalsys_t {

  /* @name General information
   * @{
   *
   * \var n_dofs
   * Total number of degrees of freedom for this system
   */

  cs_lnum_t                      n_dofs;

  /*!
   * @}
   * @name Build stage
   * Additional members which corresponds to function pointers
   * @{
   */

  cs_cdovb_scalsys_build_t          *build;

  /*!
   * @}
   * @name Assembly stage
   * Additional members which may be used to assemble the system
   * @{
   */

  /* \var assemble
   * Function pointer to manage the assembly process for the Navier-Stokes
   * system of equation
   */

  cs_cdovb_scalsys_asb_t            *assemble;

  /*!
   * @}
   * @name Solve stage
   * Additional members which may be used to solve the system
   * @{
   *
   * \var solve
   * Function dedicated to the resolution of the linear system
   */

  cs_cdovb_scalsys_solve_t          *solve;

  /*! @} */
};


/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_mesh_t              *cs_shared_mesh;
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fill the dof_val array with values collected from field values
 *         If the array is not allocated then one allocates the array (one
 *         takes into account the size needed for parallel synchronizations)
 *
 * \param[in]      n_eqs     number of equations constituting the system
 * \param[in]      sh        pointer to the system helper structure
 * \param[in]      fields    array of field pointers (one for each eq.)
 * \param[in, out] p_dof_vals  double pointer to the array to fill
 */
/*----------------------------------------------------------------------------*/

static void
_set_dof_vals(int                                 n_eqs,
              const cs_cdo_system_helper_t       *sh,
              cs_field_t                       **fields,
              cs_real_t                        **p_dof_vals)
{
  assert(sh != NULL);
  assert(sh->n_blocks == 1);
  assert(sh->blocks[0]->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);

  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  assert(n_vertices * n_eqs <= n_cols);

  /* Initialize the solution array */

  cs_real_t  *dof_vals = *p_dof_vals;

  if (dof_vals == NULL)
    BFT_MALLOC(dof_vals, n_cols, cs_real_t);

  for (int i = 0; i < n_eqs; i++) {

    const cs_field_t  *f = fields[i];
    assert(f != NULL);
    assert(f->val != NULL);

    memcpy(dof_vals + i*n_vertices, f->val, sizeof(cs_real_t)*n_vertices);

  }

  *p_dof_vals = dof_vals;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fill the dof_val array with values collected from field values
 *         If the array is not allocated then one allocates the array (one
 *         takes into account the size needed for parallel synchronizations)
 *
 * \param[in]      n_eqs     number of equations constituting the system
 * \param[in]      dof_vals  pointer to the array with values to copy
 * \param[in, out] fields    array of field pointers (one for each eq.) to fill
 */
/*----------------------------------------------------------------------------*/

static void
_set_field_vals(int                                n_eqs,
                const cs_real_t                   *dof_vals,
                cs_field_t                       **fields)
{
  if (n_eqs < 1)
    return;

  assert(dof_vals != NULL && fields != NULL);

  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  for (int i = 0; i < n_eqs; i++) {

    cs_field_t  *f = fields[i];

    memcpy(f->val, dof_vals + i*n_vertices, sizeof(cs_real_t)*n_vertices);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step. Add the block (row_id, col_id) in the
 *          full (coupled) system
 *
 * \param[in]      row_id     id of the row in the coupled system
 * \param[in]      csys       pointer to a cellwise view of the system
 * \param[in, out] sh         pointer to the system helper of the coupled sys.
 * \param[in, out] eqc        context for this kind of discretization
 * \param[in, out] asb        pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_svb_one_dblock_assemble(int                           row_id,
                         const cs_cell_sys_t          *csys,
                         cs_cdo_system_helper_t       *sh,
                         cs_equation_builder_t        *eqb,
                         cs_cdo_assembly_t            *asb)
{
  CS_UNUSED(eqb);

  /* Members of the coupled system helper structure */

  assert(sh->n_blocks == 1);
  cs_cdo_system_block_t  *block = sh->blocks[0];
  assert(block->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  cs_cdo_system_dblock_t  *db = block->block_pointer;
  cs_cdo_system_block_info_t  bi = block->info;
  assert(bi.stride > row_id);
  assert(bi.interlaced == false);

  /* Matrix assembly for the cellwise system inside the full system */

  db->slave_assembly_func(csys->mat, csys->dof_ids, db->range_set,
                          asb, db->mav);

  /* RHS assembly */

  cs_real_t  *rhs = sh->rhs + bi.n_elements * row_id;

#if CS_CDO_OMP_SYNC_SECTIONS > 0
# pragma omp critical
  {
    for (int v = 0; v < csys->n_dofs; v++)
      rhs[csys->dof_ids[v]] += csys->rhs[v];
  }
#else  /* Use atomic barrier */
  for (int v = 0; v < csys->n_dofs; v++)
#   pragma omp atomic
    rhs[csys->dof_ids[v]] += csys->rhs[v];
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with scalar-valued
 *         degrees of freedom
 *
 * \param[in]      n_eqs     number of equations constituting the system
 * \param[in]      n_dofs    local number of DoFs (may be != n_gather_elts)
 * \param[in]      sysp      parameter settings
 * \param[in, out] sh        pointer to the system helper structure
 * \param[in, out] fields    array of field pointers (one for each eq.)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

static int
_solve_mumps(int                                 n_eqs,
             cs_lnum_t                           n_dofs,
             const cs_equation_system_param_t   *sysp,
             cs_cdo_system_helper_t             *sh,
             cs_field_t                        **fields)
{
  assert(sh != NULL);
  assert(sh->n_blocks == 1);
  assert(sh->blocks[0]->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  assert(n_dofs == sh->full_rhs_size);

  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

  /* n_cols could be greater than n_dofs = n_equations*n_vertices in case of a
     parallel computation */

  assert(n_dofs <= cs_matrix_get_n_columns(matrix));

  /* Initialize the solution array */

  cs_real_t  *dof_vals = NULL;

  _set_dof_vals(n_eqs, sh, fields, &dof_vals);

  if (rset != NULL) { /* parallel/periodic operations */

    /* Move to a gathered view of the solution and rhs arrays */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE, 1, /* type, stride */
                        dof_vals,     /* in: size = n_dofs (scatter) */
                        dof_vals);    /* out: size = n_gather_dofs */

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks contribute also to a vertex
       owned by the local rank */

    if (rset->ifs != NULL) /* parallel or periodic computations */
      cs_interface_set_sum(rset->ifs,
                           n_dofs, 1, false, /* size, stride, interlaced */
                           CS_REAL_TYPE,
                           sh->rhs);

    cs_range_set_gather(rset,
                        CS_REAL_TYPE, 1, /* type, stride */
                        sh->rhs,         /* in: size = n_dofs */
                        sh->rhs);        /* out: size = n_gather_dofs */

  } /* scatter --> gather view + parallel sync. */

  /* Retrieve the SLES structure */

  cs_sles_t  *sles = cs_sles_find_or_add(-1, sysp->name);

  /* Set the input monitoring state */

  cs_solving_info_t  sinfo = {.n_it = 0, .rhs_norm = 1, .res_norm = 1e16};
  cs_real_t  eps = 1e-6;        /* useless in case of a direct solver */

  /* Solve the system as a scalar-valued system of size n_dofs */

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    sh->rhs,
                                                    dof_vals,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */

  if (sysp->linear_solver.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d>"
                  " n_iter %3d | res.norm % -8.4e | rhs.norm % -8.4e\n",
                  sysp->name, code, sinfo.n_it, sinfo.res_norm, sinfo.rhs_norm);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       dof_vals, dof_vals);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       sh->rhs, sh->rhs);

  /* dof_vals --> fields */

  _set_field_vals(n_eqs, dof_vals, fields);

  BFT_FREE(dof_vals);
  cs_sles_free(sles);

  return sinfo.n_it;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme.
 *
 * \param[in]      cur2prev  do a "current to previous" operation ?
 * \param[in]      n_eqs     number of equations
 * \param[in, out] blocks    array of the core members for an equation
 * \param[in, out] scalsys   pointer to a structure cast on-the-fly
 * \param[in, out] fields    array of pointers to the associated fields
 * \param[in, out] sh        pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdovb_scalsys_build_implicit(bool                           cur2prev,
                              int                            n_equations,
                              cs_equation_core_t           **blocks,
                              cs_cdovb_scalsys_t            *scalsys,
                              cs_field_t                   **fields,
                              cs_cdo_system_helper_t        *sh)
{
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  double  rhs_norm = 0.;

  /* Default initialization of properties associated to each block of the
     system */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    for (int j_eq = 0; j_eq < n_equations; j_eq++) {

      int ij = i_eq*n_equations + j_eq;

      cs_equation_core_t  *block_ij = blocks[ij];

      const cs_equation_param_t  *eqp = block_ij->param;
      const cs_real_t  *f_val =
        cur2prev ? fields[j_eq]->val : fields[j_eq]->val_pre;

      cs_equation_builder_t  *eqb = block_ij->builder;
      cs_cdovb_scaleq_t  *eqc = block_ij->scheme_context;

      /* Setup stage: Set useful arrays:
       * -----------
       * -> the Dirichlet values at vertices
       * -> the translation of the enforcement values at vertices if needed
       */

      if (eqb->init_step) {

        cs_cdovb_scaleq_setup(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag);

        eqb->init_step = false;

      }

#     pragma omp parallel if (quant->n_cells > CS_THR_MIN)
      {
        /* Set variables and structures inside the OMP section so that each
           thread has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
        int  t_id = omp_get_thread_num();
#else
        int  t_id = 0;
#endif

        cs_cell_builder_t  *cb = NULL;
        cs_cell_sys_t *csys = NULL;

        cs_cdovb_scaleq_get(&csys, &cb);

        cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);

        cs_cdo_assembly_set_shift(asb,
                                  i_eq * n_vertices,  /* row shift */
                                  j_eq * n_vertices); /* col shift */

        cs_cdovb_scaleq_init_properties(t_id, time_eval, eqp, eqb, eqc);

        /* --------------------------------------------- */
        /* Main loop on cells to build the linear system */
        /* --------------------------------------------- */

#       pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

          /* Build the cellwise system */

          rhs_norm += cs_cdovb_scaleq_build_block_implicit(t_id, c_id, f_val,
                                                           i_eq, j_eq,
                                                           eqp,
                                                           eqb,
                                                           eqc,
                                                           cb, csys);

          /* Assembly process
           * ================ */

          scalsys->assemble(i_eq, csys, sh, eqb, asb);

        } /* Main loop on cells */

      } /* OPENMP Block */

    } /* j_eq */
  } /* i_eq */

  cs_cdo_system_helper_finalize_assembly(sh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a cs_cdovb_scalsys_t structure
 *
 * \param[in]      n_eqs            number of equations
 * \param[in]      sysp             set of parameters to specify a system of eqs
 *
 * \return a pointer to a new allocated system context structure
 */
/*----------------------------------------------------------------------------*/

static cs_cdovb_scalsys_t *
_create_scalsys(int                                 n_eqs,
                const cs_equation_system_param_t   *sysp)
{
  if (n_eqs == 0)
    return NULL;

  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  cs_cdovb_scalsys_t  *scalsys = NULL;

  BFT_MALLOC(scalsys, 1, cs_cdovb_scalsys_t);

  scalsys->n_dofs = n_vertices * n_eqs;

  /* Set pointers to function */

  scalsys->build = _cdovb_scalsys_build_implicit;
  scalsys->assemble = NULL;
  scalsys->solve = NULL;

  switch (sysp->sles_strategy) {

  case CS_EQUATION_SYSTEM_SLES_MUMPS:
    scalsys->assemble = _svb_one_dblock_assemble;
    scalsys->solve = _solve_mumps;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy to solve the system.\n",
              __func__);

  } /* End of switch */

  return scalsys;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to the main shared structures
 *
 * \param[in]  mesh        basic mesh structure
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_init_sharing(const cs_mesh_t              *mesh,
                              const cs_cdo_connect_t       *connect,
                              const cs_cdo_quantities_t    *quant,
                              const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */

  cs_shared_mesh = mesh;
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize factories for extra-diagonal blocks
 *        Build equation builders and scheme context structures for each
 *        equation which are in the extra-diagonal blocks related to a system
 *        of equations. Structures associated to diagonal blocks should be
 *        already initialized during the treatment of the classical full
 *        equations.
 *
 *        Case of scalar-valued CDO-Vb scheme in each block
 *
 * \param[in]      n_eqs            number of equations
 * \param[in]      sysp             set of parameters to specify a system of eqs
 * \param[in, out] block_factories  array of the core members for an equation
 * \param[out]     sh               system helper structure to initialize
 *
 * \return a pointer to a new allocated system context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scalsys_define(int                                 n_eqs,
                        const cs_equation_system_param_t   *sysp,
                        cs_equation_core_t                **block_factories,
                        cs_cdo_system_helper_t            **p_sh)
{
  if (n_eqs == 0)
    return NULL;

  if (block_factories == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Array of structures is not allocated.\n", __func__);
  if (*p_sh != NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: System helper should be not allocated.\n", __func__);

  assert(sysp != NULL);

  /* Create and initialize the system helper for the system of equations
   * A similar initialization is done for each block during the
   * initialization of the context structure */

  cs_cdo_system_helper_t  *sh = NULL;

  /* Up to now, only one full block with all equations is built. More complex
   * cases are possible according to the type of linear solvers. Work in
   * progress. */

  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  switch (sysp->sles_strategy) {

  case CS_EQUATION_SYSTEM_SLES_MUMPS:
    {
      cs_lnum_t  col_block_sizes = n_eqs * n_vertices;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_COUPLED,
                                       1,   /* n_col_blocks */
                                       &col_block_sizes,
                                       1);  /* n_blocks */

      cs_cdo_system_add_dblock(sh, 0,
                               CS_CDO_SYSTEM_MATRIX_CS,
                               cs_flag_primal_vtx,
                               n_vertices,
                               n_eqs,   /* stride */
                               false,   /* interlaced */
                               true);   /* unrolled */

      cs_cdo_system_build_block(sh, 0); /* build/set structures */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid SLES strategy for system \"%s\"\n",
              __func__, sysp->name);
    break;

  } /* Switch on the SLES strategy */

#if 0
  cs_lnum_t  *col_block_sizes = NULL;
  BFT_MALLOC(col_block_sizes, n_eqs, cs_lnum_t);
  for (int i = 0; i < n_eqs; i++)
    col_block_sizes[i] = n_vertices;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_COUPLED,
                                   n_eqs,            /* n_col_blocks */
                                   &col_block_sizes,
                                   n_eqs*n_eqs);     /* n_blocks */
#endif

  /* Initialize the builder and the context structure for each extra-diagonal
     block */

  for (int j = 0; j < n_eqs; j++) { /* Loop on columns */

    /* Remark: Diagonal blocks should be already assigned */

    cs_equation_core_t  *block_jj = block_factories[j*n_eqs+j];

    assert(block_jj != NULL);

    for (int i = 0; i < n_eqs; i++) { /* Loop on rows */

      if (i != j) { /* Extra-diagonal block */

        int  ij = i*n_eqs+j;
        cs_equation_core_t  *block_ij = block_factories[ij];

        /* Copy the boundary conditions of the variable associated to the
           diagonal block j */

        cs_equation_param_copy_bc(block_jj->param, block_ij->param);

        /* Define the equation builder */

        cs_equation_builder_t  *eqb =
          cs_equation_builder_create(block_ij->param, cs_shared_mesh);

        cs_cdovb_scaleq_t  *eqc =
          cs_cdovb_scaleq_init_context(block_ij->param,
                                       -1,              /* No field */
                                       -1,              /* No field */
                                       eqb);

        block_ij->builder = eqb;
        block_ij->scheme_context = eqc;

      } /* i != j */

    } /* row i */

  } /* column j */

  /* Create and initialize the context for the coupled system */

  cs_cdovb_scalsys_t  *scalsys = _create_scalsys(n_eqs, sysp);

  /* Return the allocated and initialized system helper */

  *p_sh = sh;

  return scalsys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free an array of structures (equation parameters, equation builders
 *        or scheme context) for each equation which are in the extra-diagonal
 *        blocks related to a system of equations. Structures associated to
 *        diagonal blocks are freed during the treatment of the classical full
 *        equations.
 *
 *        Case of scalar-valued CDO-Vb scheme in each block
 *
 * \param[in]      n_eqs        number of equations
 * \param[in, out] blocks       array of the core structures for an equation
 * \param[in, out] sys_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scalsys_free(int                     n_eqs,
                      cs_equation_core_t    **blocks,
                      void                   *sys_context)
{
  if (n_eqs == 0)
    return NULL;

  if (blocks == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Structure not allocated\n", __func__);

  cs_cdovb_scalsys_t  *scalsys = sys_context;

  /* Free the extra-diagonal cs_equation_param_t structures, builders and
     scheme context structures */

  for (int i = 0; i < n_eqs; i++) {
    for (int j = 0; j < n_eqs; j++) {

      cs_equation_core_t  *block_ij = blocks[i*n_eqs+j];

      if (i != j) {

        block_ij->param = cs_equation_param_free(block_ij->param);

        cs_equation_builder_free(&(block_ij->builder));

        block_ij->scheme_context =
          cs_cdovb_scaleq_free_context(block_ij->scheme_context);

      }

      BFT_FREE(block_ij);

    } /* Loop on equations (j) */
  } /* Loop on equations (i) */

  /* Free the system context */
  /* ----------------------- */

  BFT_FREE(scalsys);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build and solve the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme using an implicit time scheme.
 *
 * \param[in]      cur2prev     do a "current to previous" operation ?
 * \param[in]      n_eqs        number of equations
 * \param[in]      sysp         set of paremeters for the system of equations
 * \param[in, out] blocks       array of the core members for an equation
 * \param[in, out] sys_context  pointer to a structure cast on-the-fly
 * \param[in, out] sh           pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_solve_implicit(bool                           cur2prev,
                                int                            n_equations,
                                cs_equation_system_param_t    *sysp,
                                cs_equation_core_t           **blocks,
                                void                          *sys_context,
                                cs_cdo_system_helper_t        *sh)
{
  assert(sysp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(sysp->block_var_dim == 1);

  cs_cdovb_scalsys_t  *scalsys = sys_context;

  const cs_lnum_t  n_dofs = scalsys->n_dofs;

  /* Retrieve the field associated to each diagonal block */

  cs_field_t  **fields = NULL;
  BFT_MALLOC(fields, n_equations, cs_field_t *);

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    cs_equation_core_t  *block_ii = blocks[i_eq*n_equations + i_eq];
    cs_cdovb_scaleq_t  *eqc = block_ii->scheme_context;

    fields[i_eq] = cs_field_by_id(eqc->var_field_id);

  }

  /* Initialize the algebraic structures
   * ->  rhs, matrix and assembler values */

  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Build the coupled system of equations */
  /* ------------------------------------- */

  scalsys->build(cur2prev, n_equations, blocks, scalsys, fields, sh);

  /* Reset builder structures and operate a current to previous op. if needed */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {
    for (int j_eq = 0; j_eq < n_equations; j_eq++) {

      cs_equation_core_t  *block_ij = blocks[i_eq*n_equations + j_eq];

      cs_equation_builder_reset(block_ij->builder);

    } /* Loop on blocks corresponding to the column of the system  */
  } /* Loop on blocks corresponding to the row of the system  */

  /* Copy current field values to previous values */

  if (cur2prev)
    for (int i_eq = 0; i_eq < n_equations; i_eq++)
      cs_field_current_to_previous(fields[i_eq]);

  /* Solve the linear system */
  /* ----------------------- */

  scalsys->solve(n_equations, n_dofs, sysp, sh, fields);

  /* Free temporary buffers and structures */

  BFT_FREE(fields);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
