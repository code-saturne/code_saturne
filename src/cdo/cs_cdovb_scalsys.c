/*============================================================================
 * Build an algebraic CDO vertex-based system of equations. These equations
 * corresponds to scalar-valued unsteady convection diffusion reaction
 * equations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

#include "cs_cdovb_priv.h"
#include "cs_cdovb_scaleq.h"
#include "cs_equation_assemble.h"
#include "cs_matrix.h"

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

/* Algebraic system for CDO vertex-based discretization */

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
 * \brief Define a new matrix structure and retrieve (or build) a range set
 *        structure
 *
 * \param[in]      n_vtx_dofs  number of DoFs per vertex
 * \param[in, out] p_rs        double pointer to a range set to define
 * \param[in, out] p_ms        double pointer to a matrix structure to build
 */
/*----------------------------------------------------------------------------*/

static void
_init_algebraic_structures(int                      n_vtx_dofs,
                           cs_range_set_t         **p_rs,
                           cs_matrix_structure_t  **p_ms)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *v2v = connect->v2v;
  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  assert(n_vtx_dofs > 1);

  /* 1. Build a range set */

  if (n_vtx_dofs == 3) { /* Perhaps already defined */

    if (connect->range_sets[CS_DOF_VTX_VECT] != NULL)
      *p_rs = connect->range_sets[CS_DOF_VTX_VECT];

  }
  else {

    cs_interface_set_t  *ifs = NULL; /* temporary */

    cs_cdo_connect_assign_vtx_ifs_rs(cs_shared_mesh,
                                     n_vtx_dofs,
                                     &ifs,
                                     p_rs);

  }

  /* 2. Intermediate structure */

  cs_matrix_assembler_t  *ma =
    cs_equation_build_matrix_assembler(n_vertices, n_vtx_dofs, v2v, *p_rs);

  /* 3. Build the matrix structure */

  cs_matrix_structure_t  *ms =
    cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

  /* Set the pointer to the structure to return */

  *p_ms = ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step for scalar-valued CDO Vb systems
 *
 * \param[in]      row_shift  shift to apply to local row numbering
 * \param[in]      eqc        context for this kind of discretization
 * \param[in]      csys       pointer to a cellwise view of the system
 * \param[in]      rs         pointer to a cs_range_set_t structure
 * \param[in, out] eqa        pointer to a cs_equation_assemble_t structure
 * \param[in, out] mav        pointer to a cs_matrix_assembler_values_t struct.
 * \param[in, out] rhs        right-hand side array
 */
/*----------------------------------------------------------------------------*/

static void
_svb_sys_assemble(cs_lnum_t                          row_shift,
                  const cs_cdovb_scaleq_t           *eqc,
                  const cs_cell_sys_t               *csys,
                  const cs_range_set_t              *rs,
                  cs_equation_assemble_t            *eqa,
                  cs_matrix_assembler_values_t      *mav,
                  cs_real_t                         *rhs)
{
  /* Matrix assembly */

  eqc->assemble(csys->mat, csys->dof_ids, rs, eqa, mav);

  /* RHS assembly */

#if CS_CDO_OMP_SYNC_SECTIONS > 0
# pragma omp critical
  {
    for (int v = 0; v < csys->n_dofs; v++)
      rhs[csys->dof_ids[v] + row_shift] += csys->rhs[v];
  }
#else  /* Use atomic barrier */
  for (int v = 0; v < csys->n_dofs; v++)
#   pragma omp atomic
    rhs[csys->dof_ids[v] + row_shift] += csys->rhs[v];
#endif
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
cs_cdovb_scalsys_init_common(const cs_mesh_t              *mesh,
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
 * \param[in, out] block_factories  array of the core members for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_init_structures(int                        n_eqs,
                                 cs_equation_core_t       **block_factories)
{
  if (n_eqs == 0)
    return;

  if (block_factories == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Array of structures is not allocated.\n", __func__);

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

        cs_equation_builder_t  *eqb = cs_equation_builder_init(block_ij->param,
                                                               cs_shared_mesh);

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
 * \param[in]      n_eqs    number of equations
 * \param[in, out] blocks   array of the core structures for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_free_structures(int                        n_eqs,
                                 cs_equation_core_t       **blocks)
{
  if (n_eqs == 0)
    return;

  if (blocks == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Structure not allocated\n", __func__);

  /* Free the extra-diagonal cs_equation_param_t structures, builders and
     scheme context structures */

  for (int i = 0; i < n_eqs; i++) {
    for (int j = 0; j < n_eqs; j++) {
      if (i != j) {

        cs_equation_core_t  *block_ij = blocks[i*n_eqs+j];

        block_ij->param = cs_equation_param_free(block_ij->param);

        cs_equation_builder_free(&(block_ij->builder));

        block_ij->scheme_context =
          cs_cdovb_scaleq_free_context(block_ij->scheme_context);

      }
    } /* Loop on equations (j) */
  } /* Loop on equations (i) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build and solve the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme using an implicit time scheme.
 *
 * \param[in]      cur2prev  true="current to previous" operation is performed
 * \param[in]      n_eqs     number of equations
 * \param[in]      sysp      set of paremeters for the system of equations
 * \param[in, out] blocks    array of the core members for an equation
 * \param[in, out] p_ms      double pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_solve_implicit(bool                           cur2prev,
                                int                            n_equations,
                                cs_equation_system_param_t    *sysp,
                                cs_equation_core_t           **blocks,
                                cs_matrix_structure_t        **p_ms)
{
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_lnum_t  n_dofs = n_equations * n_vertices;

  assert(sysp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(sysp->block_var_dim == 1);

  /* Retrieve the field associated to each diagonal block */

  cs_field_t  **fields = NULL;
  BFT_MALLOC(fields, n_equations, cs_field_t *);

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    cs_equation_core_t  *block_ii = blocks[i_eq*n_equations + i_eq];
    cs_cdovb_scaleq_t  *eqc = block_ii->scheme_context;

    fields[i_eq] = cs_field_by_id(eqc->var_field_id);

  }

  /* Set the algebraic structures */
  /* ---------------------------- */

  /* 1. Build the matrix structure if needed */

  cs_matrix_structure_t  *ms = (p_ms == NULL) ? NULL : *p_ms;
  cs_range_set_t  *rs = NULL;

  if (ms == NULL) {
    _init_algebraic_structures(n_equations, &rs, &ms);
    *p_ms = ms;
  }

  /* 3. Create the matrix and its rhs */

  cs_matrix_t  *matrix = cs_matrix_create(ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_dofs, cs_real_t);
# pragma omp parallel for if  (n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_dofs; i++) rhs[i] = 0.0;

  /* 4. Initialize the structure to assemble values */

  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* Setup stage: Set useful arrays:
   * -----------
   * -> the Dirichlet values at vertices
   * -> the translation of the enforcement values at vertices if needed
   */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    cs_equation_core_t  *block_ii = blocks[i_eq*n_equations + i_eq];

    const cs_equation_param_t  *eqp = block_ii->param;;
    cs_equation_builder_t  *eqb = block_ii->builder;
    cs_cdovb_scaleq_t  *eqc = block_ii->scheme_context;

    cs_cdovb_scaleq_setup(ts->t_cur + ts->dt[0],
                          mesh, eqp, eqb, eqc->vtx_bc_flag);

    if (eqb->init_step)
      eqb->init_step = false;

  } /* Loop on equations (diagonal blocks) */

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

  double  rhs_norm = 0.;

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_cell_builder_t  *cb = NULL;
    cs_cell_sys_t *csys = NULL;

    cs_cdovb_scaleq_get(&csys, &cb);

    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

    const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

    /* Default initialization of properties associated to each block of the
       system */

    for (int i_eq = 0; i_eq < n_equations; i_eq++) {

      for (int j_eq = 0; j_eq < n_equations; j_eq++) {

        int ij = i_eq*n_equations + j_eq;

        cs_equation_core_t  *block_ij = blocks[ij];

        cs_equation_assemble_set_shift(eqa,
                                       i_eq * n_vertices,  /* row shift */
                                       j_eq * n_vertices); /* col shift */

        const cs_equation_param_t  *eqp = block_ij->param;
        const cs_real_t  *f_val = fields[j_eq]->val;

        cs_equation_builder_t  *eqb = block_ij->builder;
        cs_cdovb_scaleq_t  *eqc = block_ij->scheme_context;

        cs_cdovb_scaleq_init_properties(t_id, time_eval, eqp, eqb, eqc);

        /* --------------------------------------------- */
        /* Main loop on cells to build the linear system */
        /* --------------------------------------------- */

#       pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

          /* Build the cellwise system */

          rhs_norm += cs_cdovb_scaleq_cw_build_implicit(t_id, c_id, f_val,
                                                        eqp,
                                                        eqb,
                                                        eqc,
                                                        cb, csys);

          /* Assembly process
           * ================ */

          _svb_sys_assemble(i_eq*n_vertices, /* row shift (for RHS assembly) */
                            eqc, csys, rs, eqa, mav, rhs);

        } /* Main loop on cells */

      } /* j_eq */
    } /* i_eq */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    for (int j_eq = 0; j_eq < n_equations; j_eq++) {

      cs_equation_core_t  *block_ij = blocks[i_eq*n_equations + j_eq];

      cs_equation_builder_reset(block_ij->builder);

    }

    /* Copy current field values to previous values */

    if (cur2prev)
      cs_field_current_to_previous(fields[i_eq]);

  }

  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building
   * Begin the solving step
   * -------------------------- */

  /* TODO */

  /* Free temporary buffers and structures */

  BFT_FREE(rhs);
  BFT_FREE(fields);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
