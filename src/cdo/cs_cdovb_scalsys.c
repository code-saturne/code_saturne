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
 * \brief  Define a new matrix structure
 *
 * \param[in]  n_vtx_dofs  number of DoFs per vertex
 *
 * \return a pointer to the newly allocated matrix structure
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_structure_t *
_init_matrix_structure(int       n_vtx_dofs)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *v2v = connect->v2v;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  assert(rs != NULL);
  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  cs_matrix_assembler_t  *ma =
    cs_equation_build_matrix_assembler(n_vertices, n_vtx_dofs, v2v, rs);

  cs_matrix_structure_t  *ms =
    cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

  return ms;
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
 * \brief Build and solve the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme using an implicit time scheme.
 *
 * \param[in]      cur2prev  true="current to previous" operation is performed
 * \param[in]      n_eqs     number of equations
 * \param[in]      eqps      pointer to a list of cs_equation_param_t struct.
 * \param[in, out] eqbs      pointer to a list of cs_equation_builder_t struct.
 * \param[in, out] eqcs      pointer to a list of cs_cdovb_scaleq_t struct.
 * \param[in, out] p_ms      double pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_solve_implicit(bool                          cur2prev,
                                int                           n_equations,
                                cs_equation_param_t         **eqps,
                                cs_equation_builder_t       **eqbs,
                                void                        **eqcs,
                                cs_matrix_structure_t       **p_ms)
{
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_lnum_t  n_dofs = n_equations * n_vertices;

  /* Retrieve the field associated to each diagonal block */

  cs_field_t  **fields = NULL;
  BFT_MALLOC(fields, n_equations, cs_field_t *);

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    cs_cdovb_scaleq_t  *eqc = eqcs[i_eq*n_equations + i_eq];

    fields[i_eq] = cs_field_by_id(eqc->var_field_id);

  }

  /* Build the matrix structure if needed and then the matrix and its rhs */

  cs_matrix_structure_t  *ms = (p_ms == NULL) ? NULL : *p_ms;

  if (ms == NULL) {
    ms = _init_matrix_structure(n_equations);
    *p_ms = ms;
  }

  cs_matrix_t  *matrix = cs_matrix_create(ms);
  cs_real_t  *rhs = NULL;
  double  rhs_norm = 0.;

  BFT_MALLOC(rhs, n_dofs, cs_real_t);
# pragma omp parallel for if  (n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_dofs; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */

  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* Build an array storing the Dirichlet values at vertices and another one
     to detect vertices with an enforcement */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    int  ii = i_eq*n_equations + i_eq;

    const cs_equation_param_t  *eqp = eqps[ii];

    cs_equation_builder_t  *eqb = eqbs[ii];
    cs_cdovb_scaleq_t  *eqc = eqcs[ii];

    cs_cdovb_scaleq_setup(ts->t_cur + ts->dt[0],
                          mesh, eqp, eqb, eqc->vtx_bc_flag);

    if (eqb->init_step)
      eqb->init_step = false;

  } /* Loop on equations (diagonal blocks) */

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

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

    const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

    /* Default initialization of properties associated to each block of the
       system */

    for (int i_eq = 0; i_eq < n_equations; i_eq++) {

      for (int j_eq = 0; j_eq < n_equations; j_eq++) {

        int ij = i_eq*n_equations + j_eq;

        const cs_equation_param_t  *eqp = eqps[ij];
        const cs_real_t  *f_val = fields[j_eq]->val;

        cs_equation_builder_t  *eqb = eqbs[ij];
        cs_cdovb_scaleq_t  *eqc = eqcs[ij];

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

          /* TODO */

        } /* Main loop on cells */

      } /* j_eq */
    } /* i_eq */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  for (int i_eq = 0; i_eq < n_equations; i_eq++) {

    cs_equation_builder_reset(eqbs[i_eq]);

    /* Copy current field values to previous values */

    if (cur2prev)
      cs_field_current_to_previous(fields[i_eq]);

  }

  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building
   * Begin the solving step
  /* -------------------------- */

  /* TODO */

  /* Free temporary buffers and structures */

  BFT_FREE(rhs);
  BFT_FREE(fields);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
