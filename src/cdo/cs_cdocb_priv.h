#ifndef __CS_CDOCB_PRIV_H__
#define __CS_CDOCB_PRIV_H__

/*============================================================================
 * Structure and functions common to all CDO cell-based schemes but not public
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_assembly.h"
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_equation_param.h"
#include "cs_hodge.h"
#include "cs_saddle_solver.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_cdocb_t  cs_cdocb_scaleq_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the assembly stage for a vector-valued system obtained with
 *        CDO-Fb schemes
 *
 * \param[in]      csys  pointer to a cs_cell_sys_t structure
 * \param[in]      cm    pointer to a cs_cell_mesh_t structure
 * \param[in, out] eqc   context structure for a scalar-valued Cb equation
 * \param[in, out] asb   pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdocb_scaleq_assemble_t)(const cs_cell_sys_t   *csys,
                             const cs_cell_mesh_t  *cm,
                             cs_cdocb_scaleq_t     *eqc,
                             cs_cdo_assembly_t     *asb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function prototype to solve the saddle-point linear system
 *        arising from the discretization of the scalar-valued CDO cell-based
 *        scheme.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

typedef int
(cs_cdocb_scaleq_solve_t)(cs_saddle_solver_t  *saddle,
                          cs_real_t           *flux,
                          cs_real_t           *pot);

/* Main structure */
/* ============== */

/* Context structure for CDO cell-based discretizations */
/* ---------------------------------------------------- */

struct _cs_cdocb_t {

  /* Ids related to the variable field and to the boundary flux field */

  int          var_field_id;
  int          bflux_field_id;

  /* System size (n_faces + n_cells) */

  cs_lnum_t    n_faces;
  cs_lnum_t    n_cells;
  cs_lnum_t    n_dofs;

  /* Solution of the algebraic system DoF unknowns (x) + BCs */

  cs_real_t   *flux;       /* At the last iteration */
  cs_real_t   *flux_pre;   /* At the previous iteration */

  /* Reconstructred variable face values */

  cs_real_t   *face_values;
  cs_real_t   *face_values_pre;

  /* Members used for the cell-wise building of the linear system */
  /* ------------------------------------------------------------ */

  /* Array storing the local cell-wise divergence operator (one by thread) */

  cs_real_t                  **div_op_cw;

  /* Pointer of function to build the diffusion term */

  cs_hodge_t                 **diff_hodge;
  cs_hodge_compute_t          *compute_diff_hodge;

  /* Boundary conditions */
  /* ------------------- */

  cs_cdo_enforce_bc_t         *enforce_dirichlet;
  cs_cdo_enforce_bc_t         *enforce_neumann;
  cs_cdo_enforce_bc_t         *enforce_robin_bc;

  /* Linear system */
  /* ------------- */

  /* \var block21_op
   * Unassembled (2,1)-block (related to the divergence). Not always
   * allocated. It depends on the solver for the saddle-point system.
   */

  cs_real_t                   *block21_op;

  /* \var system_helper
   * Set of structure to handle the saddle-point matrix and its rhs
   */

  cs_cdo_system_helper_t      *system_helper;

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

  cs_cdocb_scaleq_assemble_t  *assemble;

  /*!
   * @}
   * @name Solve stage
   * Additional members which may be used to solve the system
   * @{
   *
   * \var solve
   * Function dedicated to the resolution of the saddle-point system
   */

  cs_cdocb_scaleq_solve_t     *solve;

  /* \var saddle_solver
   * Set of pointers to enable the resolution of saddle-point system
   * with various algorithms. This structure allows us to unify the prototype
   * of "solve" functions
   */

  cs_saddle_solver_t          *saddle_solver;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the default settings for a scalar-valued CDO cell-based scheme
 *
 * \param[in, out] eqp  set of equation parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_init_default_param(cs_equation_param_t  *eqp);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOCB_PRIV_H__ */
