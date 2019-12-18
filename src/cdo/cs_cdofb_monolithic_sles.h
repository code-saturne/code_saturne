#ifndef __CS_CDOFB_MONOLITHIC_SLES_H__
#define __CS_CDOFB_MONOLITHIC_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_navsto_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Context related to the resolution of a saddle point problem */
typedef struct {

  /* Block matrices: The gradient operator is the -transpose of div_op */

  /*
   * If n_row_blocks = 1 and div_op == NULL then all the velocity-pressure
   * components are gathered in one block. In this case, the full matrix is a
   * 1x1 matrix
   *
   * If n_row_blocks = 1 and div_op != NULL then all the velocity components are
   * gathered in one block. In this case, the full matrix is a 2x2 matrix

   * If n_row_blocks = 3 and div_op == NULL then there is one block dedicated to
   * each velocity component. In this case, the full matrix is a 4x4 matrix
   *
   * Please notice that div_op if not NULL is stored in a non-assembled way.
   */

  int            n_row_blocks;

  /* Blocks related to the velocity momentum */
  cs_matrix_t  **block_matrices;

  cs_real_t     *div_op;    /* Block related to the -divergence (block A_{10} */

  /* Arrays split according to the block shape. U is interlaced or not
   * according to the SLES strategy */

  cs_lnum_t      n_faces;       /* local number of DoFs for each component
                                 * of the velocity */
  cs_lnum_t      n_cells;       /* local number of DoFs for the pressure */

  cs_real_t     *u_f;           /* velocity values at faces */
  cs_real_t     *p_c;           /* pressure values at cells */

  cs_real_t     *b_f;           /* RHS for the momentum (size = 3*n_faces) */
  cs_real_t     *b_c;           /* RHS for the mass equation (size = n_cells) */

  cs_sles_t     *sles;          /* main SLES structure */

  cs_real_t      graddiv_coef;  /* value of the grad-div coefficient in case
                                 * of augmented system */

} cs_cdofb_monolithic_sles_t;

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
cs_cdofb_monolithic_sles_create(void);

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
                              cs_cdofb_monolithic_sles_t   *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset to zero rhs and clean the cs_sles_t structure
 *
 * \param[in, out]  msles   pointer to the structure to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_reset(cs_cdofb_monolithic_sles_t   *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a part of the structure
 *
 * \param[in, out]  msles   pointer to the structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_clean(cs_cdofb_monolithic_sles_t   *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdofb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_free(cs_cdofb_monolithic_sles_t   **p_msles);

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
                                    const cs_range_set_t          *rset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(const cs_navsto_param_t    *nsp,
                             void                       *context);

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
                          cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         Navier-Stokes equation with a CDO face-based approach.
 *         The system is split into blocks to enable more efficient
 *         preconditioning techniques
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_by_blocks_solve(const cs_navsto_param_t       *nsp,
                                    const cs_equation_param_t     *eqp,
                                    cs_cdofb_monolithic_sles_t    *msles);

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
                              cs_cdofb_monolithic_sles_t    *msles);

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
                                   cs_cdofb_monolithic_sles_t    *msles);

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
                                        cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_MONOLITHIC_SLES_H__ */
