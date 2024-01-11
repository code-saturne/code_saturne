#ifndef __CS_CDOCB_MONOLITHIC_SLES_H__
#define __CS_CDOCB_MONOLITHIC_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO cell-based schemes with a monolithic velocity-pressure coupling
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdocb_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdocb_monolithic_sles_t structure
 *
 * \param[in] n_faces     number of faces (interior + border)
 * \param[in] n_cells     number of cells
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdocb_monolithic_sles_t *
cs_cdocb_monolithic_sles_create(cs_lnum_t    n_faces,
                                cs_lnum_t    n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a part of the structure
 *
 * \param[in, out]  msles   pointer to the structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_clean(cs_cdocb_monolithic_sles_t   *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdocb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_free(cs_cdocb_monolithic_sles_t   **p_msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  connect  pointer to cdo connectivities
 * \param[in]  quant    pointer to additional mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_init_sharing(const cs_cdo_connect_t        *connect,
                                      const cs_cdo_quantities_t     *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free if needed structure(s) associated CDO cell-based schemes with
 *         a monolithic velocity-pressure coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_set_sles(cs_equation_param_t   *eqp,
                             void                  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         Navier-Stokes equation with a CDO cell-based approach.
 *         The full system is treated as one block and solved as it is.
 *         In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] msles    pointer to a cs_cdocb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_monolithic_solve(const cs_equation_param_t     *eqp,
                          const cs_cdo_system_helper_t  *sh,
                          cs_cdocb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOCB_MONOLITHIC_SLES_H__ */
