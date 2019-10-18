#ifndef __CS_CDOFB_MONOLITHIC_SLES_H__
#define __CS_CDOFB_MONOLITHIC_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in]      gamma    value of the grad-div coefficient
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] u_f      initial velocity on faces
 * \param[in, out] p_c      initial pressure in cells
 * \param[in, out] b_f      right-hand side (scatter/gather if needed) on faces
 * \param[in, out] b_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_gkb_solve(const cs_matrix_t             *matrix,
                   const cs_navsto_param_t       *nsp,
                   const cs_equation_param_t     *eqp,
                   const cs_real_t               *div_op,
                   cs_real_t                      gamma,
                   cs_sles_t                     *sles,
                   cs_real_t                     *u_f,
                   cs_real_t                     *p_c,
                   cs_real_t                     *b_f,
                   cs_real_t                     *b_c);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_MONOLITHIC_SLES_H__ */
