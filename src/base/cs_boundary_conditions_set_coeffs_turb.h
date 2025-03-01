#ifndef __CS_BOUNDARY_CONDITIONS_SET_COEFFS_TURB__
#define __CS_BOUNDARY_CONDITIONS_SET_COEFFS_TURB__

/*============================================================================
 * Wall boundary condition management.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Boundary conditions for smooth walls (icodcl = 5).
 *
 * The wall functions may change the value of the diffusive flux.

 * The values at a boundary face \f$ \fib \f$ stored in the face center
 * \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
 * are written as:
 * \f[
 * P_{\face} = A_P^g + B_P^g P_{\centi}
 * \f]
 * and
 * \f[
 * Q_{\face} = A_P^f + B_P^f P_{\centi}
 * \f]
 * where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
 * neighboring cell.

 * Warning:

 * - For a vector field such as the velocity \f$ \vect{u} \f$ the boundary
 *   conditions may read:
 *   \f[
 *   \vect{u}_{\face} = \vect{A}_u^g + \tens{B}_u^g \vect{u}_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \vect{Q}_{\face} = \vect{A}_u^f + \tens{B}_u^f \vect{u}_{\centi}
 *   \f]
 *   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
 *   which coupled velocity components next to a boundary.

 * Please refer to the
 * <a href="../../theory.pdf#wallboundary"><b>wall boundary conditions</b></a>
 * section of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#clptur"><b>clptur</b></a> section.

 * \param[in]     isvhb         id of field whose exchange coeffient should be
 *                               saved at the walls, or -1.
 * \param[in]     velipb        value of the velocity at \f$ \centip \f$
 *                              of boundary cells
 * \param[in]     rijipb        value of \f$ R_{ij} \f$ at \f$ \centip \f$
 *                              of boundary cells
 * \param[out]    visvdr        dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]    hbord         exchange coefficient at boundary
 * \param[in]     theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_turb(int        isvhb,
                                       cs_real_t  velipb[][3],
                                       cs_real_t  rijipb[][6],
                                       cs_real_t  visvdr[],
                                       cs_real_t  hbord[],
                                       cs_real_t  theipb[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_SET_COEFFS_TURB__ */
