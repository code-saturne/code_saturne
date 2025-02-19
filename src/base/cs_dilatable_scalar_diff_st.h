#ifndef __CS_DILATABLE_SCALAR_DIFF_ST_H__
#define __CS_DILATABLE_SCALAR_DIFF_ST_H__

/*============================================================================
 * Weakly compressible algorithm (semi-analytic):
   Computation of scalar diffusion terms.
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

/*----------------------------------------------------------------------------
 * Standard library headers
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
 * \brief Weakly compressible algorithm (semi-analytic):
          Computation of scalar diffusion terms.
 *
 * \param[in]  iterns      Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_dilatable_scalar_diff_st(int iterns);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DILATABLE_SCALAR_DIFF_ST_H__ */
