#ifndef __CS_BAD_CELLS_REGULARISATION_H__
#define __CS_BAD_CELLS_REGULARISATION_H__

/*============================================================================
 * Divergence operators.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for scalars
 *
 * \param[in, out]  var  variable to regularize.
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_scalar(cs_real_t *var);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for vectors
 *
 * \param[in, out]  var  variable to regularize.
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_vector(cs_real_3_t *var,
                                   int          boundary_projection);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for symmetric tensors.
 *
 * \param[in, out]  var  variable to regularize.
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_sym_tensor(cs_real_6_t *var,
                                       int          boundary_projection);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for tensors
 *
 * \param[in, out]  var  variable to regularize.
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_tensor(cs_real_9_t *var,
                                   int          boundary_projection);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BAD_CELLS_REGULARISATION_H__ */
