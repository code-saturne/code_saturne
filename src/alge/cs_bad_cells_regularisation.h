#ifndef __CS_BAD_CELLS_REGULARISATION_H__
#define __CS_BAD_CELLS_REGULARISATION_H__

/*============================================================================
 * Divergence operators.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_base.h"
#include "cs_halo.h"

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
 * \brief Add comments
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_scalar(cs_real_t *var);

void
cs_bad_cells_regularisation_vector(cs_real_3_t *var,
                                   int         boundary_projection);

void
cs_bad_cells_regularisation_sym_tensor(cs_real_6_t *var,
                                       int         boundary_projection);

void
cs_bad_cells_regularisation_tensor(cs_real_9_t *var,
                                   int         boundary_projection);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BAD_CELLS_REGULARISATION_H__ */
