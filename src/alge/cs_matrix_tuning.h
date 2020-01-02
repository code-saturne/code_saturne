#ifndef __CS_MATRIX_TUNING_H__
#define __CS_MATRIX_TUNING_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_numbering.h"
#include "cs_halo_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a matrix variant tuned matrix.vector product operations.
 *
 * The variant may later be applied to matrices of the same type and fill type.
 *
 * \param[in]  m           associated matrix
 * \param[in]  verbosity   verbosity level
 * \param[in]  n_measure   minimum number of measuring runs
 * \param[in]  t_measure   minimum measure time
 *
 * \returns  pointer to tuning results structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_tuned(const cs_matrix_t  *m,
                        int                 verbosity,
                        int                 n_measure,
                        double              t_measure);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_TUNING_H__ */
