#ifndef __CS_CELL_TO_VERTEX_H__
#define __CS_CELL_TO_VERTEX_H__

/*============================================================================
 * Cell to vertex interpolation.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cell to vertex computation method
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_CELL_TO_VERTEX_UNWEIGHTED,          /*!< Uniform (constant) weights */
  CS_CELL_TO_VERTEX_SHEPARD,             /*!< Shepard interpolation
                                           (weights by inverse distance) */
  CS_CELL_TO_VERTEX_LR                   /*!< Linear regression
                                           (least-squares) */

} cs_cell_to_vertex_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for cell to vertex methods */

extern const char *cs_cell_to_vertex_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free cell to vertex interpolation weights.
 *
 * This will force subsequent calls to rebuild those weights if needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_to_vertex_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate cell values to vertex values for a scalar array.
 *
 * \param[in]       method      interpolation method
 * \param[in]       verbosity   verbosity level
 * \param[in]       tr_dim      2 for tensor with periodicity of rotation,
 *                              0 otherwise
 * \param[in]       c_weight    cell weight, or NULL
 * \param[in]       c_var       base cell-based variable
 * \param[in]       b_var       base boundary-face values, or NULL
 * \param[out]      v_var       vertex-based variable
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_to_vertex_scalar(cs_cell_to_vertex_type_t   method,
                         int                        verbosity,
                         int                        tr_dim,
                         const cs_real_t            c_weight[restrict],
                         const cs_real_t            c_var[restrict],
                         const cs_real_t            b_var[restrict],
                         cs_real_t                  v_var[restrict]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CELL_TO_VERTEX__ */
