#ifndef __CS_VERTEX_TO_CELL_H__
#define __CS_VERTEX_TO_CELL_H__

/*============================================================================
 * Vertex to cell interpolation.
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

  CS_VERTEX_TO_CELL_UNWEIGHTED,          /*!< Uniform (constant) weights */
  CS_VERTEX_TO_CELL_SHEPARD,             /*!< Shepard interpolation
                                           (weights by inverse distance) */
  CS_VERTEX_TO_CELL_LR                   /*!< Linear regression
                                           (least-squares) */

} cs_vertex_to_cell_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for cell to vertex methods */

extern const char *cs_vertex_to_cell_type_name[];

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
cs_vertex_to_cell_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate vertex values to cell values.
 *
 * \param[in]       method            interpolation method
 * \param[in]       verbosity         verbosity level
 * \param[in]       var_dim           variable dimension
 * \param[in]       v_weight          vertex weight, or NULL
 * \param[in]       v_var             base vertex-based variable
 * \param[out]      c_var             cell-based variable
 */
/*----------------------------------------------------------------------------*/

void
cs_vertex_to_cell(cs_vertex_to_cell_type_t   method,
                  int                        verbosity,
                  cs_lnum_t                  var_dim,
                  const cs_real_t            v_weight[restrict],
                  const cs_real_t            v_var[restrict],
                  cs_real_t                  c_var[restrict]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VERTEX_TO_CELL__ */
