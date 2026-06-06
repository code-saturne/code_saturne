#ifndef CS_GRID_H
#define CS_GRID_H

/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Aggregation algorithm */

typedef enum {

  CS_GRID_COARSENING_DEFAULT,        /*!< default among following choices */
  CS_GRID_COARSENING_SPD_DX,         /*!< SPD, diag/extradiag ratio based */
  CS_GRID_COARSENING_SPD_MX,         /*!< SPD, max extradiag ratio based */
  CS_GRID_COARSENING_SPD_PW,         /*!< SPD, pairwise aggregation */
  CS_GRID_COARSENING_CONV_DIFF_DX    /*!< convection+diffusion,
                                          diag/extradiag ratio based */

} cs_grid_coarsening_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names for coarsening options */

extern const char *cs_grid_coarsening_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Set factor to ensure diagonal dominance.
 *
 * \param[in]  factor  clip margin factor (ignored if < 0).
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_set_diag_dom_clip_factor(double  factor);

/*----------------------------------------------------------------------------
 * Set matrix tuning behavior for multigrid coarse meshes.
 *
 * The finest mesh (level 0) is handled by the default tuning options,
 * so only coarser meshes are considered here.
 *
 * parameters:
 *   fill_type <-- associated matrix fill type
 *   max_level <-- maximum level for which tuning is active
 *----------------------------------------------------------------------------*/

void
cs_grid_set_matrix_tuning(cs_matrix_fill_type_t  fill_type,
                          int                    max_level);

/*----------------------------------------------------------------------------*/

#endif /* CS_GRID_H */
