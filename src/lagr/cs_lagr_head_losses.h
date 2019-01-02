#ifndef __CS_LAGR_HEAD_LOSSES_H__
#define __CS_LAGR_HEAD_LOSSES_H__

/*============================================================================
 * Functions and types for head losses and porosity
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Define Head losses to take into account deposit in the flow
 *
 * \param[in]   n_hl_cells  number of cells on which to apply head losses
 * \param[in]   cell_ids    ids of cells on which to apply head losses
 * \param[in]   bc_type     boundary face type
 * \param[out]  cku         head loss coefficients at matchin cells
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_head_losses(cs_lnum_t        n_hl_cells,
                    const cs_lnum_t  cell_ids[],
                    const cs_lnum_t  bc_type[],
                    cs_real_t        cku[][6]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_HEAD_LOSSES_H__ */
