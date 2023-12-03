#ifndef __CS_SOLID_ZONE_H__
#define __CS_SOLID_ZONE_H__

/*============================================================================
 * Solid zones handling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_base.h"
#include "cs_mesh_location.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief build solid flag for mesh cells.
 *
 * If no solid cells are present, NULL is returned.

 * If non-empty, the caller is responsible for freeing the flag
 *
 * \param[in]   m         pointer to mesh
 *
 * \return  solid cell flag array (0 is fluid, 1 if solid), or NULL
 */
/*----------------------------------------------------------------------------*/

int *
cs_solid_zone_flag(const cs_mesh_t   *m);

/*----------------------------------------------------------------------------*/
/*
 * \brief Zero an array on cells of a solid zone
 *
 * \param[in]   stride  array stride
 * \param[out]  a       array of cell values
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_zone_set_zero_on_cells(int         stride,
                                cs_real_t  *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant scalar value to cells of a solid zone.
 *
 * \param[in]   ref_val  reference value
 * \param[out]  a        array of cell values
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_zone_set_scalar_on_cells(cs_real_t  ref_val,
                                  cs_real_t  a[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLID_ZONE_H__ */
