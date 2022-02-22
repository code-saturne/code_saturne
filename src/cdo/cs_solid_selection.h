#ifndef __CS_SOLID_CELLS_H__
#define __CS_SOLID_CELLS_H__

/*============================================================================
 * Manage the list of solid cells and associated helper functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_cdo_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Structure and type definitions
 *============================================================================*/

typedef struct {

  cs_lnum_t     n_cells;        /*!< local number of solid cells */
  cs_gnum_t     n_g_cells;      /*!< global number of solid cells  */
  cs_lnum_t    *cell_ids;       /*!< local list of solid cells */

  bool         *cell_is_solid;  /*!< true if this is a solid cell */
  bool         *face_is_solid;  /*!< true if the face belongs to a solid cell */

} cs_solid_selection_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the information related to the list of solid cells
 *         If this structure does not exist, there is an initialization.
 *
 * \return a pointer to a cs_solid_selection_t structure
 */
/*----------------------------------------------------------------------------*/

cs_solid_selection_t *
cs_solid_selection_get(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the solid selection
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_selection_sync(const cs_cdo_connect_t   *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure storing the information related to the list of
 *         solid cells.
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_selection_free(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLID_CELLS_H__ */
