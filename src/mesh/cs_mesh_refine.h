#ifndef __CS_MESH_REFINE_H__
#define __CS_MESH_REFINE_H__

/*============================================================================
 * Mesh refinement.
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

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine flagged mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       conforming  if true, propagate refinement to ensure
 *                              subdivision is conforming
 * \param[in]       cell_flag   subdivision type for each cell
 *                              (0: none; 1: isotropic)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_refine_simple(cs_mesh_t  *m,
                      bool        conforming,
                      const int   cell_flag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine selected mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       conforming  if true, propagate refinement to ensure
 *                              subdivision is conforming
 * \param[in]       n_cells     number of selected cells
 * \param[in]       cells       list of selected cells (0 to n-1)
 *                              or NULL if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_refine_simple_selected(cs_mesh_t        *m,
                               bool              conforming,
                               cs_lnum_t         n_cells,
                               const cs_lnum_t   cells[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_REFINE_H__ */
