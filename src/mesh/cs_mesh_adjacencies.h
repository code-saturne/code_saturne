#ifndef __CS_MESH_ADJACENCIES_H__
#define __CS_MESH_ADJACENCIES_H__

/*============================================================================
 * Additional mesh adjacencies.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

typedef struct {

  /* metadata */

  bool        single_faces_to_cells;   /* true if a single face is adjacent
                                          to 2 given cells */

  /* cells -> cells connectivity (standard) */

  cs_lnum_t  *cell_cells_idx;          /* indexes (shared) */
  cs_lnum_t  *cell_cells;              /* adjacency (shared) */

  /* cells -> cells connectivity (extended) */

  const cs_lnum_t  *cell_cells_e_idx;  /* indexes (shared) */
  const cs_lnum_t  *cell_cells_e;      /* adjacency (shared) */

  /* cells -> boundary faces connectivity */

  cs_lnum_t        *cell_b_faces_idx;
  cs_lnum_t        *cell_b_faces;

} cs_mesh_adjacencies_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Read-only pointer to global mesh additional adjacencies structure */

extern const cs_mesh_adjacencies_t  *cs_glob_mesh_adjacencies;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mesh adjacencies helper API.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize mesh adjacencies helper API.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update mesh adjacencies helper API relative to mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_mesh(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update extended cell -> cell connectivites in
 *         mesh adjacencies helper API relative to mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_cell_cells_e(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_ADJACENCIES__ */
