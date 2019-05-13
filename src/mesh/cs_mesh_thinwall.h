#ifndef __CS_MESH_THINWALL_H__
#define __CS_MESH_THINWALL_H__

/*============================================================================
 * Insert thin walls into the mesh.
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

#include "fvm_defs.h"

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 * The created faces share vertices.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * \param[in]       mesh           pointer to mesh structure to modify
 * \param[in, out]  face_list      list of selected (interior) faces (0 to n-1)
 * \param[in]       face_list_size number of selected (interior) faces
 */
/*----------------------------------------------------------------------------*/

void
cs_create_thinwall(cs_mesh_t  *mesh,
                   cs_lnum_t  *face_list,
                   cs_lnum_t   face_list_size);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_THINWALL_H__ */
