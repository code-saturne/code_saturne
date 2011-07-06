/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_MESH_THINWALL_H__
#define __CS_MESH_THINWALL_H__

/*============================================================================
 * Insert thin walls into the mesh.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * parameters:
 *   mesh           <->  pointer to mesh structure to modify
 *   face_list      <-> list of selected (interior) faces (1 to n)
 *   face_list_size <-> number of selected (interior) faces
 *----------------------------------------------------------------------------*/

void
cs_create_thinwall(cs_mesh_t   *mesh,
                   fvm_lnum_t  *face_list,
                   int          face_list_size);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_THINWALL_H__ */
