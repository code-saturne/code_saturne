/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_MESH_SELECT_H__
#define __CS_MESH_SELECT_H__

/*============================================================================
 * Functions dealing with the selection of cs_mesh_t entities
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"

/*============================================================================
 * Type definition
 *============================================================================*/

/* Structure for selection management */

typedef struct _cs_mesh_select_t  cs_mesh_select_t;

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation and initialization of a cs_mesh_select_t structure.
 *
 * parameters:
 *   n_colors    --> number of colors
 *   n_groups    --> number of groups
 *   colors      --> color list
 *   groups      --> group name list
 *   invsel      --> invert selection if CS_TRUE
 *
 * returns:
 *   pointer to created cs_mesh_select_t structure
 *----------------------------------------------------------------------------*/

cs_mesh_select_t *
cs_mesh_select_create(cs_int_t     n_colors,
                      cs_int_t     n_groups,
                      cs_int_t    *colors,
                      char       **groups,
                      cs_bool_t    invsel);

/*----------------------------------------------------------------------------
 * Destroy a cs_mesh_select_t structure
 *
 * parameters:
 *   selection --> pointer to selection structure that should be destroyed
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

cs_mesh_select_t *
cs_mesh_select_destroy(cs_mesh_select_t  *selection);

/*----------------------------------------------------------------------------
 * Get the number of colors in a cs_mesh_select_t structure.
 *
 * parameters:
 *   selection --> pointer to selection structure
 *
 * returns:
 *   Number of colors in the cs_mesh_select_t struture
 *----------------------------------------------------------------------------*/

cs_int_t
cs_mesh_select_get_n_colors(const cs_mesh_select_t *selection);

/*----------------------------------------------------------------------------
 * Get the number of groups in a cs_mesh_select_t struture.
 *
 * parameters:
 *   selection --> pointer to selection structure
 *
 * returns:
 *   Number of groups in the cs_mesh_select_t struture
 *----------------------------------------------------------------------------*/

cs_int_t
cs_mesh_select_get_n_groups(const cs_mesh_select_t *selection);

/*----------------------------------------------------------------------------
 * Extract border faces from criteria included in a cs_mesh_select_t
 * structure.
 *
 * parameters:
 *   mesh               --> pointer to a mesh structure
 *   selection          --> pointer to a selection structure
 *   n_selected_b_faces <-- number of border faces selected
 *   b_face_select_lst  <-- list of the selected border faces
 *----------------------------------------------------------------------------*/

void
cs_mesh_select_extract_b_faces(const cs_mesh_t        *const mesh,
                               const cs_mesh_select_t       *selection,
                               cs_int_t                     *n_selected_b_faces,
                               cs_int_t                    **b_face_select_lst);

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_select_t structure.
 *
 * parameters:
 *   selection --> pointer to a selection structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_select_dump(cs_mesh_select_t  *selection);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_MESH_SELECT_H__ */
