#ifndef __CS_MESH_GROUP_H__
#define __CS_MESH_GROUP_H__

/*============================================================================
 * Management of mesh groups.
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"
#include "cs_mesh.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*============================================================================
 *  Public function header for Fortran API
 *===========================================================================*/

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clean mesh group definitions.
 *
 * \param[in]  mesh  pointer to mesh structure to modify
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_clean(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Combine mesh group classes.
 *
 * \param[in, out]  mesh          pointer to mesh structure to modify
 * \param[in]       n_elts        number of local elements
 * \param[in]       gc_id_idx     element group class index (size: n_elts +1)
 * \param[in]       gc_id         input element group classes
 *                                (size: gc_id_idx[n_elts])
 * \param[in]       gc_id_merged  output element group classes (size: n_elts)
 *
 * \return  array of new element group class ids
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_combine_classes(cs_mesh_t   *mesh,
                              cs_lnum_t    n_elts,
                              cs_lnum_t    gc_id_idx[],
                              int          gc_id[],
                              int          gc_id_merged[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to cells, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected cells
 * \param[in]       n_selected_cells  number of selected cells
 * \param[in]       selected_cell_id  selected cell ids (size: n_selected_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_cells_set(cs_mesh_t        *mesh,
                        const char       *name,
                        cs_lnum_t         n_selected_cells,
                        const cs_lnum_t   selected_cell_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to interior faces, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_i_faces_set(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to boundary faces, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_b_faces_set(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected cells to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected cells
 * \param[in]       n_selected_cells  number of selected cells
 * \param[in]       selected_cell_id  selected cell ids (size: n_selected_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_cells_add(cs_mesh_t        *mesh,
                        const char       *name,
                        cs_lnum_t         n_selected_cells,
                        const cs_lnum_t   selected_cell_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected interior faces to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_i_faces_add(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected boundary faces to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_b_faces_add(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_GROUP_H__ */
