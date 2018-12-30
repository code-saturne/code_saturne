#ifndef __CS_SELECTOR_H__
#define __CS_SELECTOR_H__

/*============================================================================
 * Build selection lists for faces or cells
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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fill a list of boundary face numbers verifying a given selection criteria.
 *
 * parameters:
 *   criteria        <-- selection criteria string
 *   n_b_faces       --> number of selected interior faces
 *   b_face_num_list --> list of selected boundary face numbers
 *                       (1 to n, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_num_list(const char  *criteria,
                                cs_lnum_t   *n_b_faces,
                                cs_lnum_t    b_face_num_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of interior faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria        <-- selection criteria string
 *   n_i_faces       --> number of selected interior faces
 *   i_face_num_list --> list of selected interior face numbers
 *                       (1 to n, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_i_face_num_list(const char  *criteria,
                                cs_lnum_t   *n_i_faces,
                                cs_lnum_t    i_face_num_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of cells verifying a given selection criteria.
 *
 * parameters:
 *   criteria      <-- selection criteria string
 *   n_cells       --> number of selected cells
 *   cell_num_list --> list of selected cell numbers
 *                     (1 to n, preallocated to cs_glob_mesh->n_cells)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cell_num_list(const char  *criteria,
                              cs_lnum_t   *n_cells,
                              cs_lnum_t    cell_num_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of boundary faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_b_faces   --> number of selected interior faces
 *   b_face_list --> list of selected boundary faces
 *                   (0 to n-1, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_list(const char  *criteria,
                            cs_lnum_t   *n_b_faces,
                            cs_lnum_t    b_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of interior faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_i_faces   --> number of selected interior faces
 *   i_face_list --> list of selected interior faces
 *                   (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_i_face_list(const char  *criteria,
                            cs_lnum_t   *n_i_faces,
                            cs_lnum_t    i_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of cells verifying a given selection criteria.
 *
 * parameters:
 *   criteria  <-- selection criteria string
 *   n_cells   --> number of selected cells
 *   cell_list --> list of selected cells
 *                 (0 to n-1, preallocated to cs_glob_mesh->n_cells)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cell_list(const char  *criteria,
                          cs_lnum_t   *n_cells,
                          cs_lnum_t    cell_list[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of cells verifying a given selection criteria.
 *
 * \param[in]   criteria    selection criteria string
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cell_vertices_list(const char  *criteria,
                                   cs_lnum_t   *n_vertices,
                                   cs_lnum_t    vtx_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of vertices belonging to a given list of cells.
 *
 * \param[in]   criteria    selection criteria string
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cell_vertices_list_by_ids(cs_lnum_t         n_cells,
                                          const cs_lnum_t   cell_ids[],
                                          cs_lnum_t        *n_vertices,
                                          cs_lnum_t         vtx_ids[]);

/*----------------------------------------------------------------------------
 * Fill lists of faces at the boundary of a set of cells verifying a given
 * selection criteria.
 *
 * parameters:
 *   criteria  <-- selection criteria string
 *   n_i_faces --> number of selected interior faces
 *   n_b_faces --> number of selected interior faces
 *   i_face_id --> list of selected interior faces
 *                 (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 *   b_face_id --> list of selected boundary faces
 *                 (0 to n-1, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cells_boundary(const char  *criteria,
                               cs_lnum_t   *n_i_faces,
                               cs_lnum_t   *n_b_faces,
                               cs_lnum_t    i_face_id[],
                               cs_lnum_t    b_face_id[]);

/*----------------------------------------------------------------------------
 * Fill a list of interior faces belonging to a given periodicity.
 *
 * parameters:
 *   perio_num <-- periodicity number
 *   n_i_faces --> number of selected interior faces
 *   i_face_id --> list of selected interior faces
 *                 (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_perio_face_list(int         perio_num,
                                cs_lnum_t  *n_i_faces,
                                cs_lnum_t   i_face_id[]);

/*----------------------------------------------------------------------------
 * Fill a list of families verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_families  --> number of selected families
 *   family_list --> list of selected families faces
 *                   (0 to n-1, preallocated to cs_glob_mesh->n_families + 1)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_family_list(const char  *criteria,
                            cs_lnum_t   *n_families,
                            cs_int_t     family_list[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SELECTOR_H__ */
