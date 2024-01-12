#ifndef __CS_MESH_ADJACENCIES_H__
#define __CS_MESH_ADJACENCIES_H__

/*============================================================================
 * Additional mesh adjacencies.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include "cs_base_accel.h"
#include "cs_halo.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @defgroup adjacency flags specifying metadata related to a cs_adjacency_t
 *           structure
 * @{
 */

/*!  1: Members of the structure are shared with a parent structure */
#define  CS_ADJACENCY_SHARED    (1 << 0)
/*!  2: Ajacency relies on a stride. No index are used to access ids */
#define  CS_ADJACENCY_STRIDE    (1 << 1)
/*!  4: Ajacency has a sgn member to known how is oriented each element */
#define  CS_ADJACENCY_SIGNED    (1 << 2)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! General adjacency structure */

typedef struct {

  cs_flag_t    flag;    /*!< compact way to store metadata */
  int          stride;  /*!< stride if strided, 0 if indexed */

  cs_lnum_t    n_elts;
  cs_lnum_t   *idx;     /*!< index, or NULL if strided (size = n_elts + 1) */
  cs_lnum_t   *ids;     /*!< ids from 0 to n-1 (there is no multifold entry) */
  short int   *sgn;     /*!< -/+1 according to the orientation of the element */

} cs_adjacency_t;

/*! Additional mesh adjacencies build from mesh structure */

typedef struct {

  /* metadata */

  bool        single_faces_to_cells;   /*!< true if a single face is adjacent
                                         to 2 given cells */

  /* cells -> cells connectivity (standard) */

  cs_lnum_t  *cell_cells_idx;          /*!< cells to cells indexes */
  cs_lnum_t  *cell_cells;              /*!< cells to cells adjacency */

  /* cells -> cells connectivity (extended) */

  const cs_lnum_t  *cell_cells_e_idx;  /*!< indexes (shared) */
  const cs_lnum_t  *cell_cells_e;      /*!< adjacency (shared) */

  /* cells->interior faces connectivity
     (same index as cells->cells connectivity) */

  cs_lnum_t  *cell_i_faces;            /*!< cells to interior faces adjacency */
  short int  *cell_i_faces_sgn;        /*!< cells to interior faces orientation */

  /* cells -> boundary faces connectivity */

  cs_lnum_t        *cell_b_faces_idx;  /*!< cells to boundary faces index */
  cs_lnum_t        *cell_b_faces;      /*!< cells to boundary faces adjacency */

  /* cells -> hidden boundary faces connectivity, if present;
     "hidden" boundary faces are those whose ids are in the meshe's
     ] n_b_faces, n_b_faces_all ] range */

  cs_lnum_t        *cell_hb_faces_idx;  /*!< cells to hidden boundary faces
                                          index, if present */
  cs_lnum_t        *cell_hb_faces;      /*!< cells to hidden boundary faces
                                          adjacency, if present */

  /* cells -> faces connectivity */

  const cs_adjacency_t  *c2f;          /*!< cells to faces adjacency */

  cs_adjacency_t        *_c2f;         /*!< cells to faces adjacency if owner,
                                         NULL otherwise */

  /* cells -> vertices connectivity */

  const cs_adjacency_t  *c2v;          /*!< cells to vertices adjacency */

  cs_adjacency_t        *_c2v;         /*!< cells to vertices adjacency if owner,
                                         NULL otherwise */

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
 * \brief  Get non-const pointer to cs_glob_mesh_adacencies.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_adjacencies_t *
cs_mesh_adjacencies_get_global(void);

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
/*!
 * \brief  Ensure presence of cell -> interior face connectivites in
 *         mesh adjacencies helper API relative to mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_cell_i_faces(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map some global mesh adjacency arrays for use on device.
 *
 * More elements may be mapped dependin on which arrays are used in
 * accelerated algorithms.
 *
 * \param[in]  alloc_mode  chosen allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_device(cs_alloc_mode_t  alloc_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return cell -> face connectivites in
 *         mesh adjacencies helper API relative to mesh.
 *
 * Boundary faces appear first, interior faces second.
 *
 * This connectivity is built only when first requested, then updated later if
 * needed.
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t  *
cs_mesh_adjacencies_cell_faces(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return cell -> vertex connectivites in
 *         mesh adjacencies helper API relative to mesh.
 *
 * This connectivity is built only when first requested, then updated later if
 * needed.
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t  *
cs_mesh_adjacencies_cell_vertices(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_adjacency_t structure of size n_elts
 *
 * \param[in]  flag       metadata related to the new cs_adjacency to create
 * \param[in]  stride     > 0 if useful otherwise ignored
 * \param[in]  n_elts     number of entries of the indexed list
 *
 * \return  a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_adjacency_create(cs_flag_t    flag,
                    int          stride,
                    cs_lnum_t    n_elts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_adjacency_t structure sharing arrays scanned with a
 *          stride
 *
 * \param[in]  n_elts    number of elements
 * \param[in]  stride    value of the stride
 * \param[in]  ids       array of element ids (size = stride * n_elts)
 * \param[in]  sgn       array storing the orientation (may be NULL)
 *
 * \return  a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_adjacency_create_from_s_arrays(cs_lnum_t    n_elts,
                                  int          stride,
                                  cs_lnum_t   *ids,
                                  short int   *sgn);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_adjacency_t structure sharing arrays scanned with an
 *          index
 *
 * \param[in]  n_elts    number of elements
 * \param[in]  idx       array of size n_elts + 1
 * \param[in]  ids       array of element ids (size = idx[n_elts])
 * \param[in]  sgn       array storing the orientation (may be NULL)
 *
 * \return  a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_adjacency_create_from_i_arrays(cs_lnum_t     n_elts,
                                  cs_lnum_t    *idx,
                                  cs_lnum_t    *ids,
                                  short int    *sgn);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a cs_adjacency_t structure.
 *
 * \param[in, out]  p_adj   pointer of pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_adjacency_destroy(cs_adjacency_t   **p_adj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a new cs_adjacency_t structure from the composition of
 *          two cs_adjacency_t structures: (1) A -> B and (2) B -> C
 *          The resulting structure describes A -> C. It does not rely on a
 *          stride and has no sgn member.
 *
 * \param[in]  n_c_elts  number of elements in C set
 * \param[in]  a2b       adjacency A -> B
 * \param[in]  b2c       adjacency B -> C
 *
 *\return  a pointer to the cs_adjacency_t structure A -> C
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_adjacency_compose(int                      n_c_elts,
                     const cs_adjacency_t    *a2b,
                     const cs_adjacency_t    *b2c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a new cs_adjacency_t structure from a one corresponding to
 *          A -> B. The resulting structure deals with B -> A
 *
 * \param[in]  n_b_elts    size of the set of B elements
 * \param[in]  a2b         pointer to the A -> B cs_adjacency_t structure
 *
 * \return  a new pointer to the cs_adjacency_t structure B -> A
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_adjacency_transpose(int                     n_b_elts,
                       const cs_adjacency_t   *a2b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each sub-list related to an entry in a cs_adjacency_t
 *          structure
 *
 * \param[in]  adj     pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_adjacency_sort(cs_adjacency_t   *adj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   In an indexed list, remove id(s) corresponding to the current
 *          index. Useful for instance in order to prepare a matrix structure
 *          in MSR storage
 *
 * \param[in, out] adj     pointer to the cs_adjacency_t structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_adjacency_remove_self_entries(cs_adjacency_t   *adj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_adjacency_t structure to a file or into the
 *          standard output
 *
 * \param[in]  name    name of the dump file. Can be set to NULL
 * \param[in]  _f      pointer to a FILE structure. Can be set to NULL.
 * \param[in]  adj     pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_adjacency_dump(const char           *name,
                  FILE                 *_f,
                  cs_adjacency_t       *adj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a cells to faces adjacency structure.
 *
 * With the boundary_order option set to 0, boundary faces come first, so
 * interior face ids are shifted by the number of boundary faces.
 * With boundary_order set to 1, boundary faces come last, so face ids are
 * shifted by the number of interior faces.
 *
 * \param[in]  m               pointer to a cs_mesh_t structure
 * \param[in]  boundary_order  boundaries first (0) or last (1)
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_adjacency_c2f(const cs_mesh_t  *m,
                      int               boundary_order);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a partial cells to faces adjacency structure.
 *
 * With the boundary_order option set to 0, boundary faces come first, so
 * interior face ids are shifted by the number of boundary faces.
 * With boundary_order set to 1, boundary faces come last, so face ids are
 * shifted by the number of interior faces.
 * With boundary order set to -1, boundary faces are ignored.
 *
 * Only internal faces as seen from the cell with lowest cell id (out of 2
 * cells adjacent to that face) are added to the structure, so each face
 * appears only once, not twice. In other words, only the lower part
 * of the corresponding cell-cells adjacency matrix is built.
 *
 * By construction, face ids adjacent to each cell are ordered.
 *
 * \param[in]  m               pointer to a cs_mesh_t structure
 * \param[in]  boundary_order  boundaries first (0), last (1), or ignored (-1)
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_adjacency_c2f_lower(const cs_mesh_t  *m,
                            int               boundary_order);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a cells to boundary faces adjacency structure.
 *
 * By construction, face ids adjacent to each cell are ordered.
 *
 * \param[in]  m  pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_adjacency_c2f_boundary(const cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a cs_adjacency_t structure related to the
 *        connectivity vertex to vertices through edges.
 *
 * Adjacent vertices are accessed based on the vertex with lowest id.
 * Another v2v connectivity through cells is possible. Please read the
 * \file cs_cdo_connect.c source code if interested
 *
 * \param[in]  m  pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_adjacency_v2v(const cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_ADJACENCIES__ */
