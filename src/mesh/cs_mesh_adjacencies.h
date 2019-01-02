#ifndef __CS_MESH_ADJACENCIES_H__
#define __CS_MESH_ADJACENCIES_H__

/*============================================================================
 * Additional mesh adjacencies.
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

#include "cs_base.h"
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


typedef struct {

  cs_flag_t    flag;    /* Compact way to store metadata */
  int          stride;

  cs_lnum_t    n_elts;
  cs_lnum_t   *idx;     /* size = n_elts + 1 */
  cs_lnum_t   *ids;     /* ids from 0 to n-1 (there is no multifold entry) */
  short int   *sgn;     /* +/- 1 according to the orientation of the element */

} cs_adjacency_t;

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
 * \brief Allocate and define a cs_adjacency_t structure related to vertices.
 *
 * Adjacent vertices are accessed based on the vertex with lowest id.
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
