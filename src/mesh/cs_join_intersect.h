#ifndef __CS_JOIN_INTERSECT_H__
#define __CS_JOIN_INTERSECT_H__

/*============================================================================
 * Set of subroutines for finding intersections between bounding boxes
 *  - on faces
 *  - on edges
 *===========================================================================*/

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
 * Standard C library headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local library headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Definition of a structure defining an intersection */
/* -------------------------------------------------- */

typedef struct {

  cs_lnum_t   edge_id;    /* id of the edge implied in this intersection */
  cs_lnum_t   vtx_id;     /* id of the vertex resulting of the intersection */
  cs_coord_t  curv_abs;   /* curvilinear abscissa of the intersection */

} cs_join_inter_t;

/* Definition of a structure defining a set of intersections */
/* --------------------------------------------------------- */

typedef struct {

  cs_lnum_t         n_max_inter;  /* max. number of intersections allocated */
  cs_lnum_t         n_inter;      /* number of intersections */

  cs_join_inter_t  *inter_lst;    /* size = 2 * n_intersections
                                     one inter_t structure for each edge
                                     implied in an intersection */

} cs_join_inter_set_t;

/* Definition of a structure defining a set of intersections on edges */
/* ------------------------------------------------------------------ */

typedef struct {

  cs_lnum_t    n_edges;    /* Number of edges implied in an intersection */

  cs_gnum_t   *edge_gnum;  /* Global number of the related edges */
  cs_lnum_t   *index;      /* Indexed list of vertex num describing
                              intersections on a given edge without
                              vertices at the extremity and ordered
                              by curvilinear abscissa */

  cs_lnum_t   *vtx_lst;    /* List of new vertex num */
  cs_gnum_t   *vtx_glst;   /* List of new vertex global num */
  cs_coord_t  *abs_lst;    /* List of curvilinear abscissa */

  cs_lnum_t    max_sub_size;

} cs_join_inter_edges_t;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create a new cs_join_inter_set_t structure.
 *
 * parameters:
 *  init_size   <-- number of init. cs_join_inter_t structure to allocate
 *
 * returns:
 *  a pointer to a new inter_set_t structure.
 *---------------------------------------------------------------------------*/

cs_join_inter_set_t *
cs_join_inter_set_create(cs_lnum_t init_size);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_inter_set_t structure.
 *
 * parameters:
 *   inter_set <-> a pointer to the inter_set_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_inter_set_destroy(cs_join_inter_set_t  **inter_set);

/*----------------------------------------------------------------------------
 * Dump a cs_join_inter_set_t structure.
 *
 * parameters:
 *   f     <-- handle to output file
 *   i_set <-- cs_join_inter_set_t structure to dump
 *   edges <-- associated cs_join_edge_t structure
 *   mesh  <-- associated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_inter_set_dump(FILE                       *f,
                       const cs_join_inter_set_t  *i_set,
                       const cs_join_edges_t      *edges,
                       const cs_join_mesh_t       *mesh);

/*----------------------------------------------------------------------------
 * Allocate and initialize a new cs_join_inter_edges_t structure.
 *
 * parameters:
 *   n_edges <-- number of edges
 *
 * returns:
 *   a pointer to the created cs_join_inter_edges_t structure.
 *---------------------------------------------------------------------------*/

cs_join_inter_edges_t *
cs_join_inter_edges_create(cs_lnum_t  n_edges);

/*----------------------------------------------------------------------------
 * Build a cs_join_inter_edges_t structure (useful to find equivalence on
 * edges and to apply vertex merge to a cs_join_mesh_t structure).
 *
 * parameters:
 *   edges     <-- cs_join_edges_t structure
 *   inter_set <-- structure storing data on edge intersections
 *
 * returns:
 *  a pointer to the created cs_join_inter_edges_t structure
 *---------------------------------------------------------------------------*/

cs_join_inter_edges_t *
cs_join_inter_edges_define(const cs_join_edges_t      *edges,
                           const cs_join_inter_set_t  *inter_set);

/*----------------------------------------------------------------------------
 * Destroy an cs_join_inter_edges_t structure.
 *
 * parameters:
 *   inter_edges <-> pointer to cs_join_inter_edges_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_inter_edges_destroy(cs_join_inter_edges_t  **inter_edges);

/*----------------------------------------------------------------------------
 * Find non-trivial equivalences between vertices sharing the same edges.
 *
 * For instance an equivalence between a vertex from the extremity of an edge
 * and a vertex created by an edge-edge intersection.
 *
 * parameters:
 *  param       <-- set of user-defined parameter
 *  mesh        <-- pointer to the local cs_join_mesh_t structure
 *                  which has initial vertex data
 *  edges       <-- list of edges
 *  inter_edges <-- structure including data on edge intersections
 *  vtx_equiv   <-> structure dealing with vertex equivalences
 *---------------------------------------------------------------------------*/

void
cs_join_add_equiv_from_edges(cs_join_param_t               param,
                             cs_join_mesh_t               *mesh,
                             const cs_join_edges_t        *edges,
                             const cs_join_inter_edges_t  *inter_edges,
                             cs_join_eset_t               *vtx_equiv);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Synchronize the definition of intersections on each edge by block.
 *
 * parameters:
 *   edges <-- cs_join_edges_t structure
 *   mesh  <-- cs_join_mesh_t structure
 *   part  <-- structure storing data on edge intersections by partition
 *
 * returns:
 *   newly allocated cs_join_inter_edges_t, synchronized and defined on
 *   a block
 *---------------------------------------------------------------------------*/

cs_join_inter_edges_t *
cs_join_inter_edges_part_to_block(const cs_join_mesh_t         *mesh,
                                  const cs_join_edges_t        *edges,
                                  const cs_join_inter_edges_t  *part);

/*----------------------------------------------------------------------------
 * Synchronize the definition of intersections on each edge from a
 * cs_join_inter_edges_t structure defined on a block.
 *
 * parameters:
 *   n_g_egdes <-- global number of edges
 *   block     <-- synchronized cs_join_inter_edges_t struct. by block
 *   part      <-> cs_join_inter_edges_t to synchronized on partition
 *---------------------------------------------------------------------------*/

void
cs_join_inter_edges_block_to_part(cs_gnum_t                     n_g_edges,
                                  const cs_join_inter_edges_t  *block,
                                  cs_join_inter_edges_t        *part);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Redefine a cs_join_inter_edges_t structure to be consistent with the local
 * numbering of a given couple of cs_join_mesh_t structure and
 * cs_join_edges_t structure.
 * Add future new vertices for the face definition in cs_join_mesh_t
 *
 * parameters:
 *   verbosity   <-- verbosity level
 *   edges       <-- cs_join_edges_t structure
 *   mesh        <-> cs_join_mesh_t structure
 *   inter_edges <-> current cs_join_inter_edges_t struct. to work with
 *---------------------------------------------------------------------------*/

void
cs_join_intersect_update_struct(int                      verbosity,
                                const cs_join_edges_t   *edges,
                                cs_join_mesh_t          *mesh,
                                cs_join_inter_edges_t  **inter_edges);

/*----------------------------------------------------------------------------
 * Get all real edges intersections among possible edges intersections.
 *
 * parameters:
 *   param         <-- set of user-defined parameters for the joining
 *   edge_edge_vis <-- a pointer to a cs_join_gset_t structure
 *   edges         <-- pointer to a structure defining edges
 *   mesh          <-- pointer to the cs_join_mesh_t structure
 *                     which has the face connectivity
 *   inter_set     <-> pointer to a structure including data on edge
 *                     intersections
 *   vtx_eset      <-> pointer to a structure dealing with vertex
 *                     equivalences
 *
 * returns:
 *   the type of joining encountered (conforming or not)
 *---------------------------------------------------------------------------*/

cs_join_type_t
cs_join_intersect_edges(cs_join_param_t         param,
                        const cs_join_gset_t   *edge_edge_vis,
                        const cs_join_edges_t  *edges,
                        const cs_join_mesh_t   *mesh,
                        cs_join_eset_t        **vtx_eset,
                        cs_join_inter_set_t   **inter_set);

/*----------------------------------------------------------------------------
 * Build a tree structure on which we associate leaves and face bounding boxes.
 * Create a cs_join_gset_t structure (indexed list on global numbering)
 * storing potential intersections between face bounding boxes.
 *
 * parameters:
 *   param     <-- set of user-defined parameters
 *   join_mesh <-- cs_join_mesh_t structure where faces are defined
 *   stats     <-> joining statistics
 *
 * returns:
 *   a new allocated pointer to a cs_join_gset_t structure storing the
 *   face - face visibility.
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_intersect_faces(const cs_join_param_t   param,
                        const cs_join_mesh_t   *join_mesh,
                        cs_join_stats_t        *stats);

/*----------------------------------------------------------------------------
 * Transform face visibility into edge visibility.
 *
 * parameters:
 *   mesh       <-- pointer to a cs_join_mesh_t structure
 *   edges      <-- pointer to a cs_join_edges_t structure
 *   face_visib <-- pointer to a cs_join_gset_t structure
 *
 * returns:
 *   a new allocated cs_join_gset_t structure holding edge visibility
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_intersect_face_to_edge(const cs_join_mesh_t   *mesh,
                               const cs_join_edges_t  *edges,
                               const cs_join_gset_t   *face_visib);

/*----------------------------------------------------------------------------
 * Dump a cs_join_inter_edges_t structure.
 *
 * parameters:
 *   f           <-- handle to output file
 *   inter_edges <-- cs_join_inter_edges_t structure to dump
 *   edges       <-- list of edges
 *   mesh        <-- associated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_inter_edges_dump(FILE                         *f,
                         const cs_join_inter_edges_t  *inter_edges,
                         const cs_join_edges_t        *edges,
                         const cs_join_mesh_t         *mesh);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_INTERSECT_H__ */
