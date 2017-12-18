#ifndef __CS_JOIN_MESH_H__
#define __CS_JOIN_MESH_H__

/*============================================================================
 * Subroutines useful to manipulate a cs_join_mesh_t structure
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local library headers
 *---------------------------------------------------------------------------*/

#include "fvm_defs.h"

#include "cs_base.h"
#include "cs_join_util.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

typedef enum {

  CS_JOIN_FACE_UNDEFINED,
  CS_JOIN_FACE_BORDER,
  CS_JOIN_FACE_MULTIPLE_BORDER,
  CS_JOIN_FACE_INTERIOR

} cs_join_face_type_t;

typedef struct {

  cs_join_state_t  state;      /* State of the vertices (perio/origin/...) */
  cs_gnum_t        gnum;       /* Global vertex number */
  double           tolerance;  /* Tolerance = radius of the sphere in which
                                  intersection and merge is possible. */
  double           coord[3];   /* Coordinates of the vertex */

} cs_join_vertex_t;

typedef struct {

  /* Edge numbering is defined by the ordering of the couples of vertices
     in their global numbering */

  cs_lnum_t    n_edges;    /* Local number of edges */
  cs_gnum_t    n_g_edges;  /* Global number of edges */
  cs_lnum_t   *def;        /* Definition of each edge by a couple of vertex
                              numbers */
  cs_gnum_t   *gnum;       /* Global numbering of edges */

  /*
    Edge definition through the relation between vertices :

    vtx_idx: index on vertices : -> define first vertex :
    V1
    adj_vtx_lst: list of coupled vertices with the first vertex V1:
    V1a, V1b, ...
    edge_lst: list of edge numbers relative to the defined couple:
    (V1, V1a), (V1, V1b), ...
  */

  cs_lnum_t   n_vertices;    /* Number of vertices in index */
  cs_lnum_t  *vtx_idx;       /* Index on first vertices */
  cs_lnum_t  *adj_vtx_lst;   /* List of adjacent vertices */
  cs_lnum_t  *edge_lst;      /* List of corresponding edge ids */

} cs_join_edges_t;

/* Structure defining a mesh on selected faces for the joining operation */

typedef struct {

  char         *name;   /* For post-processing and dump purpose */

  /* Face connectivity */

  cs_lnum_t    n_faces;
  cs_gnum_t    n_g_faces;
  cs_gnum_t   *face_gnum;
  cs_lnum_t   *face_vtx_idx;
  cs_lnum_t   *face_vtx_lst;

  /* Vertex data */

  cs_lnum_t          n_vertices;
  cs_gnum_t          n_g_vertices;
  cs_join_vertex_t  *vertices;

} cs_join_mesh_t;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype for the cs_join_vertex_t structure.
 *
 * returns:
 *   an MPI_Datatype associated to the cs_join_vertex_t structure.
 *---------------------------------------------------------------------------*/

MPI_Datatype
cs_join_mesh_create_vtx_datatype(void);

/*----------------------------------------------------------------------------
 * Create a function to define an operator for MPI reduction operation
 *
 * parameters:
 *   in        <--  input vertices
 *   inout     <->  in/out vertices (vertex with the min. toelrance)
 *   len       <--  size of input array
 *   datatype  <--  MPI_datatype associated to cs_join_vertex_t
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_mpi_vertex_min(cs_join_vertex_t   *in,
                            cs_join_vertex_t   *inout,
                            int                *len,
                            MPI_Datatype       *datatype);

/*----------------------------------------------------------------------------
 * Create a function to define an operator for MPI reduction operation
 *
 * parameters:
 *   in        <--  input vertices
 *   inout     <->  in/out vertices (vertex with the max. toelrance)
 *   len       <--  size of input array
 *   datatype  <--  MPI_datatype associated to cs_join_vertex_t
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_mpi_vertex_max(cs_join_vertex_t   *in,
                            cs_join_vertex_t   *inout,
                            int                *len,
                            MPI_Datatype       *datatype);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Allocate and initialize a new cs_join_mesh_t structure.
 *
 * parameters:
 *   name <-- name of the mesh
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure.
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create(const char  *name);

/*----------------------------------------------------------------------------
 * Get a cs_join_mesh_t structure with the given list of global faces inside.
 *
 * Exchange between ranks to get the connectivity associated to each
 * face of the global numbering list.
 *
 * parameters:
 *   mesh_name       <-- name of the created mesh
 *   n_elts          <-- number of elements in the global list
 *   glob_sel        <-- list of global elements (ordered)
 *   gnum_rank_index <-- index on ranks for the global elements
 *   local_mesh      <-- pointer to the local part of the distributed
 *                       cs_join_mesh_t structure on selected elements
 *
 * returns:
 *   a pointer to a new allocated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_glob_sel(const char            *mesh_name,
                                  cs_lnum_t              n_elts,
                                  const cs_gnum_t        glob_sel[],
                                  const cs_gnum_t        gnum_rank_index[],
                                  const cs_join_mesh_t  *local_mesh);

/*----------------------------------------------------------------------------
 * Allocate and define a cs_join_mesh_t structure relative to an extraction
 * of selected faces.
 *
 * The selection must be ordered.
 *
 * parameters:
 *   mesh_name   <-- name of the name to create
 *   subset_size <-- number of selected faces in the subset
 *   selection   <-- list of selected faces. Numbering in parent mesh
 *   parent_mesh <-- parent cs_join_mesh_t structure
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_subset(const char            *mesh_name,
                                cs_lnum_t              subset_size,
                                const cs_lnum_t        selection[],
                                const cs_join_mesh_t  *parent_mesh);

/*----------------------------------------------------------------------------
 * Define a cs_join_mesh_t structure from a selection of faces and its
 * related vertices.
 *
 * parameters:
 *   name       <-- mesh name of the resulting cs_join_mesh_t structure
 *   param      <-- set of user-defined parameters for the joining
 *   selection  <-> selected entities
 *   b_f2v_idx  <-- border "face -> vertex" connectivity index
 *   b_f2v_lst  <-- border "face -> vertex" connectivity
 *   i_f2v_idx  <-- interior "face -> vertex" connectivity index
 *   i_f2v_lst  <-- interior "face -> vertex" connectivity
 *   n_vertices <-- number of vertices in the parent mesh
 *   vtx_coord  <-- vertex coordinates in parent mesh
 *   vtx_gnum   <-- global numbering of vertices
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_select(const char              *name,
                                const cs_join_param_t    param,
                                cs_join_select_t        *selection,
                                const cs_lnum_t          b_f2v_idx[],
                                const cs_lnum_t          b_f2v_lst[],
                                const cs_lnum_t          i_f2v_idx[],
                                const cs_lnum_t          i_f2v_lst[],
                                const cs_lnum_t          n_vertices,
                                const cs_real_t          vtx_coord[],
                                const cs_gnum_t          vtx_gnum[]);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_mesh_t structure.
 *
 * parameters:
 *  mesh <->  pointer to pointer to cs_join_mesh_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_destroy(cs_join_mesh_t  **mesh);

/*----------------------------------------------------------------------------
 * Re-initialize an existing cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh <-> pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_reset(cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Copy a cs_join_mesh_t structure into another.
 *
 * parameters:
 *   mesh     <-> pointer to a cs_join_mesh_t structure to fill
 *   ref_mesh <-- pointer to the reference
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_copy(cs_join_mesh_t        **mesh,
                  const cs_join_mesh_t   *ref_mesh);

/*----------------------------------------------------------------------------
 * Compute the global min/max tolerance defined on vertices and display it
 *
 * parameters:
 *   param <-- user-defined parameters for the joining algorithm
 *   mesh  <-- pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_minmax_tol(cs_join_param_t    param,
                        cs_join_mesh_t    *mesh);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the connectivity of a list of global elements distributed over the
 * ranks.
 *
 * parameters:
 *   n_send     <-- number of face/rank couples to send
 *   send_rank  <-- index on ranks for the face distribution
 *   send_faces <-- list of face ids to send
 *   send_mesh  <-- pointer to the sending cs_join_mesh_t structure
 *   recv_mesh  <-> pointer to the receiving cs_join_mesh_t structure
 *   comm       <-- mpi communicator on which take places comm.
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_exchange(cs_lnum_t              n_send,
                      const cs_lnum_t        send_rank[],
                      const cs_lnum_t        send_faces[],
                      const cs_join_mesh_t  *send_mesh,
                      cs_join_mesh_t        *recv_mesh,
                      MPI_Comm               comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Destroy a cs_join_edges_t structure.
 *
 * parameters:
 *   edges <-> pointer to pointer to cs_join_edges_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_destroy_edges(cs_join_edges_t  **edges);

/*----------------------------------------------------------------------------
 * Order a cs_join_mesh_t structure according to its global face numbering
 *
 * Delete redundancies.
 *
 * parameters:
 *   mesh <-> pointer to a cs_join_mesh_t structure to order
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_face_order(cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Synchronize vertices definition over the rank. For a vertex with the same
 * global number but a not equal tolerance, we keep the minimal tolerance.
 *
 * parameters:
 *  mesh <->  pointer to the cs_join_mesh_t structure to synchronize
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_sync_vertices(cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Delete vertices which appear several times (same global number) and
 * vertices which are not used in face definition.
 *
 * parameters:
 *   mesh <-> pointer to cs_join_mesh_t structure to clean
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_vertex_clean(cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Clean the given cs_join_mesh_t structure, removing degenerate edges.
 *
 * parameters:
 *   mesh      <-> pointer to the cs_join_mesh_t structure to clean
 *   verbosity <-- level of display
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_clean(cs_join_mesh_t  *mesh,
                   int              verbosity);

/*----------------------------------------------------------------------------
 * Define a list of edges associated to a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh <-- pointer to a cs_join_mesh_t structure
 *
 * returns:
 *   a pointer to the new defined cs_join_edges_t structure.
 *---------------------------------------------------------------------------*/

cs_join_edges_t *
cs_join_mesh_define_edges(const cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Get the edge number relative to a couple of vertex numbers.
 *
 * edge_num > 0 if couple is in the same order as the edge->def
 * edge_num < 0 otherwise
 *
 * parameters:
 *   v1_num <-- vertex number for the first vertex
 *   v2_num <-- vertex number for the second vertex
 *   edges  <-- pointer to a cs_join_edges_t structure
 *
 * returns:
 *   an edge number relative to the couple of vertices
 *---------------------------------------------------------------------------*/

cs_int_t
cs_join_mesh_get_edge(cs_lnum_t               v1_num,
                      cs_lnum_t               v2_num,
                      const cs_join_edges_t  *edges);

/*----------------------------------------------------------------------------
 * Re-organize the cs_join_mesh_t structure after a renumbering of
 * the vertices following the merge operation + a new description of each
 * face.
 *
 * parameters:
 *   mesh             <-> pointer to the cs_join_mesh_t structure to update
 *   edges            <-- pointer to a cs_join_edges_t structure
 *   edge_index       <-- index on edges for the new vertices
 *   edge_new_vtx_lst <-- list of new vertices for each edge
 *   n_new_vertices   <-- new local number of vertices after merge
 *   old2new          <-- array storing the relation between old/new vertex id
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_update(cs_join_mesh_t         *mesh,
                    const cs_join_edges_t  *edges,
                    const cs_lnum_t         edge_index[],
                    const cs_lnum_t         edge_new_vtx_lst[],
                    cs_lnum_t               n_new_vertices,
                    const cs_lnum_t         old2new[]);

/*----------------------------------------------------------------------------
 * Compute for each face of the cs_join_mesh_t structure the face normal.
 * || face_normal || = 1 (divided by the area of the face)
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   mesh <-- pointer to a cs_join_mesh_t structure
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *
 *
 * returns:
 *   an array with the face normal for each face of the mesh
 *---------------------------------------------------------------------------*/

cs_real_t *
cs_join_mesh_get_face_normal(const cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Allocate and define an "edge -> face" connectivity
 *
 * parameters:
 *   mesh          <-- pointer to a cs_join_mesh_t structure
 *   edges         <-- pointer to a cs_join_edges_t structure
 *   edge_face_idx --> pointer to the edge -> face connect. index
 *   edge_face_lst --> pointer to the edge -> face connect. list
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_get_edge_face_adj(const cs_join_mesh_t   *mesh,
                               const cs_join_edges_t  *edges,
                               cs_lnum_t              *edge_face_idx[],
                               cs_lnum_t              *edge_face_lst[]);

/*----------------------------------------------------------------------------
 * Dump a cs_join_vertex_t structure into a file.
 *
 * parameters:
 *   f      <-- handle to output file
 *   vertex <-- cs_join_vertex_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_vertex(FILE                   *f,
                         const cs_join_vertex_t  vertex);

/*----------------------------------------------------------------------------
 * Dump a cs_join_mesh_t structure into a file.
 *
 * parameters:
 *   f    <-- handle to output file
 *   mesh <-- pointer to cs_join_mesh_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump(FILE                  *f,
                  const cs_join_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Dump a list of cs_join_edge_t structures.
 *
 * parameters:
 *   f     <-- handle to output file
 *   edges <-- cs_join_edges_t structure to dump
 *   mesh  <-- associated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_edges(FILE                   *f,
                        const cs_join_edges_t  *edges,
                        const cs_join_mesh_t   *mesh);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_MESH_H__ */
