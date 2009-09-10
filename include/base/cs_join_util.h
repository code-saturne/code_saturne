/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2009 EDF S.A., France
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
 *===========================================================================*/

#ifndef __CS_JOIN_UTIL_H__
#define __CS_JOIN_UTIL_H__

/*============================================================================
 * Manipulation of low-level structures for joining operations
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local library headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_selector.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef enum {

  CS_JOIN_TYPE_CONFORM,
  CS_JOIN_TYPE_NO_CONFORM

} cs_join_type_t;

/*----------------------------------------------------------------------------
 * Set of user parameters to control the join operation
 *----------------------------------------------------------------------------*/

typedef struct {

  int  num;        /* number associated to the current join operation */

  /* Octree - Quadtree search algorithm */
  /* ---------------------------------- */

  int   tree_max_level;     /* Deepest level reachable during tree building */
  int   tree_n_max_boxes;   /* Max. number of boxes which can be related to a
                               a leaf of the tree if level != tree_max_level */

  int   tree_max_box_ratio; /* Stop tree building if:
                               n_linked_boxes > tree_max_box_ratio*n_init_boxes
                            */

  /* Geometric parameters */
  /* -------------------- */

  /* parameter used to compute the tolerance associated to each vertex.
     Also used for finding equivalent vertices during edge intersections */

  double  fraction;

  /* maximum angle between normals of two faces considered to
     be in the same plane (for face split) */

  double  plane;

  /* Coef. used to reduce the tolerance: new tol. = tol * coef
     Values between [0.0, 1.0[ */

  double  reduce_tol_factor;

  /* Coef. used to modify the tolerance associated to each vertex before the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  double  merge_tol_coef;

  /* Coef. used to modify locally the tolerance associated to each vertex
     BEFORE adding equivalences between vertices after edge intersections.
     If coef = 0.0 => add no equivalence
     If coef < 1.0 => reduce the number of equivalences between vertices
                      sharing the same edge
     If coef = 1.0 => no change
     If coef > 1.0 => increase the number of equivalences between vertices
                      sharing the same edge. NOT ADVISED. */

  double  edge_equiv_tol_coef;

  /* Parameter to switch on/off the influence of adjacent faces in the
     computation of tolerance */

  cs_bool_t  include_adj_faces;

  /* Parameter used to define if we get vertex equivalence trough the
     comparison of vertex tolerance or through the difference of curvilinear
     abscissa of vertices on edges.
     If include_adj_faces = false => this parameter should not have any
     effect. (Not a user-defined parameter) */

  cs_bool_t  edge_equiv_by_tolerance;

  int  max_sub_faces;  /* Maximum number of sub-faces when splitting a face */

  int  verbosity;  /* Level of display:
                      O : no information printed
                      1 : general information printed
                      2 : more information printed
                      5 and beyond : highest level (DEBUG LEVEL) */


} cs_join_param_t;

/*----------------------------------------------------------------------------
 * Structure used to store the result of the extraction of entities
 * implied in the joining operation
 *----------------------------------------------------------------------------*/

typedef struct { /* Structure used to synchronize single elements */

  int      n_elts;
  int      n_ranks;
  int     *ranks;
  int     *index;
  int     *array;

} cs_join_sync_t;

typedef struct {

  cs_int_t      n_faces;     /* Number of border faces selected
                                for the joining operation */
  fvm_gnum_t    n_g_faces;   /* Global number of border faces selected
                                for the joining operation */
  cs_int_t     *faces;       /* List of selected border faces */

  fvm_gnum_t   *compact_face_gnum;    /* Global face numbering defined
                                         on the selected faces */
  fvm_gnum_t   *compact_rank_index;   /* Distribution of the selected faces
                                         over the ranks */

  cs_int_t     *cell_filter;     /* Size: n_cells
                                    value = -1 if not implied
                                    value >= 0 else */

  cs_real_t    *cell_cen;        /* Cell center for cells implied */
  fvm_gnum_t   *cell_gnum;       /* Global cell numbering of the cells
                                         holding the selected face */

  cs_int_t      n_vertices;      /* Number of vertices selected
                                    for the joining operation */
  fvm_gnum_t    n_g_vertices;    /* Global number of selected vertices */
  cs_int_t     *vertices;        /* List of selected vertices */

  /* Adjacent faces of the current face selection: border and interior */

  cs_int_t      n_b_adj_faces;
  cs_int_t      n_i_adj_faces;

  cs_int_t     *b_adj_faces;
  cs_int_t     *i_adj_faces;

  /*
     Single elements (Only possible in parallel. It appears
     when the domain splitting has a poor quality and elements
     on the joining interface are prisms or tetraedrals)
     s = single / c = coupled
  */

  cs_bool_t    do_single_sync;

  cs_join_sync_t  *s_vertices;
  cs_join_sync_t  *c_vertices;
  cs_join_sync_t  *s_edges;
  cs_join_sync_t  *c_edges;

} cs_join_select_t;

/*----------------------------------------------------------------------------
 * Structure used to store information about a block distribution
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_gnum_t  n_g_elts;    /* Global number of elements to distribute */
  fvm_gnum_t  first_gnum;  /* Global number of the element in the local block */

  int  n_ranks;            /* Number of processes in the communicator used
                              to define the current distribution */
  int  local_rank;         /* Id of the current process in the communicator
                              used to define the current distribution */

  size_t  size;            /* Size of block for the given set of parameters */
  size_t  local_size;      /* Number of global elements to treat on the local
                              rank according to the distribution */

} cs_join_block_info_t;

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a set of parameters to control a contiguous distribution by block.
 *
 * parameters:
 *   n_g_elts   <-- global number of elements to treat
 *   n_ranks    <-- number of ranks in the MPI communicator related to the
 *                  cs_join_block_info_t structure to create
 *   local_rank <-- rank in the MPI communicator related to the
 *                  cs_join_block_info_t structure to create
 *
 * returns:
 *   a new defined cs_join_block_info_t structure
 *----------------------------------------------------------------------------*/

cs_join_block_info_t
cs_join_get_block_info(fvm_gnum_t  n_g_elts,
                       int         n_ranks,
                       int         local_rank);

/*----------------------------------------------------------------------------
 * Initialize a cs_join_param_t structure.
 *
 * parameters:
 *   join_id       <-- id of the current joining operation
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   rtf           <-- value of the "reduction tolerance factor" parameter
 *   ftf           <-- value of the "merge tolerance factor" parameter
 *   etf           <-- value of the "edge equiv. tolerance factor" parameter
 *   max_sub_faces <-- maximum number of sub-faces allowed during splitting
 *   tml           <-- value of the "tree max level" parameter
 *   tmb           <-- value of the "tree max boxes" parameter
 *   tmr           <-- value of the "tree max ratio" parameter
 *   verbosity     <-- level of verbosity required
 *
 * returns:
 *   a pointer to a cs_join_param_t structure
 *---------------------------------------------------------------------------*/

cs_join_param_t
cs_join_param_define(int     join_id,
                     double  fraction,
                     double  plane,
                     double  rtf,
                     double  ftf,
                     double  etf,
                     int     max_sub_faces,
                     int     tml,
                     int     tmb,
                     double  tmr,
                     int     verbosity);

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_select_t structure.
 *
 * parameters:
 *   selection_criteria <-- pointer to a cs_mesh_select_t structure
 *   verbosity          <-- level of verbosity required
 *
 * returns:
 *   pointer to a newly created cs_join_select_t structure
 *---------------------------------------------------------------------------*/

cs_join_select_t *
cs_join_select_create(const char  *selection_criteria,
                      int          verbosity);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_select_t structure.
 *
 * parameters:
 *   join_select <-- pointer to pointer to structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_select_destroy(cs_join_select_t  **join_select);

/*----------------------------------------------------------------------------
 * Build vertex -> vertex index for a selection of faces.
 *
 * "v2v_idx" is already allocated to the number of vertices in the mesh.
 * At this stage, it is just a counter.
 *
 * parameters:
 *   n_faces <-- number of selected faces
 *   faces   <-- list of selected faces
 *   f2v_idx <-- face -> vertex connectivity index
 *   f2v_lst <-- face -> vertex connectivity list
 *   v2v_idx <-> index to build (already allocated and may be used again)
 *---------------------------------------------------------------------------*/

void
cs_join_build_edges_idx(cs_int_t        n_faces,
                        const cs_int_t  faces[],
                        const cs_int_t  f2v_idx[],
                        const cs_int_t  f2v_lst[],
                        cs_int_t        v2v_idx[]);

/*----------------------------------------------------------------------------
 * Build vertex -> vertex list for a selection of faces.
 * "count" and "v2v_lst" are already allocated to the number of vertices in
 * the mesh.
 *
 * parameters:
 *   n_faces <-- number of selected faces
 *   faces   <-- list of selected faces
 *   f2v_idx <-- face -> vertex connectivity index
 *   f2v_lst <-- face -> vertex connectivity list
 *   count   <-> array used to count the number of values already added
 *   v2v_idx <-- vertex -> vertex connect. index
 *   v2v_lst <-> vertex -> vertex connect. list to build (can be used again)
 *---------------------------------------------------------------------------*/

void
cs_join_build_edges_lst(cs_int_t        n_faces,
                        const cs_int_t  faces[],
                        const cs_int_t  f2v_idx[],
                        const cs_int_t  f2v_lst[],
                        cs_int_t        count[],
                        const cs_int_t  v2v_idx[],
                        cs_int_t        v2v_lst[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_UTIL_H__ */
