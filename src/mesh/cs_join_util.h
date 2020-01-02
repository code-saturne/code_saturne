#ifndef __CS_JOIN_UTIL_H__
#define __CS_JOIN_UTIL_H__

/*============================================================================
 * Manipulation of low-level structures for joining operations
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

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_selector.h"
#include "cs_timer.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef enum {

  CS_JOIN_TYPE_NULL,
  CS_JOIN_TYPE_CONFORMING,
  CS_JOIN_TYPE_NON_CONFORMING

} cs_join_type_t;

typedef enum {

  CS_JOIN_STATE_UNDEF,
  CS_JOIN_STATE_NEW,
  CS_JOIN_STATE_ORIGIN,
  CS_JOIN_STATE_PERIO,
  CS_JOIN_STATE_MERGE,
  CS_JOIN_STATE_PERIO_MERGE,
  CS_JOIN_STATE_SPLIT

} cs_join_state_t;

/*----------------------------------------------------------------------------
 * Joining statistics
 *----------------------------------------------------------------------------*/

typedef struct {

  int        n_calls;               /* number of calls */

  /* Intersection determination info */

  int        bbox_layout;            /* bounding box layout */
  cs_gnum_t  bbox_depth[3];          /* box tree depth */
  cs_gnum_t  n_leaves[3];            /* number of leaves */
  cs_gnum_t  n_boxes[3];             /* number of boxes */
  cs_gnum_t  n_th_leaves[3];         /* number leaves over threshold */
  cs_gnum_t  n_leaf_boxes[3];        /* number of boxes per leaf */
  cs_gnum_t  box_mem_final[3];       /* final box memory required */
  cs_gnum_t  box_mem_required[3];    /* memory required */

  cs_timer_counter_t  t_box_build;   /* box build times */
  cs_timer_counter_t  t_box_query;   /* box query times */
  cs_timer_counter_t  t_inter_sort;   /* sort intersections times */

  /* other info */

  cs_timer_counter_t  t_l_join_mesh;  /* build local joining mesh times */
  cs_timer_counter_t  t_edge_inter;   /* edge intersection times */
  cs_timer_counter_t  t_new_vtx;      /* new vertices times */
  cs_timer_counter_t  t_merge_vtx;    /* merge vertices times */
  cs_timer_counter_t  t_u_merge_vtx;  /* update after merge vertices times */
  cs_timer_counter_t  t_split_faces;  /* split faces times */

  cs_timer_counter_t  t_total;        /* total time */

}  cs_join_stats_t;

/*----------------------------------------------------------------------------
 * Set of user parameters to control the join operation
 *----------------------------------------------------------------------------*/

typedef struct {

  int  num;         /* number associated to the current join operation */
  int  perio_type;  /* FVM_PERIODICITY_NULL for non-periodic joinings,
                       periodicity type for periodic joinings. */

  double perio_matrix[3][4];  /* Periodicity matrix for periodic joinings */

  /* Octree - Quadtree search algorithm */
  /* ---------------------------------- */

  int    tree_max_level;     /* Deepest level reachable during tree building */
  int    tree_n_max_boxes;   /* Max. number of boxes which can be related to
                               a leaf of the tree if level != tree_max_level */

  float  tree_max_box_ratio; /* Stop building tree when:
                                (  n_linked_boxes
                                 > tree_max_box_ratio * n_init_boxes) */
  float  tree_max_box_ratio_distrib; /* In parallel, tree_max_box_ratio for
                                        initial coarse tree used to
                                        determine load-distribution */

  /* Geometric parameters */
  /* -------------------- */

  /* parameter used to compute the tolerance associated to each vertex.
     Also used for finding equivalent vertices during edge intersections */

  float  fraction;

  /* maximum angle between normals of two faces considered to
     be in the same plane (for face split) */

  float  plane; /* in degree */
  double plane_criteria; /* cos(plane in rad)*cos(plane in rad) */

  /* Coef. used to modify the tolerance associated to each vertex before the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  float  merge_tol_coef;

  /* Coef. used to compute a limit on staightfoward merge between
     two vertices before the merge step. It should be a small value. */

  float  pre_merge_factor;

  /* Maximum number of equivalence breaks */

  int  n_max_equiv_breaks;

   /* Tolerance computation mode: tcm
      1: (default) tol = min. edge length related to a vertex * fraction
      2: tolerance is computed like in mode 1 with in addition, the
         multiplication by a coef. which is equal to the max sin(e1, e2)
         where e1 and e2 are two edges sharing the same vertex V for which
         we want to compute the tolerance
     11: like 1 but only in taking into account only the selected faces
     12: like 2 but only in taking into account only the selected faces */

  int  tcm;

   /* Intersection computation mode: icm
      1: (default) Original algorithm. Try to clip intersection on extremity
      2: New intersection algorithm. Avoid to clip intersection on extremity
   */

  int  icm;

  /* Maximum number of sub-faces when splitting a face */

  int  max_sub_faces;

  /* Verbosity:
       O : no information printed
       1 : general information printed
       2 : more information printed
       5 and beyond : highest level (DEBUG LEVEL) */

  int  verbosity;

  /* Visualization level:
       O : no visualization output
       1 : visualization output of joined faces
       2 : faces modified by joining
  */

  int  visualization;

  /* Preprocessing flag:
     true if this joining is part of preprocessing, false otherwise */

  bool preprocessing;

} cs_join_param_t;

/*----------------------------------------------------------------------------
 * Set of variables to synchronize single elements
 *---------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

typedef struct {

  cs_lnum_t   n_elts;
  int         n_ranks;
  int        *ranks;
  cs_lnum_t  *index;
  cs_lnum_t  *array;

} cs_join_sync_t;

/*----------------------------------------------------------------------------
 * Structure used to store the result of the extraction of entities
 * implied in the joining operation
 *---------------------------------------------------------------------------*/

typedef struct {

  cs_lnum_t     n_init_b_faces;  /* Number of border faces before joining */
  cs_lnum_t     n_init_i_faces;  /* Number of interior faces before joining */
  cs_lnum_t     n_init_vertices; /* Number of vertices before joining */

  cs_lnum_t     n_faces;     /* Number of border faces selected
                                for the joining operation */
  cs_gnum_t     n_g_faces;   /* Global number of border faces selected
                                for the joining operation */
  cs_lnum_t    *faces;       /* List of selected border faces */

  cs_gnum_t    *compact_face_gnum;    /* Global face numbering defined
                                         on the selected faces */
  cs_gnum_t    *compact_rank_index;   /* Distribution of the selected faces
                                         over the ranks */

  cs_lnum_t     n_vertices;      /* Number of vertices selected
                                    for the joining operation */
  cs_gnum_t     n_g_vertices;    /* Global number of selected vertices */
  cs_lnum_t     *vertices;        /* List of selected vertices */

  /* Adjacent faces of the current face selection: border and interior */

  cs_lnum_t     n_b_adj_faces;
  cs_lnum_t     n_i_adj_faces;

  cs_lnum_t    *b_adj_faces;
  cs_lnum_t    *i_adj_faces;

  /* Keep the status of all faces of the related cs_mesh_t */

  cs_join_state_t   *b_face_state;
  cs_join_state_t   *i_face_state;

  /* For periodicity handling: list of periodic vertex couples */

  cs_lnum_t    n_couples;
  cs_gnum_t   *per_v_couples;

  /*
     Single elements (Only possible in parallel). Appear mainly
     when the domain splitting has a poor quality and elements
     on the joining interface are prisms or tetrahedra
     s = single (receiver) / c = coupled (owner).
  */

  bool         do_single_sync;

  cs_join_sync_t  *s_vertices;
  cs_join_sync_t  *c_vertices;
  cs_join_sync_t  *s_edges;
  cs_join_sync_t  *c_edges;

} cs_join_select_t;

/*----------------------------------------------------------------------------
 * Highest level structure to manage the joining algorithm
 *---------------------------------------------------------------------------*/

typedef struct {

  cs_join_param_t   param;       /* Set of parameters used to control
                                    the joining operations */

  cs_join_stats_t   stats;       /* Performance statistics */

  cs_join_select_t  *selection;  /* Store entities implied in the joining
                                    operation */

  char              *criteria;   /* Criteria used to select border faces
                                    implied in the joining operation */

  char              *log_name;   /* Optional log file name */

} cs_join_t;

/*=============================================================================
 * Global variables
 *===========================================================================*/

extern int  cs_glob_join_count;
extern int  cs_glob_n_joinings;
extern cs_join_t  **cs_glob_join_array;

extern FILE  *cs_glob_join_log;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_t structure.
 *
 * parameters:
 *   join_number   <-- number related to the joining operation
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   perio_type    <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   perio_matrix  <-- periodicity transformation matrix
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   preprocessing <-- is joining part of the preprocessing stage ?
 *
 * returns:
 *   a pointer to a new allocated cs_join_t structure
 *---------------------------------------------------------------------------*/

cs_join_t *
cs_join_create(int                      join_number,
               const char              *sel_criteria,
               float                    fraction,
               float                    plane,
               fvm_periodicity_type_t   perio_type,
               double                   perio_matrix[3][4],
               int                      verbosity,
               int                      visualization,
               bool                     preprocessing);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_t structure.
 *
 * parameters:
 *  join <-> pointer to the cs_join_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_destroy(cs_join_t  **join);

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_select_t structure.
 *
 * parameters:
 *   selection_criteria <-- pointer to a cs_mesh_select_t structure
 *   perio_type         <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   verbosity          <-- level of verbosity required
 *
 * returns:
 *   pointer to a newly created cs_join_select_t structure
 *---------------------------------------------------------------------------*/

cs_join_select_t *
cs_join_select_create(const char              *selection_criteria,
                      fvm_periodicity_type_t   perio_type,
                      int                      verbosity);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_select_t structure.
 *
 * parameters:
 *   param       <-- user-defined joining parameters
 *   join_select <-- pointer to pointer to structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_select_destroy(cs_join_param_t     param,
                       cs_join_select_t  **join_select);

/*----------------------------------------------------------------------------
 * Extract vertices from a selection of faces.
 *
 * parameters:
 *   n_select_faces <-- number of selected faces
 *   select_faces   <-- list of faces selected
 *   f2v_idx        <-- "face -> vertex" connect. index
 *   f2v_lst        <-- "face -> vertex" connect. list
 *   n_vertices     <-- number of vertices
 *   n_sel_vertices <-> pointer to the number of selected vertices
 *   sel_vertices   <-> pointer to the list of selected vertices
 *---------------------------------------------------------------------------*/

void
cs_join_extract_vertices(cs_lnum_t         n_select_faces,
                         const cs_lnum_t  *select_faces,
                         const cs_lnum_t  *f2v_idx,
                         const cs_lnum_t  *f2v_lst,
                         cs_lnum_t         n_vertices,
                         cs_lnum_t        *n_select_vertices,
                         cs_lnum_t        *select_vertices[]);

/*----------------------------------------------------------------------------
 * Eliminate redundancies found between two lists of elements.
 * Delete elements in elts[] and keep elements in the reference list.
 *
 * parameters:
 *  n_elts     <-> number of elements in the list to clean
 *  elts       <-> list of elements in the list to clean
 *  n_ref_elts <-- number of elements in the reference list
 *  ref_elts   <-- list of reference elements
 *---------------------------------------------------------------------------*/

void
cs_join_clean_selection(cs_lnum_t  *n_elts,
                        cs_lnum_t  *elts[],
                        cs_lnum_t   n_ref_elts,
                        cs_lnum_t   ref_elts[]);

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
cs_join_build_edges_idx(cs_lnum_t        n_faces,
                        const cs_lnum_t  faces[],
                        const cs_lnum_t  f2v_idx[],
                        const cs_lnum_t  f2v_lst[],
                        cs_lnum_t        v2v_idx[]);

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
cs_join_build_edges_lst(cs_lnum_t        n_faces,
                        const cs_lnum_t  faces[],
                        const cs_lnum_t  f2v_idx[],
                        const cs_lnum_t  f2v_lst[],
                        cs_lnum_t        count[],
                        const cs_lnum_t  v2v_idx[],
                        cs_lnum_t        v2v_lst[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_UTIL_H__ */
