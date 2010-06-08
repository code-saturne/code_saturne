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
 * Set of user parameters to control the join operation
 *----------------------------------------------------------------------------*/

typedef struct {

  int  num;        /* number associated to the current join operation */
  int  perio_num;  /* periodicity number associated to the joining op. */

  /* Octree - Quadtree search algorithm */
  /* ---------------------------------- */

  int    tree_max_level;     /* Deepest level reachable during tree building */
  int    tree_n_max_boxes;   /* Max. number of boxes which can be related to
                               a leaf of the tree if level != tree_max_level */

  float  tree_max_box_ratio; /* Stop tree building if:
                                n_linked_boxes > tree_max_box_ratio*n_init_boxes */

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

  /* Level of display:
       O : no information printed
       1 : general information printed
       2 : more information printed
       5 and beyond : highest level (DEBUG LEVEL) */

  int  verbosity;

} cs_join_param_t;

/*----------------------------------------------------------------------------
 * Set of variables to synchronize single elements
 *---------------------------------------------------------------------------*/

typedef struct {

  int      n_elts;
  int      n_ranks;
  int     *ranks;
  int     *index;
  int     *array;

} cs_join_sync_t;

/*----------------------------------------------------------------------------
 * Structure used to store the result of the extraction of entities
 * implied in the joining operation
 *---------------------------------------------------------------------------*/

typedef struct {

  cs_int_t      n_init_b_faces;  /* Number of border faces before joining */
  cs_int_t      n_init_i_faces;  /* Number of interior faces before joining */
  cs_int_t      n_init_vertices; /* Number of vertices before joining */

  cs_int_t      n_faces;     /* Number of border faces selected
                                for the joining operation */
  fvm_gnum_t    n_g_faces;   /* Global number of border faces selected
                                for the joining operation */
  cs_int_t     *faces;       /* List of selected border faces */

  fvm_gnum_t   *compact_face_gnum;    /* Global face numbering defined
                                         on the selected faces */
  fvm_gnum_t   *compact_rank_index;   /* Distribution of the selected faces
                                         over the ranks */

  fvm_gnum_t   *cell_gnum;       /* Global cell numbering of the cells
                                    bearing the selected face */

  cs_int_t      n_vertices;      /* Number of vertices selected
                                    for the joining operation */
  fvm_gnum_t    n_g_vertices;    /* Global number of selected vertices */
  cs_int_t     *vertices;        /* List of selected vertices */

  /* Adjacent faces of the current face selection: border and interior */

  cs_int_t      n_b_adj_faces;
  cs_int_t      n_i_adj_faces;

  cs_int_t     *b_adj_faces;
  cs_int_t     *i_adj_faces;

  /* Keep the status of all faces of the related cs_mesh_t */

  cs_join_state_t   *b_face_state;
  cs_join_state_t   *i_face_state;

  /* For periodicity handling: list of periodic vertex couples */

  cs_int_t     n_couples;
  fvm_gnum_t  *per_v_couples;

  /*
     Single elements (Only possible in parallel). It appears
     when the domain splitting has a poor quality and elements
     on the joining interface are prisms or tetraedrals
     s = single / c = coupled
  */

  cs_bool_t    do_single_sync;

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

  cs_join_select_t  *selection;  /* Store entities implied in the joining
                                    operation */

  char              *criteria;   /* Criteria used to select border faces
                                    implied in the joining operation */

} cs_join_t;

/*----------------------------------------------------------------------------
 * Structure used to store information about a block distribution
 *---------------------------------------------------------------------------*/

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

/*=============================================================================
 * Global variables
 *===========================================================================*/

extern int  cs_glob_join_count;
extern int  cs_glob_n_joinings;
extern cs_join_t  **cs_glob_join_array;

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_t structure.
 *
 * parameters:
 *   join_number  <-- number related to the joining operation
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- edge fraction tolerance parameter
 *   plane        <-- plane normal angle tolerance
 *   perio_num    <-- periodicity number (0 if not a periodic joining)
 *   verbosity    <-- verbosity level
 *
 * returns:
 *   pointer to a newly allocated cs_join_t structure
 *---------------------------------------------------------------------------*/

cs_join_t *
cs_join_create(int          join_number,
               const char  *sel_criteria,
               float        fraction,
               float        plane,
               int          perio_num,
               int          verbosity);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_t structure.
 *
 * parameters:
 *  join <-> pointer to the cs_join_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_destroy(cs_join_t  **join);

/*----------------------------------------------------------------------------
 * Set advanced parameters to user-defined values.
 *
 * parameters:
 *   join           <-> pointer a to cs_join_t struct. to update
 *   mtf            <-- merge tolerance coefficient
 *   pmf            <-- pre-merge factor
 *   tcm            <-- tolerance computation mode
 *   icm            <-- intersection computation mode
 *   maxbrk         <-- max number of equivalences to break (merge step)
 *   max_sub_faces  <-- max. possible number of sub-faces by splitting a face
 *   tml            <-- tree max level
 *   tmb            <-- tree max boxes
 *   tmr            <-- tree max ratio
 *---------------------------------------------------------------------------*/

void
cs_join_set_advanced_param(cs_join_t   *join,
                           cs_real_t    mtf,
                           cs_real_t    pmf,
                           cs_int_t     tcm,
                           cs_int_t     icm,
                           cs_int_t     maxbrk,
                           cs_int_t     max_sub_faces,
                           cs_int_t     tml,
                           cs_int_t     tmb,
                           cs_real_t    tmr);

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
 *   a newly defined cs_join_block_info_t structure
 *---------------------------------------------------------------------------*/

cs_join_block_info_t
cs_join_get_block_info(fvm_gnum_t  n_g_elts,
                       int         n_ranks,
                       int         local_rank);

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
cs_join_extract_vertices(cs_int_t         n_select_faces,
                         const cs_int_t  *select_faces,
                         const cs_int_t  *f2v_idx,
                         const cs_int_t  *f2v_lst,
                         cs_int_t         n_vertices,
                         cs_int_t        *n_select_vertices,
                         cs_int_t        *select_vertices[]);

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
cs_join_clean_selection(cs_int_t  *n_elts,
                        cs_int_t  *elts[],
                        cs_int_t   n_ref_elts,
                        cs_int_t   ref_elts[]);

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
