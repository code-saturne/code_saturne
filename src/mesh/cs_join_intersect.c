/*============================================================================
 * Set of subroutines useful for intersecting entities during the joining
 * operations:
 *  - Intersections of bounding boxes thanks to a box-tree structure,
 *  - Creation of new vertices from edge intersections,
 *  - Synchronizing intersections in parallel mode.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_neighborhood.h"
#include "fvm_io_num.h"

#include "cs_all_to_all.h"
#include "cs_block_dist.h"
#include "cs_join_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_search.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_intersect.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *===========================================================================*/

/*============================================================================
 * Structure and type definitions
 *===========================================================================*/

/* -------------------------------------------------------------
 * Definition of a structure useful for exchanging intersections
 * between ranks (_sync_edge_inter())
 * ------------------------------------------------------------- */

typedef struct {

  cs_gnum_t   vtx_gnum;   /* vertex global number relative to the inter. */
  cs_coord_t  curv_abs;   /* curvilinear abscissa of the intersection */

} exch_inter_t;

/*============================================================================
 * Global variable definitions
 *===========================================================================*/

static int  _n_inter_tolerance_warnings;

static const double  _cs_join_tol_eps_coef = 1.0001;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (local
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *   b <-> associated array
 *---------------------------------------------------------------------------*/

inline static void
_adapted_lshellsort(cs_lnum_t   l,
                    cs_lnum_t   r,
                    cs_coord_t  a[],
                    cs_lnum_t   b[])
{
  int  i, j, h;
  cs_lnum_t  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_coord_t  va = a[i];
      cs_lnum_t   vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */
}

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (global
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *   b <-> associated array
 *---------------------------------------------------------------------------*/

inline static void
_adapted_gshellsort(cs_lnum_t   l,
                    cs_lnum_t   r,
                    cs_coord_t  a[],
                    cs_gnum_t   b[])
{
  int  i, j, h;
  cs_lnum_t  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_coord_t  va = a[i];
      cs_gnum_t   vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */
}

/*----------------------------------------------------------------------------
 * Compute the dot product of two vectors.
 *
 * parameters:
 *   v1 <-- first vector
 *   v2 <-- second vector
 *
 * returns:
 *   the resulting dot product (v1.v2)
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(const double  v1[],
             const double  v2[])
{
  int  i;
  double  result = 0.0;

  for (i = 0; i < 3; i++)
    result += v1[i] * v2[i];

  return result;
}

/*----------------------------------------------------------------------------
 * Compute the min/max coordinates for the current face in taking into account
 * the fraction parameter.
 *
 * parameters:
 *   face_start    <--  start index in face_vtx_lst
 *   face_end      <--  end index in face_vtx_lst
 *   face_vtx_lst  <--  list of vertices for "face->vertices" connectivity
 *   vertices      <--  array on data associated to each selected vertex
 *   extents       -->  min. and max coordinates of the bounding box
 *---------------------------------------------------------------------------*/

static void
_get_face_extents(const cs_lnum_t          face_start,
                  const cs_lnum_t          face_end,
                  const cs_lnum_t          face_vtx_lst[],
                  const cs_join_vertex_t  *vertices,
                  cs_coord_t               extents[6])
{
  cs_lnum_t  i, j;

  /* Initalization */

  for (i = 0; i < 3; i++) {
    extents[i] =  DBL_MAX;
    extents[3 + i] = -DBL_MAX;
  }

  /* Loop on vertices */

  for (i = face_start; i < face_end; i++) {

    cs_lnum_t  vtx_id = face_vtx_lst[i];
    cs_join_vertex_t  vtx = vertices[vtx_id];

    for (j = 0; j < 3; j++) {
      extents[j] = CS_MIN(extents[j], vtx.coord[j] - vtx.tolerance);
      extents[3 + j] = CS_MAX(extents[3 + j], vtx.coord[j] + vtx.tolerance);
    }

  } /* End of loop on face->vertices connectivity */

}

/*----------------------------------------------------------------------------
 * Compute the length of a segment between two vertices.
 *
 * parameters:
 *   v1 <-- cs_join_vertex_t structure for the first vertex of the segment
 *   v2 <-- cs_join_vertex_t structure for the second vertex of the segment
 *
 * returns:
 *   length of the segment
 *---------------------------------------------------------------------------*/

inline static cs_real_t
_compute_length(cs_join_vertex_t  v1,
                cs_join_vertex_t  v2)
{
  cs_lnum_t  k;
  cs_real_t  len = 0.0, d2 = 0.0;

  for (k = 0; k < 3; k++) {
    cs_real_t  d = v1.coord[k] - v2.coord[k];
    d2 += d * d;
  }
  len = sqrt(d2);

  return len;
}

/*----------------------------------------------------------------------------
 * Compute a new cs_join_vertex_t structure.
 *
 * parameters:
 *   curv_abs   <-- curvilinear abscissa of the intersection
 *   gnum       <-- global number associated to the new
 *                  cs_join_vertex_t structure
 *   vtx_couple <-- couple of vertex numbers defining the current edge
 *   work       <-- local cs_join_mesh_t structure
 *
 * returns:
 *   a new cs_join_vertex_t structure
 *---------------------------------------------------------------------------*/

static cs_join_vertex_t
_get_new_vertex(cs_coord_t             curv_abs,
                cs_gnum_t              gnum,
                const cs_lnum_t       *vtx_couple,
                const cs_join_mesh_t  *work)
{
  cs_lnum_t  k;
  cs_join_vertex_t  new_vtx_data;

#if defined(DEBUG) && !defined(NDEBUG)
  /* Avoid Valgrind warnings in byte copies due to padding */
  memset(&new_vtx_data, 0, sizeof(cs_join_vertex_t));
#endif

  cs_join_vertex_t  v1 = work->vertices[vtx_couple[0]-1];
  cs_join_vertex_t  v2 = work->vertices[vtx_couple[1]-1];

  assert(curv_abs >= 0.0);
  assert(curv_abs <= 1.0);

  /* New global number */

  new_vtx_data.gnum = gnum;
  new_vtx_data.state = CS_JOIN_STATE_NEW;

  /* New tolerance */

  new_vtx_data.tolerance = (1-curv_abs)*v1.tolerance + curv_abs*v2.tolerance;

  /* New coordinates */

  for (k = 0; k < 3; k++)
    new_vtx_data.coord[k] = (1-curv_abs)*v1.coord[k] + curv_abs*v2.coord[k];

  return new_vtx_data;
}

/*----------------------------------------------------------------------------
 * Check the coherency of an equivalence. Return true if vertices are each
 * other under their tolerance.
 *
 * parameters:
 *   edges      <--  cs_join_edges_t structure
 *   mesh       <--  cs_join_mesh_t structure associated
 *   e1_id      <--  first edge implied in the equivalence
 *   curv_abs1  <--  curvilinear abscissa of the intersection on edge 1
 *   e2_id      <--  second edge implied in the equivalence
 *   curv_abs2  <--  curvilinear abscissa of the intersection on edge 2
 *   verbosity  <--  level of accuracy in information display
 *   logfile    <--  handle to log file if verbosity > 1
 *
 * returns:
 *  true if the check is ok, false otherwise
 *---------------------------------------------------------------------------*/

static bool
_check_equiv(const cs_join_edges_t  *edges,
             const cs_join_mesh_t   *mesh,
             cs_lnum_t               e1_id,
             cs_coord_t              curv_abs1,
             cs_lnum_t               e2_id,
             cs_coord_t              curv_abs2,
             int                     verbosity,
             FILE                   *logfile)
{
  cs_join_vertex_t  p1 = _get_new_vertex(curv_abs1, 1,
                                         &(edges->def[2*e1_id]), mesh);
  cs_join_vertex_t  p2 = _get_new_vertex(curv_abs2, 2,
                                          &(edges->def[2*e2_id]), mesh);
  cs_real_t  d12 = _compute_length(p1, p2);

  if (p1.tolerance < d12 || p2.tolerance < d12) {

    cs_lnum_t  v1e1_id = edges->def[2*e1_id]-1;
    cs_lnum_t  v2e1_id = edges->def[2*e1_id+1]-1;
    cs_lnum_t  v1e2_id = edges->def[2*e2_id]-1;
    cs_lnum_t  v2e2_id = edges->def[2*e2_id+1]-1;

    _n_inter_tolerance_warnings++;

    if (verbosity > 3) {

      fprintf(logfile,
              "\n"
              "  Edge - Edge intersection warning between:\n"
              "    edge 1: %d (%llu) [%d (%llu), %d (%llu)]\n"
              "    edge 2: %d (%llu) [%d (%llu), %d (%llu)]\n"
              "  Intersection found for curv. abs. %f (e1) - %f (e2)"
              " will be ignored.\n",
              e1_id+1, (unsigned long long)edges->gnum[e1_id],
              v1e1_id+1, (unsigned long long)mesh->vertices[v1e1_id].gnum,
              v2e1_id+1, (unsigned long long)mesh->vertices[v2e1_id].gnum,
              e2_id+1, (unsigned long long)edges->gnum[e2_id],
              v1e2_id+1, (unsigned long long)mesh->vertices[v1e2_id].gnum,
              v2e2_id+1, (unsigned long long)mesh->vertices[v2e2_id].gnum,
              curv_abs1, curv_abs2);

      if (p1.tolerance < d12 && verbosity > 4)
        fprintf(logfile,
                " Failure for edge 1: "
                " Distance [v_inter1, v_inter2]: %e > v_inter1.tol: %e\n",
                d12, p1.tolerance);

      if (p2.tolerance < d12 && verbosity > 4)
        fprintf(logfile,
                " Failure for edge 2: "
                " Distance [v_inter1, v_inter2]: %e > v_inter2.tol: %e\n",
                d12, p1.tolerance);
    }
    return false;
  }
  else
    return true;
}

/*----------------------------------------------------------------------------
 * Add new equivalences between vertices in cs_join_eset_t structure.
 *
 * parameters:
 *   e1_id     <-- data relative to edge E1
 *   e2_id     <-- data relative to edge E2
 *   abs_e1    <-- curvilinear abscissa of intersection for edge E1
 *   abs_e2    <-- curvilinear abscissa of intersection for edge E2
 *   edges     <-- list of edges
 *   vtx_equiv <-> pointer to a structure dealing with vertex equivalences
 *---------------------------------------------------------------------------*/

static void
_add_trivial_equiv(cs_lnum_t               e1_id,
                   cs_lnum_t               e2_id,
                   cs_coord_t              abs_e1,
                   cs_coord_t              abs_e2,
                   const cs_join_edges_t  *edges,
                   cs_join_eset_t         *vtx_equiv)
{
  cs_lnum_t  v1_num, v2_num;

  cs_lnum_t  equiv_id = vtx_equiv->n_equiv;

  cs_join_eset_check_size(equiv_id, &vtx_equiv);

  /* abs_ei is either 0 or 1, but we avoid exact floating-point comparisons */

  if (abs_e1 < 0.5)
    v1_num = edges->def[2*e1_id];
  else
    v1_num = edges->def[2*e1_id+1];

  if (abs_e2 < 0.5)
    v2_num = edges->def[2*e2_id];
  else
    v2_num = edges->def[2*e2_id+1];

  if (v1_num < v2_num) {
    vtx_equiv->equiv_couple[2*equiv_id] = v1_num;
    vtx_equiv->equiv_couple[2*equiv_id+1] = v2_num;
  }
  else {
    vtx_equiv->equiv_couple[2*equiv_id] = v2_num;
    vtx_equiv->equiv_couple[2*equiv_id+1] = v1_num;
  }

  vtx_equiv->n_equiv += 1;
}

/*----------------------------------------------------------------------------
 * Add a no trivial intersection in a cs_join_inter_set_t structure.
 *
 * parameters:
 *  e1_id     <-- data relative to edge E1
 *  e2_id     <-- data relative to edge E2
 *  abs_e1    <-- curvilinear abscissa of intersection(s) for edge E1
 *  abs_e2    <-- curvilinear abscissa of intersection(s) for edge E2
 *  inter_set <-> pointer to the structure maintaining data about edge-edge
 *                intersections
 *---------------------------------------------------------------------------*/

static void
_add_inter(cs_lnum_t             e1_id,
           cs_lnum_t             e2_id,
           cs_coord_t            abs_e1,
           cs_coord_t            abs_e2,
           cs_join_inter_set_t  *inter_set)
{
  cs_join_inter_t  new_inter_e1, new_inter_e2;

  cs_lnum_t  inter_id = inter_set->n_inter;

  if (inter_id + 1 > inter_set->n_max_inter) {

    inter_set->n_max_inter *= 2;

    BFT_REALLOC(inter_set->inter_lst,
                2*inter_set->n_max_inter,
                cs_join_inter_t);

  }

  new_inter_e1.edge_id = e1_id;
  new_inter_e1.vtx_id = -1;
  new_inter_e1.curv_abs = abs_e1;

  new_inter_e2.edge_id = e2_id;
  new_inter_e2.vtx_id = -1;
  new_inter_e2.curv_abs = abs_e2;

  if (e1_id < e2_id) {
    inter_set->inter_lst[2*inter_id] = new_inter_e1;
    inter_set->inter_lst[2*inter_id+1] = new_inter_e2;
  }
  else {
    inter_set->inter_lst[2*inter_id] = new_inter_e2;
    inter_set->inter_lst[2*inter_id+1] = new_inter_e1;
  }

  inter_set->n_inter += 1;
}

/*----------------------------------------------------------------------------
 * Reduce tolerance for the weakest equivalence in order to break this
 * equivalence
 *
 * parameters:
 *  n_elts      <-- size of vtx_lst (n_inter on edges + the two extremities)
 *  edge_length <-- length of the edge
 *  equiv_lst   <-- true if equiv between two elements else false
 *  abs_lst     <-- curvilinear abscissa of each element
 *  tol_lst     <-> tolerance of each element
 *---------------------------------------------------------------------------*/

static void
_break_equivalence(cs_lnum_t         n_elts,
                   double            edge_length,
                   bool              equiv_lst[],
                   const cs_coord_t  abs_lst[],
                   const double      tol_lst[])
{
  int  i1, i2;
  double  range, _rtf, rtf12, rtf21;

  int  k = 0, i1_save = 0;
  double rtf = -1.0;

  for (i1 = 0; i1 < n_elts - 1; i1++) {

    i2 = i1 + 1;

    if (i1 == 0)
      k = 0;
    else
      k += n_elts - i1;

    if (equiv_lst[k] == true) {

      range = (abs_lst[i2] - abs_lst[i1]) * edge_length;
      rtf12 = range/tol_lst[i1];
      rtf21 = range/tol_lst[i2];

      assert(range >= 0.0);
      assert(rtf12 < 1.0);
      assert(rtf21 < 1.0);

      _rtf = CS_MAX(rtf12, rtf21);
      if (_rtf > rtf) {
        i1_save = i1;
        rtf = _rtf;
      }

    }

  } /* End of loop on i1 */

  if (rtf > 0.0) { /* We have find an equivalence to break */

    /* Break equivalence between i1_save and i2_save which the weakest
       equivalence between two consecutive vertices in the set */

    for (i1 = 0; i1 < n_elts - 1; i1++) {

      i2 = i1 + 1;

      if (i1 == 0)
        k = 0;
      else
        k += n_elts - i1;

      if (i1 == i1_save) {
#if 0 && defined(DEBUG) && !defined(NDEBUG)
        fprintf(cs_glob_join_log,
                " Break equivalence between [%d, %d]\n", i1, i1_save);
#endif
        for (i2 = i1 + 1; i2 < n_elts; i2++, k++)
          equiv_lst[k] = false;
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Fill equiv_lst[] according to the tolerance related to each vertex
 *
 * parameters:
 *  param       <-- set of user-defined parameters
 *  n_elts      <-- size of vtx_lst (n_inter on edges + the two extremities)
 *  abs_lst     <-- curvilinear abscissa of each element
 *  tol_lst     <-> tolerance of each element
 *  equiv_lst   <-> true if equiv between two elements else false
 *  tag         <-- array used to mark equivalent vertices by the same tag
 *  edge_length <-- length of the edge
 *
 * returns:
 *  number of tolerance reductions applied
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_find_edge_equiv(cs_join_param_t  param,
                 cs_lnum_t        n_elts,
                 cs_coord_t       abs_lst[],
                 double           tol_lst[],
                 bool             equiv_lst[],
                 int              tag[],
                 double           edge_length)
{
  int  i, i1, i2, k;
  double  dist;

  cs_lnum_t  n_loops = 0;
  bool  do_reduction = true;

  /* Find equivalence between vertices sharing the same edge */

  for (i1 = 0, k = 0; i1 < n_elts - 1; i1++) {
    for (i2 = i1 + 1; i2 < n_elts; i2++, k++) {

      dist = (abs_lst[i2] - abs_lst[i1]) * edge_length;
      assert(dist >= 0.0);

      /* Tag equivalence if test is true */

      if (tol_lst[i1] < dist || tol_lst[i2] < dist)
        equiv_lst[k] = false;
      else
        equiv_lst[k] = true;

    }
  }

  while (do_reduction == true && n_loops <= param.n_max_equiv_breaks) {

    /* Reduce tolerance after the first loop if necessary */

    if (do_reduction == true && n_loops > 0)
      _break_equivalence(n_elts,
                         edge_length,
                         equiv_lst,
                         abs_lst,
                         tol_lst);

    /* Tag vertices in equivalence with the same even if they are
       in equivalence thanks to transitivity */

    for (i = 0; i < n_elts; i++)  /* Re-initialize array */
      tag[i] = i+1;

    for (i1 = 0; i1 < n_elts-1; i1++) {

      if (i1 == 0)
        k = 0;
      else
        k += n_elts - i1;

      if (equiv_lst[k] == true)
        tag[i1+1] = tag[i1];

    }

    /* Check tolerance consistency: avoid excessive transitivity
       on an edge */

    do_reduction = false;
    for (i1 = 0, k = 0; i1 < n_elts - 1; i1++) {
      for (i2 = i1 + 1; i2 < n_elts; i2++, k++) {

        if (equiv_lst[k] == false)
          if (tag[i1] == tag[i2])
            do_reduction = true;

      } /* End of loop on i2 */
    } /* End of loop on i1 */

    n_loops++;

  } /* End of while reduc_tol == false */

  return n_loops - 1;
}

/*----------------------------------------------------------------------------
 * Get 3D intersection location(s) between two edges in curvilinear abscissa.
 *
 * We successively try to find vertices-vertices matching,
 *                             vertex-vertex matching,
 *                             a "real" intersection.
 *
 * Each vertex owns a tolerance which drives to a sphere in which intersection
 * is possible.
 * When an intersection is found, we store the related curvilinear abscissa
 * for each edge.
 * If curvilinear abscissa : 0 => -0.01 / 1 => 1.01 in order to test an
 * inequality rather than an equality.
 *
 * Let be two edges E1 (P1E1, P2E1) and E2 (P1E2, P2E2) :
 * a point on edge E1 is defined by :
 *   P1E1 + s * (P2E1 - P1E1)   with 0 <= s <= 1
 * a point on edge B is defined by :
 *   P1E2 + t * (P2E2 - P1E2)   with 0 <= t <= 1
 *
 * The squared distance between a point from A and a point from B is :
 *  len(s, t) = || P1E2 - P1E1 + t * (P2E2 - P1E1) - s * (P2E1 - P1E1) ||^2
 * equivalent to :
 *   len(s, t) = a.s^2 + 2b.s.t + c.t^2 + 2d.s + 2e.t + f = d2_e1e2(s,t)
 *
 * We try to find (s, t) such as len(s,t) is minimal
 * with 0 <= s <= 1 and 0 <= t <= 1
 *
 *       ->   ->
 * a =   v0 . v0                  ->
 *                                v0
 *       ->   ->     P1E1  x--------------> P2E1
 * b = - v0 . v1            \
 *                           \               P2E2
 *       ->   ->              \             ^
 * c =   v1 . v1               \           /
 *                           -> \         /
 *       ->   ->             v2  \       / ->
 * d = - v0 . v2                  \     /  v1
 *                                 \   /
 *       ->   ->                    \ /
 * e =   v1 . v2                     x
 *                                P1E2
 *
 *    ->
 *    v0 = vector (P1E1, P2E1)
 *    ->
 *    v1 = vector (P1E2, P2E2)
 *    ->
 *    v2 = vector (P1E1, P1E2)
 *
 * parameters:
 *   mesh       <-- pointer to joining mesh
 *   edges      <-- pointer to edges
 *   fraction   <--  global tolerance for geometrical intersection.
 *   e1_id      <--  id of edge 1
 *   abs_e1     <--  intersection location on E1 (curvilinear abscissa)
 *   e2_id      <--  id of edge 2
 *   abs_e2     <--  intersection location on E2 (curvilinear abscissa)
 *   verbosity  <--  level of accuracy in information display
 *   logfile    <--  handle to optional log file
 *   n_inter    <->  number of intersections detected.
 *---------------------------------------------------------------------------*/

static void
_new_edge_edge_3d_inter(const cs_join_mesh_t   *mesh,
                        const cs_join_edges_t  *edges,
                        double                  fraction,
                        cs_lnum_t               e1_id,
                        double                  abs_e1[2],
                        cs_lnum_t               e2_id,
                        double                  abs_e2[2],
                        double                  parall_eps2,
                        int                     verbosity,
                        FILE                   *logfile,
                        int                    *n_inter)
{
  int  k;
  double  a, b, c, d, e, f, s, t;
  double  d2_limit_e1, d_limit_e1, d2_limit_e2, d_limit_e2, d2_e1e2;
  double  inv_cross_norm2, cross_norm2;
  double  int_inter[2] = {0., 0.}, ext_inter[4], v0[3], v1[3], v2[3];

  int  n_int_inter = 0, n_ext_inter = 0;
  bool  int_p1e2 = false, int_p2e2 = false;

  cs_lnum_t  p1e1_id = edges->def[2*e1_id] - 1;
  cs_lnum_t  p2e1_id = edges->def[2*e1_id+1] - 1;
  cs_lnum_t  p1e2_id = edges->def[2*e2_id] - 1;
  cs_lnum_t  p2e2_id = edges->def[2*e2_id+1] - 1;

  const cs_join_vertex_t  p1e1 = mesh->vertices[p1e1_id];
  const cs_join_vertex_t  p2e1 = mesh->vertices[p2e1_id];
  const cs_join_vertex_t  p1e2 = mesh->vertices[p1e2_id];
  const cs_join_vertex_t  p2e2 = mesh->vertices[p2e2_id];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bool  tst_dbg = (verbosity > 6 &&
                   (p1e1.gnum == 716852 || p2e1.gnum == 716852 ||
                    p1e2.gnum == 716852 || p2e2.gnum == 716852) ?
                   true : false);
#endif

  /* Initialize parameters */

  *n_inter = 0;

  assert(p1e1.gnum != p2e1.gnum);
  assert(p1e2.gnum != p2e2.gnum);

  if (p1e1.gnum == p1e2.gnum || p2e1.gnum == p2e2.gnum)
    return;

  if (p1e1.gnum == p2e2.gnum || p2e1.gnum == p1e2.gnum)
    return;

  /* Compute local vectors and parameters */

  for (k = 0; k < 3; k++) {
    v0[k] = p2e1.coord[k] - p1e1.coord[k];
    v1[k] = p2e2.coord[k] - p1e2.coord[k];
    v2[k] = p1e2.coord[k] - p1e1.coord[k];
  }

  a =   _dot_product(v0, v0);
  b = - _dot_product(v0, v1);
  c =   _dot_product(v1, v1);
  d = - _dot_product(v0, v2);
  e =   _dot_product(v1, v2);
  f =   _dot_product(v2, v2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL) {
    fprintf(logfile,
            "\n\np1e1 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p1e1.gnum,
            p1e1.coord[0], p1e1.coord[1], p1e1.coord[2], p1e1.tolerance);
    fprintf(logfile,
            "p2e1 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p2e1.gnum,
            p2e1.coord[0], p2e1.coord[1], p2e1.coord[2], p2e1.tolerance);
    fprintf(logfile,
            "p1e2 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p1e2.gnum,
            p1e2.coord[0], p1e2.coord[1], p1e2.coord[2], p1e2.tolerance);
    fprintf(logfile,
            "p2e2 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n\n",
            (unsigned long long)p2e2.gnum,
            p2e2.coord[0], p2e2.coord[1], p2e2.coord[2], p2e2.tolerance);
    fprintf(logfile,
            "v0 : [ %10.8e %10.8e %10.8e]\n", v0[0], v0[1], v0[2]);
    fprintf(logfile,
            "v1 : [ %10.8e %10.8e %10.8e]\n", v1[0], v1[1], v1[2]);
    fprintf(logfile, "v2 : [ %10.8e %10.8e %10.8e]\n\n", v2[0], v2[1], v2[2]);
    fprintf(logfile, "a : %10.8e, b : %10.8e, c : %10.8e\n", a, b, c);
    fprintf(logfile, "d : %10.8e, e : %10.8e, f : %10.8e\n", d, e, f);
  }
#endif

  /* Check size of each edge is not equal to 0 */

  assert(a > 0);
  assert(c > 0);

  /* Check computation of the tolerance */

  assert(sqrt(a) * fraction * 1.00001 >= p1e1.tolerance);
  assert(sqrt(a) * fraction * 1.00001 >= p2e1.tolerance);
  assert(sqrt(c) * fraction * 1.00001 >= p1e2.tolerance);
  assert(sqrt(c) * fraction * 1.00001 >= p2e2.tolerance);

  /*------------------------------------------------------------------
   * First part : Search for interior/interior intersection
   * The only possibility is an intersection with 0 < s < 1 and
   * 0 < t < 1.
   *------------------------------------------------------------------*/

  /* Edges are parallel if a.c - b^2 = 0.
   * However, this term is o(length^4). So, we will prefer work with angle
   * If two edges are not parallel => a.c - b^2 > 0.
   *
   * cos(theta) = v0.v1 / ||v0||.||v1||
   * cos^2 (theta) = (v0.v1)^2 / ((v0.v0).(v1.v1))^2 = (-b)^2 / (a.c)
   * sin^2 (theta) = 1 - b^2 / (a.c)
   *
   * || v0 ^ v1 ||^2 = ||v0||^2 . ||v1||^2 . sin^2(thetha)
   *
   */

  cross_norm2 = CS_ABS(a * c - b * b);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile,
            " [I]: Test inter. E1 - E2: "
            "\t cross_norm2: %12.8e - parall_limit: %12.8e\n",
            cross_norm2, parall_eps2*a*c);
#endif

  if (cross_norm2 > parall_eps2 * a * c) {

    /* <=> 1 - b^2/(a.c) < epsilon^2 <=> sin^2(theta)  < epsilon^2

       Edges are not parallel.
       Define s and t if there is an intersection on interior points.

       ->     ->     ->
       v2 + s.v0 = t.v1 and cross_norm2 = a.c - b^2 != 0.0

       => s = (b.e - c.d)/cross_norm2
       => t = (e.a - b.d)/cross_norm2
    */

    s = b * e - c * d;
    t = b * d - a * e;
    inv_cross_norm2 = 1.0 / cross_norm2;
    s *= inv_cross_norm2;
    t *= inv_cross_norm2;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile,
            " [I]-1: Test inter. E1 - E2: "
            "\t s: %12.8e - t: %12.8e\n", s, t);
#endif

    if (s >= 0. && s <= 1.0) {
      if (t >= 0. && t <= 1.0)  {

        /* If tests are OK, we are on an interior point for E1 and E2 */

        d2_e1e2 = CS_ABS(s*(a*s + 2.0*(b*t + d)) + t*(c*t + 2.0*e) + f);
        d_limit_e1 = (1.0-s) * p1e1.tolerance + s * p2e1.tolerance;
        d_limit_e2 = (1.0-t) * p1e2.tolerance + t * p2e2.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;
        d2_limit_e2 = d_limit_e2 * d_limit_e2;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[I]-2: s = %10.8e, t = %10.8e\n"
                  "\t  d2_e1e2: %10.8e, d2_limit_e1: %10.8e, "
                  "d2_limit_e2: %10.8e\n",
                  s, t, d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

          /* Under tolerance for edges E1 and E2 => intersection is possible */

          assert(0.0 <= s && s <= 1.0);
          assert(0.0 <= t && t <= 1.0);

          n_int_inter = 1;
          int_inter[0] = s;
          int_inter[1] = t;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[I]-3: Add int. inter. Edge-Edge (%10.8e,%10.8e)\n",
                    s, t);
#endif

        } /* If we are under tolerance */

      } /* If it is an interior point for edge E2 */
    } /* If it is an interior point for edge E1 */

  } /* End of last test : inter. E1 - E2 */

  /*-----------------------------------------------------------------
   * Second part : we try to find under the current tolerance if
   * "parallel edges", "extremity-extremity" or "extremity-interior"
   * are possible.
   *-----------------------------------------------------------------*/

  /* Distance between P1E1 (s = 0) and a point on E2 */
  /* ----------------------------------------------- */

  t = -e / c; /* Nearest point on E2 if s = 0 => t = -e/c */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile, " [II] Test inter. P1E1 - E2: t = %10.8e\n", t);
#endif

  if (t > -fraction && t < 1.0 + fraction) { /* filter */

    d2_e1e2 = t*(c*t + 2.0*e) + f;
    assert(d2_e1e2 >= -1.e-8);

    /* Tolerance for vertex p1e1 (s = 0) */

    d2_limit_e1 = p1e1.tolerance * p1e1.tolerance;

    if (d2_e1e2 <= d2_limit_e1) { /* Under the tolerance for vertex P1E1 */

      /* Tolerance for a vertex on E2 located at "t" */

      d_limit_e2 = (1-t)*p1e2.tolerance + t*p2e2.tolerance;
      d2_limit_e2 = d_limit_e2 * d_limit_e2;

      if (d2_e1e2 <= d2_limit_e2) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[II]-3: Under tol. for P1E1 and a point of E2\n"
                  "\td2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
                  d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        /* We are under tolerance for edge E2. We try to find if it is
           a vertex-vertex intersection */

        if (t < 0.0) { /* Test intersection P1E1-P1E2 */

          d2_e1e2 = f; /* s = t = 0.0 */
          d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 0.0, e2_id, 0.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = t = 0.0
                 under the tolerance is possible */

              int_p1e2 = true;
              n_ext_inter = 1;
              ext_inter[0] = 0.0;
              ext_inter[1] = 0.0;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile, "\t[II]-4: Add inter. Vtx-Vtx (0.00, 0.00)\n");
#endif
            }
          }
        }
        else if (t > 1.0) { /* Test intersection P1E1-P2E2 */

          d2_e1e2 = c + 2.0*e + f; /* s = 0.0 and t = 1.0 */
          d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 0.0, e2_id, 1.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = 0.0, t = 1.0
                 under the tolerance is possible */

              int_p2e2 = true;
              n_ext_inter = 1;
              ext_inter[0] = 0.0;
              ext_inter[1] = 1.0;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[II]-5: Add inter. Vtx-Vtx (0.00, 1.00)\n");
#endif

            }
          }
        }
        else {

          assert(0.0 <= t && t <= 1.0);

          /* It's an "extremity-interior" intersection */

          n_ext_inter = 1;
          ext_inter[0] = 0.0;
          ext_inter[1] = t;

          if (t <= 0.0)
            int_p1e2 = true;
          if (t >= 1.0)
            int_p2e2 = true;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[II]-6: Add inter. Vtx-Edge (0.00, %10.8e)\n", t);
#endif

        }

      } /* End if d2_e1e2 <= d2_limit_e2 */

    } /* End if d2_e1e2 <= d2_limit_e1 */

  } /* End if -fraction < t < 1 + fraction */

  /* Distance between P2E1 (s = 1) and a point on edge E2 */
  /* ---------------------------------------------------- */

  t = -(b + e) / c; /* Minimum for d2_e1e2 is reached for t = -(b+e) / c
                       when s = 1. */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile,
            " [III] Test inter. P2E1 - E2: t = %10.8e\n", t);
#endif

  if (t > -fraction && t < 1.0 + fraction) { /* filter */

    d2_e1e2 = CS_ABS((a + 2.0*(b*t + d)) + t*(c*t + 2.0*e) + f);

    assert(d2_e1e2 >= -1.e-8);

    /* Tolerance for vertex p2e1 (s = 1) */

    d2_limit_e1 = p2e1.tolerance * p2e1.tolerance;

    if (d2_e1e2 <= d2_limit_e1) { /* Under the tolerance for vertex P2E1 */

      d_limit_e2 = (1.0-t)*p1e2.tolerance + t*p2e2.tolerance;
      d2_limit_e2 = d_limit_e2 * d_limit_e2;

      if (d2_e1e2 <= d2_limit_e2) { /* Under tolerance for edge 2 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[III]-3: Under tol. for P2E1 and a point of E2\n"
                  "\td2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
                  d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        /* We are under tolerance for edge E2. We try to find if it is
           a vertex-vertex intersection */

        if (t < 0.0) { /* Test intersection P2E1-P1E2 */

          d2_e1e2 = CS_ABS(a + 2.0*d + f);
          d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 1.0, e2_id, 0.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = 1.0, t = 0.0
                 under the tolerance is possible */

              int_p1e2 = true;
              ext_inter[2*n_ext_inter] = 1.0;
              ext_inter[2*n_ext_inter+1] = 0.0;
              n_ext_inter += 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[III]-4: Add inter. Vtx-Vtx (1.00, 0.00)\n");
#endif

            }
          }
        }
        else if (t > 1.0) { /* Test intersection P2E1-P2E2 */

          d2_e1e2 = CS_ABS(a + 2.0*(b + d + e) + c + f);
          d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 1.0, e2_id, 1.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = t = 1.0
                 under the tolerance is possible */

              int_p2e2 = true;
              ext_inter[2*n_ext_inter] = 1.0;
              ext_inter[2*n_ext_inter+1] = 1.0;
              n_ext_inter += 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[III]-5: Add inter. Vtx-Vtx (1.00, 1.00)\n");
#endif

            }
          }
        }
        else {

          assert(0.0 <= t && t <= 1.0);

          /* It's an "extremity-interior" intersection */

          ext_inter[2*n_ext_inter] = 1.0;
          ext_inter[2*n_ext_inter+1] = t;
          n_ext_inter += 1;

          if (t <= 0.0)
            int_p1e2 = true;
          if (t >= 1.0)
            int_p2e2 = true;

          assert(n_ext_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[III]-6: Add inter. Vtx-Edge (1.00, %10.8e)\n", t);
#endif
        }

      } /* End if d2_e1e2 <= d2_limit_e2 */

    } /* End if d2_e1e2 <= d2_limit_e1 */

  } /* End if -fraction < t < 1 + fraction */

  /* If vertices from edge E2 are not implied in an intersection
     we do a test for each vertex */

  if (int_p1e2 == false && n_ext_inter < 2) {

    /* Distance between P1E2 (t = 0.0) on edge E2 and a point on edge E1 */
    /* ----------------------------------------------------------------- */

    s = -d / a; /* t = 0.0 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              " [IV] Test inter. P1E2 - E1: s = %10.8e\n", s);
#endif

    if (s > 0.0 && s < 1.0) { /* s = 0.0 and s = 1.0 are already tested */

      d2_e1e2 = CS_ABS(s*(a*s + 2.0*d) + f);
      d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              "\t [IV-1a] d2_e1e2: %10.8e - d2_lim_e2: %10.8e\n",
              d2_e1e2, d2_limit_e2);
#endif

      if (d2_e1e2 <= d2_limit_e2) { /* Under the tolerance for vertex P1E2 */

        d_limit_e1 = (1.0-s)*p1e1.tolerance + s*p2e1.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t [IV-1b] d2_e1e2: %10.8e - d2_lim_e1: %10.8e\n",
                  d2_e1e2, d2_limit_e1);
#endif

        if (d2_e1e2 <= d2_limit_e1) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf
              (logfile,
               "\t[IV]-2: Under tol. for P1E2 and a point of E1\n"
               "\t d2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
               d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

          /* We are under tolerance for edge E1. There is an intersection
             between P1E2 and a point on edge E1 */

          ext_inter[2*n_ext_inter] = s;
          ext_inter[2*n_ext_inter+1] = 0.0;
          n_ext_inter += 1;
          assert(n_ext_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[IV]-3: Add inter. Edge-Vtx (%10.8e, 0.00)\n", s);
#endif

        } /* End if d2_e1e2 < d2_limit_e1 */

      } /* End if d2_e1e2 < d2_limit_e2 */

    } /* 0.0 < s < 1.0 */

  } /* If int_p1e2 == false && n_ext_inter < 2 */

  if (int_p2e2 == false && n_ext_inter < 2) {

    /* Distance between P2E2 (t = 1.0) on edge E2 and a point on edge E1 */
    /* ----------------------------------------------------------------- */

    s = -(b + d) / a; /* t = 1 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              " [V] Test inter. P2E2 - E1: s = %10.8e\n", s);
#endif

    if (s > 0.0 && s < 1.0) { /* s = 0.0 and s = 1.0 are already tested */

      d2_e1e2 = CS_ABS(s*(a*s + 2.0*(b + d)) + c + 2.0*e + f);
      d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              "\t [V-1a] d2_e1e2: %10.8e - d2_lim_e2: %10.8e\n",
              d2_e1e2, d2_limit_e2);
#endif

      if (d2_e1e2 <= d2_limit_e2) { /* Under the tolerance for vertex P2E2 */

        d_limit_e1 = (1.0-s)*p1e1.tolerance + s*p2e1.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t [V-1b] d2_e1e2: %10.8e - d2_lim_e1: %10.8e\n",
                  d2_e1e2, d2_limit_e1);
#endif

        if (d2_e1e2 <= d2_limit_e1) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf
              (logfile,
               "\t[V]-2: Under tol. for P2E2 and a point of E1\n"
               "\t d2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
               d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

          /* We are under tolerance for edge E1. There is an intersection
             between P2E2 and a point on edge E1 */

          ext_inter[2*n_ext_inter] = s;
          ext_inter[2*n_ext_inter+1] = 1.0;
          n_ext_inter += 1;
          assert(n_ext_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[V]-3: Add inter (%10.8e, 1.01)\n", s);
#endif

        } /* End if d2_e1e2 < d2_limit_e1 */

      } /* End if d2_e1e2 < d2_limit_e2 */

    } /* 0.0 < s < 1.0 */

  } /* If int_p2e2 == false && n_ext_inter < 2 */

  /* Define intersection(s) to return */

  if (n_int_inter == 1) {

    if (n_ext_inter < 2) { /* Interior intersection is kept */
      *n_inter = 1;
      abs_e1[0] = int_inter[0];
      abs_e2[0] = int_inter[1];
    }
    else { /* A choice between interior or extremity intersection
              has to be done */

      double  ds = CS_ABS(ext_inter[2] - ext_inter[0]);
      double  dt = CS_ABS(ext_inter[3] - ext_inter[1]);

      assert(n_ext_inter == 2);

      if (ds > fraction && dt > fraction) { /* Keep extremity inter. */
        *n_inter = 2;
        for (k = 0; k < 2; k++) {
          abs_e1[k] = ext_inter[2*k];
          abs_e2[k] = ext_inter[2*k+1];
        }
      }
      else { /* Keep interior inter. */
        *n_inter = 1;
        abs_e1[0] = int_inter[0];
        abs_e2[0] = int_inter[1];
      }

    }

  }
  else { /* n_int_inter == 0 */

    *n_inter = n_ext_inter;

    for (k = 0; k < n_ext_inter; k++) {
      abs_e1[k] = ext_inter[2*k];
      abs_e2[k] = ext_inter[2*k+1];
    }

  }

}

/*----------------------------------------------------------------------------
 * Get 3D intersection location(s) between two edges in curvilinear abscissa.
 *
 * We successively try to find vertices-vertices matching,
 *                             vertex-vertex matching,
 *                             a "real" intersection.
 *
 * Each vertex owns a tolerance which drives to a sphere in which intersection
 * is possible.
 * When an intersection is found, we store the related curvilinear abscissa
 * for each edge.
 * If curvilinear abscissa : 0 => -0.01 / 1 => 1.01 in order to test an
 * inequality rather than an equality.
 *
 * Let be two edges E1 (P1E1, P2E1) and E2 (P1E2, P2E2) :
 * a point on edge E1 is defined by :
 *   P1E1 + s * (P2E1 - P1E1)   with 0 <= s <= 1
 * a point on edge B is defined by :
 *   P1E2 + t * (P2E2 - P1E2)   with 0 <= t <= 1
 *
 * The squared distance between a point from A and a point from B is :
 *  len(s, t) = || P1E2 - P1E1 + t * (P2E2 - P1E1) - s * (P2E1 - P1E1) ||^2
 * equivalent to :
 *   len(s, t) = a.s^2 + 2b.s.t + c.t^2 + 2d.s + 2e.t + f = d2_e1e2(s,t)
 *
 * We try to find (s, t) such as len(s,t) is minimal
 * with 0 <= s <= 1 and 0 <= t <= 1
 *
 *       ->   ->
 * a =   v0 . v0                  ->
 *                                v0
 *       ->   ->     P1E1  x--------------> P2E1
 * b = - v0 . v1            \
 *                           \               P2E2
 *       ->   ->              \             ^
 * c =   v1 . v1               \           /
 *                           -> \         /
 *       ->   ->             v2  \       / ->
 * d = - v0 . v2                  \     /  v1
 *                                 \   /
 *       ->   ->                    \ /
 * e =   v1 . v2                     x
 *                                P1E2
 *
 *    ->
 *    v0 = vector (P1E1, P2E1)
 *    ->
 *    v1 = vector (P1E2, P2E2)
 *    ->
 *    v2 = vector (P1E1, P1E2)
 *
 * parameters:
 *   mesh       <-- pointer to joining mesh
 *   edges      <-- pointer to edges
 *   fraction   <--  global tolerance for geometrical intersection.
 *   e1_id      <--  id of edge 1
 *   abs_e1     <--  intersection location on E1 (curvilinear abscissa)
 *   e2_id      <--  id of edge 2
 *   abs_e2     <--  intersection location on E2 (curvilinear abscissa)
 *   verbosity  <--  level of detail in information display
 *   logfile    <--  handle to log file if verbosity > 1
 *   n_inter    <->  number of intersections detected.
 *---------------------------------------------------------------------------*/

static void
_edge_edge_3d_inter(const cs_join_mesh_t   *mesh,
                    const cs_join_edges_t  *edges,
                    double                  fraction,
                    cs_lnum_t               e1_id,
                    double                  abs_e1[2],
                    cs_lnum_t               e2_id,
                    double                  abs_e2[2],
                    double                  parall_eps2,
                    int                     verbosity,
                    FILE                   *logfile,
                    int                    *n_inter)
{
  int  k;
  double  a, b, c, d, e, f, s, t;
  double  d2_limit_e1, d_limit_e1, d2_limit_e2, d_limit_e2, d2_e1e2;
  cs_real_t v0[3], v1[3], v2[3];

  bool  int_p1e2 = false, int_p2e2 = false;

  cs_lnum_t  p1e1_id = edges->def[2*e1_id] - 1;
  cs_lnum_t  p2e1_id = edges->def[2*e1_id+1] - 1;
  cs_lnum_t  p1e2_id = edges->def[2*e2_id] - 1;
  cs_lnum_t  p2e2_id = edges->def[2*e2_id+1] - 1;

  const cs_join_vertex_t  p1e1 = mesh->vertices[p1e1_id];
  const cs_join_vertex_t  p2e1 = mesh->vertices[p2e1_id];
  const cs_join_vertex_t  p1e2 = mesh->vertices[p1e2_id];
  const cs_join_vertex_t  p2e2 = mesh->vertices[p2e2_id];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bool  tst_dbg = (verbosity > 6 &&
                   (p1e1.gnum == 716852 || p2e1.gnum == 716852 ||
                    p1e2.gnum == 716852 || p2e2.gnum == 716852) ?
                   true : false);
#endif

  /* Initialize parameters */

  *n_inter = 0 ;

  assert(p1e1.gnum != p2e1.gnum);
  assert(p1e2.gnum != p2e2.gnum);

  if (p1e1.gnum == p1e2.gnum || p2e1.gnum == p2e2.gnum)
    return;

  if (p1e1.gnum == p2e2.gnum || p2e1.gnum == p1e2.gnum)
    return;

  /* Compute local vectors and parameters */

  for (k = 0; k < 3; k++) {
    v0[k] = p2e1.coord[k] - p1e1.coord[k];
    v1[k] = p2e2.coord[k] - p1e2.coord[k];
    v2[k] = p1e2.coord[k] - p1e1.coord[k];
  }

  a =   _dot_product(v0, v0);
  b = - _dot_product(v0, v1);
  c =   _dot_product(v1, v1);
  d = - _dot_product(v0, v2);
  e =   _dot_product(v1, v2);
  f =   _dot_product(v2, v2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL) {
    fprintf(logfile,
            "\n\np1e1 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p1e1.gnum,
            p1e1.coord[0], p1e1.coord[1], p1e1.coord[2], p1e1.tolerance);
    fprintf(logfile,
            "p2e1 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p2e1.gnum,
            p2e1.coord[0], p2e1.coord[1], p2e1.coord[2], p2e1.tolerance);
    fprintf(logfile,
            "p1e2 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n",
            (unsigned long long)p1e2.gnum,
            p1e2.coord[0], p1e2.coord[1], p1e2.coord[2], p1e2.tolerance);
    fprintf(logfile,
            "p2e2 : %10llu - [%10.8e %10.8e %10.8e] - tol: %10.8g\n\n",
            (unsigned long long)p2e2.gnum,
            p2e2.coord[0], p2e2.coord[1], p2e2.coord[2], p2e2.tolerance);
    fprintf(logfile, "v0 : [ %10.8e %10.8e %10.8e]\n", v0[0], v0[1], v0[2]);
    fprintf(logfile, "v1 : [ %10.8e %10.8e %10.8e]\n", v1[0], v1[1], v1[2]);
    fprintf(logfile, "v2 : [ %10.8e %10.8e %10.8e]\n\n", v2[0], v2[1], v2[2]);
    fprintf(logfile, "a : %10.8e, b : %10.8e, c : %10.8e\n", a, b, c);
    fprintf(logfile, "d : %10.8e, e : %10.8e, f : %10.8e\n", d, e, f);
  }
#endif

  /* Check size of each edge is not equal to 0 */

  assert(a > 0);
  assert(c > 0);

  /* Check computation of the tolerance */

  assert(sqrt(a) * fraction * 1.00001 >= p1e1.tolerance);
  assert(sqrt(a) * fraction * 1.00001 >= p2e1.tolerance);
  assert(sqrt(c) * fraction * 1.00001 >= p1e2.tolerance);
  assert(sqrt(c) * fraction * 1.00001 >= p2e2.tolerance);

  /*-----------------------------------------------------------------
   * First part : we try to find under the current tolerance if
   * "parallel edges", "extremity-extremity" or "extremity-interior"
   * are possible.
   *-----------------------------------------------------------------*/

  /* Distance between P1E1 (s = 0) and a point on E2 */
  /* ----------------------------------------------- */

  t = -e / c; /* Nearest point on E2 if s = 0 => t = -e/c */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile, " [I] Test inter. P1E1 - E2: t = %10.8e\n", t);
#endif

  if (t > -fraction && t < 1.0 + fraction) { /* filter */

    d2_e1e2 = t*(c*t + 2.0*e) + f;
    assert(d2_e1e2 >= -1.e-8);

    /* Tolerance for vertex p1e1 (s = 0) */

    d2_limit_e1 = p1e1.tolerance * p1e1.tolerance;

    if (d2_e1e2 <= d2_limit_e1) { /* Under the tolerance for vertex P1E1 */

      /* Tolerance for a vertex on E2 located at "t" */

      d_limit_e2 = (1-t)*p1e2.tolerance + t*p2e2.tolerance;
      d2_limit_e2 = d_limit_e2 * d_limit_e2;

      if (d2_e1e2 <= d2_limit_e2) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[I]-3: Under tol. for P1E1 and a point of E2\n"
                  "\td2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
                  d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        /* We are under tolerance for edge E2. We try to find if it is
           a vertex-vertex intersection */

        if (t < 0.0) { /* Test intersection P1E1-P1E2 */

          d2_e1e2 = f; /* s = t = 0.0 */
          d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 0.0, e2_id, 0.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = t = 0.0
                 under the tolerance is possible */

              int_p1e2 = true;
              abs_e1[*n_inter] = 0.0;
              abs_e2[*n_inter] = 0.0;

              *n_inter += 1;
              assert(*n_inter == 1);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[I]-4: Add inter. Vtx-Vtx (0.00, 0.00)\n");
#endif
            }
          }
        }
        else if (t > 1.0) { /* Test intersection P1E1-P2E2 */

          d2_e1e2 = c + 2.0*e + f; /* s = 0.0 and t = 1.0 */
          d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 0.0, e2_id, 1.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = 0.0, t = 1.0
                 under the tolerance is possible */

              int_p2e2 = true;
              abs_e1[*n_inter] =  0.0;
              abs_e2[*n_inter] =  1.0;

              *n_inter += 1;
              assert(*n_inter == 1);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[I]-5: Add inter. Vtx-Vtx (0.00, 1.00)\n");
#endif

            }
          }
        }
        else {

          assert(0.0 <= t && t <= 1.0);

          /* It's an "extremity-interior" intersection */

          abs_e1[*n_inter] = 0.00;
          abs_e2[*n_inter] = t;

          if (abs_e2[*n_inter] <= 0.0)
            int_p1e2 = true;
          if (abs_e2[*n_inter] >= 1.0)
            int_p2e2 = true;

          *n_inter += 1;
          assert(*n_inter == 1);


#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[I]-6: Add inter. Vtx-Edge (0.00, %10.8e)\n", t);
#endif

        }

      } /* End if d2_e1e2 <= d2_limit_e2 */

    } /* End if d2_e1e2 <= d2_limit_e1 */

  } /* End if -fraction < t < 1 + fraction */

  /* Distance between P2E1 (s = 1) and a point on edge E2 */
  /* ---------------------------------------------------- */

  t = -(b + e) / c; /* Minimum for d2_e1e2 is reached for t = -(b+e) / c
                       when s = 1. */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && logfile != NULL)
    fprintf(logfile, " [II] Test inter. P2E1 - E2: t = %10.8e\n", t);
#endif

  if (t > -fraction && t < 1.0 + fraction) { /* filter */

    d2_e1e2 = CS_ABS((a + 2.0*(b*t + d)) + t*(c*t + 2.0*e) + f);

    assert(d2_e1e2 >= -1.e-8);

    /* Tolerance for vertex p2e1 (s = 1) */

    d2_limit_e1 = p2e1.tolerance * p2e1.tolerance;

    if (d2_e1e2 <= d2_limit_e1) { /* Under the tolerance for vertex P2E1 */

      d_limit_e2 = (1.0-t)*p1e2.tolerance + t*p2e2.tolerance;
      d2_limit_e2 = d_limit_e2 * d_limit_e2;

      if (d2_e1e2 <= d2_limit_e2) { /* Under tolerance for edge 2 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[II]-3: Under tol. for P2E1 and a point of E2\n"
                  "\td2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
                  d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        /* We are under tolerance for edge E2. We try to find if it is
           a vertex-vertex intersection */

        if (t < 0.0) { /* Test intersection P2E1-P1E2 */

          d2_e1e2 = CS_ABS(a + 2.0*d + f);
          d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 1.0, e2_id, 0.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = 1.0, t = 0.0
                 under the tolerance is possible */

              int_p1e2 = true;
              abs_e1[*n_inter] = 1.0;
              abs_e2[*n_inter] = 0.0;

              *n_inter += 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[II]-4: Add inter. Vtx-Vtx (1.00, 0.00)\n");
#endif

            }
          }
        }
        else if (t > 1.0) { /* Test intersection P2E1-P2E2 */

          d2_e1e2 = CS_ABS(a + 2.0*(b + d + e) + c + f);
          d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

          if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

            if (_check_equiv(edges, mesh,
                             e1_id, 1.0, e2_id, 1.0,
                             verbosity, logfile)) {

              /* "extremity-extremity" intersection with s = t = 1.0
                 under the tolerance is possible */

              int_p2e2 = true;
              abs_e1[*n_inter] = 1.00;
              abs_e2[*n_inter] = 1.00;

              *n_inter += 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && logfile != NULL)
                fprintf(logfile,
                        "\t[II]-5: Add inter. Vtx-Vtx (1.00, 1.00)\n");
#endif

            }
          }
        }
        else {

          assert(0.0 <= t && t <= 1.0);

          /* It's an "extremity-interior" intersection */

          abs_e1[*n_inter] = 1.00;
          abs_e2[*n_inter] = t;

          if (abs_e2[*n_inter] <= 0.0)
            int_p1e2 = true;
          if (abs_e2[*n_inter] >= 1.0)
            int_p2e2 = true;

          *n_inter += 1;
          assert(*n_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[II]-6: Add inter. Vtx-Edge (1.00, %10.8e)\n", t);
#endif
        }

      } /* End if d2_e1e2 <= d2_limit_e2 */

    } /* End if d2_e1e2 <= d2_limit_e1 */

  } /* End if -fraction < t < 1 + fraction */

  /* If two intersections are already detected, we stop here */

  if (*n_inter == 2)
    return;

  /* If vertices from edge E2 are not implied in an intersection
     we do a test for each vertex */

  if (int_p1e2 == false) {

    /* Distance between P1E2 (t = 0.0) on edge E2 and a point on edge E1 */
    /* ----------------------------------------------------------------- */

    s = -d / a; /* t = 0.0 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              " [III] Test inter. P1E2 - E1: s = %10.8e\n", s);
#endif

    if (s > 0.0 && s < 1.0) { /* s = 0.0 and s = 1.0 are already tested */

      d2_e1e2 = CS_ABS(s*(a*s + 2.0*d) + f);
      d2_limit_e2 = p1e2.tolerance * p1e2.tolerance;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              "\t [III-1a] d2_e1e2: %10.8e - d2_lim_e2: %10.8e\n",
              d2_e1e2, d2_limit_e2);
#endif

      if (d2_e1e2 <= d2_limit_e2) { /* Under the tolerance for vertex P1E2 */

        d_limit_e1 = (1.0-s)*p1e1.tolerance + s*p2e1.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t [IV-1b] d2_e1e2: %10.8e - d2_lim_e1: %10.8e\n",
                  d2_e1e2, d2_limit_e1);
#endif

        if (d2_e1e2 <= d2_limit_e1) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf
              (logfile,
               "\t[III]-2: Under tol. for P1E2 and a point of E1\n"
               "\t d2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
               d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

          /* We are under tolerance for edge E1. There is an intersection
             between P1E2 and a point on edge E1 */

          abs_e1[*n_inter] = s;
          abs_e2[*n_inter] = 0.0;

          *n_inter += 1;
          assert(*n_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[III]-3: Add inter. Edge-Vtx (%10.8e, 0.00)\n", s);
#endif

          if (*n_inter == 2)
            return;

        } /* End if d2_e1e2 < d2_limit_e1 */

      } /* End if d2_e1e2 < d2_limit_e2 */

    } /* 0.0 < s < 1.0 */

  } /* If int_p1e2 == false */

  if (int_p2e2 == false) {

    /* Distance between P2E2 (t = 1.0) on edge E2 and a point on edge E1 */
    /* ----------------------------------------------------------------- */

    s = -(b + d) / a; /* t = 1 */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile, " [IV] Test inter. P2E2 - E1: s = %10.8e\n", s);
#endif

    if (s > 0.0 && s < 1.0) { /* s = 0.0 and s = 1.0 are already tested */

      d2_e1e2 = CS_ABS(s*(a*s + 2.0*(b + d)) + c + 2.0*e + f);
      d2_limit_e2 = p2e2.tolerance * p2e2.tolerance;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile, "\t [IV-1a] d2_e1e2: %10.8e - d2_lim_e2: %10.8e\n",
              d2_e1e2, d2_limit_e2);
#endif

      if (d2_e1e2 <= d2_limit_e2) { /* Under the tolerance for vertex P2E2 */

        d_limit_e1 = (1.0-s)*p1e1.tolerance + s*p2e1.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile, "\t [IV-1b] d2_e1e2: %10.8e - d2_lim_e1: %10.8e\n",
                  d2_e1e2, d2_limit_e1);
#endif

        if (d2_e1e2 <= d2_limit_e1) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf
              (logfile,
               "\t[IV]-2: Under tol. for P2E2 and a point of E1\n"
               "\t d2_e1e2 = %10.8e, d2_lim_e1: %10.8e, d2_lim_e2: %10.8e\n",
               d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

          /* We are under tolerance for edge E1. There is an intersection
             between P2E2 and a point on edge E1 */

          abs_e1[*n_inter] = s;
          abs_e2[*n_inter] = 1.0;

          *n_inter += 1;
          assert(*n_inter <= 2);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile, "\t[IV]-3: Add inter (%10.8e, 1.01)\n", s);
#endif

          if (*n_inter == 2)
            return;

        } /* End if d2_e1e2 < d2_limit_e1 */

      } /* End if d2_e1e2 < d2_limit_e2 */

    } /* 0.0 < s < 1.0 */

  } /* If int_p2e2 == false */

  /* If there is at least one intersection, we stop here. */

  if (*n_inter > 0)
    return;

  {
    /*------------------------------------------------------------------
     * Second part : no intersection has been found ("parallel edges",
     * "extremity-extremity" or "extremity-interior"). The only
     * remaining possibility is an intersection with 0 < s < 1 and
     * 0 < t < 1.
     *------------------------------------------------------------------*/

    /* Edges are parallel if a.c - b^2 = 0.
     * However, this term is o(length^4). So, we will prefer work with angle
     * If two edges are not parallel => a.c - b^2 > 0.
     *
     * cos(theta) = v0.v1 / ||v0||.||v1||
     * cos^2 (theta) = (v0.v1)^2 / ((v0.v0).(v1.v1))^2 = (-b)^2 / (a.c)
     * sin^2 (theta) = 1 - b^2 / (a.c)
     *
     * || v0 ^ v1 ||^2 = ||v0||^2 . ||v1||^2 . sin^2(thetha)
     *
     */

    double  inv_cross_norm2;
    double  cross_norm2 = CS_ABS(a * c - b * b);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && logfile != NULL)
      fprintf(logfile,
              " [V]: Test inter. E1 - E2: cross_norm2: %10.8e\n"
              "\ta = %10.8e, b = %10.8e, c = %10.8e, a*c = %10.8e\n",
              cross_norm2, a, b, c, a*c);
#endif

    if (cross_norm2 < parall_eps2 * a * c) /* <=> 1 - b^2/(a.c) < epsilon^2
                                              <=> sin^2(theta)  < epsilon^2 */
      return;

    /* Edges are not parallel.
       Define s and t if there is an intersection on interior points.

       ->     ->     ->
       v2 + s.v0 = t.v1 and cross_norm2 = a.c - b^2 != 0.0

       => s = (b.e - c.d)/cross_norm2
       => t = (e.a - b.d)/cross_norm2
    */

    s = b * e - c * d;
    t = b * d - a * e;

    if (s >= 0. && s <= cross_norm2) {
      if (t >= 0. && t <= cross_norm2)  {

        /* If tests are OK, we are on an interior point for E1 and E2 */

        inv_cross_norm2 = 1.0 / cross_norm2;
        s *= inv_cross_norm2;
        t *= inv_cross_norm2;
        d2_e1e2 = CS_ABS(s*(a*s + 2.0*(b*t + d)) + t*(c*t + 2.0*e) + f);
        d_limit_e1 = (1.0-s) * p1e1.tolerance + s * p2e1.tolerance;
        d_limit_e2 = (1.0-t) * p1e2.tolerance + t * p2e2.tolerance;
        d2_limit_e1 = d_limit_e1 * d_limit_e1;
        d2_limit_e2 = d_limit_e2 * d_limit_e2;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && logfile != NULL)
          fprintf(logfile,
                  "\t[V]-2: s = %10.8e, t = %10.8e\n"
                  "\t  d2_e1e2: %10.8e, d2_limit_e1: %10.8e, "
                  "d2_limit_e2: %10.8e\n",
                  s, t, d2_e1e2, d2_limit_e1, d2_limit_e2);
#endif

        if (d2_e1e2 <= d2_limit_e1 && d2_e1e2 <= d2_limit_e2) {

          /* Under tolerance for edges E1 and E2 => intersection is possible */

          assert(0.0 <= s && s <= 1.0);
          assert(0.0 <= t && t <= 1.0);

          abs_e1[*n_inter] = s;
          abs_e2[*n_inter] = t;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg && logfile != NULL)
            fprintf(logfile,
                    "\t[V]-3: Add inter. Edge-Edge (%10.8e,%10.8e)\n", s, t);
#endif

          *n_inter += 1;
          assert(*n_inter == 1);

        } /* If we are under tolerance */

      } /* If it is an interior point for edge E2 */
    } /* If it is an interior point for edge E1 */

  } /* End of last test : inter. E1 - E2 */

}

/*----------------------------------------------------------------------------
 * Check if we have to add the current intersection into the definition
 * of the reference edge description.
 *
 * parameters:
 *   inter     <-- intersection to test
 *   ref       <-- list of reference edges
 *   edge_id   <-- id of the edge to treat
 *   cur_shift <-- shift on edges (current used positions)
 *
 * returns:
 *   true if we should add this intersection to the cs_join_edge_inter_t
 *   structure, false otherwise.
 *---------------------------------------------------------------------------*/

inline static bool
_need_to_add_exch_inter(exch_inter_t                  inter,
                        const cs_join_inter_edges_t  *ref,
                        cs_lnum_t                     edge_id,
                        const cs_lnum_t              *cur_shift)
{
  cs_lnum_t  i;

  bool ret = true;

  for (i = ref->index[edge_id]; i < cur_shift[edge_id]; i++) {

    if (ref->vtx_glst[i] == inter.vtx_gnum) {
      if (fabs(ref->abs_lst[i] - inter.curv_abs) < 1e-30) {
        ret = false;
        break;
      }
    }

  } /* End of loop on sub-elements of the current edge */

  return ret;
}

/*----------------------------------------------------------------------------
 * Update statistics and timings for face bounding-box intersection search.
 *
 * parameters:
 *   face_neighborhood <-- pointer to face neighborhood management structure.
 *   box_wtime         <-- bounding box construction time
 *   box_dim           --> dimension of bounding box layout
 *   stats             <-> joining statistics
 *---------------------------------------------------------------------------*/

static void
_face_bbox_search_stats(const fvm_neighborhood_t  *face_neighborhood,
                        cs_timer_counter_t         box_time,
                        int                       *box_dim,
                        cs_join_stats_t           *stats)
{
  int i;
  int  depth[3];
  cs_lnum_t  _n_leaves[3], _n_boxes[3];
  cs_lnum_t  _n_threshold_leaves[3], _n_leaf_boxes[3];
  size_t  _mem_final[3], _mem_required[3];
  double  build_wtime, build_cpu_time, query_wtime, query_cpu_time;

  int dim = fvm_neighborhood_get_box_stats(face_neighborhood,
                                           depth,
                                           _n_leaves,
                                           _n_boxes,
                                           _n_threshold_leaves,
                                           _n_leaf_boxes,
                                           _mem_final,
                                           _mem_required);

  fvm_neighborhood_get_times(face_neighborhood,
                             &build_wtime,
                             &build_cpu_time,
                             &query_wtime,
                             &query_cpu_time);

  cs_timer_counter_t  build_time, query_time;

  build_time.wall_nsec = build_wtime*1e9;
  build_time.cpu_nsec = build_cpu_time*1e9;
  query_time.wall_nsec = query_wtime*1e9;
  query_time.cpu_nsec = query_cpu_time*1e9;

  for (i = 0; i < 3; i++) {
    _mem_final[i] /= 1024;
    _mem_required[i] /= 1024;
  }

  *box_dim = dim;

  stats->bbox_layout = CS_MAX(stats->bbox_layout, dim);

  if (stats->n_calls < 1) {
    stats->bbox_depth[1] = depth[1];
    stats->n_leaves[1] = _n_leaves[1];
    stats->n_boxes[1] = _n_boxes[1];
    stats->n_th_leaves[1] = _n_threshold_leaves[1];
    stats->n_leaf_boxes[1] = _n_leaf_boxes[1];
    stats->box_mem_final[1] = _mem_final[1];
    stats->box_mem_required[1] = _mem_required[1];
  }

  stats->bbox_depth[0] += depth[0];
  stats->bbox_depth[1] = CS_MIN(stats->bbox_depth[1],
                                (cs_gnum_t)depth[1]);
  stats->bbox_depth[2] = CS_MAX(stats->bbox_depth[2],
                                (cs_gnum_t)depth[2]);

  stats->n_leaves[0] += _n_leaves[0];
  stats->n_leaves[1] = CS_MIN(stats->n_leaves[1],
                               (cs_gnum_t)_n_leaves[1]);
  stats->n_leaves[2] = CS_MAX(stats->n_leaves[2],
                              (cs_gnum_t)_n_leaves[2]);

  stats->n_boxes[0] += _n_boxes[0];
  stats->n_boxes[1] = CS_MIN(stats->n_boxes[1],
                             (cs_gnum_t)_n_boxes[1]);
  stats->n_boxes[2] = CS_MAX(stats->n_boxes[2],
                             (cs_gnum_t)_n_boxes[2]);

  stats->n_th_leaves[0] += _n_threshold_leaves[0];
  stats->n_th_leaves[1] = CS_MIN(stats->n_th_leaves[1],
                                 (cs_gnum_t)_n_threshold_leaves[1]);
  stats->n_th_leaves[2] = CS_MAX(stats->n_th_leaves[2],
                                 (cs_gnum_t)_n_threshold_leaves[2]);

  stats->n_leaf_boxes[0] += _n_leaf_boxes[0];
  stats->n_leaf_boxes[1] = CS_MIN(stats->n_leaf_boxes[1],
                                  (cs_gnum_t)_n_leaf_boxes[1]);
  stats->n_leaf_boxes[2] = CS_MAX(stats->n_leaf_boxes[2],
                                  (cs_gnum_t)_n_leaf_boxes[2]);

  stats->box_mem_final[0] += _mem_final[0];
  stats->box_mem_final[1] = CS_MIN(stats->box_mem_final[1], _mem_final[1]);
  stats->box_mem_final[2] = CS_MAX(stats->box_mem_final[2], _mem_final[2]);

  stats->box_mem_required[0] += _mem_required[0];
  stats->box_mem_required[1] = CS_MIN(stats->box_mem_required[1],
                                       (cs_gnum_t)_mem_required[1]);
  stats->box_mem_required[2] = CS_MAX(stats->box_mem_required[2],
                                       (cs_gnum_t)_mem_required[2]);

  CS_TIMER_COUNTER_ADD(stats->t_box_build, stats->t_box_build, box_time);
  CS_TIMER_COUNTER_ADD(stats->t_box_build, stats->t_box_build, build_time);

  CS_TIMER_COUNTER_ADD(stats->t_box_query, stats->t_box_query, query_time);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
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
cs_join_inter_set_create(cs_lnum_t  init_size)
{
  cs_join_inter_set_t  *new_set = NULL;

  BFT_MALLOC(new_set, 1, cs_join_inter_set_t);

  new_set->n_max_inter = init_size; /* default value */
  new_set->n_inter = 0;

  BFT_MALLOC(new_set->inter_lst, 2*new_set->n_max_inter, cs_join_inter_t);

  return new_set;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_join_inter_set_t structure.
 *
 * parameter:
 *   inter_set <-> a pointer to the inter_set_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_inter_set_destroy(cs_join_inter_set_t  **inter_set)
{
  if (inter_set != NULL) {
    if (*inter_set != NULL) {
      BFT_FREE((*inter_set)->inter_lst);
      BFT_FREE(*inter_set);
    }
  }
}

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
                       const cs_join_mesh_t       *mesh)
{
  int  i;

  fprintf(f, "\n  Dump an inter_set_t structure (%p)\n", (const void *)i_set);

  fprintf(f, "  n_max_inter: %10d\n", i_set->n_max_inter);
  fprintf(f, "  n_inter    : %10d\n\n", i_set->n_inter);

  for (i = 0; i < i_set->n_inter; i++) {

    cs_join_inter_t  inter1 = i_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = i_set->inter_lst[2*i+1];
    cs_lnum_t  v1e1_id = edges->def[2*inter1.edge_id] - 1;
    cs_lnum_t  v2e1_id = edges->def[2*inter1.edge_id+1] - 1;
    cs_lnum_t  v1e2_id = edges->def[2*inter2.edge_id] - 1;
    cs_lnum_t  v2e2_id = edges->def[2*inter2.edge_id+1] - 1;
    cs_gnum_t  v1e1 = (mesh->vertices[v1e1_id]).gnum;
    cs_gnum_t  v2e1 = (mesh->vertices[v2e1_id]).gnum;
    cs_gnum_t  v1e2 = (mesh->vertices[v1e2_id]).gnum;
    cs_gnum_t  v2e2 = (mesh->vertices[v2e2_id]).gnum;

    fprintf(f, "\n%5d - (%9llu - %9llu)\n",
            i, (unsigned long long)edges->gnum[inter1.edge_id],
            (unsigned long long)edges->gnum[inter2.edge_id]);
    fprintf(f, "E1 [%5llu %5llu]  (%6.3f)\n",
            (unsigned long long)v1e1, (unsigned long long)v2e1,
            inter1.curv_abs);
    fprintf(f, "E2 [%5llu %5llu]  (%6.3f)\n",
            (unsigned long long)v1e2, (unsigned long long)v2e2,
            inter2.curv_abs);

  } /* End of loop on edge intersections */

  fflush(f);
}

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
cs_join_inter_edges_create(cs_lnum_t  n_edges)
{
  cs_lnum_t  i;

  cs_join_inter_edges_t  *inter_edges = NULL;


  /* Allocate and initialize structure */

  BFT_MALLOC(inter_edges, 1, cs_join_inter_edges_t);

  inter_edges->n_edges = n_edges;

  /* Build inter_edges_idx */

  BFT_MALLOC(inter_edges->index, n_edges + 1, cs_lnum_t);

  for (i = 0; i < n_edges + 1; i++)
    inter_edges->index[i] = 0;

  BFT_MALLOC(inter_edges->edge_gnum, n_edges, cs_gnum_t);

  for (i = 0; i < n_edges; i++)
    inter_edges->edge_gnum[i] = 0;

  inter_edges->max_sub_size = 0;

  inter_edges->vtx_lst = NULL;
  inter_edges->vtx_glst = NULL;
  inter_edges->abs_lst = NULL;

  return inter_edges;
}

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
                           const cs_join_inter_set_t  *inter_set)
{
  cs_lnum_t  i, lst_size, n_edge_inter;

  cs_lnum_t  max_n_sub_inter = 0;
  cs_lnum_t  *counter = NULL;
  cs_join_inter_edges_t  *inter_edges = NULL;

  /* Sanity checks */

  assert(edges != NULL);
  assert(inter_set != NULL);

  /* Allocate and initialize structure */

  inter_edges = cs_join_inter_edges_create(edges->n_edges);

  for (i = 0; i < edges->n_edges; i++)
    inter_edges->edge_gnum[i] = edges->gnum[i];

  n_edge_inter = 2*inter_set->n_inter;

  if (n_edge_inter == 0)
    return inter_edges;

  for (i = 0; i < n_edge_inter; i++) {

    cs_join_inter_t  inter = inter_set->inter_lst[i];
    cs_lnum_t  edge_id = inter.edge_id;

    assert(edge_id < edges->n_edges);

    if (inter.curv_abs > 0.0 && inter.curv_abs < 1.0)
      inter_edges->index[edge_id+1] += 1;

  }

  for (i = 0; i < edges->n_edges; i++) {

    cs_lnum_t  n_sub_inter = inter_edges->index[i+1];

    max_n_sub_inter = CS_MAX(max_n_sub_inter, n_sub_inter);
    inter_edges->index[i+1] += inter_edges->index[i];

  }

  inter_edges->max_sub_size = max_n_sub_inter;
  lst_size = inter_edges->index[edges->n_edges];

  /* Fill structures */

  BFT_MALLOC(inter_edges->vtx_lst, lst_size, cs_lnum_t);
  BFT_MALLOC(inter_edges->abs_lst, lst_size, cs_coord_t);

  BFT_MALLOC(counter, edges->n_edges, cs_lnum_t);

  for (i = 0; i < edges->n_edges; i++)
    counter[i] = 0;

  for (i = 0; i < n_edge_inter; i++) {

    cs_join_inter_t  inter = inter_set->inter_lst[i];
    cs_lnum_t  edge_id = inter.edge_id;

    if (inter.curv_abs > 0.0 && inter.curv_abs < 1.0) {

      cs_lnum_t  shift = inter_edges->index[edge_id] + counter[edge_id];

      inter_edges->vtx_lst[shift] = inter.vtx_id+1;
      inter_edges->abs_lst[shift] = inter.curv_abs;
      counter[edge_id] += 1;

    }

  } /* End of loop on intersections */

  /* Order lists */

  for (i = 0; i < edges->n_edges; i++) {

    cs_lnum_t  start = inter_edges->index[i];
    cs_lnum_t  end = inter_edges->index[i+1];

    if (end - start > 1)
      _adapted_lshellsort(start,
                          end,
                          inter_edges->abs_lst,
                          inter_edges->vtx_lst);

  } /* End of loop on edges */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check */

  for (i = 0; i < edges->n_edges; i++) {

    cs_lnum_t  j;
    cs_lnum_t  start = inter_edges->index[i];
    cs_lnum_t  end = inter_edges->index[i+1];

    if (end - start > 0) {

      assert(inter_edges->abs_lst[start] < 1.0);
      assert(inter_edges->abs_lst[start] > 0.0);

      for (j = start; j < end - 1; j++) {

        assert(inter_edges->abs_lst[j+1] < 1.0);
        assert(inter_edges->abs_lst[j+1] > 0.0);

        if (inter_edges->abs_lst[j] > inter_edges->abs_lst[j+1])
          bft_error(__FILE__, __LINE__, 0,
                    _("\n  Incoherency found in inter_edges structure for"
                      " edge %d (%u):\n"
                      "  Bad ordering of curvilinear abscissa.\n"
                      "   Vertex %d (%u) of abscissa: %f is before vertex %d"
                      " (%u) of abscissa: %f\n"),
                    i+1, edges->gnum[i], inter_edges->vtx_lst[j],
                    inter_edges->abs_lst[j], inter_edges->vtx_lst[j+1],
                    inter_edges->abs_lst[j+1]);

        if (inter_edges->vtx_lst[j] == inter_edges->vtx_lst[j+1])
          bft_error(__FILE__, __LINE__, 0,
                    _("\n  Incoherency found in inter_edges structure.\n"
                      "  Redundancy for edge %d (%u) :\n"
                      "   v1_num :  %d\n"
                      "   v2_num :  %d\n"
                      "  Vertex %d appears twice.\n"), i+1, edges->gnum[i],
                    edges->def[2*i], edges->def[2*i+1],
                    inter_edges->vtx_lst[j]);

      }

    } /* n_sub_inter > 1 */

  } /* End of loop on edges */
#endif

  /* Free memory */

  BFT_FREE(counter);

  /* Return pointer */

  return inter_edges;
}

/*----------------------------------------------------------------------------
 * Destroy an cs_join_inter_edges_t structure.
 *
 * parameters:
 *   inter_edges <-> pointer to cs_join_inter_edges_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_inter_edges_destroy(cs_join_inter_edges_t  **inter_edges)
{
  if (inter_edges != NULL) {
    cs_join_inter_edges_t  *ie = *inter_edges;
    if (ie != NULL) {
      BFT_FREE(ie->index);
      BFT_FREE(ie->edge_gnum);
      BFT_FREE(ie->vtx_lst);
      BFT_FREE(ie->vtx_glst);
      BFT_FREE(ie->abs_lst);
      BFT_FREE(*inter_edges);
    }
  }
}

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
                             cs_join_eset_t               *vtx_equiv)
{
  cs_lnum_t  i, j, k, i1, i2, size, esize, n_breaks;

  bool  *equiv_lst = NULL;
  cs_lnum_t   *vtx_lst = NULL, *tag_lst = NULL;
  cs_coord_t  *abs_lst = NULL;
  double  *tol_lst = NULL;
  FILE  *logfile = cs_glob_join_log;

  assert(mesh != NULL);
  assert(edges != NULL);
  assert(inter_edges != NULL);
  assert(vtx_equiv != NULL);

  int  n_break_counter = 0, n_max_breaks = 0;

  if (inter_edges != NULL) {
    if (inter_edges->index[inter_edges->n_edges] > 0) {

      assert(inter_edges->vtx_lst != NULL);
      assert(inter_edges->abs_lst != NULL);

      size = inter_edges->max_sub_size + 2;
      BFT_MALLOC(vtx_lst, size, cs_lnum_t);
      BFT_MALLOC(tag_lst, size, cs_lnum_t);
      BFT_MALLOC(abs_lst, size, cs_coord_t);
      BFT_MALLOC(tol_lst, size, double);
      esize = size*(size-1)/2;
      BFT_MALLOC(equiv_lst, esize, bool);

      /* Main loop */

      for (i = 0; i < inter_edges->n_edges; i++) {

        cs_lnum_t  v1_num = edges->def[2*i];
        cs_lnum_t  v2_num = edges->def[2*i+1];
        cs_lnum_t  v1_id = v1_num - 1;
        cs_lnum_t  v2_id = v2_num - 1;
        cs_lnum_t  start = inter_edges->index[i];
        cs_lnum_t  end = inter_edges->index[i+1];
        cs_lnum_t  n_sub_elts = 2 + end - start;
        double  edge_length = _compute_length(mesh->vertices[v1_id],
                                              mesh->vertices[v2_id]);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (param.verbosity > 4) {
          int  vid;
          double  ds_tol;
          cs_gnum_t  v1_gnum = (mesh->vertices[v1_num-1]).gnum;
          cs_gnum_t  v2_gnum = (mesh->vertices[v2_num-1]).gnum;

          fprintf(logfile,
                  "\n%6d: [%9llu] = (%7d [%9llu] - %7d [%9llu] - len: %8.6e)\n",
                  i, (unsigned long long)edges->gnum[i],
                  v1_num, (unsigned long long)v1_gnum,
                  v2_num, (unsigned long long)v2_gnum,
                  edge_length);
          fflush(logfile);

          if (inter_edges->vtx_glst == NULL) {

            for (j = start, k = 0; j < end; j++, k++) {
              vid = inter_edges->vtx_lst[j] - 1;
              if (vid > mesh->n_vertices) {
                cs_join_mesh_dump(logfile, mesh);
                fprintf(logfile, "vid: %d - n_vertices: %d\n",
                        vid, mesh->n_vertices);
                bft_error(__FILE__, __LINE__, 0,
                          _("  Vertex number out of bounds.\n"));
              }
              ds_tol = mesh->vertices[vid].tolerance/edge_length;
              fprintf(logfile,
                      "    %7d (%9llu) - (%7d, s = %8.6e, ds_tol = %8.6e)\n",
                      vid+1, (unsigned long long)mesh->vertices[vid].gnum,
                      k+1, inter_edges->abs_lst[j], ds_tol);
            }
          }
          else {

            for (j = start, k = 0; j < end; j++, k++) {
              vid = inter_edges->vtx_lst[j] - 1;
              if (vid > mesh->n_vertices) {
                cs_join_mesh_dump(logfile, mesh);
                fprintf(logfile, "vid: %d - n_vertices: %d\n",
                        vid, mesh->n_vertices);
                bft_error(__FILE__, __LINE__, 0,
                          _("  Vertex number out of bounds.\n"));
              }
            ds_tol = mesh->vertices[vid].tolerance/edge_length;
            fprintf(logfile,
                    "   %9llu - (%7d, s = %8.6e, ds_tol = %8.6e)\n",
                    (unsigned long long)inter_edges->vtx_glst[j],
                    k+1, inter_edges->abs_lst[j], ds_tol);
            }
          }
        } /* param.verbosity > 4 */
#endif

        /* Build temporary lists */

        vtx_lst[0] = v1_num;
        abs_lst[0] = 0.0;
        tol_lst[0] = (mesh->vertices[v1_id]).tolerance * _cs_join_tol_eps_coef;

        for (j = start, k = 1; j < end; j++, k++) {
          vtx_lst[k] = inter_edges->vtx_lst[j];
          abs_lst[k] = inter_edges->abs_lst[j];
          tol_lst[k] =  (mesh->vertices[vtx_lst[k]-1]).tolerance
                       * _cs_join_tol_eps_coef;
        }

        vtx_lst[k] = v2_num;
        abs_lst[k] = 1.0;
        tol_lst[k] =  (mesh->vertices[v2_id]).tolerance
                    * _cs_join_tol_eps_coef;

        /* Loop on couples of vertices to find if two vertices are equivalent
           Apply a tolerance reduction if necessary. */

        n_breaks = _find_edge_equiv(param,
                                    n_sub_elts,
                                    abs_lst,
                                    tol_lst,
                                    equiv_lst,
                                    tag_lst,
                                    edge_length);

        if (n_breaks > 0) {
          n_break_counter += 1;
          if (param.verbosity > 3)
            fprintf(logfile,
                    " Edge %8d: n_equiv. broken: %d\n", i+1, n_breaks);
        }

        n_max_breaks = CS_MAX(n_max_breaks, n_breaks);

        /* Add new equivalences */

        for (i1 = 0; i1 < n_sub_elts - 1; i1++) {
          for (i2 = i1 + 1; i2 < n_sub_elts; i2++) {

            if (tag_lst[i1] == tag_lst[i2]) {
              if (vtx_lst[i1] != vtx_lst[i2]) {

                cs_lnum_t  equiv_id = vtx_equiv->n_equiv;
                cs_join_eset_check_size(equiv_id, &vtx_equiv);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
                if (logfile != NULL)
                  fprintf
                    (logfile, "  Add equiv %d between [%llu, %llu]\n",
                     equiv_id+1,
                     (unsigned long long)mesh->vertices[vtx_lst[i1]-1].gnum,
                     (unsigned long long)mesh->vertices[vtx_lst[i2]-1].gnum);
#endif

                if (vtx_lst[i1] < vtx_lst[i2]) {
                  vtx_equiv->equiv_couple[2*equiv_id] = vtx_lst[i1];
                  vtx_equiv->equiv_couple[2*equiv_id+1] = vtx_lst[i2];
                }
                else {
                  vtx_equiv->equiv_couple[2*equiv_id] = vtx_lst[i2];
                  vtx_equiv->equiv_couple[2*equiv_id+1] = vtx_lst[i1];
                }
                vtx_equiv->n_equiv += 1;

              }
            } /* Equivalence found */

          } /* End of loop on i2 */
        } /* End of loop on i1 */

      } /* End of loop on edge intersections */

      /* Free memory */

      BFT_FREE(vtx_lst);
      BFT_FREE(tag_lst);
      BFT_FREE(abs_lst);
      BFT_FREE(tol_lst);
      BFT_FREE(equiv_lst);

    } /* inter_edges->index[inter_edges->n_edges] > 0 */
  } /* inter_edges != NULL */

  if (param.verbosity > 0) {

    cs_gnum_t n_g_break_counter = n_break_counter;
    cs_parall_counter(&n_g_break_counter, 1);

    bft_printf(_("\n  Equivalences broken for %llu edges.\n"),
               (unsigned long long)n_g_break_counter);

    if (param.verbosity > 1) {
      cs_lnum_t g_n_max_breaks = n_max_breaks;
      cs_parall_counter_max(&g_n_max_breaks, 1);
      bft_printf(_("\n  Max. number of equiv. breaks: %llu\n"),
                 (unsigned long long)g_n_max_breaks);
    }
  }

}

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
                                  const cs_join_inter_edges_t  *part)
{
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = cs_glob_rank_id;
  const cs_lnum_t  n_edges = edges->n_edges;

  /* Sanity check */

  assert(mesh != NULL);
  assert(part != NULL);
  assert(part->n_edges == n_edges);
  assert(edges != NULL);

  cs_block_dist_info_t  bi
    = cs_block_dist_compute_sizes(local_rank,
                                  n_ranks,
                                  1,
                                  0,
                                  edges->n_g_edges);

  cs_all_to_all_t
    *d = cs_all_to_all_create_from_block(n_edges,
                                         0, /* flags */
                                         part->edge_gnum,
                                         bi,
                                         mpi_comm);

  /* Send global numbers and index */

  cs_gnum_t *orig_gnum = cs_all_to_all_copy_array(d,
                                                  CS_GNUM_TYPE,
                                                  1,
                                                  false, /* reverse */
                                                  part->edge_gnum,
                                                  NULL);

  cs_lnum_t n_recv = cs_all_to_all_n_elts_dest(d);

  cs_lnum_t *orig_index
    = cs_all_to_all_copy_index(d,
                               false, /* reverse */
                               part->index,
                               NULL);

  /* Build send_inter_list and exchange information */

  cs_lnum_t send_inter_list_size = part->index[n_edges];
  exch_inter_t  *send_inter_list;
  BFT_MALLOC(send_inter_list, send_inter_list_size, exch_inter_t);

  for (cs_lnum_t i = 0; i < n_edges; i++) {

    cs_lnum_t  s_id = part->index[i];
    cs_lnum_t  e_id = part->index[i+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t  vtx_id = part->vtx_lst[j] - 1;
      (send_inter_list[j]).vtx_gnum = (mesh->vertices[vtx_id]).gnum;
      (send_inter_list[j]).curv_abs = part->abs_lst[j];
    }

  }

  /* Exchange buffers;
     index size if adapted temporarily for the datatype (a better
     solution would be to allow user datatypes in cs_def.c) */

  for (cs_lnum_t i = 0; i < n_recv+1; i++) {
    orig_index[i] *= sizeof(exch_inter_t);
  }

  cs_lnum_t *part_index;
  BFT_MALLOC(part_index, n_edges+1, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_edges+1; i++) {
    part_index[i] = part->index[i] * sizeof(exch_inter_t);
  }

  exch_inter_t *recv_inter_list
    = cs_all_to_all_copy_indexed(d,
                                 CS_CHAR,
                                 false, /* reverse */
                                 part_index,
                                 send_inter_list,
                                 orig_index,
                                 NULL);

  BFT_FREE(part_index);

  for (cs_lnum_t i = 0; i < n_recv+1; i++) {
    orig_index[i] /= sizeof(exch_inter_t);
  }

  BFT_FREE(send_inter_list);

  cs_all_to_all_destroy(&d);

  /* Synchronize the definition of each edge.
     Define a new cs_join_inter_edges_t struct. */

  cs_lnum_t block_size = 0;
  if (bi.gnum_range[1] > bi.gnum_range[0])
    block_size = bi.gnum_range[1] - bi.gnum_range[0];

  cs_join_inter_edges_t  *block
    = cs_join_inter_edges_create(block_size);

  for (cs_lnum_t i = 0; i < block_size; i++) {
    cs_gnum_t g_id = i;
    block->edge_gnum[i] = bi.gnum_range[0] + g_id;
  }

  /* First scan fill index_ref */

  for (cs_lnum_t i = 0; i < n_recv; i++) {

    cs_gnum_t  num = orig_gnum[i] - bi.gnum_range[0] + 1;
    cs_lnum_t  n_sub_elts = orig_index[i+1] - orig_index[i];

    assert(num <= (cs_gnum_t)block_size);
    block->index[num] += n_sub_elts;

  }

  cs_lnum_t *shift_ref;
  BFT_MALLOC(shift_ref, block_size, cs_lnum_t);

  for (cs_lnum_t i = 0; i < block_size; i++) {
    block->index[i+1] += block->index[i];
    shift_ref[i] = block->index[i];
  }

  BFT_MALLOC(block->vtx_glst,
             block->index[block_size],
             cs_gnum_t);

  BFT_MALLOC(block->abs_lst,
             block->index[block_size],
             cs_coord_t);

  /* Second scan: fill buffers */

  for (cs_lnum_t i = 0; i < n_recv; i++) {

    cs_gnum_t block_id = orig_gnum[i] - bi.gnum_range[0];

    assert(block_id < (cs_gnum_t)block_size);

    cs_lnum_t s_id = orig_index[i];
    cs_lnum_t e_id = orig_index[i+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {

      exch_inter_t  exch_inter = recv_inter_list[j];

      if (_need_to_add_exch_inter(exch_inter,
                                  block,
                                  block_id,
                                  shift_ref)) {

        cs_lnum_t  _shift = shift_ref[block_id];

        assert(_shift < block->index[block_id+1]);

        block->vtx_glst[_shift] = exch_inter.vtx_gnum;
        block->abs_lst[_shift] = exch_inter.curv_abs;
        shift_ref[block_id] += 1;

      } /* End of adding a new intersection in the edge definition */

    } /* End of loop on sub_elts */

  } /* End of loop on edge descriptions */

  BFT_FREE(recv_inter_list);

  /* Compact block */

  BFT_FREE(orig_gnum);
  BFT_FREE(orig_index);

  cs_lnum_t shift = 0;
  for (cs_lnum_t i = 0; i < block_size; i++) {

    for (cs_lnum_t j = block->index[i]; j < shift_ref[i]; j++) {

      block->vtx_glst[shift] = block->vtx_glst[j];
      block->abs_lst[shift] = block->abs_lst[j];
      shift++;

    }

  }

  cs_lnum_t *new_index;
  BFT_MALLOC(new_index, block_size + 1, cs_lnum_t);

  new_index[0] = 0;
  for (cs_lnum_t i = 0; i < block_size; i++)
    new_index[i+1] = new_index[i] + shift_ref[i] - block->index[i];

  BFT_FREE(shift_ref);
  BFT_FREE(block->index);

  block->index = new_index;

  BFT_REALLOC(block->vtx_glst, block->index[block_size], cs_gnum_t);
  BFT_REALLOC(block->abs_lst, block->index[block_size], cs_coord_t);

  /* Sort intersection by increasing curvilinear abscissa for each edge */

  cs_lnum_t  _max = 0;

  for (cs_lnum_t i = 0; i < block_size; i++)
    _max = CS_MAX(_max, new_index[i+1] - new_index[i]);

  block->max_sub_size = _max;

  for (cs_lnum_t i = 0; i < block_size; i++) {

    cs_lnum_t  start = block->index[i];
    cs_lnum_t  end = block->index[i+1];

    if (end - start > 0)
      _adapted_gshellsort(start, end, block->abs_lst, block->vtx_glst);

  } /* End of loop on edges of block */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_join_inter_edges_dump(cs_glob_join_log, block, edges, mesh);
#endif

  return block;
}

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
                                  cs_join_inter_edges_t        *part)
{
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = cs_glob_rank_id;
  const cs_block_dist_info_t  bi = cs_block_dist_compute_sizes(local_rank,
                                                               n_ranks,
                                                               1,
                                                               0,
                                                               n_g_edges);

  /* Sanity check */

  assert(block != NULL);
  assert(part != NULL);

#if defined(DEBUG) && !defined(NDEBUG)
  {
    cs_lnum_t  block_size = 0;
    if (bi.gnum_range[1] > bi.gnum_range[0])
      block_size = bi.gnum_range[1] - bi.gnum_range[0];
    assert(block->n_edges == block_size);
  }
#endif

  cs_all_to_all_t
    *d = cs_all_to_all_create_from_block(part->n_edges,
                                         CS_ALL_TO_ALL_USE_DEST_ID, /* flags */
                                         part->edge_gnum,
                                         bi,
                                         mpi_comm);

  /* Exchange global numbers of edges which should be returned
     from block to partition */

  cs_gnum_t *orig_gnum = cs_all_to_all_copy_array(d,
                                                  CS_GNUM_TYPE,
                                                  1,
                                                  false, /* reverse */
                                                  part->edge_gnum,
                                                  NULL);

  cs_lnum_t send_list_size = cs_all_to_all_n_elts_dest(d);

  assert(send_list_size == block->n_edges);

  /* Send the number of sub_elements for each requested edge */

  cs_lnum_t *block_index;
  BFT_MALLOC(block_index, send_list_size+1, cs_lnum_t);

  block_index[0] = 0;
  for (cs_lnum_t i = 0; i < send_list_size; i++) {
    cs_gnum_t  s_id = orig_gnum[i] - bi.gnum_range[0];
    block_index[i+1] =   block_index[i]
                       + block->index[s_id+1] - block->index[s_id];
  }

  cs_all_to_all_copy_index(d,
                           true, /* reverse */
                           block_index,
                           part->index);

  /* block struct. send the requested global edges to all the part struct
     on the distant ranks */

  cs_lnum_t send_inter_list_size = block_index[send_list_size];

  exch_inter_t *send_inter_list = NULL;
  BFT_MALLOC(send_inter_list, send_inter_list_size, exch_inter_t);

  for (cs_lnum_t i = 0; i < send_list_size; i++) {
    cs_gnum_t  b_id = orig_gnum[i] - bi.gnum_range[0];

    cs_lnum_t  s_id = block->index[b_id];
    cs_lnum_t  e_id = block->index[b_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      (send_inter_list[j]).vtx_gnum = block->vtx_glst[j];
      (send_inter_list[j]).curv_abs = block->abs_lst[j];
    }

  } /* End of loop on elements to send */

  BFT_FREE(orig_gnum);

  /* Exchange buffers */

  for (cs_lnum_t i = 0; i < send_list_size+1; i++) {
    block_index[i] *= sizeof(exch_inter_t);
  }
  for (cs_lnum_t i = 0; i < part->n_edges+1; i++) {
    part->index[i] *= sizeof(exch_inter_t);
  }

  exch_inter_t *recv_inter_list
    = cs_all_to_all_copy_indexed(d,
                                 CS_CHAR,
                                 true, /* reverse */
                                 block_index,
                                 send_inter_list,
                                 part->index,
                                 NULL);

  for (cs_lnum_t i = 0; i < part->n_edges+1; i++) {
    part->index[i] /= sizeof(exch_inter_t);
  }

  /* Partial free memory */

  BFT_FREE(send_inter_list);
  BFT_FREE(block_index);

  cs_all_to_all_destroy(&d);

  /* Update part definitions */

  cs_lnum_t  max_sub_size = 0;

  for (cs_lnum_t i = 0; i < part->n_edges; i++) {
    cs_lnum_t nsub = part->index[i+1] - part->index[i];
    max_sub_size = CS_MAX(max_sub_size, nsub);
  }

  part->max_sub_size = max_sub_size;

  cs_lnum_t part_indexed_size = part->index[part->n_edges];

  BFT_FREE(part->vtx_lst);
  part->vtx_lst = NULL;
  BFT_REALLOC(part->vtx_glst, part_indexed_size, cs_gnum_t);
  BFT_REALLOC(part->abs_lst, part_indexed_size, cs_coord_t);

  for (cs_lnum_t i = 0; i < part_indexed_size; i++) {
    part->vtx_glst[i] = (recv_inter_list[i]).vtx_gnum;
    part->abs_lst[i] = (recv_inter_list[i]).curv_abs;
  }

  /* Free memory */

  BFT_FREE(recv_inter_list);
}

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
                                cs_join_inter_edges_t  **inter_edges)
{
  cs_lnum_t  i, j, shift, o_id, max_size;

  cs_lnum_t  n_new_vertices = 0;
  cs_gnum_t  *vtx_gnum = NULL, *edge_gnum = NULL;
  cs_lnum_t  *edge_order = NULL, *vtx_order = NULL;
  cs_join_inter_edges_t  *_inter_edges = *inter_edges;
  cs_join_inter_edges_t  *new_inter_edges = NULL;
  cs_join_vertex_t  *new_vertices = NULL;

  assert(edges != NULL);
  assert(mesh != NULL);

  const cs_lnum_t  n_edges = edges->n_edges;
  const cs_lnum_t  n_init_vertices = mesh->n_vertices;

  assert(inter_edges != NULL);
  assert(n_edges == _inter_edges->n_edges);

  /* Check if we have to re-order inter_edges */

  for (i = 0; i < n_edges; i++)
    if (_inter_edges->edge_gnum[i] != edges->gnum[i])
      break;

  if (i != n_edges) {

    new_inter_edges = cs_join_inter_edges_create(n_edges);

    /* Order global edge numbering */

    BFT_MALLOC(edge_order, n_edges, cs_lnum_t);
    BFT_MALLOC(edge_gnum, n_edges, cs_gnum_t);

    cs_order_gnum_allocated(NULL, edges->gnum, edge_order, n_edges);

    for (i = 0; i < n_edges; i++)
      edge_gnum[i] = edges->gnum[edge_order[i]];

    /* edge_gnum -> edge_id */

    for (i = 0; i < n_edges; i++) {

      cs_gnum_t  e_gnum = _inter_edges->edge_gnum[i];
      cs_lnum_t  e_id = cs_search_g_binary(n_edges, e_gnum, edge_gnum);

      if (e_id == -1)
        bft_error(__FILE__, __LINE__, 0,
                  _("  The received edge global number (%llu) is unknown"
                    " on the current rank.\n"),
                  (unsigned long long)_inter_edges->edge_gnum[i]);

      o_id = edge_order[e_id];
      new_inter_edges->edge_gnum[o_id] = e_gnum;

      new_inter_edges->index[o_id+1] =
        _inter_edges->index[i+1] - _inter_edges->index[i];

    }

    for (i = 0; i < n_edges; i++)
      new_inter_edges->index[i+1] += new_inter_edges->index[i];

    BFT_MALLOC(new_inter_edges->vtx_glst,
               new_inter_edges->index[n_edges], cs_gnum_t);
    BFT_MALLOC(new_inter_edges->abs_lst,
               new_inter_edges->index[n_edges], cs_coord_t);

    for (i = 0; i < n_edges; i++) {

      cs_gnum_t  e_gnum = _inter_edges->edge_gnum[i];
      cs_lnum_t  e_id = cs_search_g_binary(n_edges, e_gnum, edge_gnum);

      o_id = edge_order[e_id];
      shift = new_inter_edges->index[o_id];

      for (j = _inter_edges->index[i]; j < _inter_edges->index[i+1]; j++) {
        new_inter_edges->vtx_glst[shift] = _inter_edges->vtx_glst[j];
        new_inter_edges->abs_lst[shift] = _inter_edges->abs_lst[j];
        shift++;
      }

    }

    /* Partial memory free */

    BFT_FREE(edge_gnum);
    BFT_FREE(edge_order);

    cs_join_inter_edges_destroy(&_inter_edges);

  } /* End if re-order necessary */

  else
    new_inter_edges = _inter_edges;

  if (new_inter_edges->vtx_lst == NULL)
    BFT_MALLOC(new_inter_edges->vtx_lst,
               new_inter_edges->index[n_edges], cs_lnum_t);

  /* Order global vertex numbering */

  BFT_MALLOC(vtx_gnum, n_init_vertices, cs_gnum_t);
  BFT_MALLOC(vtx_order, n_init_vertices, cs_lnum_t);

  for (i = 0; i < n_init_vertices; i++)
    vtx_gnum[i] = mesh->vertices[i].gnum;

  cs_order_gnum_allocated(NULL, vtx_gnum, vtx_order, n_init_vertices);

  for (i = 0; i < n_init_vertices; i++)
    vtx_gnum[i] = mesh->vertices[vtx_order[i]].gnum;

  /* Pre-allocate a buffer to store data on possible new vertices */

  max_size = 100;
  BFT_MALLOC(new_vertices, max_size, cs_join_vertex_t);

  /* Fill vtx_lst array of the cs_join_inter_edges_t structure */

  for (i = 0; i < n_edges; i++) {

    cs_lnum_t  start = new_inter_edges->index[i];
    cs_lnum_t  end = new_inter_edges->index[i+1];

    for (j = start; j < end; j++) {

      cs_lnum_t  id = cs_search_g_binary(n_init_vertices,
                                        new_inter_edges->vtx_glst[j],
                                        vtx_gnum);

      if (id == -1) { /* New vertex to add in the mesh structure */

        if (n_new_vertices >= max_size) {
          max_size *= 2;
          BFT_REALLOC(new_vertices, max_size, cs_join_vertex_t);
        }

        new_vertices[n_new_vertices]
          = _get_new_vertex(new_inter_edges->abs_lst[j],
                            new_inter_edges->vtx_glst[j],
                            &(edges->def[2*i]),
                            mesh);

        /* update new_inter_edges */

        n_new_vertices++;
        new_inter_edges->vtx_lst[j] = n_init_vertices + n_new_vertices;

      }
      else
        new_inter_edges->vtx_lst[j] = vtx_order[id] + 1;

    } /* End of loop on intersection description */

  } /* End of loop on edges */

  if (n_new_vertices > 0) {

    if (verbosity > 2)
      fprintf(cs_glob_join_log,
              "\n  Add %d new vertices in the %s mesh definition"
              " during update of the edge definition.\n",
              n_new_vertices, mesh->name);

    BFT_REALLOC(mesh->vertices,
                n_init_vertices + n_new_vertices,
                cs_join_vertex_t);

    for (i = 0; i < n_new_vertices; i++)
      mesh->vertices[n_init_vertices + i] = new_vertices[i];

    mesh->n_vertices = n_init_vertices + n_new_vertices;

  }

  /* Free memory */

  BFT_FREE(vtx_gnum);
  BFT_FREE(vtx_order);
  BFT_FREE(new_vertices);

  /* Returns pointer */

  *inter_edges = new_inter_edges;
}

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
                        cs_join_inter_set_t   **inter_set)
{
  cs_lnum_t  i, j, k;
  double  abs_e1[2], abs_e2[2];

  cs_join_type_t  join_type = CS_JOIN_TYPE_CONFORMING;
  cs_lnum_t  n_inter = 0;
  cs_lnum_t  n_inter_detected = 0, n_real_inter = 0, n_trivial_inter = 0;
  cs_gnum_t  n_g_inter[3] = {0, 0, 0};
  cs_join_inter_set_t  *_inter_set = NULL;
  cs_join_eset_t  *_vtx_eset = NULL;
  FILE  *logfile = cs_glob_join_log;

  const double  merge_limit = param.fraction * param.pre_merge_factor;
  const double  parall_eps2 = 1e-6;

  /* Sanity checks */

  assert(mesh != NULL);
  assert(edges != NULL);
  assert(edge_edge_vis != NULL);

  assert(vtx_eset != NULL);
  assert(inter_set != NULL);

  if (param.verbosity > 3)
    fprintf(logfile, "  Parallel intersection criterion: %8.5e\n",
            parall_eps2);

  /* Initialization of structures */

  _n_inter_tolerance_warnings = 0;

  _inter_set = cs_join_inter_set_create(50);
  _vtx_eset = cs_join_eset_create(30);

  /* Loop on edges */

  for (i = 0; i < edge_edge_vis->n_elts; i++) {

    int  e1 = edge_edge_vis->g_elts[i]; /* This is a local number */

    for (j = edge_edge_vis->index[i]; j < edge_edge_vis->index[i+1]; j++) {

      int  e2 = edge_edge_vis->g_list[j]; /* This is a local number */
      int  e1_id = (e1 < e2 ? e1 - 1 : e2 - 1);
      int  e2_id = (e1 < e2 ? e2 - 1 : e1 - 1);

      assert(e1 != e2);

      /* Get edge-edge intersection */

      if (param.icm == 1)
        _edge_edge_3d_inter(mesh,
                            edges,
                            param.fraction,
                            e1_id, abs_e1,
                            e2_id, abs_e2,
                            parall_eps2,
                            param.verbosity,
                            logfile,
                            &n_inter);

      else if (param.icm == 2)
        _new_edge_edge_3d_inter(mesh,
                                edges,
                                param.fraction,
                                e1_id, abs_e1,
                                e2_id, abs_e2,
                                parall_eps2,
                                param.verbosity,
                                logfile,
                                &n_inter);

      n_inter_detected += n_inter;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (param.verbosity > 3 && n_inter > 0) {

        cs_lnum_t  v1e1 = edges->def[2*e1_id] - 1;
        cs_lnum_t  v2e1 = edges->def[2*e1_id+1] - 1;
        cs_lnum_t  v1e2 = edges->def[2*e2_id] - 1;
        cs_lnum_t  v2e2 = edges->def[2*e2_id+1] - 1;

        fprintf(logfile,
                "\n Intersection: "
                "E1 (%llu) [%llu - %llu] / E2 (%llu) [%llu - %llu]\n",
                (unsigned long long)edges->gnum[e1_id],
                (unsigned long long)mesh->vertices[v1e1].gnum,
                (unsigned long long)mesh->vertices[v2e1].gnum,
                (unsigned long long)edges->gnum[e2_id],
                (unsigned long long)mesh->vertices[v1e2].gnum,
                (unsigned long long)mesh->vertices[v2e2].gnum);
        fprintf(logfile, "  n_inter: %d ", n_inter);
        for (k = 0; k < n_inter; k++)
          fprintf(logfile,
                  " (%d) - s_e1 = %g, s_e2 = %g", k, abs_e1[k], abs_e2[k]);
        fflush(logfile);
      }
#endif

      for (k = 0; k < n_inter; k++) {

        bool  trivial = false;

        if (abs_e1[k] <= merge_limit || abs_e1[k] >= 1.0 - merge_limit)
          if (abs_e2[k] <= merge_limit || abs_e2[k] >= 1.0 - merge_limit)
            trivial = true;

        if (trivial) {

          _add_trivial_equiv(e1_id,
                             e2_id,
                             abs_e1[k],
                             abs_e2[k],
                             edges,
                             _vtx_eset);

          n_trivial_inter += 1;

        }
        else {

          if (join_type == CS_JOIN_TYPE_CONFORMING)
            join_type = CS_JOIN_TYPE_NON_CONFORMING;

          _add_inter(e1_id, e2_id, abs_e1[k], abs_e2[k], _inter_set);

        }

      } /* End of loop on detected edge_edge_vis */

    } /* End of loop on entities intersecting elements */

  } /* End of loop on elements in intersection list */

  n_real_inter = n_inter_detected - n_trivial_inter;

  if (n_inter_detected == 0)
    join_type = CS_JOIN_TYPE_NULL;

  /* Order and delete redundant equivalences */

  cs_join_eset_clean(&_vtx_eset);

  /* Synchronize join_type and counts */

  n_g_inter[0] = n_inter_detected;
  n_g_inter[1] = n_trivial_inter;
  n_g_inter[2] = n_real_inter;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    int  tag = (int)join_type;
    int  sync_tag = tag;
    cs_gnum_t tmp_inter[3];

    MPI_Allreduce(&tag, &sync_tag, 1, MPI_INT, MPI_MAX, cs_glob_mpi_comm);

    join_type = sync_tag;

    MPI_Allreduce(n_g_inter, tmp_inter, 3, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0; i < 3; i++)
      n_g_inter[i] = tmp_inter[i];
  }
#endif

  if (param.verbosity > 0) {

    bft_printf(_("\n"
                 "  Global number of intersections detected: %12llu\n"
                 "    Vertex-Vertex intersections:    %12llu\n"
                 "    Other intersections:            %12llu\n"),
               (unsigned long long)n_g_inter[0],
               (unsigned long long)n_g_inter[1],
               (unsigned long long)n_g_inter[2]);

    if (param.verbosity > 2) {
      fprintf(logfile,
              "\n"
              "    Local number of intersections detected: %10d\n"
              "      Vertex-Vertex intersections:          %10d\n"
              "      Other intersections:                  %10d\n",
              (int)n_inter_detected, (int)n_trivial_inter,
              (int)n_real_inter);
      fprintf(logfile,
              "\n  Local number of edge-edge intersection warnings: %9d\n",
              _n_inter_tolerance_warnings);
      fprintf(logfile,
              "\n  Local number of equivalences between vertices: %9d\n",
              _vtx_eset->n_equiv);
    }

  }

  bft_printf_flush();

  /* Return pointer */

  *inter_set = _inter_set;
  *vtx_eset = _vtx_eset;

  return join_type;
}

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
                        cs_join_stats_t        *stats)
{
  cs_lnum_t  i;

  int  box_dim = 0;
  cs_coord_t  *f_extents = NULL;
  fvm_neighborhood_t  *face_neighborhood = NULL;
  cs_join_gset_t  *face_visibility = NULL;

  assert(join_mesh != NULL);

  cs_timer_t  t0 = cs_timer_time();

#if defined HAVE_MPI
  face_neighborhood = fvm_neighborhood_create(cs_glob_mpi_comm);
#else
  face_neighborhood = fvm_neighborhood_create();
#endif

  fvm_neighborhood_set_options(face_neighborhood,
                               param.tree_max_level,
                               param.tree_n_max_boxes,
                               param.tree_max_box_ratio,
                               param.tree_max_box_ratio_distrib);

  /* Allocate temporary extent arrays */

  BFT_MALLOC(f_extents, join_mesh->n_faces*6, cs_coord_t);

  /* Define each bounding box for the selected faces */

  for (i = 0; i < join_mesh->n_faces; i++)
    _get_face_extents(join_mesh->face_vtx_idx[i],
                      join_mesh->face_vtx_idx[i+1],
                      join_mesh->face_vtx_lst,
                      join_mesh->vertices,
                      f_extents + i*6);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  extents_time = cs_timer_diff(&t0, &t1);

  fvm_neighborhood_by_boxes(face_neighborhood,
                            3, /* spatial dimension */
                            join_mesh->n_faces,
                            join_mesh->face_gnum,
                            NULL,
                            NULL,
                            &f_extents);

  _face_bbox_search_stats(face_neighborhood,
                          extents_time,
                          &box_dim,
                          stats);

  if (param.verbosity > 0) {
    bft_printf(_("  Determination of possible face intersections:\n\n"
                 "    bounding-box tree layout: %dD\n"), box_dim);
    bft_printf_flush();
  }

  /* Retrieve face -> face visibility */

  BFT_MALLOC(face_visibility, 1, cs_join_gset_t);

  assert(sizeof(cs_lnum_t) == sizeof(cs_lnum_t));

  fvm_neighborhood_transfer_data(face_neighborhood,
                                 &(face_visibility->n_elts),
                                 &(face_visibility->g_elts),
                                 &(face_visibility->index),
                                 &(face_visibility->g_list));

  fvm_neighborhood_destroy(&face_neighborhood);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_FaceVis.dat")+1+2+4;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Join%02dDBG_FaceVis%04d.dat",
            param.num, CS_MAX(cs_glob_rank_id, 0));
    dbg_file = fopen(filename, "w");

    cs_join_gset_dump(dbg_file, face_visibility);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);
  }
#endif /* defined(DEBUG) && !defined(NDEBUG) */

  return face_visibility;
}

/*----------------------------------------------------------------------------
 * Transform face visibility into edge visibility (mesh->face_gnum must be
 * ordered).
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
                               const cs_join_gset_t   *face_visib)
{
  cs_lnum_t  i, j, k, edge_num, edge_id, shift;

  assert(mesh != NULL);
  assert(edges != NULL);
  assert(face_visib != NULL);

  cs_lnum_t  size = 0, size_max = 0;
  cs_lnum_t  *count = NULL, *face2edge_idx = NULL, *face2edge_lst = NULL;
  cs_gnum_t  *tmp = NULL;
  cs_join_gset_t  *edge_visib = NULL;

  /* Create a local "face -> edge" connectivity
     First, count number of edges for each face */

  BFT_MALLOC(face2edge_idx, mesh->n_faces + 1, cs_lnum_t);

  face2edge_idx[0] = 0;
  for (i = 0; i < mesh->n_faces; i++)
    face2edge_idx[i+1] = mesh->face_vtx_idx[i+1] - mesh->face_vtx_idx[i];

  for (i = 0; i < mesh->n_faces; i++)
    face2edge_idx[i+1] += face2edge_idx[i];

  /* Build face2edge_lst */

  BFT_MALLOC(face2edge_lst, face2edge_idx[mesh->n_faces], cs_lnum_t);
  BFT_MALLOC(count, mesh->n_faces, cs_lnum_t);

  for (i = 0; i < mesh->n_faces; i++)
    count[i] = 0;

  for (i = 0; i < mesh->n_faces; i++) {

    cs_lnum_t  start = mesh->face_vtx_idx[i];
    cs_lnum_t  end = mesh->face_vtx_idx[i+1];

    for (j = start; j < end - 1; j++) {

      edge_num = cs_join_mesh_get_edge(mesh->face_vtx_lst[j] + 1,
                                       mesh->face_vtx_lst[j+1] + 1,
                                       edges);

      shift = face2edge_idx[i] + count[i];
      count[i] += 1;
      face2edge_lst[shift] = CS_ABS(edge_num);

    }

    edge_num = cs_join_mesh_get_edge(mesh->face_vtx_lst[end-1] + 1,
                                     mesh->face_vtx_lst[start] + 1,
                                     edges);

    shift = face2edge_idx[i] + count[i];
    count[i] += 1;
    face2edge_lst[shift] = CS_ABS(edge_num);

  } /* End of loop on faces */

  /* Transform numbering in face_visib to match numbering
     in the current cs_join_mesh_t structure */

  for (i = 0; i < face_visib->n_elts; i++) {

    cs_lnum_t  face_id = cs_search_g_binary(mesh->n_faces,
                                           face_visib->g_elts[i],
                                           mesh->face_gnum);

    face_visib->g_elts[i] = face_id;

    for (j = face_visib->index[i]; j < face_visib->index[i+1]; j++) {

      cs_lnum_t  adj_id = cs_search_g_binary(mesh->n_faces,
                                            face_visib->g_list[j],
                                            mesh->face_gnum);

      face_visib->g_list[j] = adj_id;

    }

  } /* End of loop on bounding boxes */

  /* Create edge_visib. Unfold face_visib */

  for (i = 0; i < face_visib->n_elts; i++) {
    j = face_visib->g_elts[i];
    size += face2edge_idx[j+1] - face2edge_idx[j];
  }

  edge_visib = cs_join_gset_create(size);

  edge_id = 0;

  for (i = 0; i < face_visib->n_elts; i++) {

    cs_lnum_t  face_id = face_visib->g_elts[i];
    cs_lnum_t  start = face2edge_idx[face_id];
    cs_lnum_t  end = face2edge_idx[face_id+1];

    size = 0;
    for (j = face_visib->index[i]; j < face_visib->index[i+1]; j++) {

      cs_lnum_t  adj_id = face_visib->g_list[j];

      size += face2edge_idx[adj_id+1] - face2edge_idx[adj_id];

    }

    size_max = CS_MAX(size, size_max);

    for (j = start; j < end; j++) {

      edge_visib->g_elts[edge_id] = face2edge_lst[j];
      edge_visib->index[edge_id+1] = size;
      edge_id++;

    }

  } /* End of loop on bounding boxes */

  assert(edge_id == edge_visib->n_elts);

  /* Build index */

  for (i = 0; i < edge_visib->n_elts; i++)
    edge_visib->index[i+1] += edge_visib->index[i];

  BFT_MALLOC(edge_visib->g_list,
             edge_visib->index[edge_visib->n_elts],
             cs_gnum_t);

  BFT_MALLOC(tmp, size_max, cs_gnum_t);

  /* Build list */

  edge_id = 0;

  for (i = 0; i < face_visib->n_elts; i++) {

    cs_lnum_t  _count = 0;
    cs_lnum_t  face_id = face_visib->g_elts[i];
    cs_lnum_t  n_edges = face2edge_idx[face_id+1] - face2edge_idx[face_id];
    cs_lnum_t  b_start = face_visib->index[i];
    cs_lnum_t  b_end =  face_visib->index[i+1];

    /* Unfold face->edge connectivity for the current list of bounding boxes */

    for (j = b_start; j < b_end; j++) {

      cs_lnum_t  adj_face_id = face_visib->g_list[j];
      cs_lnum_t  n_adj_edges =  face2edge_idx[adj_face_id+1]
                              - face2edge_idx[adj_face_id];

      shift = face2edge_idx[adj_face_id];

      for (k = 0; k < n_adj_edges; k++)
        tmp[_count + k] = face2edge_lst[shift + k];
      _count += n_adj_edges;

    }

    for (j = 0; j < n_edges; j++, edge_id++) {

      assert(_count == edge_visib->index[edge_id+1]-edge_visib->index[edge_id]);

      shift = edge_visib->index[edge_id];

      for (k = 0; k < _count; k++)
        edge_visib->g_list[shift + k] = tmp[k];

    } /* End of loop on edges */

  } /* End of loop on bounding boxes */

  /* Free memory */

  BFT_FREE(face2edge_idx);
  BFT_FREE(face2edge_lst);
  BFT_FREE(count);
  BFT_FREE(tmp);

  /* Delete redundancies in g_elts, order g_elts and compact data */

  cs_join_gset_merge_elts(edge_visib, 0); /* 0 = g_elts is not ordered */

  /* Delete redundancies in g_list */

  cs_join_gset_clean(edge_visib);

  cs_join_gset_compress(edge_visib);

  /* Return pointers */

  return edge_visib;
}

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
                         const cs_join_mesh_t         *mesh)
{
  int  i, j, k;

  fprintf(f, "\n  Dump of a cs_join_inter_edges_t structure (%p)\n",
          (const void *)inter_edges);

  if (inter_edges == NULL)
    return;

  fprintf(f, "  n_edges:      %10d\n", inter_edges->n_edges);
  fprintf(f, "  max_sub_size: %10d\n\n", inter_edges->max_sub_size);

  for (i = 0; i < inter_edges->n_edges; i++) {

    assert(edges != NULL);
    assert(mesh != NULL);

    cs_lnum_t  v1_num = edges->def[2*i];
    cs_lnum_t  v2_num = edges->def[2*i+1];
    cs_gnum_t  v1_gnum = (mesh->vertices[v1_num-1]).gnum;
    cs_gnum_t  v2_gnum = (mesh->vertices[v2_num-1]).gnum;
    cs_lnum_t  start = inter_edges->index[i];
    cs_lnum_t  end = inter_edges->index[i+1];

    fprintf(f, "\n%6d: [%9llu] = (%7d [%9llu] - %7d [%9llu])\n",
            i, (unsigned long long)edges->gnum[i],
            v1_num, (unsigned long long)v1_gnum,
            v2_num, (unsigned long long)v2_gnum);

    fprintf(f, "    n_sub_inter: %4d - index : %7d <-- %7d\n",
            end-start, start, end);

    if (inter_edges->vtx_glst == NULL) {

      for (j = start, k = 0; j < end; j++, k++)
        fprintf
          (f, "       %7d (%9d) - (%7llu, %8.6e)\n",
           k, inter_edges->vtx_lst[j],
           (unsigned long long)mesh->vertices[inter_edges->vtx_lst[j]-1].gnum,
           inter_edges->abs_lst[j]);

    }
    else {

      for (j = start, k = 0; j < end; j++, k++)
        fprintf(f, "       %9d - (%7llu, %8.6e)\n",
                k, (unsigned long long)inter_edges->vtx_glst[j],
                inter_edges->abs_lst[j]);

    }

  } /* End of loop on edge intersections */

  fflush(f);
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
