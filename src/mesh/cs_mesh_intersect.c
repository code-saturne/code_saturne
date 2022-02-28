/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_geom.h"
#include "fvm_nodal.h"
#include "fvm_point_location.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_intersect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a post-processing mesh for probes
 *
 * parameters:
 *   n_points <-- number of points to locate
 *   coords   <-- point coordinates
 *   m        <-- mesh handle
 *   cell_id  --> cell id containing point, or -1
 *----------------------------------------------------------------------------*/

static void
_locate_points_in_mesh(cs_lnum_t         n_points,
                       const cs_real_t  *coords,
                       const cs_mesh_t  *m,
                       cs_lnum_t         cell_id[])
{
  fvm_nodal_t  *location_mesh = cs_mesh_connect_cells_to_nodal(m,
                                                               "location",
                                                               false,
                                                               m->n_cells,
                                                               NULL);
  float *distance;
  BFT_MALLOC(distance, n_points, float);

  for (cs_lnum_t i = 0; i < n_points; i++) {
    cell_id[i] = -1;
    distance[i] = -1;
  }

  fvm_point_location_nodal(location_mesh,
                           0.,   /* tolerance_base */
                           0.05, /* relative tolerance */
                           1,    /* locate_on_parents */
                           n_points,
                           NULL, /* point_tag */
                           (const cs_coord_t *)coords,
                           cell_id,
                           distance);

  BFT_FREE(distance);

  location_mesh = fvm_nodal_destroy(location_mesh);

  /* Switch to 0-based location */

  for (cs_lnum_t i = 0; i < n_points; i++) {
    if (cell_id[i] > 0)
      cell_id[i] -= 1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check intersection of segment with mesh boundary of parallel
 *        boundary
 *
 * \param[in]   sx0   segment start coordinates
 * \param[in]   sx1   segment end coordinates
 * \param[out]  f_id  id of intersecting face, or -1
 *                    + n_i_faces for boundary faces
 * \param[out]  s     curvilinear coordinates of segment intersection
 */
/*----------------------------------------------------------------------------*/

static void
_check_segment_intersect_boundary(const cs_real_t  *sx0,
                                  const cs_real_t  *sx1,
                                  cs_lnum_t        *f_id,
                                  cs_real_t        *s)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_3_t *vtx_coord= (const cs_real_3_t *)m->vtx_coord;

  cs_real_t t_min    = HUGE_VAL;
  cs_lnum_t f_id_min = -1;

  /* Contribution from interior (parallel boundary) faces. */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t  c_id0 = m->i_face_cells[face_id][0];
    cs_lnum_t  c_id1 = m->i_face_cells[face_id][1];

    /* No need to check intersection for interior faces of the same
       domain: if point is on the exterior, a segment should first
       cross through a boundary or parallel boundary face. */

    if (c_id0 < n_cells && c_id1 < n_cells)
      continue;

    int n_inout[2] = {0, 0};

    cs_lnum_t vtx_start = m->i_face_vtx_idx[face_id];
    cs_lnum_t vtx_end = m->i_face_vtx_idx[face_id+1];
    cs_lnum_t n_vertices = vtx_end - vtx_start;
    const cs_lnum_t *vertex_ids = m->i_face_vtx_lst + vtx_start;

    const cs_real_t *face_center = mq->i_face_cog + (3*face_id);

    double t = cs_geom_segment_intersect_face(0,
                                              n_vertices,
                                              vertex_ids,
                                              vtx_coord,
                                              face_center,
                                              sx0,
                                              sx1,
                                              n_inout,
                                              NULL);

    if (t >= 0 && t <= 1) {
      if (t < t_min) {
        t_min = t;
        f_id_min = face_id;
      }
    }

  }

  /* Contribution from boundary faces*/

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    int n_inout[2] = {0, 0};

    cs_lnum_t vtx_start = m->b_face_vtx_idx[face_id];
    cs_lnum_t vtx_end = m->b_face_vtx_idx[face_id+1];
    cs_lnum_t n_vertices = vtx_end - vtx_start;
    const cs_lnum_t *vertex_ids = m->b_face_vtx_lst + vtx_start;

    const cs_real_t *face_center = mq->b_face_cog + (3*face_id);

    double t = cs_geom_segment_intersect_face(0,
                                              n_vertices,
                                              vertex_ids,
                                              vtx_coord,
                                              face_center,
                                              sx0,
                                              sx1,
                                              n_inout,
                                              NULL);

    if (t >= 0 && t <= 1) {
      if (t < t_min) {
        t_min = t;
        f_id_min = face_id + n_i_faces;
      }
    }

  }

  *f_id = f_id_min;
  *s = t_min;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a given segment
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \param[in]   input     pointer to segment start and end:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_intersect_segment_cell_select(void        *input,
                                      cs_lnum_t   *n_cells,
                                      cs_lnum_t  **cell_ids)
{
  cs_real_t *sx = (cs_real_t *)input;

  const cs_real_t sx0[3] = {sx[0], sx[1], sx[2]};
  const cs_real_t sx1[3] = {sx[3], sx[4], sx[5]};

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t _n_cells = m->n_cells;
  cs_lnum_t *_cell_ids = NULL;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  BFT_MALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Allocate selection list */

  /* Mark for each cell */
  /*--------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < _n_cells; cell_id++) {
    _cell_ids[cell_id] = -1;
  }

  const cs_real_3_t *vtx_coord= (const cs_real_3_t *)m->vtx_coord;

  /* Contribution from interior faces;
     note the to mark cells, we could use a simple loop,
     as thread races would not lead to a incorrect result, but
     even if is slightly slower, we prefer to have a clean
     behavior under thread debuggers. */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        int n_inout[2] = {0, 0};

        cs_lnum_t vtx_start = m->i_face_vtx_idx[face_id];
        cs_lnum_t vtx_end = m->i_face_vtx_idx[face_id+1];
        cs_lnum_t n_vertices = vtx_end - vtx_start;
        const cs_lnum_t *vertex_ids = m->i_face_vtx_lst + vtx_start;

        const cs_real_t *face_center = fvq->i_face_cog + (3*face_id);

        double t = cs_geom_segment_intersect_face(0,
                                                  n_vertices,
                                                  vertex_ids,
                                                  vtx_coord,
                                                  face_center,
                                                  sx0,
                                                  sx1,
                                                  n_inout,
                                                  NULL);

        if (t >= 0 && t <= 1) {
          cs_lnum_t  c_id0 = m->i_face_cells[face_id][0];
          cs_lnum_t  c_id1 = m->i_face_cells[face_id][1];
          if (c_id0 < _n_cells)
            _cell_ids[c_id0] = 1;
          if (c_id1 < _n_cells)
            _cell_ids[c_id1] = 1;
        }

      }

    }

  }

  /* Contribution from boundary faces*/

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        int n_inout[2] = {0, 0};

        cs_lnum_t vtx_start = m->b_face_vtx_idx[face_id];
        cs_lnum_t vtx_end = m->b_face_vtx_idx[face_id+1];
        cs_lnum_t n_vertices = vtx_end - vtx_start;
        const cs_lnum_t *vertex_ids = m->b_face_vtx_lst + vtx_start;

        const cs_real_t *face_center = fvq->b_face_cog + (3*face_id);

        double t = cs_geom_segment_intersect_face(0,
                                                  n_vertices,
                                                  vertex_ids,
                                                  vtx_coord,
                                                  face_center,
                                                  sx0,
                                                  sx1,
                                                  n_inout,
                                                  NULL);

        if (t >= 0 && t <= 1) {
          cs_lnum_t  c_id = m->b_face_cells[face_id];
          _cell_ids[c_id] = 1;
        }

      }

    }

  }

  /* Now check marked cells */

  _n_cells = 0;
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    if (_cell_ids[cell_id] >= 0)
      _cell_ids[_n_cells++] = cell_id;
  }

  BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                  but not required) */

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a line composed of segments
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \param[in]   input     pointer to segments starts and ends:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[in]   n_points  number of vertices in the polyline
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 * \param[out]  seg_c_len array of length of the segment in the selected cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_intersect_polyline_cell_select(void        *input,
                                       cs_lnum_t    n_points,
                                       cs_lnum_t   *n_cells,
                                       cs_lnum_t  **cell_ids,
                                       cs_real_t  **seg_c_len)
{
  cs_real_t *sx = (cs_real_t *)input;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t _n_cells = m->n_cells;
  cs_lnum_t *_cell_ids = NULL;
  cs_lnum_t *_in = NULL;
  cs_lnum_t *_out = NULL;
  cs_real_t *_seg_c_len = NULL;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  BFT_MALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Allocate selection list */
  BFT_MALLOC(_seg_c_len, _n_cells, cs_real_t); /* Allocate selection list length */
  BFT_MALLOC(_in, _n_cells, cs_lnum_t);
  BFT_MALLOC(_out, _n_cells, cs_lnum_t);

  /* Mark for each cell */
  /*--------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < _n_cells; cell_id++) {
    _cell_ids[cell_id] = -1;
    _seg_c_len[cell_id] = 0.;
  }

  const cs_real_3_t *vtx_coord= (const cs_real_3_t *)m->vtx_coord;

  /* Loop over the vertices of the polyline */
  for (cs_lnum_t s_id = 0; s_id < n_points - 1; s_id++) {

    const cs_real_t *sx0 = &(sx[3*s_id]);
    const cs_real_t *sx1 = &(sx[3*(s_id+1)]);

    cs_real_t length =  cs_math_3_distance(sx0, sx1);

    /* Count the number of ingoing and outgoing intersections
     * to check if the segment is inside the cell */
    for (cs_lnum_t cell_id = 0; cell_id < _n_cells; cell_id++) {
      _in[cell_id] = 0;
      _out[cell_id] = 0;
    }

    /* Contribution from interior faces;
       note the to mark cells, we could use a simple loop,
       as thread races would not lead to a incorrect result, but
       even if is slightly slower, we prefer to have a clean
       behavior under thread debuggers. */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
            face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
            face_id++) {

          cs_lnum_t vtx_start = m->i_face_vtx_idx[face_id];
          cs_lnum_t vtx_end = m->i_face_vtx_idx[face_id+1];
          cs_lnum_t n_vertices = vtx_end - vtx_start;
          const cs_lnum_t *vertex_ids = m->i_face_vtx_lst + vtx_start;

          const cs_real_t *face_center = fvq->i_face_cog + (3*face_id);

          cs_lnum_t c_id0 = m->i_face_cells[face_id][0];
          cs_lnum_t c_id1 = m->i_face_cells[face_id][1];

          /* The line (OD) goes in (n_inout[0]++)
           * or goes out (n_inout[1]++) the cell */
          int n_inout[2] = {0, 0};

          double t = cs_geom_segment_intersect_face(0,
                                                    n_vertices,
                                                    vertex_ids,
                                                    vtx_coord,
                                                    face_center,
                                                    sx0,
                                                    sx1,
                                                    n_inout,
                                                    NULL);

          /* Segment is inside cell i if
           *  n_inout[0] > 0 and t < 0 for a face
           *  and
           *  n_inout[1] > 0 and t > 0 for an other face
           */
          if (c_id0 < _n_cells) {
            /* Intersection of (OD) with the face
             * may be on [OD)
             * It may leave c_id0 */
            if (t >= 0.)
              _out[c_id0] += n_inout[1];

            /* Intersection of (OD) with the face
             * may be on (OD] (because t < 0)
             * It may enter c_id0 */
            if (t < 0)
              _in[c_id0] += n_inout[0];
          }
          /* Segment is inside cell j if
           *  n_inout[0] > 0 and t < 0 for a face
           *  and
           *  n_inout[1] > 0 and t > 0 for an other face
           */
          if (c_id1 < _n_cells) {
            /* Intersection of (OD) with the face
             * may be on [OD) (because t > 0)
             * It may leave c_id1 (if OD.n < 0) */
            if (t >= 0.)
              _out[c_id1] += n_inout[0];

            /* Intersection of (OD) with the face
             * may be on (OD]
             * It may enter c_id1 (if 0D.n > 0) */
            if (t < 0.)
              _in[c_id1] += n_inout[1];
          }

          /* Segment crosses the face */
          if (t >= 0 && t <= 1) {
            /* length upwind the face*/
            cs_real_t length_up =  t * length;
            /* length downwind the face*/
            cs_real_t length_down =  (1.-t) * length;
            if (c_id0 < _n_cells) {

              /* Mark cell by segment id (the cell may already be marked by another
               * segment */
              _cell_ids[c_id0] = s_id;

              /* OD enters cell i from cell j */
              if (n_inout[0] > 0)
                _seg_c_len[c_id0] -= length_up;

              /* OD leaves cell i to cell j */
              if (n_inout[1] > 0)
                _seg_c_len[c_id0] -= length_down;

            }
            if (c_id1 < _n_cells) {

              /* Mark cell by segment id (the cell may already be marked by another
               * segment */
              _cell_ids[c_id1] = s_id;

              /* OD enters cell i from cell j
               * so leaves cell j */
              if (n_inout[0] > 0)
                _seg_c_len[c_id1] -= length_down;

              /* OD leaves cell i to cell j
               * so enters cell j */
              if (n_inout[1] > 0)
                _seg_c_len[c_id1] -= length_up;

            }
          }
        }

      }

    }

    /* Contribution from boundary faces*/

    for (int g_id = 0; g_id < n_b_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_b_threads; t_id++) {

        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
            face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
            face_id++) {

          cs_lnum_t vtx_start = m->b_face_vtx_idx[face_id];
          cs_lnum_t vtx_end = m->b_face_vtx_idx[face_id+1];
          cs_lnum_t n_vertices = vtx_end - vtx_start;
          const cs_lnum_t *vertex_ids = m->b_face_vtx_lst + vtx_start;

          const cs_real_t *face_center = fvq->b_face_cog + (3*face_id);
          cs_lnum_t  c_id = m->b_face_cells[face_id];

          int n_inout[2] = {0, 0};

          double t = cs_geom_segment_intersect_face(0,
                                                    n_vertices,
                                                    vertex_ids,
                                                    vtx_coord,
                                                    face_center,
                                                    sx0,
                                                    sx1,
                                                    n_inout,
                                                    NULL);
          /* Segment is inside cell i if
           *  n_inout[0] > 0 and t < 0 for a face
           *  and
           *  n_inout[1] > 0 and t > 0 for an other face
           */
          if (c_id < _n_cells) {
            /* Intersection of (OD) with the face
             * may be on [OD)
             * It may leave c_id */
            if (t >= 0.)
              _out[c_id] += n_inout[1];

            /* Intersection of (OD) with the face
             * may be on (OD]
             * It may enter c_id */
            if (t < 0)
              _in[c_id] += n_inout[0];
          }

          /* Segment crosses the face */
          if (t >= 0 && t <= 1) {

            /* length upwind the face*/
            cs_real_t length_up =  t * length;
            /* length downwind the face*/
            cs_real_t length_down =  (1.-t) * length;

            /* Mark cell by segment id (the cell may already be marked by another
             * segment */
            _cell_ids[c_id] = s_id;

            /* OD enters cell i */
            if (n_inout[0] > 0)
              _seg_c_len[c_id] -= length_up;

            /* OD leaves cell i */
            if (n_inout[1] > 0)
              _seg_c_len[c_id] -= length_down;

          }

        }
      }

    }

    /* Finalize the length computation to deal with cases where the segment
     * is inside the cell */
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
      /* There is one intersection on the left of [OD)
       * and one on the right of [OD) which means that
       * O is inside the cell */
      if ((_in[cell_id] > 0 && _out[cell_id] > 0)
          || (_cell_ids[cell_id] == s_id)) {
        _cell_ids[cell_id] = s_id;
        _seg_c_len[cell_id] += length;
      }
    }

  } /* End loop over the segments */

  BFT_FREE(_in);
  BFT_FREE(_out);

  /* Now check marked cells and renumber */
  _n_cells = 0;
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    if (_cell_ids[cell_id] >= 0) {
      _cell_ids[_n_cells] = cell_id;
      _seg_c_len[_n_cells] = _seg_c_len[cell_id];
      _n_cells++;
    }
  }

  BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                  but not required) */
  BFT_REALLOC(_seg_c_len, _n_cells, cs_real_t);

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
  *seg_c_len = _seg_c_len;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map a polyline of segments to the local mesh.
 *
 * The caller is responsible for freeing the returned seg_cell* arrays.
 *
 * \param[in]   n_points           number of vertices in the polyline
 * \param[in]   point_coords       pointer to segments starts and ends:
 *                                 [x0, y0, z0, x1, y1, z1, ...]
 * \param[in]   min_fraction       minimum fraction of each edge under which
 *                                 fraction is ignored
 * \param[out]  seg_cell_idx       index of start an end cell ids per segment
 *                                 (size: n_points == n_segments-1)
 * \param[out]  seg_cell           cells intersected by each segment i are given
 *                                 by seg_cell[], with
 *                                 seg_cell_idx[i] <= j < with seg_cell_idx[i+1]
 * \param[out]  seg_cell_fraction  fraction of each segment in intersected
 *                                 cells (same index as seg_cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_intersect_polyline_map(cs_lnum_t          n_points,
                               const cs_real_t    point_coords[],
                               cs_real_t          min_fraction,
                               cs_lnum_t        **seg_cell_idx,
                               cs_lnum_t        **seg_cell,
                               cs_real_t        **seg_cell_fraction)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;

  cs_lnum_t n_cells = m->n_cells;
  cs_lnum_t *_seg_cell_idx = NULL;
  cs_lnum_t *_seg_cell = NULL;
  cs_real_t *_seg_cell_fraction = NULL;

  *seg_cell_idx = _seg_cell_idx;
  *seg_cell = _seg_cell;
  *seg_cell_fraction = _seg_cell_fraction;

  if (n_points < 1 || n_cells < 1)
    return;

  /* Locate cells in mesh, so as to robustly handle various cases,
     including when polyline enters/exists mesh multiple times. */

  cs_lnum_t *point_cell_id;
  BFT_MALLOC(point_cell_id, n_points, cs_lnum_t);

  _locate_points_in_mesh(n_points,
                         point_coords,
                         m,
                         point_cell_id);

  /* Base return array allocation;

     Size is unknown at start, so we resize when needed; We use a small initial
     estimate, toallows checking that realloction is handled correctly. */

  cs_lnum_t alloc_size = 10;
  BFT_MALLOC(_seg_cell_idx, n_points, cs_lnum_t);
  BFT_MALLOC(_seg_cell, alloc_size, cs_lnum_t);
  BFT_MALLOC(_seg_cell_fraction, alloc_size, cs_real_t);

  /* Build cell to face mapping
     (in the future, this may be avaiable as a mes ajacecy optional feature) */

  cs_lnum_t *cell_face_idx, *cell_face_lst;
  cs_lnum_t *counter = NULL;

  BFT_MALLOC(counter, n_cells, cs_lnum_t);
  BFT_MALLOC(cell_face_idx, n_cells + 1, cs_lnum_t);

  cell_face_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cell_face_idx[i+1] = 0;
    counter[i] = 0;
  }

  /* Counting pass */

  for (cs_lnum_t i = 0; i < m->n_i_faces; i++)
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t iel = m->i_face_cells[i][j] + 1;
      if (iel <= m->n_cells)
        cell_face_idx[iel] += 1;
    }

  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    cell_face_idx[m->b_face_cells[i] + 1] += 1;

  /* Build index */

  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    cell_face_idx[i+1] += cell_face_idx[i];

  BFT_MALLOC(cell_face_lst,
             cell_face_idx[m->n_cells], cs_lnum_t);

  /* Build pass: border faces are < 0 and interior faces > 0 */

  for (cs_lnum_t i = 0; i < m->n_i_faces; i++) {
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t c_id = m->i_face_cells[i][j];
      if (c_id < m->n_cells) {
        cs_lnum_t  shift = cell_face_idx[c_id] + counter[c_id];
        cell_face_lst[shift] = i+1;
        counter[c_id] += 1;
      }
    }
  }

  for (cs_lnum_t i = 0; i < m->n_b_faces; i++) {
    cs_lnum_t  c_id = m->b_face_cells[i];
    cs_lnum_t  shift = cell_face_idx[c_id] + counter[c_id];
    cell_face_lst[shift] = -(i+1);
    counter[c_id] += 1;
  }

  BFT_FREE(counter);

  /* Mark for each cell */
  /*--------------------*/

  _seg_cell_idx[0] = 0;
  _seg_cell_idx[1] = 0;

  cs_lnum_t idx = 0;
  cs_real_t frac_truncated = 0;

  /* Loop over the vertices of the polyline */
  for (cs_lnum_t s_id = 0; s_id < n_points - 1; s_id++) {

    _seg_cell_idx[s_id] = idx;

    if (s_id > 0 && frac_truncated > 0) {
      cs_real_t mult = 1. + frac_truncated;
      for (cs_lnum_t i = _seg_cell_idx[s_id-1]; i < _seg_cell_idx[s_id]; i++)
        _seg_cell_fraction[i] *= mult;
    }

    cs_real_t s_start[3] = {point_coords[3*s_id],
                            point_coords[3*s_id + 1],
                            point_coords[3*s_id + 2]};

    const cs_real_t *s_end = point_coords + 3*(s_id+1);
    const cs_real_t length_ini = cs_math_3_distance(s_start, s_end);

    cs_lnum_t cur_cell_id = point_cell_id[s_id];
    cs_real_t length = length_ini;

    while (length > 1e-12) {

      /* If a point starts outside the local mesh, check it segment
         intersects the boundary (we cannot simply check if the next point
         is inside the mesh, as two points could be outsie the mesh
         but their joining segment cross that mesh). */

      if (point_cell_id[s_id] < 0) {
        cs_lnum_t f_id = -1;
        cs_real_t s = HUGE_VAL;
        _check_segment_intersect_boundary(s_start,
                                          s_end,
                                          &f_id,
                                          &s);

        if (f_id >= 0) {
          for (cs_lnum_t k = 0; k < 3; k++) {
            s_start[k] += (s_end[k] - s_start[k])*s;
          }
          length = cs_math_3_distance(s_start, s_end);
          if (f_id >= m->n_i_faces) {
            f_id -= m->n_i_faces;
            cur_cell_id = m->b_face_cells[f_id];
          }
          else {
            cs_lnum_t c_id0 = m->i_face_cells[f_id][0];
            cs_lnum_t c_id1 = m->i_face_cells[f_id][1];
            if (c_id0 < m->n_cells)
              cur_cell_id = m->b_face_cells[c_id0];
            else
              cur_cell_id = m->b_face_cells[c_id1];
          }
        }
        else {
          /* No intersection of this segment with mesh face
             (i.e. segment outside local mesh) */
          break;
        }
      }

      /* Case where segment is contained in cell (ignore potential
         crossings if cell is non-convex, and handle it as convex). */

      if (cur_cell_id == point_cell_id[s_id+1]) {
        if (idx >= alloc_size) {
          alloc_size *= 2;
          BFT_REALLOC(_seg_cell, alloc_size, cs_lnum_t);
          BFT_REALLOC(_seg_cell_fraction, alloc_size, cs_real_t);
        }

        _seg_cell[idx] = cur_cell_id;
        _seg_cell_fraction[idx] = length/length_ini;
        length = 0;

        if (_seg_cell_fraction[idx] >= min_fraction)
          idx += 1;
        else
          frac_truncated += _seg_cell_fraction[idx];

        break;  /* switch to next segment */
      }

      cs_lnum_t exit_face = 0; /* > 0 for interior faces,
                                  < 0 for boundary faces */

      double adist_min = 2.;

      double t_intersect = -1;

      int n_in = 0;
      int n_out = 0;

      /* Loop on faces to see if the segment crosses it */
      for (cs_lnum_t i = cell_face_idx[cur_cell_id];
           i < cell_face_idx[cur_cell_id+1];
           i++) {

        cs_lnum_t face_id, vtx_start, vtx_end, n_vertices;
        const cs_lnum_t *face_connect;
        const cs_real_t *face_cog;

        /* Outward normal: always well oriented for external faces,
           depends on the connectivity for internal faces */
        int reorient_face = 1;

        cs_lnum_t face_num = cell_face_lst[i];

        if (face_num > 0) {

          /* Interior face */

          face_id = face_num - 1;
          if (cur_cell_id == m->i_face_cells[face_id][1])
            reorient_face = -1;
          vtx_start = m->i_face_vtx_idx[face_id];
          vtx_end = m->i_face_vtx_idx[face_id+1];
          n_vertices = vtx_end - vtx_start;

          face_connect = m->i_face_vtx_lst + vtx_start;
          face_cog = mq->i_face_cog + 3*face_id;

        }
        else {

          assert(face_num < 0);

          /* Boundary faces */

          face_id = -face_num - 1;
          vtx_start = m->b_face_vtx_idx[face_id];
          vtx_end = m->b_face_vtx_idx[face_id+1];
          n_vertices = vtx_end - vtx_start;

          face_connect = m->b_face_vtx_lst + vtx_start;
          face_cog = mq->b_face_cog + 3*face_id;

        }

        /*
          adimensional distance estimation of face intersection
          (2 if no chance of intersection)
        */

        int n_crossings[2] = {0, 0};

        double t = cs_geom_segment_intersect_face
                     (reorient_face,
                      n_vertices,
                      face_connect,
                      (const cs_real_3_t *)m->vtx_coord,
                      face_cog,
                      s_start,
                      s_end,
                      n_crossings,
                      NULL);

        n_in += n_crossings[0];
        n_out += n_crossings[1];

        /* Store the nearest intesection from the O point...*/
        if (t < adist_min && t >= 0) {
          exit_face = face_num;
          t_intersect = t;
          adist_min = t_intersect;
        }

      }

      if (idx >= alloc_size) {
        alloc_size *= 2;
        BFT_REALLOC(_seg_cell, alloc_size, cs_lnum_t);
        BFT_REALLOC(_seg_cell_fraction, alloc_size, cs_real_t);
      }

      _seg_cell[idx] = cur_cell_id;

      /* If segment is inside cell (should not occur if tests are consistent,
         but could occur dur to location tests not being based on segment/face
         intersections: in this case, we assume the end-point is very close
         to a face, but inside that face) */

      if (exit_face == 0) {
        _seg_cell_fraction[idx] = length;
        length = 0;
      }

      /* If segment exits cell */
      else {

        _seg_cell_fraction[idx] = length/length_ini*t_intersect;

        for (cs_lnum_t k = 0; k < 3; k++) {
          s_start[k] += (s_end[k] - s_start[k])*t_intersect;
        }
        length = cs_math_3_distance(s_start, s_end);

        if (exit_face > 0) { /* Segment crosses to the neighbor cell
                                through the current face "face_num" */
          cs_lnum_t face_id = exit_face - 1;

          cs_lnum_t  c_id0 = m->i_face_cells[face_id][0];
          cs_lnum_t  c_id1 = m->i_face_cells[face_id][1];

          if (c_id0 == cur_cell_id)
            cur_cell_id = c_id1;
          else
            cur_cell_id = c_id0;

          if (cur_cell_id >= n_cells)
            cur_cell_id = -1;
        }
        else { /* (exit_face < 0), segment crosses the boundary
                  through the current face "face_num" */
          cur_cell_id = -1;
        }
      }

      if (_seg_cell_fraction[idx] >= min_fraction)
        idx += 1;
      else
        frac_truncated += _seg_cell_fraction[idx];

    }  /* End for this segment (length > 0) */

  } /* End loop over the segments */

  _seg_cell_idx[n_points-1] = idx;

  BFT_FREE(cell_face_lst);
  BFT_FREE(cell_face_idx);
  BFT_FREE(point_cell_id);

  /* Set return values */

  alloc_size = _seg_cell_idx[n_points-1];
  BFT_REALLOC(_seg_cell, alloc_size, cs_lnum_t);
  BFT_REALLOC(_seg_cell_fraction, alloc_size, cs_real_t);

  *seg_cell_idx = _seg_cell_idx;
  *seg_cell = _seg_cell;
  *seg_cell_fraction = _seg_cell_fraction;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
