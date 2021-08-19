/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include "cs_math.h"
#include "cs_mesh.h"
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

END_C_DECLS
