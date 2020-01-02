/*============================================================================
 * Mesh smoothing.
 *============================================================================*/

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
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_mesh_quality.h"
#include "cs_all_to_all.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_smoother.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/
#undef _DOT_PRODUCT_3D

#define UNWARPING_MAX_LOOPS 50
#define UNWARPING_MVT 0.1
#define _PI_ atan(1.0)*4.0

#define _DOT_PRODUCT_3D(v1, v2) ( \
 v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the distance between two vertices.
 *
 * parameters:
 *   a <-- coordinates of the first vertex.
 *   b <-- coordinates of the second vertex.
 *
 * returns:
 *   distance between a and b.
 *---------------------------------------------------------------------------*/

inline static double
_compute_distance(const double  a[3],
                  const double  b[3])
{
  double  distance;

  distance = sqrt(  (b[0] - a[0])*(b[0] - a[0])
                  + (b[1] - a[1])*(b[1] - a[1])
                  + (b[2] - a[2])*(b[2] - a[2]));

  return distance;
}

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_minmax(cs_int_t            n_vals,
                const cs_real_t     var[],
                cs_real_t          *min,
                cs_real_t          *max)
{
  cs_int_t  i;
  cs_real_t  _min = DBL_MAX, _max = -DBL_MAX;

  for (i = 0; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(&_min, min, 1, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);

    MPI_Allreduce(&_max, max, 1, CS_MPI_REAL, MPI_MAX,
                  cs_glob_mpi_comm);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    *min = _min;
    *max = _max;
  }
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector.
 *
 * parameters:
 *   n_steps <-- number of histogram steps
 *   var_min <-- minimum variable value (histogram scaling)
 *   var_max <-- maximum variable value (histogram scaling)
 *   min     <-- minimum variable value to display
 *   max     <-- maximum variable value to display
 *   count   <-> count for each histogram slice (size: n_steps)
 *               local values in, global values out
 *----------------------------------------------------------------------------*/

static void
_display_histograms(int           n_steps,
                    cs_real_t     var_min,
                    cs_real_t     var_max,
                    cs_real_t     min,
                    cs_real_t     max,
                    cs_gnum_t     count[])
{
  int  i, j;
  double var_step;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t _g_count[10];
    cs_gnum_t *g_count = _g_count;

    MPI_Allreduce(count, g_count, n_steps, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0; i < n_steps; i++)
      count[i] = g_count[i];
  }

#endif

  /* Print base min, max, and increment */

  bft_printf(_("    minimum value =         %10.5e\n"), (double)min);
  bft_printf(_("    maximum value =         %10.5e\n\n"), (double)max);

  var_step = CS_ABS(var_max - var_min) / n_steps;

  if (CS_ABS(var_max - var_min) > 0.) {

    /* Number of elements in each subdivision */

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10.5e ; %10.5e [ = %10llu\n",
                 i+1, var_min + i*var_step, var_min + j*var_step,
                 (unsigned long long)(count[i]));

    bft_printf("    %3d : [ %10.5e ; %10.5e ] = %10llu\n",
               n_steps,
               var_min + (n_steps - 1)*var_step,
               var_max,
               (unsigned long long)(count[n_steps - 1]));

  }
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector on cells or
 * boundary faces.
 *
 * parameters:
 *   n_vals     <-- number of values
 *   var        <-- pointer to vector (size: n_vals)
 *   _min       <-- min variable value for histogram
 *   _max       <-- max variable value for histogram
 *   n_min      <-- min variable value for display
 *   n_max      <-- max variable value for display
 *----------------------------------------------------------------------------*/

static void
_histogram(cs_lnum_t             n_vals,
           const cs_real_t       var[],
           cs_real_t             min,
           cs_real_t             max,
           cs_real_t             n_min,
           cs_real_t             n_max)
{
  cs_lnum_t  i;
  int        j, k;

  cs_real_t  step;

  cs_gnum_t count[10];
  const int  n_steps = 10;

  assert (sizeof(double) == sizeof(cs_real_t));

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (CS_ABS(max - min) > 0.) {

    step = CS_ABS(max - min) / n_steps;

    /* Loop on values */

    for (i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < min + k*step)
          break;
      }
      count[j] += 1;

    }

  }

  _display_histograms(n_steps, min, max, n_min, n_max, count);
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector on interior faces.
 *
 * parameters:
 *   mesh   <-- pointer to mesh structure
 *   var    <-- pointer to vector (size: mesh->n_i_faces)
 *   _min   <-- min variable value for histogram
 *   _max   <-- max variable value for histogram
 *   n_min  <-- min variable value for display
 *   n_max  <-- max variable value for display
 *----------------------------------------------------------------------------*/

static void
_int_face_histogram(const cs_mesh_t      *mesh,
                    const cs_real_t       var[],
                    cs_real_t             min,
                    cs_real_t             max,
                    cs_real_t             n_min,
                    cs_real_t             n_max)
{
  cs_lnum_t  i;
  int        j, k;

  cs_real_t  step;

  cs_gnum_t count[8];
  const int  n_steps = 8;

  assert(sizeof(double) == sizeof(cs_real_t));


  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (CS_ABS(max - min) > 0.) {

    step = CS_ABS(max - min) / n_steps;

    /* Loop on faces */

    for (i = 0; i < mesh->n_i_faces; i++) {

      if (mesh->i_face_cells[i][0] >= mesh->n_cells)
        continue;

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < min + k*step)
          break;
      }
      count[j] += 1;

    }

  }

  _display_histograms(n_steps, min, max, n_min, n_max, count);
}

/*----------------------------------------------------------------------------
 * Move the vertices from local displacement
 *
 * parameters:
 *   mesh         <-> pointer to a cs_mesh_t structure
 *   vtx_mvt      <-- local displacement
 *   vtx_is_fixed <-- array to mark fixed vertices (1 : fixed, 0 : free)
 *---------------------------------------------------------------------------*/

static void
_move_vertices(cs_mesh_t  *mesh,
               cs_real_t  *vtx_mvt,
               const int   vtx_is_fixed[])
{
  int i, k;
  for (i = 0; i < mesh->n_vertices; i++) {
    if (vtx_is_fixed[i] == 0) {
      for (k = 0; k < 3; k++)
        mesh->vtx_coord[3*i + k] += vtx_mvt[3*i + k];
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute tolerance
 * tolerance = shortest edge length * fraction
 *
 * parameters:
 *   vertex_coords    <--  coordinates of vertices.
 *   vertex_tolerance <->  local tolerance affected to each vertex and
 *                         to be updated
 *   n_faces          <--  number of selected faces
 *   face_vtx_idx     <--  "face -> vertex" connect. index
 *   face_vtx_lst     <--  "face -> vertex" connect. list
 *   fraction         <--  parameter used to compute the tolerance
 *---------------------------------------------------------------------------*/

static void
_get_local_tolerance(const cs_real_t   vtx_coords[],
                     double            vtx_tolerance[],
                     const cs_int_t    n_faces,
                     const cs_int_t    face_vtx_idx[],
                     const cs_int_t    face_vtx_lst[],
                     double            fraction)
{
  cs_lnum_t  j, k, start, end, face_id, vtx_id1, vtx_id2;
  cs_real_t  length, tolerance;
  cs_real_t  a[3], b[3];

  for (face_id = 0; face_id < n_faces; face_id++) {

    start = face_vtx_idx[face_id];
    end = face_vtx_idx[face_id + 1];

    /* Loop on the vertices of the face */

    for (j = start; j < end - 1; j++) {

      vtx_id1 = face_vtx_lst[j];
      vtx_id2 = face_vtx_lst[j+1];

      for (k = 0; k < 3; k++) {
        a[k] = vtx_coords[3*vtx_id1 + k];
        b[k] = vtx_coords[3*vtx_id2 + k];
      }

      length = _compute_distance(a, b);
      tolerance = length * fraction;
      vtx_tolerance[vtx_id1] = CS_MIN(vtx_tolerance[vtx_id1], tolerance);
      vtx_tolerance[vtx_id2] = CS_MIN(vtx_tolerance[vtx_id2], tolerance);

    }

    /* Case end - start */

    vtx_id1 = face_vtx_lst[end-1];
    vtx_id2 = face_vtx_lst[start];

    for (k = 0; k < 3; k++) {
      a[k] = vtx_coords[3*vtx_id1 + k];
      b[k] = vtx_coords[3*vtx_id2 + k];
    }

    length = _compute_distance(a, b);
    tolerance = length * fraction;
    vtx_tolerance[vtx_id1] = CS_MIN(vtx_tolerance[vtx_id1], tolerance);
    vtx_tolerance[vtx_id2] = CS_MIN(vtx_tolerance[vtx_id2], tolerance);

  } /* End of loop on faces */

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Exchange local vertex tolerances to get a global vertex tolerance.
 *
 * parameters:
 *   mesh              <--  pointer to a cs_mesh_t structure
 *   vtx_tolerance     <->  local tolerance affected to each vertex and
 *                          to be updated
 *---------------------------------------------------------------------------*/

static void
_get_global_tolerance(cs_mesh_t            *mesh,
                      cs_real_t            *vtx_tolerance)
{
  cs_lnum_t  i, vtx_id;
  cs_gnum_t  first_vtx_gnum;

  cs_lnum_t n_vertices = mesh->n_vertices;
  double  *g_vtx_tolerance = NULL, *send_list = NULL;
  cs_gnum_t  *send_glist = NULL;
  cs_gnum_t  n_g_vertices = mesh->n_g_vertices;
  const cs_gnum_t  *io_gnum = mesh->global_vtx_num;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;

  /* Define a fvm_io_num_t structure on vertices */

  cs_block_dist_info_t
    bi = cs_block_dist_compute_sizes(local_rank,
                                     n_ranks,
                                     1,
                                     0,
                                     n_g_vertices);

  cs_int_t block_size = bi.block_size;

  cs_all_to_all_t
    *d = cs_all_to_all_create_from_block(n_vertices,
                                         0, /* flags */
                                         io_gnum,
                                         bi,
                                         mpi_comm);

  /* Send the global numbering for each vertex */

  cs_gnum_t *recv_glist = cs_all_to_all_copy_array(d,
                                                   CS_GNUM_TYPE,
                                                   1,
                                                   false, /* reverse */
                                                   io_gnum,
                                                   NULL);

  /* Send the vertex tolerance for each vertex */

  double *recv_list = cs_all_to_all_copy_array(d,
                                               CS_REAL_TYPE,
                                               1,
                                               false, /* reverse */
                                               vtx_tolerance,
                                               NULL);

  cs_lnum_t n_recv = cs_all_to_all_n_elts_dest(d);

  /* Define the global tolerance array */

  BFT_MALLOC(g_vtx_tolerance, block_size, double);

  for (i = 0; i < block_size; i++)
    g_vtx_tolerance[i] = DBL_MAX;

  first_vtx_gnum = block_size * local_rank + 1;

  for (i = 0; i < n_recv; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    g_vtx_tolerance[vtx_id] = CS_MIN(g_vtx_tolerance[vtx_id], recv_list[i]);
  }

  /* Replace local vertex tolerance by the new computed global tolerance */

  for (i = 0; i < n_recv; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    recv_list[i] = g_vtx_tolerance[vtx_id];
  }

  cs_all_to_all_copy_array(d,
                           CS_REAL_TYPE,
                           1,
                           true, /* reverse */
                           recv_list,
                           vtx_tolerance);

  /* Free memory */

  cs_all_to_all_destroy(&d);

  BFT_FREE(recv_glist);
  BFT_FREE(send_glist);
  BFT_FREE(send_list);
  BFT_FREE(recv_list);
  BFT_FREE(g_vtx_tolerance);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Define for each vertex a tolerance which is the radius of the
 * sphere in which the vertex can be fused with another vertex.
 *
 * parameters:
 *   mesh             <--  pointer to a cs_mesh_t structure
 *   vtx_tolerance    -->  tolerance affected to each vertex
 *   fraction         <--  parameter used to compute the tolerance
 *---------------------------------------------------------------------------*/

static void
_get_tolerance(cs_mesh_t   *mesh,
               cs_real_t   *vtx_tolerance,
               double       fraction)
{
  cs_lnum_t i;

  for (i = 0; i < mesh->n_vertices; i++)
    vtx_tolerance[i] = DBL_MAX;

  /* Define local tolerance */

  _get_local_tolerance(mesh->vtx_coord,
                       vtx_tolerance,
                       mesh->n_b_faces,
                       mesh->b_face_vtx_idx,
                       mesh->b_face_vtx_lst,
                       fraction);

  _get_local_tolerance(mesh->vtx_coord,
                       vtx_tolerance,
                       mesh->n_i_faces,
                       mesh->i_face_vtx_idx,
                       mesh->i_face_vtx_lst,
                       fraction);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    /* Global number of selected vertices and associated
       fvm_io_num_t structure */

    _get_global_tolerance(mesh,
                          vtx_tolerance);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Unwarping algorithm, called by _unwarping
 *
 * parameters:
 *   mesh                <--  pointer to a cs_mesh_t structure
 *   i_face_norm         <--  surface normal of interior faces
 *   b_face_norm         <--  surface normal of border faces
 *   i_face_cog          <--  centre of gravity of interior faces
 *   b_face_cog          <--  centre of gravity of border faces
 *   loc_vtx_mvt         -->  local vertices displacement
 *   i_face_warp         <--  value of interior faces warping
 *   b_face_warp         <--  value of border faces warping
 *   vtx_tolerance       <--  local mouvement tolerance
 *   frac                <--  tolerance fraction
 *
 * returns:
 *   value of the most warped face
 *----------------------------------------------------------------------------*/

static cs_real_t
_unwarping_mvt(cs_mesh_t            *mesh,
               cs_real_t            *i_face_norm,
               cs_real_t            *b_face_norm,
               cs_real_t            *i_face_cog,
               cs_real_t            *b_face_cog,
               cs_real_t            *loc_vtx_mvt,
               cs_real_t            *i_face_warp,
               cs_real_t            *b_face_warp,
               cs_real_t            *vtx_tolerance,
               double                frac)
{
  cs_lnum_t face_id, i;
  int coord_id;
  cs_lnum_t start_id, end_id, vtx;
  cs_real_t lambda;
  cs_real_t max_vtxtol = 0.;
  cs_real_t maxwarp = 0.;

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++)
    if (maxwarp < i_face_warp[face_id])
      maxwarp = i_face_warp[face_id];
  for (face_id = 0; face_id < mesh->n_b_faces; face_id++)
    if (maxwarp < b_face_warp[face_id])
      maxwarp = b_face_warp[face_id];

  for (i = 0; i < mesh->n_vertices*3; i++)
    loc_vtx_mvt[i] = 0.0;
  for (i = 0; i < mesh->n_vertices; i++)
    if (vtx_tolerance[i] > max_vtxtol)
      max_vtxtol = vtx_tolerance[i];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_real_t maxpar[2];
    cs_real_t _maxpar[2];
    maxpar[0] = maxwarp;
    maxpar[1] = max_vtxtol;

    MPI_Allreduce(maxpar, _maxpar, 2, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    maxwarp = _maxpar[0];
    max_vtxtol = _maxpar[1];
  }
#endif

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
    start_id = mesh->b_face_vtx_idx[face_id];
    end_id = mesh->b_face_vtx_idx[face_id + 1];
    for (i = start_id; i < end_id; i++) {
      vtx = mesh->b_face_vtx_lst[i];
      lambda = 0.0;
      for (coord_id = 0; coord_id < 3; coord_id++)
        lambda +=  (mesh->vtx_coord[3*vtx + coord_id]
                    - b_face_cog[3*face_id + coord_id])
                    * b_face_norm[3*face_id + coord_id];

      for (coord_id = 0; coord_id < 3; coord_id++) {
        loc_vtx_mvt[vtx*3 + coord_id] -=
          lambda * b_face_norm[3*face_id + coord_id]
                 * UNWARPING_MVT * (b_face_warp[face_id]/maxwarp)
                 * (vtx_tolerance[vtx]/(max_vtxtol*frac));
      }
    }
  }


  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    if (mesh->i_face_cells[face_id][0] < mesh->n_cells) {
      start_id = mesh->i_face_vtx_idx[face_id];
      end_id = mesh->i_face_vtx_idx[face_id + 1];
      for (i = start_id; i < end_id; i++) {
        vtx = mesh->i_face_vtx_lst[i];
        lambda = 0.0;
        for (coord_id = 0; coord_id < 3; coord_id++)
          lambda += (mesh->vtx_coord[3*vtx + coord_id]
                     - i_face_cog[3*face_id + coord_id])
                     * i_face_norm[3*face_id + coord_id];

        for (coord_id = 0; coord_id < 3; coord_id++) {
          loc_vtx_mvt[vtx*3 + coord_id] -=
            lambda * i_face_norm[3*face_id + coord_id]
                   * UNWARPING_MVT * (i_face_warp[face_id]/maxwarp)
                   * (vtx_tolerance[vtx]/(max_vtxtol*frac));
        }
      }
    }
  }

  if (mesh->vtx_interfaces != NULL) { /* Parallel or periodic treatment */
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         loc_vtx_mvt);
  }

  for (i = 0; i < mesh->n_vertices; i++)
    for (coord_id = 0; coord_id < 3; coord_id++)
      loc_vtx_mvt[3*i + coord_id] = CS_MIN(loc_vtx_mvt[3*i + coord_id],
                                           vtx_tolerance[i]);

  return maxwarp;
}

/*----------------------------------------------------------------------------
 * Compute normals for boundary vertices
 *
 * parameters:
 *   mesh         <--  pointer to a mesh structure
 *   b_face_norm  <--  normals associated with boundary faces
 *   b_vtx_norm   -->  normals associated with boundary vertices
 *----------------------------------------------------------------------------*/

static void
_compute_vtx_normals(cs_mesh_t           *mesh,
                     cs_real_t           *b_face_norm,
                     cs_real_t           *b_vtx_norm)
{
  int coord_id;
  cs_lnum_t i, j;
  cs_real_t norm;

  for (i = 0; i < mesh->n_vertices*3; i++)
    b_vtx_norm[i] = 0.;

  for (i = 0; i < mesh->n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i];
         j < mesh->b_face_vtx_idx[i+1];
         j++) {
      for (coord_id = 0; coord_id < 3; coord_id++) {
        b_vtx_norm[(mesh->b_face_vtx_lst[j])*3 + coord_id]
          += b_face_norm[i*3 + coord_id];
      }
    }
  }

  /* summing upon processors (or periodic vertices) if necessary */

  if (mesh->vtx_interfaces != NULL)
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         b_vtx_norm);

  /* normalizing */
  for (i = 0; i < mesh->n_vertices; i++) {
    norm = sqrt(  b_vtx_norm[i*3    ]*b_vtx_norm[i*3    ]
                + b_vtx_norm[i*3 + 1]*b_vtx_norm[i*3 + 1]
                + b_vtx_norm[i*3 + 2]*b_vtx_norm[i*3 + 2]);

    if (norm > DBL_MIN) {
      b_vtx_norm[i*3    ] /= norm;
      b_vtx_norm[i*3 + 1] /= norm;
      b_vtx_norm[i*3 + 2] /= norm;
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Set fixed vertices flag based on feature angle criterion.
 *
 * <a name="fix_by_feature"></a>
 *
 * Please refer to the
 * <a href="../../theory.pdf#fixbyfeature"><b>specific treatment for boundary faces</b></a>
 * section of the theory guide for more informations.
 *
 * parameters:
 *  \param[in]  mesh           pointer to a cs_mesh_t structure
 *  \param[in]  feature_angle  feature angle (bounded between 0 and 90 degrees)
 *  \param[out] vtx_is_fixed   array to define vertices mobility (1: fixed, 0: free)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_smoother_fix_by_feature(cs_mesh_t   *mesh,
                                cs_real_t    feature_angle,
                                int          vtx_is_fixed[])
{
  cs_lnum_t face, j;

  cs_real_t rnorm_b;
  cs_real_t *face_norm, *vtx_norm;
  cs_real_t *b_face_norm = NULL;
  cs_real_t *b_face_cog = NULL;
  cs_real_t *b_vtx_norm = NULL;
  cs_real_t *_vtx_is_fixed = NULL;

  BFT_MALLOC(_vtx_is_fixed, mesh->n_vertices, cs_real_t);
  BFT_MALLOC(b_vtx_norm, 3*(mesh->n_vertices), cs_real_t);

  cs_mesh_quantities_b_faces(mesh,
                             &(b_face_cog),
                             &(b_face_norm));
  BFT_FREE(b_face_cog);

  for (face = 0; face < mesh->n_b_faces; face++) {
    rnorm_b = sqrt(  b_face_norm[3*face    ]*b_face_norm[3*face    ]
                   + b_face_norm[3*face + 1]*b_face_norm[3*face + 1]
                   + b_face_norm[3*face + 2]*b_face_norm[3*face + 2]);

    b_face_norm[3*face    ] /= rnorm_b;
    b_face_norm[3*face + 1] /= rnorm_b;
    b_face_norm[3*face + 2] /= rnorm_b;
  }

  _compute_vtx_normals(mesh,
                       b_face_norm,
                       b_vtx_norm);

  for (j = 0; j < mesh->n_vertices; j++)
    _vtx_is_fixed[j] = 0;

  for (face = 0; face < mesh->n_b_faces; face++) {
    for (j = mesh->b_face_vtx_idx[face];
         j < mesh->b_face_vtx_idx[face +1];
         j++) {
      face_norm = &b_face_norm[face*3];
      vtx_norm = &b_vtx_norm[(mesh->b_face_vtx_lst[j])*3];

      if (_DOT_PRODUCT_3D(face_norm, vtx_norm) < cos(feature_angle*_PI_/180.0)
          || feature_angle < DBL_MIN)
        _vtx_is_fixed[mesh->b_face_vtx_lst[j]] += 1;
    }
  }

  if (mesh->vtx_interfaces != NULL) {
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         _vtx_is_fixed);
  }

  for (j = 0; j < mesh->n_vertices; j++) {
    if (_vtx_is_fixed[j] > 0.1)
      vtx_is_fixed[j] = 1;
    else
      vtx_is_fixed[j] = 0;
  }

  BFT_FREE(b_face_norm);
  BFT_FREE(b_vtx_norm);

  BFT_FREE(_vtx_is_fixed);
}

/*----------------------------------------------------------------------------*/
/*! \brief Unwarping smoother.
 *
 * <a name="unwarp"></a>
 *
 * Please refer to the
 * <a href="../../theory.pdf#unwarp"><b>unwarping algorithm</b></a>
 * section of the theory guide for more informations.
 *
 * parameters:
 *   \param[in]  mesh          pointer to a cs_mesh_t structure
 *   \param[out] vtx_is_fixed  array to define vertices mobility (1 : fixed, 0 : free)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_smoother_unwarp(cs_mesh_t  *mesh,
                        const int   vtx_is_fixed[])
{
  int face;
  cs_real_t maxwarp, minhist_i, minhist_b, maxhist_i, maxhist_b;
  cs_real_t rnorm_b, rnorm_i;
  bool conv = false;
  int iter = 0;
  int max_iter = UNWARPING_MAX_LOOPS;
  double frac = 0.1;
  double eps = 1.e-4;
  cs_real_t maxwarp_p = 90;
  cs_real_t *vtx_tolerance = NULL;
  cs_real_t *loc_vtx_mvt = NULL;
  cs_real_t *i_face_norm = NULL;
  cs_real_t *i_face_cog = NULL;
  cs_real_t *b_face_norm = NULL;
  cs_real_t *b_face_cog = NULL;
  cs_real_t *b_face_warp = NULL;
  cs_real_t *i_face_warp = NULL;

  if (mesh->have_rotation_perio)
    bft_error(__FILE__, __LINE__, 0,
              "Smoothing in case of periodicity of rotation not yet handled.");

  bft_printf(_("\n Start unwarping algorithm\n\n"));

  BFT_MALLOC(b_face_warp, mesh->n_b_faces, cs_real_t);
  BFT_MALLOC(i_face_warp, mesh->n_i_faces, cs_real_t);

  BFT_MALLOC(vtx_tolerance, mesh->n_vertices, cs_real_t);
  BFT_MALLOC(loc_vtx_mvt, 3*(mesh->n_vertices), cs_real_t);

  while (!conv) {

    cs_mesh_quantities_i_faces(mesh,
                               &(i_face_cog),
                               &(i_face_norm));

    cs_mesh_quantities_b_faces(mesh,
                               &(b_face_cog),
                               &(b_face_norm));

    cs_mesh_quality_compute_warping(mesh,
                                    i_face_norm,
                                    b_face_norm,
                                    i_face_warp,
                                    b_face_warp);

    _get_tolerance(mesh,
                   vtx_tolerance,
                   frac);

    for (face = 0; face < mesh->n_i_faces; face++) {
      rnorm_i = sqrt (  i_face_norm[3*face]*i_face_norm[3*face]
                      + i_face_norm[3*face + 1]*i_face_norm[3*face + 1]
                      + i_face_norm[3*face + 2]*i_face_norm[3*face + 2]);

      i_face_norm[3*face   ] /= rnorm_i;
      i_face_norm[3*face +1] /= rnorm_i;
      i_face_norm[3*face +2] /= rnorm_i;
    }

    for (face = 0; face < mesh->n_b_faces; face++) {
      rnorm_b = sqrt(  b_face_norm[3*face]*b_face_norm[3*face]
                     + b_face_norm[3*face + 1]*b_face_norm[3*face + 1]
                     + b_face_norm[3*face + 2]*b_face_norm[3*face + 2]);

      b_face_norm[3*face   ] /= rnorm_b;
      b_face_norm[3*face +1] /= rnorm_b;
      b_face_norm[3*face +2] /= rnorm_b;
    }

    maxwarp = _unwarping_mvt(mesh,
                             i_face_norm,
                             b_face_norm,
                             i_face_cog,
                             b_face_cog,
                             loc_vtx_mvt,
                             i_face_warp,
                             b_face_warp,
                             vtx_tolerance,
                             frac);

    if (iter == 0) {
     _compute_minmax(mesh->n_i_faces,
                     i_face_warp,
                     &minhist_i,
                     &maxhist_i);
     _compute_minmax(mesh->n_b_faces,
                     b_face_warp,
                     &minhist_b,
                     &maxhist_b);
     bft_printf(_("\n  Histogram of the boundary faces warping"
                  " before unwarping algorithm:\n\n"));

     _histogram(mesh->n_b_faces,
                b_face_warp,
                minhist_b,
                maxhist_b,
                minhist_b,
                maxhist_b);
     bft_printf(_("\n  Histogram of the interior faces warping"
                  " before unwarping algorithm:\n\n"));

     _int_face_histogram(mesh,
                         i_face_warp,
                         minhist_i,
                         maxhist_i,
                         minhist_i,
                         maxhist_i);
    }

    if (maxwarp/maxwarp_p > 1.005) {
      if (iter <= 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("\nUnwarping algorithm failed."));
      else {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("\nUnwarping algorithm stopped at iteration %d"
                     " because it starting to diverge.\n"), iter);
        iter = max_iter +100;
        conv = true;
      }
    }
    if (   ((1 - maxwarp/maxwarp_p) > 0 && (1 - maxwarp/maxwarp_p) < eps)
        || iter == max_iter) {
      conv = true;
      bft_printf(_("\nUnwarping algorithm converged at iteration %d \n"), iter +1);
    }
    maxwarp_p = maxwarp;

    if (iter <= max_iter)
      _move_vertices(mesh,
                     loc_vtx_mvt,
                     vtx_is_fixed);

    BFT_FREE(i_face_norm);
    BFT_FREE(b_face_norm);
    BFT_FREE(i_face_cog);
    BFT_FREE(b_face_cog);
    iter++;
  }

  /* Output quality histograms */

  {
    cs_real_t min_b, max_b, max_i, min_i;

    _compute_minmax(mesh->n_i_faces,
                    i_face_warp,
                    &min_i,
                    &max_i);
    _compute_minmax(mesh->n_b_faces,
                    b_face_warp,
                    &min_b,
                    &max_b);

    bft_printf(_("\n  Histogram of the boundary faces warping"
                 " after unwarping algorithm:\n\n"));

    _histogram(mesh->n_b_faces,
               b_face_warp,
               minhist_b,
               maxhist_b,
               min_b,
               max_b);
    bft_printf(_("\n  Histogram of the interior faces warping"
                 " after unwarping algorithm:\n\n"));

    _int_face_histogram(mesh,
                        i_face_warp,
                        minhist_i,
                        maxhist_i,
                        min_i,
                        max_i);
  }

  BFT_FREE(vtx_tolerance);
  BFT_FREE(loc_vtx_mvt);

  BFT_FREE(i_face_warp);
  BFT_FREE(b_face_warp);

  bft_printf(_("\n End unwarping algorithm\n\n"));
}

/*----------------------------------------------------------------------------*/

/* Delete local macros */

#undef _DOT_PRODUCT_3D

/*----------------------------------------------------------------------------*/

END_C_DECLS
