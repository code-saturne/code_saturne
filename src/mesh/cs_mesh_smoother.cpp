/*============================================================================
 * Mesh smoothing.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_all_to_all.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"

#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_quality.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_smoother.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define UNWARPING_MAX_LOOPS 50
#define UNWARPING_MVT 0.1
#define _PI_ atan(1.0)*4.0

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
_compute_minmax(cs_lnum_t           n_vals,
                const cs_real_t     var[],
                cs_real_t          *min,
                cs_real_t          *max)
{
  cs_real_t  _min = DBL_MAX, _max = -DBL_MAX;

  for (cs_lnum_t i = 0; i < n_vals; i++) {
    _min = cs::min(_min, var[i]);
    _max = cs::max(_max, var[i]);
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

  var_step = cs::abs(var_max - var_min) / n_steps;

  if (cs::abs(var_max - var_min) > 0.) {

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
  cs_gnum_t count[10];
  const int  n_steps = 10;

  assert (sizeof(double) == sizeof(cs_real_t));

  /* Define axis subdivisions */

  for (cs_lnum_t j = 0; j < n_steps; j++)
    count[j] = 0;

  if (cs::abs(max - min) > 0.) {

    cs_real_t step = cs::abs(max - min) / n_steps;

    /* Loop on values */

    for (cs_lnum_t i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      cs_lnum_t j, k;
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
  cs_gnum_t count[8];
  const int  n_steps = 8;

  assert(sizeof(double) == sizeof(cs_real_t));

  /* Define axis subdivisions */

  for (cs_lnum_t j = 0; j < n_steps; j++)
    count[j] = 0;

  if (cs::abs(max - min) > 0.) {

    cs_real_t step = cs::abs(max - min) / n_steps;

    /* Loop on faces */

    for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++) {

      if (mesh->i_face_cells[i][0] >= mesh->n_cells)
        continue;

      /* Associated subdivision */

      cs_lnum_t j, k;
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
  cs_real_t  *vtx_coord = mesh->vtx_coord;
  for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
    if (vtx_is_fixed[i] == 0) {
      for (cs_lnum_t k = 0; k < 3; k++)
        vtx_coord[3*i + k] += vtx_mvt[3*i + k];
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
                     const cs_lnum_t   n_faces,
                     const cs_lnum_t   face_vtx_idx[],
                     const cs_lnum_t   face_vtx_lst[],
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
      vtx_tolerance[vtx_id1] = cs::min(vtx_tolerance[vtx_id1], tolerance);
      vtx_tolerance[vtx_id2] = cs::min(vtx_tolerance[vtx_id2], tolerance);

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
    vtx_tolerance[vtx_id1] = cs::min(vtx_tolerance[vtx_id1], tolerance);
    vtx_tolerance[vtx_id2] = cs::min(vtx_tolerance[vtx_id2], tolerance);

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
  double  *g_vtx_tolerance = nullptr, *send_list = nullptr;
  cs_gnum_t  *send_glist = nullptr;
  cs_gnum_t  n_g_vertices = mesh->n_g_vertices;
  const cs_gnum_t  *io_gnum = mesh->global_vtx_num;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;
  const int  local_rank = cs::max(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;

  /* Define a fvm_io_num_t structure on vertices */

  cs_block_dist_info_t
    bi = cs_block_dist_compute_sizes(local_rank,
                                     n_ranks,
                                     1,
                                     0,
                                     n_g_vertices);

  cs_lnum_t block_size = bi.block_size;

  cs_all_to_all_t
    *d = cs_all_to_all_create_from_block(n_vertices,
                                         0, /* flags */
                                         io_gnum,
                                         bi,
                                         mpi_comm);

  /* Send the global numbering for each vertex */

  cs_gnum_t *recv_glist = cs_all_to_all_copy_array(d,
                                                   1,
                                                   false, /* reverse */
                                                   io_gnum);

  /* Send the vertex tolerance for each vertex */

  double *recv_list = cs_all_to_all_copy_array(d,
                                               1,
                                               false, /* reverse */
                                               vtx_tolerance);

  cs_lnum_t n_recv = cs_all_to_all_n_elts_dest(d);

  /* Define the global tolerance array */

  CS_MALLOC(g_vtx_tolerance, block_size, double);

  for (i = 0; i < block_size; i++)
    g_vtx_tolerance[i] = DBL_MAX;

  first_vtx_gnum = block_size * local_rank + 1;

  for (i = 0; i < n_recv; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    g_vtx_tolerance[vtx_id] = cs::min(g_vtx_tolerance[vtx_id], recv_list[i]);
  }

  /* Replace local vertex tolerance by the new computed global tolerance */

  for (i = 0; i < n_recv; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    recv_list[i] = g_vtx_tolerance[vtx_id];
  }

  cs_all_to_all_copy_array(d,
                           1,
                           true, /* reverse */
                           recv_list,
                           vtx_tolerance);

  /* Free memory */

  cs_all_to_all_destroy(&d);

  CS_FREE(recv_glist);
  CS_FREE(send_glist);
  CS_FREE(send_list);
  CS_FREE(recv_list);
  CS_FREE(g_vtx_tolerance);
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
 *   i_face_cog          <--  center of gravity of interior faces
 *   b_face_cog          <--  center of gravity of border faces
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
               cs_nreal_3_t         *i_face_u_norm,
               cs_nreal_3_t         *b_face_u_norm,
               cs_real_3_t          *i_face_cog,
               cs_real_3_t          *b_face_cog,
               cs_real_t            *loc_vtx_mvt,
               cs_real_t            *i_face_warp,
               cs_real_t            *b_face_warp,
               cs_real_t            *vtx_tolerance,
               double                frac)
{
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_vertices = mesh->n_vertices;

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_real_t *vtx_coord = mesh->vtx_coord;

  cs_real_t max_vtxtol = 0.;
  cs_real_t maxwarp = 0.;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    if (maxwarp < i_face_warp[face_id])
      maxwarp = i_face_warp[face_id];
  }
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (maxwarp < b_face_warp[face_id])
      maxwarp = b_face_warp[face_id];
  }

  for (cs_lnum_t i = 0; i < n_vertices*3; i++)
    loc_vtx_mvt[i] = 0.0;
  for (cs_lnum_t i = 0; i < n_vertices; i++)
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

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t start_id = b_face_vtx_idx[face_id];
    cs_lnum_t end_id = b_face_vtx_idx[face_id + 1];
    for (cs_lnum_t i = start_id; i < end_id; i++) {
      cs_lnum_t vtx = b_face_vtx[i];
      cs_real_t lambda = 0.0;
      for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++)
        lambda +=  (vtx_coord[3*vtx + coord_id]
                    - b_face_cog[face_id][coord_id])
                    * b_face_u_norm[face_id][coord_id];

      for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++) {
        loc_vtx_mvt[vtx*3 + coord_id] -=
          lambda * b_face_u_norm[face_id][coord_id]
                 * UNWARPING_MVT * (b_face_warp[face_id]/maxwarp)
                 * (vtx_tolerance[vtx]/(max_vtxtol*frac));
      }
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    if (i_face_cells[face_id][0] < n_cells) {
      cs_lnum_t start_id = i_face_vtx_idx[face_id];
      cs_lnum_t end_id = i_face_vtx_idx[face_id + 1];
      for (cs_lnum_t i = start_id; i < end_id; i++) {
        cs_lnum_t vtx = i_face_vtx[i];
        cs_real_t lambda = 0.0;
        for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++)
          lambda += (vtx_coord[3*vtx + coord_id]
                     - i_face_cog[face_id][coord_id])
                     * i_face_u_norm[face_id][coord_id];

        for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++) {
          loc_vtx_mvt[vtx*3 + coord_id] -=
            lambda * i_face_u_norm[face_id][coord_id]
                   * UNWARPING_MVT * (i_face_warp[face_id]/maxwarp)
                   * (vtx_tolerance[vtx]/(max_vtxtol*frac));
        }
      }
    }
  }

  if (mesh->vtx_interfaces != nullptr) { /* Parallel or periodic treatment */
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         loc_vtx_mvt);
  }

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++)
      loc_vtx_mvt[3*i + coord_id] = cs::min(loc_vtx_mvt[3*i + coord_id],
                                            vtx_tolerance[i]);

  return maxwarp;
}

/*----------------------------------------------------------------------------
 * Compute normals for boundary vertices
 *
 * parameters:
 *   mesh           <--  pointer to a mesh structure
 *   b_face_u_norm  <--  unit normals associated with boundary faces
 *   b_vtx_norm     -->  normals associated with boundary vertices
 *----------------------------------------------------------------------------*/

static void
_compute_vtx_normals(cs_mesh_t       *mesh,
                     cs_nreal_3_t    *b_face_u_norm,
                     cs_real_3_t     *b_vtx_norm)
{
  for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
    for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++)
      b_vtx_norm[i][coord_id] = 0.;
  }

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    for (cs_lnum_t j = mesh->b_face_vtx_idx[i];
         j < mesh->b_face_vtx_idx[i+1];
         j++) {
      for (cs_lnum_t coord_id = 0; coord_id < 3; coord_id++) {
        b_vtx_norm[mesh->b_face_vtx_lst[j]][coord_id]
          += b_face_u_norm[i][coord_id];
      }
    }
  }

  /* summing upon processors (or periodic vertices) if necessary */

  if (mesh->vtx_interfaces != nullptr)
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         b_vtx_norm);

  /* normalizing */
  for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
    cs_math_3_normalize(b_vtx_norm[i], b_vtx_norm[i]);
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
 * This function locks a vertex in place if one of its feature angles is
 * less than the maximum feature angle (in degrees) defined by the user.
 *
 * Vertex normals are based on the average of the normals of the adjacent
 * boundary faces.
 * The feature angle between a vertex and one of its adjacent faces is defined
 * as the angle between the vertex normal and the face normal.
 *
 * Please refer to the
 * <a href="../../theory.pdf#fixbyfeature"><b>specific treatment for boundary faces</b></a>
 * section of the theory guide for more information.
 *
 *  \param[in]  mesh           pointer to a cs_mesh_t structure
 *  \param[in]  feature_angle  feature angle (between 0 and 90 degrees)
 *  \param[out] vtx_is_fixed   array to define vertices mobility
 *                             (1: fixed, 0: free)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_smoother_fix_by_feature(cs_mesh_t   *mesh,
                                cs_real_t    feature_angle,
                                int          vtx_is_fixed[])
{
  cs_real_3_t *b_vtx_norm = nullptr;
  cs_real_t *_vtx_is_fixed = nullptr;

  CS_MALLOC(_vtx_is_fixed, mesh->n_vertices, cs_real_t);
  CS_MALLOC(b_vtx_norm, mesh->n_vertices, cs_real_3_t);

  cs_real_3_t  *b_face_cog = nullptr;
  cs_nreal_3_t *b_face_u_norm = nullptr;
  CS_MALLOC_HD(b_face_cog, mesh->n_b_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(b_face_u_norm, mesh->n_b_faces, cs_nreal_3_t, cs_alloc_mode);

  cs_mesh_quantities_compute_face_cog_un
    (mesh->n_b_faces,
     reinterpret_cast<const cs_real_3_t *>(mesh->vtx_coord),
     mesh->b_face_vtx_idx,
     mesh->b_face_vtx_lst,
     b_face_cog,
     b_face_u_norm);

  CS_FREE(b_face_cog);

  _compute_vtx_normals(mesh, b_face_u_norm, b_vtx_norm);

  for (cs_lnum_t j = 0; j < mesh->n_vertices; j++)
    _vtx_is_fixed[j] = 0;

  cs_real_t pi_o_180 = cs_math_pi / 180.;

  for (cs_lnum_t face = 0; face < mesh->n_b_faces; face++) {
    for (cs_lnum_t j = mesh->b_face_vtx_idx[face];
         j < mesh->b_face_vtx_idx[face +1];
         j++) {
      const cs_nreal_t *face_u_norm = b_face_u_norm[face];
      const cs_real_t *vtx_norm = b_vtx_norm[mesh->b_face_vtx_lst[j]];

      if (  (  cs_math_3_dot_product(face_u_norm, vtx_norm)
             < cos(feature_angle*pi_o_180))
          || feature_angle < DBL_MIN)
        _vtx_is_fixed[mesh->b_face_vtx_lst[j]] += 1;
    }
  }

  if (mesh->vtx_interfaces != nullptr) {
    cs_interface_set_sum(mesh->vtx_interfaces,
                         mesh->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         _vtx_is_fixed);
  }

  for (cs_lnum_t j = 0; j < mesh->n_vertices; j++) {
    if (_vtx_is_fixed[j] > 0.1)
      vtx_is_fixed[j] = 1;
    else
      vtx_is_fixed[j] = 0;
  }

  CS_FREE(b_face_u_norm);
  CS_FREE(b_vtx_norm);

  CS_FREE(_vtx_is_fixed);
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
 * \param[in]  mesh          pointer to a cs_mesh_t structure
 * \param[out] vtx_is_fixed  array to define vertices mobility
 *                           (1 : fixed, 0 : free)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_smoother_unwarp(cs_mesh_t  *mesh,
                        const int   vtx_is_fixed[])
{
  cs_real_t maxwarp, minhist_i, minhist_b, maxhist_i, maxhist_b;
  bool conv = false;
  int iter = 0;
  int max_iter = UNWARPING_MAX_LOOPS;
  double frac = 0.1;
  double eps = 1.e-4;
  cs_real_t maxwarp_p = 90;
  cs_real_t *vtx_tolerance = nullptr;
  cs_real_t *loc_vtx_mvt = nullptr;
  cs_real_t *b_face_warp = nullptr;
  cs_real_t *i_face_warp = nullptr;

  if (mesh->have_rotation_perio)
    bft_error(__FILE__, __LINE__, 0,
              "Smoothing in case of periodicity of rotation not yet handled.");

  bft_printf(_("\n Start unwarping algorithm\n\n"));

  CS_MALLOC(b_face_warp, mesh->n_b_faces, cs_real_t);
  CS_MALLOC(i_face_warp, mesh->n_i_faces, cs_real_t);

  CS_MALLOC(vtx_tolerance, mesh->n_vertices, cs_real_t);
  CS_MALLOC(loc_vtx_mvt, 3*(mesh->n_vertices), cs_real_t);

  cs_real_3_t  *i_face_cog = nullptr, *b_face_cog = nullptr;
  cs_nreal_3_t *i_face_u_norm = nullptr, *b_face_u_norm = nullptr;
  CS_MALLOC_HD(i_face_cog, mesh->n_i_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(i_face_u_norm, mesh->n_i_faces, cs_nreal_3_t, cs_alloc_mode);
  CS_MALLOC_HD(b_face_cog, mesh->n_b_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(b_face_u_norm, mesh->n_b_faces, cs_nreal_3_t, cs_alloc_mode);

  while (!conv) {

    cs_mesh_quantities_compute_face_cog_un
      (mesh->n_i_faces,
       reinterpret_cast<const cs_real_3_t *>(mesh->vtx_coord),
       mesh->i_face_vtx_idx,
       mesh->i_face_vtx_lst,
       i_face_cog,
       i_face_u_norm);

    cs_mesh_quantities_compute_face_cog_un
      (mesh->n_b_faces,
       reinterpret_cast<const cs_real_3_t *>(mesh->vtx_coord),
       mesh->b_face_vtx_idx,
       mesh->b_face_vtx_lst,
       b_face_cog,
       b_face_u_norm);

    cs_mesh_quality_compute_warping(mesh,
                                    i_face_u_norm,
                                    b_face_u_norm,
                                    i_face_warp,
                                    b_face_warp);

    _get_tolerance(mesh, vtx_tolerance, frac);

    maxwarp = _unwarping_mvt(mesh,
                             i_face_u_norm,
                             b_face_u_norm,
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
      bft_printf(_("\nUnwarping algorithm converged at iteration %d \n"),
                 iter +1);
    }
    maxwarp_p = maxwarp;

    if (iter <= max_iter)
      _move_vertices(mesh,
                     loc_vtx_mvt,
                     vtx_is_fixed);

    iter++;
  }

  CS_FREE(i_face_u_norm);
  CS_FREE(b_face_u_norm);
  CS_FREE(i_face_cog);
  CS_FREE(b_face_cog);

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

  CS_FREE(vtx_tolerance);
  CS_FREE(loc_vtx_mvt);

  CS_FREE(i_face_warp);
  CS_FREE(b_face_warp);

  bft_printf(_("\n End unwarping algorithm\n\n"));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
