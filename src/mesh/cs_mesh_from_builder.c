/*============================================================================
 * \file Define cs_mesh_t fields from cs_mesh_builder_t fields.
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_all_to_all.h"
#include "cs_block_dist.h"
#include "cs_block_to_part.h"
#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_order.h"
#include "cs_partition.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_from_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local Type definitions
 *============================================================================*/

typedef double  _vtx_coords_t[3];

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in parallel mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * parameters:
 *   mesh              <-> pointer to mesh structure
 *   n_faces           <-- number of local faces
 *   face_ifs          <-- parallel and periodic faces interfaces set
 *   face_cell         <-- local face -> cell connectivity
 *   face_vertices_idx <-- local face -> vertices index
 *   face_type         --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_g(cs_mesh_t                 *mesh,
             cs_lnum_t                  n_faces,
             const cs_interface_set_t  *face_ifs,
             const cs_lnum_2_t          face_cell[],
             const cs_lnum_t            face_vertices_idx[],
             char                       face_type[])
{
  cs_lnum_t i;
  int j;

  const int n_interfaces = cs_interface_set_size(face_ifs);

  /* Mark base interior faces */

  for (i = 0; i < n_faces; i++) {
    if (face_cell[i][0] > -1 && face_cell[i][1] > -1)
      face_type[i] = '\0';
    else if (face_cell[i][0] > -1)
      face_type[i] = '\1';
    else if (face_cell[i][1] > -1)
      face_type[i] = '\2';
    else {
      face_type[i] = '\3';
    }
  }

  /* Also mark parallel and periodic faces as interior */

  for (j = 0; j < n_interfaces; j++) {

    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, j);
    cs_lnum_t face_if_size = cs_interface_size(face_if);
    const cs_lnum_t *loc_id = cs_interface_get_elt_ids(face_if);

    for (i = 0; i < face_if_size; i++)
      face_type[loc_id[i]] = '\0';

  }

  /* Now count faces of each type */

  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else {
      mesh->n_b_faces += 1;
      mesh->b_face_vtx_connect_size += n_f_vertices;
    }
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in serial mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * parameters:
 *   mesh               <-> pointer to mesh structure
 *   n_faces            <-- number of local faces
 *   n_periodic_couples <-- number of periodic couples associated with
 *                          each periodic list
 *   periodic_couples   <-- array indicating periodic couples (using
 *                          global numberings) for each list
 *   face_cell          <-- local face -> cell connectivity
 *   face_vertices_idx  <-- local face -> vertices index
 *   face_type          --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_l(cs_mesh_t                  *mesh,
             cs_lnum_t                   n_faces,
             const cs_lnum_t             n_periodic_couples[],
             const cs_gnum_t      *const periodic_couples[],
             const cs_lnum_2_t           face_cell[],
             const cs_lnum_t             face_vertices_idx[],
             char                        face_type[])
{
  cs_lnum_t i;
  int j;

  /* Mark base interior faces */

  for (i = 0; i < n_faces; i++) {
    if (face_cell[i][0] > -1 && face_cell[i][1] > -1)
      face_type[i] = '\0';
    else if (face_cell[i][0] > -1)
      face_type[i] = '\1';
    else if (face_cell[i][1] > -1)
      face_type[i] = '\2';
    else
      face_type[i] = '\3';
  }

  /* Also mark parallel and periodic faces as interior */

  for (i = 0; i < mesh->n_init_perio; i++) {

    const cs_gnum_t *p_couples = periodic_couples[i];

    for (j = 0; j < n_periodic_couples[i]; j++) {
      face_type[p_couples[j*2] - 1] = '\0';
      face_type[p_couples[j*2 + 1] - 1] = '\0';
    }

  }

  /* Now count faces of each type */

  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else {
      mesh->n_b_faces += 1;
      mesh->b_face_vtx_connect_size += n_f_vertices;
    }
  }

  mesh->n_g_i_faces = mesh->n_i_faces;
  mesh->n_g_b_faces = mesh->n_b_faces;
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> cell connectivity using a common
 * face -> cell connectivity and a face type marker.
 *
 * At this stage, isolated faces, if present, are considered to be
 * boundary faces, as they may participate in future mesh joining
 * operations. Their matching cell number will be set to -1.
 * Remaining isolated faces should be removed before completing
 * the mesh structure.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh      <-> pointer to mesh structure
 *   n_faces   <-- number of local faces
 *   face_cell <-- local face -> cell connectivity
 *   face_type <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_cell(cs_mesh_t         *mesh,
                   cs_lnum_t          n_faces,
                   const cs_lnum_2_t  face_cell[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);
  BFT_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_lnum_t);

  /* Now copy face -> cell connectivity */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0') {
      mesh->i_face_cells[n_i_faces][0] = face_cell[i][0];
      mesh->i_face_cells[n_i_faces][1] = face_cell[i][1];
      n_i_faces++;
    }

    else if (face_type[i] == '\1') {
      mesh->b_face_cells[n_b_faces] = face_cell[i][0];
      n_b_faces++;
    }

    else if (face_type[i] == '\2') {
      mesh->b_face_cells[n_b_faces] = face_cell[i][1];
      n_b_faces++;
    }

    else if (face_type[i] == '\3') {
      mesh->b_face_cells[n_b_faces] = -1;
      mesh->n_g_free_faces += 1;
      n_b_faces++;
    }
  }
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> vertices connectivity using a common
 * face -> vertices connectivity and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh              <-> pointer to mesh structure
 *   n_faces           <-- number of local faces
 *   face_vertices_idx <-- local face -> vertices index
 *   face_vertices     <-- local face -> vertices connectivity
 *   face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_vertices(cs_mesh_t         *mesh,
                       cs_lnum_t          n_faces,
                       const cs_lnum_t    face_vertices_idx[],
                       const cs_lnum_t    face_vertices[],
                       const char         face_type[])
{
  cs_lnum_t i;
  size_t j;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate and initialize */

  BFT_MALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces+1, cs_int_t);
  BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);

  mesh->i_face_vtx_idx[0] = 0;

  BFT_MALLOC(mesh->b_face_vtx_idx, mesh->n_b_faces+1, cs_int_t);
  mesh->b_face_vtx_idx[0] = 0;

  if (mesh->n_b_faces > 0)
    BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);

  /* Now copy face -> vertices connectivity */

  for (i = 0; i < n_faces; i++) {

    size_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    const cs_lnum_t *_face_vtx = face_vertices + face_vertices_idx[i];

    if (face_type[i] == '\0') {
      cs_lnum_t *_i_face_vtx =   mesh->i_face_vtx_lst
                               + mesh->i_face_vtx_idx[n_i_faces];
      for (j = 0; j < n_f_vertices; j++)
        _i_face_vtx[j] = _face_vtx[j] - 1;
      mesh->i_face_vtx_idx[n_i_faces + 1] =   mesh->i_face_vtx_idx[n_i_faces]
                                            + n_f_vertices;
      n_i_faces++;
    }

    else if (face_type[i] == '\1' || face_type[i] == '\3') {
      cs_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                               + mesh->b_face_vtx_idx[n_b_faces];
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[j] - 1;
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

    else if (face_type[i] == '\2') {
      cs_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                               + mesh->b_face_vtx_idx[n_b_faces];
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[n_f_vertices - j - 1] - 1;
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> global numberings using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh            <-> pointer to mesh structure
 *   n_faces         <-- number of local faces
 *   global_face_num <-- global face numbers
 *   face_type       <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gnum(cs_mesh_t         *mesh,
                   cs_lnum_t          n_faces,
                   const cs_gnum_t    global_face_num[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  cs_lnum_t *global_i_face = NULL;
  cs_lnum_t *global_b_face = NULL;

  fvm_io_num_t *tmp_face_num = NULL;

  /* Allocate arrays (including temporary arrays) */

  BFT_MALLOC(mesh->global_i_face_num, mesh->n_i_faces, cs_gnum_t);
  BFT_MALLOC(mesh->global_b_face_num, mesh->n_b_faces, cs_gnum_t);

  BFT_MALLOC(global_i_face, mesh->n_i_faces, cs_lnum_t);
  BFT_MALLOC(global_b_face, mesh->n_b_faces, cs_lnum_t);

  /* Now build internal and boundary face lists */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0')
      global_i_face[n_i_faces++] = i+1;

    else
      global_b_face[n_b_faces++] = i+1;

  }

  /* Build an I/O numbering on internal faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_i_face,
                                   global_face_num,
                                   n_i_faces,
                                   0);

  memcpy(mesh->global_i_face_num,
         fvm_io_num_get_global_num(tmp_face_num),
         n_i_faces*sizeof(cs_gnum_t));

  mesh->n_g_i_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (cs_lnum_t)n_i_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Build an I/O numbering on boundary faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_b_face,
                                   global_face_num,
                                   n_b_faces,
                                   0);

  if (n_b_faces > 0)
    memcpy(mesh->global_b_face_num,
           fvm_io_num_get_global_num(tmp_face_num),
           n_b_faces*sizeof(cs_gnum_t));

  mesh->n_g_b_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (cs_lnum_t)n_b_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Free remaining temporary arrays */

  BFT_FREE(global_i_face);
  BFT_FREE(global_b_face);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> group class id using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh       <-> pointer to mesh structure
 *   n_faces    <-- number of local faces
 *   face_gc_id <-- local face group class id
 *   face_type  <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gc_id(cs_mesh_t        *mesh,
                   cs_lnum_t          n_faces,
                   const cs_lnum_t    face_gc_id[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_family, mesh->n_i_faces, int);
  BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, int);

  /* Now copy face group class (family) id */

  for (i = 0; i < n_faces; i++) {

    assert(face_gc_id[i] > -1 && face_gc_id[i] <= mesh->n_families);

    if (face_type[i] == '\0')
      mesh->i_face_family[n_i_faces++] = face_gc_id[i];

    else
      mesh->b_face_family[n_b_faces++] = face_gc_id[i];

  }
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> refinement generation using a common
 * face level and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh       <-> pointer to mesh structure
 *   n_faces    <-- number of local faces
 *   face_r_gen <-- local face level
 *   face_type  <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_r_gen(cs_mesh_t        *mesh,
                    cs_lnum_t         n_faces,
                    const char        face_r_gen[],
                    const char        face_type[])
{
  size_t n_i_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_r_gen, mesh->n_i_faces, char);

  /* Now copy face group class (family) id */

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      mesh->i_face_r_gen[n_i_faces++] = face_r_gen[i];
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Renumber face interface references from mixed faces to interior faces.
 *
 * parameters:
 *   face_ifs          <-> parallel and periodic faces interfaces set
 *   n_faces           <-- number of local faces
 *   face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_face_ifs_to_interior(cs_interface_set_t  *face_ifs,
                      cs_lnum_t            n_faces,
                      const char           face_type[])
{
  cs_lnum_t i;

  cs_lnum_t   i_face_count = 0;
  cs_lnum_t  *i_face_id = NULL;

  /* Build face renumbering */

  BFT_MALLOC(i_face_id, n_faces, cs_lnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_id[i] = i_face_count++;
    else
      i_face_id[i] = -1;
  }

  cs_interface_set_renumber(face_ifs, i_face_id);

  BFT_FREE(i_face_id);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information for mesh builder when
 * running in serial mode.
 *
 * Arrays are simply transferred from the mesh reader to the builder and
 * renumbered
 *
 * parameters:
 *   mb            <-> pointer to mesh builder structure
 *   n_init_perio  <-- number of initial periodicities
 *   n_faces       <-- number of local faces
 *   face_type     <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_periodic_faces_l(cs_mesh_builder_t  *mb,
                          int                 n_init_perio,
                          cs_lnum_t           n_faces,
                          const char          face_type[])
{
  int i;

  cs_gnum_t   next_face_num = 1;
  cs_gnum_t  *i_face_num = NULL;

  /* Transfer arrays from reader to builder, then renumber couples */

  assert(mb != NULL);

  mb->n_perio = n_init_perio;

  /* Build face renumbering */

  BFT_MALLOC(i_face_num, n_faces, cs_gnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_num[i] = next_face_num++;
    else
      i_face_num[i] = 0;
  }

  /* Apply new numbering */

  for (i = 0; i < n_init_perio; i++) {

    size_t j;
    cs_gnum_t *p_couples = mb->per_face_couples[i];
    const size_t n_vals = mb->n_per_face_couples[i] * 2;

    for (j = 0; j < n_vals; j++) {
      p_couples[j] = i_face_num[p_couples[j] - 1];
      assert(p_couples[j] > 0);
    }
  }

  BFT_FREE(i_face_num);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute free (isolated) face centers using minimal local data.
 *
 * parameters:
 *   n_f_faces     <-- number of faces
 *   f_face_ids    <-- list of free faces
 *   face_vtx_idx  <-- face -> vertices connectivity index
 *   face_vtx      <-- face -> vertices connectivity
 *   vtx_coord     <-- vertex coordinates
 *   f_face_center --> free face centers
 *----------------------------------------------------------------------------*/

static void
_f_face_center(cs_lnum_t         n_f_faces,
               cs_lnum_t         f_face_ids[],
               const cs_lnum_t   face_vtx_idx[],
               const cs_lnum_t   face_vtx[],
               const cs_real_t   vtx_coord[],
               cs_coord_t        f_face_center[])
{
  cs_lnum_t i, j, k;
  cs_lnum_t vtx_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3];

  cs_lnum_t n_max_face_vertices = 0;

  _vtx_coords_t *face_vtx_coord = NULL;

  const double surf_epsilon = 1e-24;

  assert(face_vtx_idx[0] == 0);

  for (i = 0; i < n_f_faces; i++) {
    for (j = 0; j < 3; j++)
      f_face_center[i*3 + j] = 0.0;
  }

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (k = 0; k < n_f_faces; k++) {
    cs_lnum_t face_id = f_face_ids[k];
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices, _vtx_coords_t);

  /* Loop on each face */

  for (k = 0; k < n_f_faces; k++) {

    cs_lnum_t tri_id;

    /* Initialization */

    cs_lnum_t face_id = f_face_ids[k];
    cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
    cs_coord_t face_surface = 0.0;
    cs_coord_t *face_center = f_face_center + (k*3);

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id];
    end_id = face_vtx_idx[face_id + 1];

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      cs_lnum_t shift = 3 * (face_vtx[vtx_id]);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    /* Compute the barycenter of the face vertices */

    for (i = 0; i < 3; i++) {
      vtx_cog[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        vtx_cog[i] += face_vtx_coord[vtx_id][i];
      vtx_cog[i] /= n_face_vertices;
    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycenter) */

    for (i = 0; i < 3; i++) {
      ref_normal[i] = 0.;
      face_center[i] = 0.0;
    }

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      cs_coord_t tri_surface;
      cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

      cs_lnum_t id0 = tri_id;
      cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

      /* Normal for each triangle */

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
        vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
      }

      tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
      tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
      tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

      if (tri_id == 0) {
        for (i = 0; i < 3; i++)
          ref_normal[i] = tri_normal[i];
      }

      /* Center of gravity for a triangle */

      for (i = 0; i < 3; i++) {
        tri_center[i] = (  vtx_cog[i]
                         + face_vtx_coord[id0][i]
                         + face_vtx_coord[id1][i]) / 3.0;
      }

      tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                         + tri_normal[1]*tri_normal[1]
                         + tri_normal[2]*tri_normal[2]) * 0.5;

      if ((  tri_normal[0]*ref_normal[0]
           + tri_normal[1]*ref_normal[1]
           + tri_normal[2]*ref_normal[2]) < 0.0)
        tri_surface *= -1.0;

      /* Now compute contribution to face center and surface */

      face_surface += tri_surface;

      for (i = 0; i < 3; i++) {
        face_center[i] += tri_surface * tri_center[i];
        unweighted_center[i] = tri_center[i];
      }

    } /* End of loop  on triangles of the face */

    if (face_surface > surf_epsilon) {
      for (i = 0; i < 3; i++)
        face_center[i] /= face_surface;
    }
    else {
      face_surface = surf_epsilon;
      for (i = 0; i < 3; i++)
        face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
    }
  } /* End of loop on faces */

  BFT_FREE(face_vtx_coord);
}

/*----------------------------------------------------------------------------
 * Compute face centers using block data.
 *
 * parameters:
 *   mb          <-- pointer to mesh builder helper structure
 *   n_f_faces   <-- local number of free faces
 *   f_face_ids  <-- free face ids
 *   face_center --> cell centers array
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_precompute_free_face_center(const cs_mesh_builder_t  *mb,
                             cs_lnum_t                 n_f_faces,
                             cs_lnum_t                 f_face_ids[],
                             cs_coord_t                f_face_center[],
                             MPI_Comm                  comm)
{
  int n_ranks = 0;

  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  cs_lnum_t _n_faces = 0;

  cs_lnum_t *_face_vertices = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  _n_faces = mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0];

  /* Distribute vertices */
  /*---------------------*/

  size_t _n_vertices = 0;
  cs_gnum_t *_vtx_num = NULL;

  cs_order_single_gnum(mb->face_vertices_idx[_n_faces],
                       1, /* base */
                       mb->face_vertices,
                       &_n_vertices,
                       &_vtx_num);

  cs_all_to_all_t *d
    = cs_all_to_all_create_from_block(_n_vertices,
                                      CS_ALL_TO_ALL_USE_DEST_ID,
                                      _vtx_num,
                                      mb->vertex_bi,
                                      comm);

  cs_real_t *_vtx_coord
    = cs_all_to_all_copy_array(d,
                               real_type,
                               3,
                               true, /* reverse */
                               mb->vertex_coords,
                               NULL);

  cs_all_to_all_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, mb->face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(mb->face_vertices_idx[_n_faces],
                                   1,
                                   _n_vertices,
                                   true,
                                   _vtx_num,
                                   mb->face_vertices,
                                   _face_vertices);

  _f_face_center(n_f_faces,
                 f_face_ids,
                 mb->face_vertices_idx,
                 _face_vertices,
                 _vtx_coord,
                 f_face_center);

  BFT_FREE(_vtx_coord);
  BFT_FREE(_vtx_num);
  BFT_FREE(_face_vertices);
}

/*----------------------------------------------------------------------------
 * Compute default face destination rank array in case of isolated faces.
 *
 * parameters:
 *   mb           <-- pointer to mesh builder helper structure
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *  default rank array for faces (>= for isolated faces)
 *----------------------------------------------------------------------------*/

static int *
_default_face_rank(const cs_mesh_builder_t  *mb,
                   MPI_Comm                  comm)
{
  cs_lnum_t i;
  cs_block_dist_info_t free_face_bi;

  int n_ranks = 0, rank_id = -1;

  cs_lnum_t _n_faces = 0, n_free_faces = 0;
  cs_gnum_t _n_g_free_faces = 0, n_g_free_faces = 0;

  cs_lnum_t *free_face_ids = NULL;
  cs_coord_t *free_face_centers = NULL;

  fvm_io_num_t *free_face_io_num = NULL;
  const cs_gnum_t *free_face_num = NULL;

  int *default_rank = NULL;

  /* Initialization */

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  /* Count number of isolated faces */

  _n_faces = mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0];
  n_free_faces = 0;

  for (i = 0; i < _n_faces; i++) {
    if (mb->face_cells[2*i] == 0 && mb->face_cells[2*i+1] == 0)
      n_free_faces += 1;
  }

  _n_g_free_faces = n_free_faces;
  MPI_Allreduce(&_n_g_free_faces, &n_g_free_faces, 1,
                CS_MPI_GNUM, MPI_SUM, comm);

  /* Return if we do not have isolated faces */

  if (n_g_free_faces == 0)
    return NULL;

  /* Initialize rank info */

  MPI_Comm_size(comm, &n_ranks);
  MPI_Comm_size(comm, &rank_id);
  free_face_bi = cs_block_dist_compute_sizes(rank_id,
                                             n_ranks,
                                             0,
                                             0,
                                             n_g_free_faces);

  /* Define distribution of isolated faces based on sfc;
   *
   *  As those faces are not connected, the main objective of this function
   *  is to ensure some measure of load balancing. */

  BFT_MALLOC(default_rank, _n_faces, int);
  for (i = 0; i < _n_faces; i++)
    default_rank[i] = -1;

  BFT_MALLOC(free_face_ids, n_free_faces, cs_lnum_t);
  BFT_MALLOC(free_face_centers, n_free_faces*3, cs_coord_t);

  n_free_faces = 0;
  for (i = 0; i < _n_faces; i++) {
    if (mb->face_cells[2*i] == 0 && mb->face_cells[2*i+1] == 0)
      free_face_ids[n_free_faces++] = i;
  }

  _precompute_free_face_center(mb,
                               n_free_faces,
                               free_face_ids,
                               free_face_centers,
                               comm);

  free_face_io_num = fvm_io_num_create_from_sfc(free_face_centers,
                                                3,
                                                n_free_faces,
                                                FVM_IO_NUM_SFC_MORTON_BOX);

  BFT_FREE(free_face_centers);

  free_face_num = fvm_io_num_get_global_num(free_face_io_num);

  /* Determine rank based on global numbering with SFC ordering */
  for (i = 0; i < n_free_faces; i++) {
    default_rank[free_face_ids[i]]
      =    ((free_face_num[i] - 1) / free_face_bi.block_size)
         * free_face_bi.rank_step;
  }

  free_face_io_num = fvm_io_num_destroy(free_face_io_num);
  BFT_FREE(free_face_ids);

  return default_rank;
}

/*----------------------------------------------------------------------------
 * Organize data read by blocks in parallel and build most mesh structures.
 *
 * parameters:
 *   mesh    <-> pointer to mesh structure
 *   mb      <-> pointer to mesh builder structure
 *   comm    <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_decompose_data_g(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb,
                  MPI_Comm            comm)
{
  int n_ranks = 0;

  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  cs_lnum_t _n_faces = 0;
  cs_gnum_t *_face_num = NULL;
  cs_gnum_t *_face_gcells = NULL;
  cs_gnum_t *_face_gvertices = NULL;

  cs_lnum_2_t *_face_cells = NULL;
  cs_lnum_t *_face_gc_id = NULL;
  char *_face_r_gen = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  int *_periodicity_num = NULL;

  int  *default_face_rank = NULL;
  char *face_type = NULL;
  cs_interface_set_t *face_ifs = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  /* Different handling of cells depending on whether decomposition
     data is available or not. */

  if (mb->have_cell_rank == true) {

    cs_lnum_t n_block_ents = 0;
    if (mb->cell_bi.gnum_range[1] > mb->cell_bi.gnum_range[0])
      n_block_ents = (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]);

    cs_all_to_all_t *d
      = cs_all_to_all_create(n_block_ents,
                             CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,   /* flags */
                             NULL,
                             mb->cell_rank,
                             comm);

    mesh->n_cells = cs_all_to_all_n_elts_dest(d);

    cs_gnum_t *b_global_num;
    BFT_MALLOC(b_global_num, n_block_ents, cs_gnum_t);

    mesh->cell_family = cs_all_to_all_copy_array(d,
                                                 CS_LNUM_TYPE,
                                                 1,
                                                 false, /* reverse */
                                                 mb->cell_gc_id,
                                                 NULL);

    BFT_FREE(mb->cell_gc_id);

    cs_gnum_t  gnum_shift = mb->cell_bi.gnum_range[0];
    for (cs_lnum_t i = 0; i < n_block_ents; i++)
      b_global_num[i] = (cs_gnum_t)i + gnum_shift;

    mesh->global_cell_num = cs_all_to_all_copy_array(d,
                                                     CS_GNUM_TYPE,
                                                     1,
                                                     false, /* reverse */
                                                     b_global_num,
                                                     NULL);

    BFT_FREE(b_global_num);

    cs_all_to_all_destroy(&d);

  }
  else {

    mesh->n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

    BFT_MALLOC(mesh->global_cell_num, mesh->n_cells, cs_gnum_t);

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      mesh->global_cell_num[i] = mb->cell_bi.gnum_range[0] + i;

    mesh->cell_family = mb->cell_gc_id;
    mb->cell_gc_id = NULL;
  }

  if (mesh->n_cells == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Number of cells on rank %d is zero.\n"
                "(number of cells / number of processes ratio too low)."),
              (int)cs_glob_rank_id);

  mesh->n_cells_with_ghosts = mesh->n_cells; /* will be increased later */

  /* Distribute faces */
  /*------------------*/

  default_face_rank = _default_face_rank(mb, comm);

  cs_block_to_part_t *d
    = cs_block_to_part_create_by_adj_s(comm,
                                       mb->face_bi,
                                       mb->cell_bi,
                                       2,
                                       mb->face_cells,
                                       mb->cell_rank,
                                       default_face_rank);

  if (default_face_rank != NULL)
    BFT_FREE(default_face_rank);

  BFT_FREE(mb->cell_rank); /* Not needed anymore */

  _n_faces = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_face_gcells, _n_faces*2, cs_gnum_t);

  /* Face -> cell connectivity */

  cs_block_to_part_copy_array(d,
                              CS_GNUM_TYPE,
                              2,
                              mb->face_cells,
                              _face_gcells);

  BFT_FREE(mb->face_cells);

  /* Now convert face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces, cs_lnum_2_t);

  cs_block_to_part_global_to_local(_n_faces*2,
                                   0,
                                   mesh->n_cells,
                                   true,
                                   mesh->global_cell_num,
                                   _face_gcells,
                                   (cs_lnum_t *)_face_cells);

  BFT_FREE(_face_gcells);

  /* Face family */

  BFT_MALLOC(_face_gc_id, _n_faces, cs_lnum_t);

  cs_block_to_part_copy_array(d,
                              CS_LNUM_TYPE,
                              1,
                              mb->face_gc_id,
                              _face_gc_id);

  BFT_FREE(mb->face_gc_id);

  /* Face level */

  BFT_MALLOC(_face_r_gen, _n_faces, char);

  if (mb->have_face_r_gen)
    cs_block_to_part_copy_array(d,
                                CS_CHAR,
                                1,
                                mb->face_r_gen,
                                _face_r_gen);
  else {
    for (cs_lnum_t i = 0; i < _n_faces; i++)
      _face_r_gen[i] = 0;
  }

  BFT_FREE(mb->face_r_gen);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  cs_block_to_part_copy_index(d,
                              mb->face_vertices_idx,
                              _face_vertices_idx);

  BFT_MALLOC(_face_gvertices, _face_vertices_idx[_n_faces], cs_gnum_t);

  cs_block_to_part_copy_indexed(d,
                                CS_GNUM_TYPE,
                                mb->face_vertices_idx,
                                mb->face_vertices,
                                _face_vertices_idx,
                                _face_gvertices);

  BFT_FREE(mb->face_vertices_idx);
  BFT_FREE(mb->face_vertices);

  _face_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Vertices */

  size_t _n_vertices = 0;

  cs_order_single_gnum(_face_vertices_idx[_n_faces],
                       1, /* base */
                       _face_gvertices,
                       &_n_vertices,
                       &(mesh->global_vtx_num));

  mesh->n_vertices = _n_vertices;

  cs_all_to_all_t *dv
    = cs_all_to_all_create_from_block(mesh->n_vertices,
                                      CS_ALL_TO_ALL_USE_DEST_ID,
                                      mesh->global_vtx_num,
                                      mb->vertex_bi,
                                      comm);

  mesh->vtx_coord = cs_all_to_all_copy_array(dv,
                                             real_type,
                                             3,
                                             true, /* reverse */
                                             mb->vertex_coords,
                                             NULL);

  BFT_FREE(mb->vertex_coords);

  cs_all_to_all_destroy(&dv);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(_face_vertices_idx[_n_faces],
                                   1,
                                   mesh->n_vertices,
                                   true,
                                   mesh->global_vtx_num,
                                   _face_gvertices,
                                   _face_vertices);

  BFT_FREE(_face_gvertices);

  /* In case of periodicity, build a cs_interface so as to obtain
     periodic face correspondants in local numbering (periodic couples
     need not be defined by the ranks owning one of the 2 members
     for the interface to be built correctly). */

  BFT_MALLOC(_periodicity_num, mb->n_perio, int);

  for (int i = 0; i < mb->n_perio; i++)
    _periodicity_num[i] = i+1;

  face_ifs
    = cs_interface_set_create(_n_faces,
                              NULL,
                              _face_num,
                              mesh->periodicity,
                              mb->n_perio,
                              _periodicity_num,
                              mb->n_per_face_couples,
                              (const cs_gnum_t *const *)mb->per_face_couples);


  if (mb->n_perio > 0) {
    BFT_FREE(_periodicity_num);
    for (int i = 0; i < mb->n_perio; i++)
      BFT_FREE(mb->per_face_couples[i]);
    BFT_FREE(mb->per_face_couples);
    BFT_FREE(mb->n_g_per_face_couples);
    BFT_FREE(mb->n_per_face_couples);
    BFT_FREE(mb->per_face_bi);
  }

  /* We may now separate interior from boundary faces */

  BFT_MALLOC(face_type, _n_faces, char);

  _face_type_g(mesh,
               _n_faces,
               face_ifs,
               (const cs_lnum_2_t *)_face_cells,
               (const cs_lnum_t *)_face_vertices_idx,
               face_type);

  _extract_face_cell(mesh,
                     _n_faces,
                     (const cs_lnum_2_t *)_face_cells,
                     face_type);

  {
    cs_gnum_t _n_g_free_faces = mesh->n_g_free_faces;
    MPI_Allreduce(&_n_g_free_faces, &(mesh->n_g_free_faces), 1,
                  CS_MPI_GNUM, MPI_SUM, comm);
  }

  BFT_FREE(_face_cells);

  if (mb->n_perio == 0)
    cs_interface_set_destroy(&face_ifs);

  _extract_face_vertices(mesh,
                         _n_faces,
                         _face_vertices_idx,
                         _face_vertices,
                         face_type);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  _extract_face_gnum(mesh,
                     _n_faces,
                     _face_num,
                     face_type);

  BFT_FREE(_face_num);

  if (mb->n_perio > 0) {
    _face_ifs_to_interior(face_ifs, _n_faces, face_type);
    cs_mesh_builder_extract_periodic_faces_g(mesh->n_init_perio,
                                             mb,
                                             mesh->periodicity,
                                             mesh->global_i_face_num,
                                             face_ifs);
    cs_interface_set_destroy(&face_ifs);
  }

  _extract_face_gc_id(mesh,
                      _n_faces,
                      _face_gc_id,
                      face_type);

  BFT_FREE(_face_gc_id);

  _extract_face_r_gen(mesh,
                      _n_faces,
                      _face_r_gen,
                      face_type);
  BFT_FREE(_face_r_gen);

  BFT_FREE(face_type);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Organize data read locally and build most mesh structures
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

static void
_decompose_data_l(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb)
{
  cs_lnum_t _n_faces = 0;

  cs_lnum_2_t *_face_cells = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  char *face_type = NULL;

  /* Initialization */

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  mesh->n_cells = mb->cell_bi.gnum_range[1] - 1;
  mesh->n_cells_with_ghosts = mesh->n_cells; /* may be increased later */

  /* Cell families are already of the correct type,
     so they can simply be moved */

  mesh->cell_family = mb->cell_gc_id;
  mb->cell_gc_id = NULL;

  /* Build faces */
  /*-------------*/

  _n_faces = mb->face_bi.gnum_range[1] - 1;

  /* Now copy face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces, cs_lnum_2_t);

  for (cs_lnum_t i = 0; i < _n_faces; i++) {
    _face_cells[i][0] = mb->face_cells[i*2] - 1;
    _face_cells[i][1] = mb->face_cells[i*2+1] - 1;
  }

  BFT_FREE(mb->face_cells);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  for (cs_lnum_t i = 0; i < _n_faces+1; i++)
    _face_vertices_idx[i] = mb->face_vertices_idx[i];

  BFT_FREE(mb->face_vertices_idx);

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  for (cs_lnum_t i = 0; i < _face_vertices_idx[_n_faces]; i++)
    _face_vertices[i] = mb->face_vertices[i];

  BFT_FREE(mb->face_vertices);

  /* Vertices */

  mesh->n_vertices = mb->vertex_bi.gnum_range[1] - 1;

  mesh->vtx_coord = mb->vertex_coords;
  mb->vertex_coords = NULL;

  /* We may now separate interior from boundary faces */

  BFT_MALLOC(face_type, _n_faces, char);

  _face_type_l(mesh,
               _n_faces,
               mb->n_per_face_couples,
               (const cs_gnum_t *const *)mb->per_face_couples,
               (const cs_lnum_2_t *)_face_cells,
               _face_vertices_idx,
               face_type);

  _extract_face_cell(mesh,
                     _n_faces,
                     (const cs_lnum_2_t *)_face_cells,
                     face_type);

  BFT_FREE(_face_cells);

  if (mb->n_perio > 0) {

    /* Transfer arrays from reader to builder, then renumber couples */

    _extract_periodic_faces_l(mb,
                              mesh->n_init_perio,
                              _n_faces,
                              face_type);

    BFT_FREE(mb->n_g_per_face_couples);
    BFT_FREE(mb->per_face_bi);

  }

  _extract_face_vertices(mesh,
                         _n_faces,
                         _face_vertices_idx,
                         _face_vertices,
                         face_type);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  _extract_face_gc_id(mesh,
                      _n_faces,
                      mb->face_gc_id,
                      face_type);

  BFT_FREE(mb->face_gc_id);

  if (mb->have_face_r_gen) {
    _extract_face_r_gen(mesh,
                        _n_faces,
                        mb->face_r_gen,
                        face_type);
    BFT_FREE(mb->face_r_gen);
  }
  else {
    BFT_MALLOC(mesh->i_face_r_gen, mesh->n_i_faces, char);
    for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++)
      mesh->i_face_r_gen[i] = 0;
  }

  BFT_FREE(face_type);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer mesh builder to mesh structure.
 *
 * \param[in, out]  mesh          pointer to mesh structure
 * \param[in, out]  mesh_builder  pointer to mesh builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_from_builder(cs_mesh_t             *mesh,
                     cs_mesh_builder_t     *mesh_builder)
{
  /* Clear previous builder data if present (periodicity done separately) */

  cs_mesh_free_rebuildable(mesh, true);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _decompose_data_g(mesh,
                      mesh_builder,
                      cs_glob_mpi_comm);

#endif

  if (cs_glob_n_ranks == 1)
    _decompose_data_l(mesh, mesh_builder);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
