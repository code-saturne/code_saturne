/*============================================================================
 * Cut warped faces in serial or parallel with/without periodicity.
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
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_triangulate.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_connect.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_warping.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global variables
 *============================================================================*/

static int cs_glob_mesh_warping_post = 0;
static double cs_glob_mesh_warping_threshold = -1.0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create the list of faces to cut in order to respect warping criterion.
 *
 * parameters:
 *   n_faces         <-- number of faces
 *   max_warp_angle  <-- criterion above which face is cut
 *   face_warping    <-- face warping angle
 *   p_n_warp_faces  <-> pointer to the number of warped faces
 *   p_warp_face_lst <-> pointer to the warped face list
 *----------------------------------------------------------------------------*/

static void
_select_warped_faces(cs_lnum_t       n_faces,
                     double          max_warp_angle,
                     double          face_warping[],
                     cs_lnum_t      *p_n_warp_faces,
                     cs_lnum_t      *p_warp_face_lst[])
{
  cs_lnum_t  face_id;

  cs_lnum_t  n_warp_faces = 0;
  cs_lnum_t  *warp_face_lst = NULL;

  if (n_faces > 0) {

    for (face_id = 0; face_id < n_faces; face_id++)
      if (face_warping[face_id] >= max_warp_angle)
        n_warp_faces++;

    BFT_MALLOC(warp_face_lst, n_warp_faces, cs_lnum_t);

    n_warp_faces = 0;

    for (face_id = 0; face_id < n_faces; face_id++)
      if (face_warping[face_id] >= max_warp_angle)
        warp_face_lst[n_warp_faces++] = face_id + 1;

  }

  *p_n_warp_faces = n_warp_faces;
  *p_warp_face_lst = warp_face_lst;
}

/*----------------------------------------------------------------------------
 * Cut faces if necessary and update connectivity without periodicity
 *
 * parameters:
 *   mesh                    <-> pointer to a mesh structure
 *   stride                  <-- 2 or 1 for internal or border faces
 *   p_n_cut_faces           <-> in:  number of faces to cut
 *                               out: number of cut faces
 *   p_cut_face_lst          <-> pointer to the cut face list
 *   p_n_sub_elt_lst         <-> pointer to the sub-elt count list
 *   p_n_faces               <-> pointer to the number of faces
 *   p_face_vtx_connect_size <-> size of the "face -> vertex" connectivity
 *   p_face_cells            <-> "face -> cells" connectivity
 *   p_face_family           <-> face family
 *   p_face_r_gen            <-> face refinement generation
 *   p_face_vtx_idx          <-> pointer on "face -> vertices" connect. index
 *   p_face_vtx_lst          <-> pointer on "face -> vertices" connect. list
 *----------------------------------------------------------------------------*/

static void
_cut_warped_faces(cs_mesh_t      *mesh,
                  int             stride,
                  cs_lnum_t      *p_n_cut_faces,
                  cs_lnum_t      *p_cut_face_lst[],
                  cs_lnum_t      *p_n_sub_elt_lst[],
                  cs_lnum_t      *p_n_faces,
                  cs_lnum_t      *p_face_vtx_connect_size,
                  cs_lnum_t      *p_face_cells[],
                  int            *p_face_family[],
                  char           *p_face_r_gen[],
                  cs_lnum_t      *p_face_vtx_idx[],
                  cs_lnum_t      *p_face_vtx_lst[])
{
  cs_lnum_t  i, j, face_id, idx_start, idx_end;
  cs_lnum_t  n_triangles;

  cs_lnum_t  n_face_vertices = 0, n_max_face_vertices = 0;
  cs_lnum_t  n_new_faces = 0, n_cut_faces = 0, connect_size = 0;

  fvm_triangulate_state_t  *triangle_state = NULL;
  cs_lnum_t  *new_face_vtx_idx = NULL, *new_face_vtx_lst = NULL;
  cs_lnum_t  *new_face_cells = NULL;
  int        *new_face_family = NULL;
  char       *new_face_r_gen = NULL;
  cs_lnum_t  *cut_face_lst = NULL;
  cs_lnum_t  *n_sub_elt_lst = NULL;
  char *cut_flag = NULL;

  const cs_lnum_t  dim = mesh->dim;
  const cs_lnum_t  n_init_faces = *p_n_faces;

  assert(stride == 1 || stride ==2);
  assert(dim == 3);

  BFT_MALLOC(n_sub_elt_lst, n_init_faces, cs_lnum_t);

  /* Build flag for each face from list of faces to cut */

  BFT_MALLOC(cut_flag, n_init_faces, char);

  for (face_id = 0; face_id < n_init_faces; face_id++)
    cut_flag[face_id] = 0;

  for (i = 0; i < *p_n_cut_faces; i++)
    cut_flag[(*p_cut_face_lst)[i] - 1] = 1;

  BFT_FREE(*p_cut_face_lst);

  /* First loop: count */

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id];
    idx_end = (*p_face_vtx_idx)[face_id + 1];

    n_face_vertices = idx_end - idx_start;
    n_max_face_vertices = CS_MAX(n_max_face_vertices, n_face_vertices);

    if (cut_flag[face_id] != 0) {

      n_triangles = n_face_vertices - 2;
      connect_size += n_triangles*3;
      n_new_faces += n_triangles;
      n_cut_faces += n_triangles;
      n_sub_elt_lst[face_id] = n_triangles;

    }
    else {

      connect_size += n_face_vertices;
      n_new_faces += 1;
      n_sub_elt_lst[face_id] = 1;

    }

  } /* End of loop on faces */

  *p_n_sub_elt_lst = n_sub_elt_lst;

  if (n_cut_faces == 0) {
    BFT_FREE(cut_flag);
    return;
  }

  BFT_MALLOC(new_face_vtx_idx, n_new_faces + 1, cs_lnum_t);
  BFT_MALLOC(new_face_vtx_lst, connect_size, cs_lnum_t);
  BFT_MALLOC(new_face_cells, n_new_faces*stride, cs_lnum_t);
  BFT_MALLOC(new_face_family, n_new_faces, int);
  if (p_face_r_gen != NULL)
    BFT_MALLOC(new_face_r_gen, n_new_faces, char);

  BFT_MALLOC(cut_face_lst, n_cut_faces, cs_lnum_t);

  triangle_state = fvm_triangulate_state_create(n_max_face_vertices);

  /* Second loop: define the new connectivity after triangulation */

  new_face_vtx_idx[0] = 0;
  connect_size = 0;
  n_new_faces = 0;
  n_cut_faces = 0;

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id];
    idx_end = (*p_face_vtx_idx)[face_id + 1];
    n_face_vertices = idx_end - idx_start;

    if (cut_flag[face_id] != 0) {

      n_triangles = fvm_triangulate_polygon(dim,
                                            0,
                                            n_face_vertices,
                                            mesh->vtx_coord,
                                            NULL,
                                            (*p_face_vtx_lst) + idx_start,
                                            FVM_TRIANGULATE_MESH_DEF,
                                            new_face_vtx_lst + connect_size,
                                            triangle_state);

      assert(n_triangles == n_face_vertices - 2);

      /* Update face -> vertex connectivity */

      for (i = 0; i < n_triangles; i++) {

        cut_face_lst[n_cut_faces++] = n_new_faces + 1;

        /* Update "face -> cells" connectivity */

        for (j = 0; j < stride; j++)
          new_face_cells[stride*n_new_faces + j]
            = (*p_face_cells)[stride*face_id + j];

        /* Update family and refinement generation for each face */

        new_face_family[n_new_faces] = (*p_face_family)[face_id];
        if (p_face_r_gen != NULL)
          new_face_r_gen[n_new_faces] = (*p_face_r_gen)[face_id];

        /* Update "face -> vertices" connectivity index
           (list has already been defined by fvm_triangulate_polygon) */

        n_new_faces++;
        connect_size += 3;
        new_face_vtx_idx[n_new_faces] = new_face_vtx_idx[n_new_faces-1] + 3;

      } /* End of loop on triangles */

    }
    else {

      /* Update "face -> cells" connectivity */

      for (j = 0; j < stride; j++)
        new_face_cells[stride*n_new_faces + j]
          = (*p_face_cells)[stride*face_id + j];

      /* Update family and refinemen generation for each face */

      new_face_family[n_new_faces] = (*p_face_family)[face_id];
      if (p_face_r_gen != NULL)
        new_face_r_gen[n_new_faces] = (*p_face_r_gen)[face_id];

      /* Update "face -> vertices" connectivity */

      for (j = 0, i = idx_start; i < idx_end; i++, j++)
        new_face_vtx_lst[connect_size + j] = (*p_face_vtx_lst)[i];

      n_new_faces++;
      connect_size += n_face_vertices;
      new_face_vtx_idx[n_new_faces] =
        new_face_vtx_idx[n_new_faces-1] + n_face_vertices;

    }

  } /* End of loop on internal faces */

  triangle_state = fvm_triangulate_state_destroy(triangle_state);

  BFT_FREE(cut_flag);

  BFT_FREE(*p_face_vtx_idx);
  BFT_FREE(*p_face_vtx_lst);
  BFT_FREE(*p_face_cells);
  BFT_FREE(*p_face_family);
  if (p_face_r_gen != NULL)
    BFT_FREE(*p_face_r_gen);

  /* Define returned pointers */

  *p_face_vtx_idx = new_face_vtx_idx;
  *p_face_vtx_lst = new_face_vtx_lst;
  *p_face_cells = new_face_cells;
  *p_face_family = new_face_family;
  *p_face_r_gen = new_face_r_gen;
  *p_face_vtx_connect_size = connect_size;
  *p_n_faces = n_new_faces;
  *p_n_cut_faces = n_cut_faces;

  *p_cut_face_lst = cut_face_lst;
}

/*----------------------------------------------------------------------------
 * Update warped faces global numbers after cutting
 *
 * parameters:
 *   mesh              <-> pointer to a mesh structure
 *   n_faces           <-- number of faces
 *   n_init_faces      <-- initial number of faces
 *   n_cut_faces       <-- number of cut faces
 *   cut_face_lst      <-- pointer to the cut face list
 *   n_sub_elt_lst     <-- sub-elt count list
 *   n_g_faces         <-> global number of faces
 *   p_global_face_num <-> pointer to the global face numbers
 *----------------------------------------------------------------------------*/

static void
_update_cut_faces_num(cs_mesh_t      *mesh,
                      cs_lnum_t       n_faces,
                      cs_lnum_t       n_init_faces,
                      cs_lnum_t       n_sub_elt_lst[],
                      cs_gnum_t      *n_g_faces,
                      cs_gnum_t     **p_global_face_num)
{
  size_t  size;

  fvm_io_num_t *new_io_num = NULL, *previous_io_num = NULL;
  const cs_gnum_t  *global_num = NULL;

  /* Simply update global number of faces in trivial case */

  *n_g_faces = n_faces;

  if (*p_global_face_num == NULL)
    return;

  /* Faces should not have been reordered */

  if (cs_order_gnum_test(NULL, *p_global_face_num, n_init_faces) == false)
    bft_error(__FILE__, __LINE__, 0,
              _("The faces have been renumbered before cutting.\n"
                "This case should not arise, because the mesh entities\n"
                "should be cut before renumbering."));

  /* Update global number of internal faces and its global numbering */

  if (mesh->n_domains > 1) {

    previous_io_num = fvm_io_num_create(NULL,
                                        *p_global_face_num,
                                        n_init_faces,
                                        0);
    new_io_num = fvm_io_num_create_from_sub(previous_io_num,
                                            n_sub_elt_lst);

    previous_io_num = fvm_io_num_destroy(previous_io_num);

    *n_g_faces = fvm_io_num_get_global_count(new_io_num);

    global_num = fvm_io_num_get_global_num(new_io_num);

    BFT_REALLOC(*p_global_face_num, n_faces, cs_gnum_t);
    size = sizeof(cs_gnum_t) * n_faces;
    memcpy(*p_global_face_num, global_num, size);

    new_io_num = fvm_io_num_destroy(new_io_num);

  }
}

/*----------------------------------------------------------------------------
 * Post-process the warped faces before cutting.
 *
 * parameters:
 *   n_b_warp_faces  <-- number of border warped faces
 *   b_warp_face_lst <-- border warped face list
 *   b_face_warping  <-- face warping angle for internal faces
 *----------------------------------------------------------------------------*/

static void
_post_before_cutting(cs_lnum_t       n_b_warp_faces,
                     cs_lnum_t       b_warp_face_lst[],
                     double          b_face_warping[])
{
  cs_lnum_t  parent_num_shift[2];

  int  n_parent_lists = 2;
  fvm_nodal_t  *fvm_mesh = NULL;
  fvm_writer_t  *writer = NULL;

  const cs_lnum_t  writer_id = -1; /* default writer */
  const void  *var_ptr[2] = {NULL, NULL};

  parent_num_shift[0] = 0;
  parent_num_shift[1] = cs_glob_mesh->n_b_faces;

  if (cs_post_writer_exists(writer_id) == false)
    return;

  assert(sizeof(double) == sizeof(cs_real_t));

  fvm_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            _("Warped faces to cut"),
                                            false,
                                            0,
                                            n_b_warp_faces,
                                            NULL,
                                            b_warp_face_lst);

  writer = cs_post_get_writer(writer_id);

  fvm_writer_set_mesh_time(writer, -1, 0.0);

  /* Write a mesh from the selected faces */

  fvm_writer_export_nodal(writer, fvm_mesh);

  /* Write the warping field */

  var_ptr[0] = ((const char *)b_face_warping);
  var_ptr[1] = (NULL);

  fvm_writer_export_field(writer,
                          fvm_mesh,
                          _("Face warping"),
                          FVM_WRITER_PER_ELEMENT,
                          1,
                          CS_INTERLACE,
                          n_parent_lists,
                          parent_num_shift,
                          CS_DOUBLE,
                          (int)-1,
                          (double)0.0,
                          (const void **)var_ptr);

  fvm_mesh = fvm_nodal_destroy(fvm_mesh);
}

/*----------------------------------------------------------------------------
 * Post-process the warped faces after cutting.
 *
 * parameters:
 *   n_b_cut_faces  <-- number of border faces generated by cutting
 *   b_cut_face_lst <-- face warping angle for internal faces
 *----------------------------------------------------------------------------*/

static void
_post_after_cutting(cs_lnum_t   n_b_cut_faces,
                    cs_lnum_t   b_cut_face_lst[])
{
  fvm_nodal_t  *fvm_mesh = NULL;
  fvm_writer_t  *writer = NULL;

  const int  writer_id = -1; /* default writer */

  if (cs_post_writer_exists(writer_id) == false)
    return;

  fvm_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            _("Warped faces after cutting"),
                                            false,
                                            0,
                                            n_b_cut_faces,
                                            NULL,
                                            b_cut_face_lst);

  writer = cs_post_get_writer(writer_id);

  fvm_writer_set_mesh_time(writer, -1, 0.0);

  /* Write a mesh from the selected faces */

  fvm_writer_export_nodal(writer, fvm_mesh);

  fvm_mesh = fvm_nodal_destroy(fvm_mesh);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut warped boundary faces.
 *
 * Update border face connectivity and associated mesh quantities.
 *
 * \param[in]  mesh            pointer to mesh structure
 * \param[in]  max_warp_angle  criterion to know which face to cut
 * \param[in]  post_flag       1 if we have to post-process cut faces,
 *                             0 otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_warping_cut_faces(cs_mesh_t  *mesh,
                          double      max_warp_angle,
                          bool        post_flag)
{
  cs_lnum_t  n_b_cut_faces = 0;
  cs_lnum_t  *b_face_lst = NULL;
  double  *b_face_warping = NULL;
  cs_lnum_t  *n_b_sub_elt_lst = NULL;
  cs_gnum_t  n_g_faces_ini = 0;
  cs_gnum_t  n_g_b_cut_faces = 0;

  const cs_lnum_t  n_init_b_faces = mesh->n_b_faces;

#if 0   /* DEBUG */
  cs_mesh_dump(mesh);
#endif

  bft_printf(_("\n\n Cutting of warped faces requested\n"
               " ---------------------------------\n\n"
               " Maximum allowed angle (deg): %7.4f\n\n"), max_warp_angle);

  /* Compute face warping */

  BFT_MALLOC(b_face_warping, n_init_b_faces, double);

  for (cs_lnum_t i = 0; i < n_init_b_faces; i++)
    b_face_warping[i] = 0.;

  cs_real_t  *b_face_normal = NULL, *b_face_cog = NULL;;

  cs_mesh_quantities_b_faces(mesh,
                             &b_face_cog,
                             &b_face_normal);

  BFT_FREE(b_face_cog);

  cs_mesh_quality_compute_b_face_warping(mesh,
                                         b_face_normal,
                                         b_face_warping);

  BFT_FREE(b_face_normal);

  _select_warped_faces(n_init_b_faces,
                       max_warp_angle,
                       b_face_warping,
                       &n_b_cut_faces,
                       &b_face_lst);

  /* Define the global number of faces which need to be cut */

  n_g_b_cut_faces = n_b_cut_faces;
  cs_parall_counter(&n_g_b_cut_faces, 1);

  /* Test if there are faces to cut to continue */

  if (n_g_b_cut_faces == 0) {

    BFT_FREE(b_face_lst);
    BFT_FREE(b_face_warping);

    bft_printf(_("\n No face to cut. Verify the criterion if necessary.\n"));
    return;
  }

  /* Post-processing management */

  if (post_flag == true)
    _post_before_cutting(n_b_cut_faces,
                         b_face_lst,
                         b_face_warping);

  BFT_FREE(b_face_warping);

  /* Border face treatment */
  /* --------------------- */

  n_g_faces_ini = mesh->n_g_b_faces;

  _cut_warped_faces(mesh,
                    1,
                    &n_b_cut_faces,
                    &b_face_lst,
                    &n_b_sub_elt_lst,
                    &mesh->n_b_faces,
                    &mesh->b_face_vtx_connect_size,
                    &mesh->b_face_cells,
                    &mesh->b_face_family,
                    NULL,
                    &mesh->b_face_vtx_idx,
                    &mesh->b_face_vtx_lst);

  /* Update global number of border faces and its global numbering */

  _update_cut_faces_num(mesh,
                        mesh->n_b_faces,
                        n_init_b_faces,
                        n_b_sub_elt_lst,
                        &(mesh->n_g_b_faces),
                        &(mesh->global_b_face_num));

  bft_printf(_(" Boundary faces:\n\n"
               "   %12llu faces before cutting\n"
               "   %12llu faces after cutting\n\n"),
             (unsigned long long)n_g_faces_ini,
             (unsigned long long)(mesh->n_g_b_faces));

  /* Partial memory free */

  BFT_FREE(n_b_sub_elt_lst);

  /* post processing of the selected faces */

  if (post_flag == true)
    _post_after_cutting(n_b_cut_faces,
                        b_face_lst);

  /* Free memory */

  BFT_FREE(b_face_lst);

  /* Set mesh modification flag */

  mesh->modified |= CS_MESH_MODIFIED;
}

/*----------------------------------------------------------------------------
 * Set defaults for cutting of warped faces.
 *
 * parameters:
 *   max_warp_angle <-- maximum warp angle (in degrees) over which faces will
 *                      be cut; negative (-1) if faces should not be cut
 *   postprocess    <-- 1 if postprocessing should be activated when cutting
 *                      warped faces, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_set_defaults(double  max_warp_angle,
                             int     postprocess)
{
  if (max_warp_angle >= 0.0 && max_warp_angle <= 180.0)
    cs_glob_mesh_warping_threshold = max_warp_angle;
  else
    cs_glob_mesh_warping_threshold = -1.0;

  if (postprocess != 0)
    cs_glob_mesh_warping_post = 1;
}

/*----------------------------------------------------------------------------
 * Get defaults for cutting of warped faces.
 *
 * parameters:
 *   max_warp_angle --> if non NULL, returns maximum warp angle (in degrees)
 *                      over which faces will be cut, or -1 if faces should
 *                      not be cut
 *   postprocess    --> if non NULL, returns 1 if postprocessing should be
 *                      activated when cutting warped faces, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_get_defaults(double  *max_warp_angle,
                             int     *postprocess)
{
  if (max_warp_angle != NULL)
    *max_warp_angle = cs_glob_mesh_warping_threshold;

  if (postprocess != NULL)
    *postprocess = cs_glob_mesh_warping_post;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
