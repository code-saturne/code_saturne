/*============================================================================
 * Definition of the calculation mesh.
 *
 * In this example, group information from free faces is projected to
 * the neighboring boundary faces that are missing such information.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh-groups_from_free_faces.c
 *
 * \brief Mesh modification example.
 *
 * See \subpage cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Transfer group information from free faces to boundary faces.
 *
 * Values are transferred only to boundary faces with no group.
 *
 * parameters :
 * - mesh : pointer to cs_mesh_t structure
 * - tolerance : tolerance for mesh location
 */
/*----------------------------------------------------------------------------*/

static void
_mesh_groups_from_free_faces(cs_mesh_t  *mesh,
                             double      tolerance)
{
  cs_lnum_t  i, j;
  cs_lnum_t  n_free_faces = 0, n_b_faces = 0, n_no_group = 0;
  cs_lnum_t  *free_faces_list = NULL, *no_group_list = NULL;
  int *family_flag = NULL;

  if (mesh->n_g_free_faces == 0)
    return;

  /* Mark families matching "no groups" */

  BFT_MALLOC(family_flag, mesh->n_families, int);
  for (i = 0; i < mesh->n_families; i++)
    family_flag[i] = 1;

  for (i = 0; i < mesh->n_families; i++) {
    for (j = 0; j <  mesh->n_max_family_items; j++) {
      if (mesh->family_item[j * mesh->n_families + i] != 0)
        family_flag[i] = 0;
    }
  }

  /* Count values to locate and number of free faces */

  for (i = 0; i < mesh->n_b_faces; i++) {

    if (mesh->b_face_cells[i] < 0)
      n_free_faces += 1;

    else {
      n_b_faces += 1;
      if (family_flag[mesh->b_face_family[i] - 1] == 1)
        n_no_group += 1;
    }

  }

  /* Build associated lists */

  BFT_MALLOC(free_faces_list, n_free_faces, cs_lnum_t);
  BFT_MALLOC(no_group_list, n_no_group, cs_lnum_t);

  n_free_faces = 0;
  n_no_group = 0;

  for (i = 0; i < mesh->n_b_faces; i++) {

    if (mesh->b_face_cells[i] < 0)
      free_faces_list[n_free_faces++] = i+1;    /* 1-based */

    else if (family_flag[mesh->b_face_family[i] - 1] == 1)
      no_group_list[n_no_group++] = i;          /* 0-based */

  }

  BFT_FREE(family_flag);

  /* Build nodal mesh associated to isolated faces */

  fvm_nodal_t *free_faces = cs_mesh_connect_faces_to_nodal(mesh,
                                                           "Free faces",
                                                           false,
                                                           0,
                                                           n_free_faces,
                                                           NULL,
                                                           free_faces_list);

  /* Associated PLE locator */

#if defined(PLE_HAVE_MPI)
  ple_locator_t *locator = ple_locator_create(cs_glob_mpi_comm,
                                              cs_glob_n_ranks,
                                              0);
#else
  ple_locator_t *locator = ple_locator_create();
#endif

  cs_real_t *b_face_cog = NULL, * b_face_normal = NULL;

  cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

  BFT_FREE(b_face_normal);

  ple_locator_set_mesh(locator,
                       free_faces,
                       NULL,      /* options */
                       0.,        /* absolute tolerance */
                       tolerance, /* relative tolerance */
                       3,
                       n_no_group,
                       no_group_list,
                       NULL,
                       b_face_cog,
                       NULL,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  BFT_FREE(b_face_cog);

  /* Log number of found and free faces */

  {
    cs_gnum_t n_g_faces[5] = {n_free_faces, n_b_faces, n_no_group, 0, 0};
    n_g_faces[3] = ple_locator_get_n_interior(locator);
    n_g_faces[4] = ple_locator_get_n_exterior(locator);

    cs_parall_counter(n_g_faces, 5);

    bft_printf
      ("\n"
       "Projecting group information from free faces to boundary faces\n"
       "--------------------------------------------------------------\n\n"
       "  number of free faces:                     %llu\n"
       "  number of boundary faces:                 %llu\n"
       "  number of boundary faces with no group:   %llu\n"
       "  number of faces with match:               %llu\n"
       "  number of faces without  match:           %llu\n\n",
       (unsigned long long)n_g_faces[0], (unsigned long long)n_g_faces[1],
       (unsigned long long)n_g_faces[2], (unsigned long long)n_g_faces[3],
       (unsigned long long)n_g_faces[4]);
  }

  /* Now transfer information */

  ple_lnum_t n_dist_points = ple_locator_get_n_dist_points(locator);

  const ple_lnum_t *dist_loc = ple_locator_get_dist_locations(locator);

  int *dist_fm_id = NULL;
  BFT_MALLOC(dist_fm_id, n_dist_points, int);

  for (i = 0; i < n_dist_points; i++)
    dist_fm_id[i] = mesh->b_face_family[dist_loc[i] - 1];

  ple_locator_exchange_point_var(locator,
                                 dist_fm_id,
                                 mesh->b_face_family,
                                 no_group_list,
                                 sizeof(int),
                                 1,
                                 0);

  locator = ple_locator_destroy(locator);

  free_faces = fvm_nodal_destroy(free_faces);

  BFT_FREE(free_faces_list);
  BFT_FREE(no_group_list);
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * In this example, group information from free faces is projected to
 * the neighboring boundary faces that are missing such information. *
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  /* Try to transfer group information from free faces to boundary faces
     with no such information */

  _mesh_groups_from_free_faces(mesh, 1);

  /* Set mesh modification flag if it should be saved for future re-use. */

  mesh->modified = 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
