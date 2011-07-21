
/* VERS */

/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
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
 *============================================================================*/

/*============================================================================
 * User function for geometry and mesh modification.
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_selector.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include "cs_mesh_thinwall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Insert thin wall into a mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_thinwall(cs_mesh_t  *mesh)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Example: modify vertex coordinates */
  /*------------------------------------*/

  cs_lnum_t   n_selected_faces = 0;
  cs_lnum_t  *selected_faces = NULL; 

  cs_real_t  *i_face_cog = NULL, *i_face_normal = NULL;

  /* example of multi-line character string */

  const char criteria[] = "plane[0, -1, 0, 0.5, epsilon = 0.0001]"
                          " or plane[-1, 0, 0, 0.5, epsilon = 0.0001]";

  cs_mesh_init_group_classes(mesh);

  cs_mesh_quantities_i_faces(mesh, &i_face_cog, &i_face_normal);

  cs_glob_mesh->select_i_faces = fvm_selector_create(mesh->dim,
                                                     mesh->n_i_faces,
                                                     mesh->class_defs,
                                                     mesh->i_face_family,
                                                     1,
                                                     i_face_cog,
                                                     i_face_normal);

  BFT_MALLOC(selected_faces, mesh->n_i_faces, cs_int_t);
                              
  cs_selector_get_i_face_list(criteria,
                              &n_selected_faces,
                              selected_faces);
  cs_create_thinwall(mesh,
                     selected_faces,
                     n_selected_faces);
  
  BFT_FREE(i_face_cog);
  BFT_FREE(i_face_normal);

  mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);
  fvm_selector_destroy(mesh->select_i_faces);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
