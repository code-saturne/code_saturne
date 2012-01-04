/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between Code_Saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smooth the mesh.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
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
#include "cs_join.h"
#include "cs_join_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_smoother.h"
#include "cs_mesh_thinwall.h"
#include "cs_preprocessor_data.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define mesh files to read and optional associated transformations.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_input(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Determine list of files to add */
  /*--------------------------------*/

  /* Read input mesh with no modification */
  {
    cs_preprocessor_data_add_file("mesh_input/mesh_01", 0, NULL, NULL);
  }

  /* Add same mesh with transformations */
  {
    const char *renames[] = {"Inlet", "Injection_2",
                             "Group_to_remove", NULL};
    const double transf_matrix[3][4] = {{1., 0., 0., 5.},
                                        {0., 1., 0., 0.},
                                        {0., 0., 1., 0.}};

    cs_preprocessor_data_add_file("mesh_input/mesh_02",
                                  2, renames,
                                  transf_matrix);
  }
}

/*----------------------------------------------------------------------------
 * Define mesh joinings.
 *
 * This is done by calling the cs_join_add() function for each joining
 * operation to add.
 *
 * The arguments to cs_join_add() are:
 *   sel_criteria <-- boundary face selection criteria string
 *   fraction     <-- value of the fraction parameter;
 *                    the initial tolerance radius associated to each vertex
 *                    is equal to the lenght of the shortest incident edge,
 *                    multiplied by this fraction.
 *   plane        <-- value of the plane parameter;
 *                    when subdividing faces, 2 faces are considered
 *                    coplanar and may be joined if angle between their
 *                    normals (in degrees) does not exceed this parameter.
 *   verbosity    <-- level of verbosity required
 *
 * The function returns a number (1 to n) associated with the
 * new joining. This number may optionnally be used to assign advanced
 * parameters to the joining.
 *----------------------------------------------------------------------------*/

void
cs_user_join(void)
{
  int    join_num;
  int    verbosity = 1, visualization;
  float  fraction = 0.10, plane = 0.25;
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  fraction = 0.10;
  plane = 25.0;
  verbosity = 1; /* debug level if >= 3 */
  visualization = 1; /* debug level if >= 3 */

  /* Add a joining operation */
  /* ----------------------- */

  join_num = cs_join_add("98 or 99",
                         fraction,
                         plane,
                         verbosity,
                         visualization);

  /*--------------------------------------------------------------------------*/

  /* Example with advanced parameters;
     Advanced parameters may be modified to solve errors during the
     joining step or to get a better mesh quality. */

  {
    /* Merge tolerance factor:
     * used to locally modify the tolerance associated to each
     * vertex BEFORE the merge step.
     *   = 0 => no vertex merge;
     *   < 1 => vertex merge is more strict. It may increase the number
     *          of tolerance reduction and so define smaller subset of
     *          vertices to merge together but it can drive to loose
     *          intersections;
     *   = 1 => no change;
     *   > 1 => vertex merge is less strict. The subset of vertices able
     *          to be merged together is greater. */

     double mtf = 1.00;

     /* Pre-merge factor:
      * this parameter is used to define a bound under which two vertices
      * are merged before the merge step.
      * Tolerance limit for the pre-merge = pmf * fraction. */

     double pmf = 0.10;

     /* Tolerance computation mode:
      *
      *   1: (default) tol = min. edge length related to a vertex * fraction
      *   2: tolerance is computed like in mode 1 with in addition, the
      *      multiplication by a coefficient equal to the max sin(e1, e2)
      *      where e1 and e2 are two edges sharing the same vertex V for which
      *      we want to compute the tolerance.
      *  11: as 1 but taking into account only the selected faces
      *  12: as 2 but taking into account only the selected faces */

     int tcm = 1;

      /* Intersection computation mode:
       *  1: (default) Original algorithm.
       *     Try to snap intersection to extremity.
       *  2: New intersection algorithm.
       *     Avoid snapping intersection to extremity. */

     int icm = 1;

     /* Maximum number of equivalence breaks which is
      * enabled during the merge step. */

     int max_break = 500;

     /* Maximum number of sub-faces when splitting a selected face. */

     int max_sub_face = 100;

     /* tml, tmb and tmr are parameters of the searching algorithm for
      * face intersections between selected faces (octree structure).
      * Useful to adjust speed vs. memory consumption. */

     /* Tree Max Level:
      * deepest level reachable when building the tree */

     int tml = 30;

     /* Tree Max Boxes:
      * max. number of bounding boxes which can be linked to a leaf of the tree
      * (not necessary true for the deepest level) */

     int tmb = 25;

     /* Tree Max. Ratio:
      * stop refining the tree structure when
      * number of bounding boxes > tmr * number of faces to locate.
      * In parallel, a separate (usually lower) value may be set for
      * the initial coarse tree used to determine distribution.
      * Reducing this will help reduce memory consumption. */

     double tmr = 5.0;
     double tmr_distrib = 2.0;

     /* Set advanced parameters */

     cs_join_set_advanced_param(join_num,
                                mtf, pmf, tcm, icm,
                                max_break, max_sub_face,
                                tml, tmb, tmr, tmr_distrib);
  }
}

/*----------------------------------------------------------------------------
 * Define periodic faces.
 *
 * This is done by calling one of the cs_join_perio_add_*() functions for
 * each periodicity to add.
 *
 * The first arguments to cs_join_perio_add_() are the same as for
 * mesh joining:
 *   sel_criteria <-- boundary face selection criteria string
 *   fraction     <-- value of the fraction parameter;
 *                    the initial tolerance radius associated to each vertex
 *                    is equal to the lenght of the shortest incident edge,
 *                    multiplied by this fraction.
 *   plane        <-- value of the plane parameter;
 *                    when subdividing faces, 2 faces are considered
 *                    coplanar and may be joined if angle between their
 *                    normals (in degrees) does not exceed this parameter.
 *   verbosity    <-- level of verbosity required
 *
 * The last arguments depend on the type of periodicity to define,
 * and are described below.
 *
 * The function returns a number (1 to n) associated with the
 * new joining. This number may optionnally be used to assign advanced
 * parameters to the joining.
 *----------------------------------------------------------------------------*/

void
cs_user_periodicity(void)
{
  int    join_num;
  int    verbosity = 1, visualization = 1;
  float  fraction = 0.10, plane = 0.25;
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  fraction = 0.10;
  plane = 25.0;
  verbosity = 1; /* processor-local files if > 1, debug level if >= 3 */
  visualization = 1; /* debug level if >= 3 */

  /* Example 1: define a periodicity of translation */
  /* ---------------------------------------------- */

  {
    const double translation[3] = {1.0, 0.0, 0.0}; /* Translation vector */

    join_num = cs_join_perio_add_translation("98 or 99",
                                             fraction,
                                             plane,
                                             verbosity,
                                             visualization,
                                             translation);
  }

  /* Example 2: define a periodicity of rotation */
  /* ------------------------------------------- */

  {
    double  theta = 20;                /* angle in degrees */
    double  axis[3] = {1.0, 0, 0};     /* axis of rotation */
    double  invariant[3] = {0, 0, 0};  /* invariant point */

    /* change default values */
    fraction = 0.2;
    verbosity = 2;

    join_num = cs_join_perio_add_rotation("3",
                                          fraction,
                                          plane,
                                          verbosity,
                                          visualization,
                                          theta,
                                          axis,
                                          invariant);

    /* restore default values */
    fraction = 0.1;
    verbosity = 1;
  }

  /* Example 3: define a general periodicity */
  /* --------------------------------------- */

  /* We define a general transformation using a homogeneous matrix:
   *
   * We define the first 3 rows of a 4x4 matrix:
   *    _               _
   *   | r11 r12 r13 tx  |  t(x,y,z) : translation vector
   *   | r21 r22 r23 ty  |  r(i,j)   : rotation matrix
   *   | r31 r32 r33 tz  |
   *   |_  0   0   0  1 _|
   *
   * Transformations may be combined using matrix multiplication,
   * so this be used for helecoidal transformations for instance. */

  {
    double matrix[3][4] = {{1., 0., 0., 0.5},
                           {0., 1., 0., 0.},
                           {0., 0., 1., 0.}};

    join_num = cs_join_perio_add_mixed("all[]",
                                       fraction,
                                       plane,
                                       verbosity,
                                       visualization,
                                       matrix);
  }

  /*--------------------------------------------------------------------------*/

  /* Example with advanced parameters;
     Advanced parameters may be modified to solve errors during the
     joining step or to get a better mesh quality. */

  {
    /* Merge tolerance factor:
     * used to locally modify the tolerance associated to each
     * vertex BEFORE the merge step.
     *   = 0 => no vertex merge;
     *   < 1 => vertex merge is more strict. It may increase the number
     *          of tolerance reduction and so define smaller subset of
     *          vertices to merge together but it can drive to loose
     *          intersections;
     *   = 1 => no change;
     *   > 1 => vertex merge is less strict. The subset of vertices able
     *          to be merged together is greater. */

     double mtf = 1.00;

     /* Pre-merge factor:
      * this parameter is used to define a bound under which two vertices
      * are merged before the merge step.
      * Tolerance limit for the pre-merge = pmf * fraction. */

     double pmf = 0.10;

     /* Tolerance computation mode:
      *
      *   1: (default) tol = min. edge length related to a vertex * fraction
      *   2: tolerance is computed like in mode 1 with in addition, the
      *      multiplication by a coefficient equal to the max sin(e1, e2)
      *      where e1 and e2 are two edges sharing the same vertex V for which
      *      we want to compute the tolerance.
      *  11: as 1 but taking into account only the selected faces
      *  12: as 2 but taking into account only the selected faces */

     int tcm = 1;

      /* Intersection computation mode:
       *  1: (default) Original algorithm.
       *     Try to snap intersection to extremity.
       *  2: New intersection algorithm.
       *     Avoid snapping intersection to extremity. */

     int icm = 1;

     /* Maximum number of equivalence breaks which is
      * enabled during the merge step. */

     int max_break = 500;

     /* Maximum number of sub-faces when splitting a selected face. */

     int max_sub_face = 100;

     /* tml, tmb and tmr are parameters of the searching algorithm for
      * face intersections between selected faces (octree structure).
      * Useful to adjust speed vs. memory consumption. */

     /* Tree Max Level:
      * deepest level reachable when building the tree */

     int tml = 30;

     /* Tree Max Boxes:
      * max. number of bounding boxes which can be linked to a leaf of the tree
      * (not necessary true for the deepest level) */

     int tmb = 25;

     /* Tree Max. Ratio:
      * stop refining the tree structure when
      * number of bounding boxes > tmr * number of faces to locate.
      * In parallel, a separate (usually lower) value may be set for
      * the initial coarse tree used to determine distribution.
      * Reducing this will help reduce memory consumption. */

     double tmr = 5.0;
     double tmr_distrib = 2.0;

     /* Set advanced parameters */

     cs_join_set_advanced_param(join_num,
                                mtf, pmf, tcm, icm,
                                max_break, max_sub_face,
                                tml, tmb, tmr, tmr_distrib);
  }
}

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

/*----------------------------------------------------------------------------
 * Modify geometry and mesh.
 *
 * The mesh structure is described in cs_mesh.h
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Example: modify vertex coordinates */
  /*------------------------------------*/

  /* Divide coordinates by 1000 (millimetres to metres).
   *
   * Warning:
   *
   *   This is incompatible with pre-processed periodicity,
   *   as the periodicity transformation is not updated.
   *
   *   With periodicity, using a coordinate transformation matrix
   *   in cs_user_mesh_input is preferred. */

  {
    cs_lnum_t  vtx_id;
    const double  coo_mult = 1. / 1000.;

    for (vtx_id = 0; vtx_id < mesh->n_vertices; vtx_id++) {
      mesh->vtx_coord[vtx_id*3]     *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 1] *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 2] *= coo_mult;
    }

    /* Set mesh modification flag it i should be saved for future re-use */

    mesh->modified = 1;
  }
}

/*----------------------------------------------------------------------------
 * Mesh smoothing.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_smoothe(cs_mesh_t  *mesh)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  double feature_angle = 25; /* bounded between 0 and 90 degrees */
  int *vtx_is_fixed = NULL;

  BFT_MALLOC(vtx_is_fixed, mesh->n_vertices, int);

  /* Get fixed boundary vertices flag */

  cs_mesh_smoother_fix_by_feature(mesh,
                                  feature_angle,
                                  vtx_is_fixed);

  /* Call unwarping smoother */

  cs_mesh_smoother_unwarp(mesh, vtx_is_fixed);

  /* Free memory */

  BFT_FREE(vtx_is_fixed);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
