/*============================================================================
 * Management of the GUI parameters file: mesh related options
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_join.h"
#include "cs_join_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_warping.h"
#include "cs_mesh_smoother.h"
#include "cs_mesh_boundary.h"
#include "cs_mesh_extrude.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the value to a face joining markup
 *
 * parameters:
 *   keyword <-- label of the markup
 *   number  <-- joining number
 *----------------------------------------------------------------------------*/

static char *
_get_face_joining(const char  *keyword,
                  int          number)
{
  char* value = NULL;
  char *path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "joining");
  cs_xpath_add_element_num(&path, "face_joining", number);
  cs_xpath_add_element(&path, keyword);
  cs_xpath_add_function_text(&path);
  value = cs_gui_get_text_value(path);
  BFT_FREE(path);
  return value;
}

/*-----------------------------------------------------------------------------
 * Return the value to a face periodicity markup
 *
 * parameters:
 *   keyword <-- label of the markup
 *   number  <-- joining number
 *----------------------------------------------------------------------------*/

static char *
_get_periodicity(const char  *keyword,
                 int          number)
{
  char* value = NULL;
  char *path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "periodicity");
  cs_xpath_add_element_num(&path, "face_periodicity", number);
  cs_xpath_add_element(&path, keyword);
  cs_xpath_add_function_text(&path);
  value = cs_gui_get_text_value(path);
  BFT_FREE(path);
  return value;
}

/*-----------------------------------------------------------------------------
 * Get transformation parameters associated with a translational periodicity
 *
 * parameter:
 *   number <-- number of the syrthes coupling
 *   trans  --> translation values
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_translation(int     number,
                             double  trans[3])
{
  size_t dir_id = 0;
  char *path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "periodicity");
  cs_xpath_add_element_num(&path, "face_periodicity", number);
  cs_xpath_add_elements(&path, 2, "translation", "translation_x");
  dir_id = strlen(path) - 1;
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, trans + 0))
    trans[0] = 0.0;

  path[dir_id] = 'y';
  if (!cs_gui_get_double(path, trans + 1))
    trans[1] = 0.0;

  path[dir_id] = 'z';
  if (!cs_gui_get_double(path, trans + 2))
    trans[2] = 0.0;

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf("==> _get_periodicity_translation\n");
  bft_printf("--translation = [%f %f %f]\n",
             trans[0], trans[1], trans[2]);
#endif

  return;
}

/*-----------------------------------------------------------------------------
 * Get transformation parameters associated with a rotational periodicity
 *
 * parameter:
 *   number    <-- number of the syrthes coupling
 *   angle     --> rotation angle
 *   axis      --> rotation axis
 *   invariant --> invariant point
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_rotation(int      number,
                          double  *angle,
                          double   axis[3],
                          double   invariant[3])
{
  size_t coo_id = 0;
  char *path = NULL;
  char *path_0 = cs_xpath_init_path();
  cs_xpath_add_elements(&path_0, 2, "solution_domain", "periodicity");
  cs_xpath_add_element_num(&path_0, "face_periodicity", number);
  cs_xpath_add_element(&path_0, "rotation");

  BFT_MALLOC(path, strlen(path_0) + 1, char);
  strcpy(path, path_0);

  /* Angle */

  cs_xpath_add_element(&path, "angle");
  cs_xpath_add_function_text(&path);
  if (!cs_gui_get_double(path, angle))
    *angle = 0.0;

  /* Axis */

  strcpy(path, path_0);
  cs_xpath_add_element(&path, "axis_x");
  coo_id = strlen(path) - 1;
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, axis + 0))
    axis[0] = 0.0;

  path[coo_id] = 'y';
  if (!cs_gui_get_double(path, axis + 1))
    axis[1] = 0.0;

  path[coo_id] = 'z';
  if (!cs_gui_get_double(path, axis + 2))
    axis[2] = 0.0;

  /* Invariant */

  strcpy(path, path_0);
  cs_xpath_add_element(&path, "invariant_x");
  coo_id = strlen(path) - 1;
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, invariant + 0))
    invariant[0] = 0.0;

  path[coo_id] = 'y';
  if (!cs_gui_get_double(path, invariant + 1))
    invariant[1] = 0.0;

  path[coo_id] = 'z';
  if (!cs_gui_get_double(path, invariant + 2))
    invariant[2] = 0.0;

  BFT_FREE(path);
  BFT_FREE(path_0);

#if _XML_DEBUG_
  bft_printf("==> _get_periodicity_translation\n");
  bft_printf("--angle = %f\n",
             *angle);
  bft_printf("--axis = [%f %f %f]\n",
             axis[0], axis[1], axis[2]);
  bft_printf("--invariant = [%f %f %f]\n",
             invariant[0], invariant[1], invariant[2]);
#endif

  return;
}

/*-----------------------------------------------------------------------------
 * Get transformation parameters associated with a mixed periodicity
 *
 * parameter:
 *   number <-- number of the syrthes coupling
 *   matrix --> translation values (m11, m12, m13, m14, ..., m33, m34)
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_mixed(int     number,
                       double  matrix[3][4])
{
  int i, j;
  size_t coeff_id = 0;
  char *path = cs_xpath_init_path();
  const char id_str[] = {'1', '2', '3','4'};

  cs_xpath_add_elements(&path, 2, "solution_domain", "periodicity");
  cs_xpath_add_element_num(&path, "face_periodicity", number);
  cs_xpath_add_elements(&path, 2, "mixed", "matrix_11");
  coeff_id = strlen(path) - 2;
  cs_xpath_add_function_text(&path);

  for (i = 0; i < 3; i++) {
    path[coeff_id] = id_str[i];

    for (j = 0; j < 4; j++) {
      path[coeff_id + 1] = id_str[j];

      if (!cs_gui_get_double(path, &(matrix[i][j]))) {
        if (i != j)
          matrix[i][j] = 0.0;
        else
          matrix[i][j] = 1.0;
      }
    }
  }

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf("==> _get_periodicity_translation\n");
  bft_printf("--matrix = [[%f %f %f %f]\n"
             "            [%f %f %f %f]\n"
             "            [%f %f %f %f]]\n",
             matrix[0][0], matrix[0][1] ,matrix[0][2], matrix[0][3],
             matrix[1][0], matrix[1][1] ,matrix[1][2], matrix[1][3],
             matrix[2][0], matrix[2][1] ,matrix[2][2], matrix[2][3]);
#endif

  return;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Determine whether warped faces should be cut.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_warping(void)
{
  char  *path = NULL;
  int cut_warped_faces = 0;
  double max_warp_angle = -1;

  if (!cs_gui_file_is_loaded())
    return;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "faces_cutting");
  cs_xpath_add_attribute(&path, "status");

  cs_gui_get_status(path, &cut_warped_faces);

  if (cut_warped_faces) {

    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3,
                          "solution_domain",
                          "faces_cutting",
                          "warp_angle_max");
    cs_xpath_add_function_text(&path);

    if (!cs_gui_get_double(path, &max_warp_angle))
      max_warp_angle = -1;

#if _XML_DEBUG_
  bft_printf("==> uicwf\n");
  bft_printf("--cut_warped_faces = %d\n"
             "--warp_angle_max   = %f\n",
             cut_warped_faces, max_warp_angle);
#endif

  }

  BFT_FREE(path);

  /* Apply warp angle options now */

  if (cut_warped_faces && max_warp_angle > 0.0)
    cs_mesh_warping_set_defaults(max_warp_angle, 0);
}

/*-----------------------------------------------------------------------------
 * Define joinings using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_joinings(void)
{
  int join_id;
  int n_join = 0;

  if (!cs_gui_file_is_loaded())
    return;

  n_join = cs_gui_get_tag_count("/solution_domain/joining/face_joining", 1);

  if (n_join == 0)
    return;

  for (join_id = 0; join_id < n_join; join_id++) {

    char *selector_s  =  _get_face_joining("selector", join_id+1);
    char *fraction_s  =  _get_face_joining("fraction", join_id+1);
    char *plane_s     =  _get_face_joining("plane", join_id+1);
    char *verbosity_s =  _get_face_joining("verbosity", join_id+1);
    char *visu_s      =  _get_face_joining("visualization", join_id+1);

    double fraction = (fraction_s != NULL) ? atof(fraction_s) : 0.1;
    double plane = (plane_s != NULL) ? atof(plane_s) : 25.0;
    int verbosity = (verbosity_s != NULL) ? atoi(verbosity_s) : 1;
    int visualization = (visu_s != NULL) ? atoi(visu_s) : 1;

    cs_join_add(selector_s,
                fraction,
                plane,
                verbosity,
                visualization);

#if _XML_DEBUG_
    bft_printf("==> cs_gui_mesh_define_joinings\n");
    bft_printf("--selector  = %s\n", selector_s);
    bft_printf("--fraction  = %s\n", fraction_s);
    bft_printf("--plane     = %s\n", plane_s);
    bft_printf("--verbosity = %s\n", verbosity_s);
    bft_printf("--visualization = %s\n", visu_s);
#endif

    BFT_FREE(selector_s);
    BFT_FREE(fraction_s);
    BFT_FREE(plane_s);
    BFT_FREE(verbosity_s);
    BFT_FREE(visu_s);
  }
}

/*-----------------------------------------------------------------------------
 * Define periodicities using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_periodicities(void)
{
  int perio_id;
  double angle, trans[3], axis[3], invariant[3], matrix[3][4];

  int  n_perio = 0;
  int  n_modes = 0;
  char  **modes = NULL;
  char  *path = NULL;

  if (!cs_gui_file_is_loaded())
    return;

  n_perio
    = cs_gui_get_tag_count("/solution_domain/periodicity/face_periodicity",
                           1);
  if (n_perio == 0)
    return;

  /* Get modes associated with each periodicity */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "solution_domain",
                        "periodicity", "face_periodicity");
  cs_xpath_add_attribute(&path, "mode");
  modes = cs_gui_get_attribute_values(path, &n_modes);

  if (n_modes != n_perio)
    bft_error(__FILE__, __LINE__, 0,
              _("Number of periodicities (%d) and modes (%d) do not match."),
              n_perio, n_modes);

  BFT_FREE(path);

  /* loop on periodicities */

  for (perio_id = 0; perio_id < n_perio; perio_id++) {

    char *selector_s  =  _get_periodicity("selector", perio_id+1);
    char *fraction_s  =  _get_periodicity("fraction", perio_id+1);
    char *plane_s     =  _get_periodicity("plane", perio_id+1);
    char *verbosity_s =  _get_periodicity("verbosity", perio_id+1);
    char *visu_s      =  _get_periodicity("visualization", perio_id+1);

    double fraction = (fraction_s != NULL) ? atof(fraction_s) : 0.1;
    double plane = (plane_s != NULL) ? atof(plane_s) : 25.0;
    int verbosity = (verbosity_s != NULL) ? atoi(verbosity_s) : 1;
    int visualization = (visu_s != NULL) ? atoi(visu_s) : 1;

    if (!strcmp(modes[perio_id], "translation")) {
      _get_periodicity_translation(perio_id+1, trans);
      cs_join_perio_add_translation(selector_s,
                                    fraction,
                                    plane,
                                    verbosity,
                                    visualization,
                                    trans);
    }

    else if (!strcmp(modes[perio_id], "rotation")) {
      _get_periodicity_rotation(perio_id+1, &angle, axis, invariant);
      cs_join_perio_add_rotation(selector_s,
                                 fraction,
                                 plane,
                                 verbosity,
                                 visualization,
                                 angle,
                                 axis,
                                 invariant);
    }

    else if (!strcmp(modes[perio_id], "mixed")) {
      _get_periodicity_mixed(perio_id+1, matrix);
      cs_join_perio_add_mixed(selector_s,
                              fraction,
                              plane,
                              verbosity,
                              visualization,
                              matrix);
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Periodicity mode \"%s\" unknown."), modes[perio_id]);

#if _XML_DEBUG_
    bft_printf("==> cs_gui_mesh_define_periodicities\n");
    bft_printf("--selector      = %s\n", selector_s);
    bft_printf("--fraction      = %s\n", fraction_s);
    bft_printf("--plane         = %s\n", plane_s);
    bft_printf("--verbosity     = %s\n", verbosity_s);
    bft_printf("--visualization = %s\n", visu_s);
#endif

    BFT_FREE(selector_s);
    BFT_FREE(fraction_s);
    BFT_FREE(plane_s);
    BFT_FREE(verbosity_s);
    BFT_FREE(visu_s);
  }

  for (perio_id = 0; perio_id < n_perio; perio_id++)
    BFT_FREE(modes[perio_id]);
  BFT_FREE(modes);
}

/*----------------------------------------------------------------------------
 * Mesh smoothing.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_smoothe(cs_mesh_t  *mesh)
{
  char  *path = NULL;
  int mesh_smooting = 0;
  double angle = 25.;

  if (!cs_gui_file_is_loaded())
    return;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "mesh_smoothing");
  cs_xpath_add_attribute(&path, "status");

  cs_gui_get_status(path, &mesh_smooting);

  if (mesh_smooting) {

    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3,
                          "solution_domain",
                          "mesh_smoothing",
                          "smooth_angle");
    cs_xpath_add_function_text(&path);

    if (!cs_gui_get_double(path, &angle))
      angle = 25.;

#if _XML_DEBUG_
  bft_printf("==> uicwf\n");
  bft_printf("--mesh_smoothing = %d\n"
             "--angle          = %f\n",
             mesh_smooting, angle);
#endif

    int *vtx_is_fixed = NULL;

    BFT_MALLOC(vtx_is_fixed, mesh->n_vertices, int);

    /* Get fixed boundary vertices flag */

    cs_mesh_smoother_fix_by_feature(mesh,
                                    angle,
                                    vtx_is_fixed);

    /* Call unwarping smoother */

    cs_mesh_smoother_unwarp(mesh, vtx_is_fixed);

    /* Free memory */

    BFT_FREE(vtx_is_fixed);
  }
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Define user thin wall through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_boundary(cs_mesh_t  *mesh)
{
  if (!cs_gui_file_is_loaded())
    return;

  int nwall = cs_gui_get_tag_count("/solution_domain/thin_walls/thin_wall", 1);

  if (nwall == 0)
    return;

  for (int iw = 0; iw < nwall; iw++) {
    char *path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "thin_walls");
    cs_xpath_add_element_num(&path, "thin_wall", iw + 1);
    cs_xpath_add_element(&path, "selector");
    cs_xpath_add_function_text(&path);
    char *value = cs_gui_get_text_value(path);
    BFT_FREE(path);

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;
    BFT_MALLOC(selected_faces, mesh->n_i_faces, cs_lnum_t);

    cs_selector_get_i_face_list(value,
                               &n_selected_faces,
                               selected_faces);

    cs_mesh_boundary_insert(mesh,
                            n_selected_faces,
                            selected_faces);

#if _XML_DEBUG_
    bft_printf("cs_gui_mesh_boundary==> \n");
    bft_printf("--selector  = %s\n", value);
#endif
    BFT_FREE(selected_faces);
    BFT_FREE(value);
  }
}

/*----------------------------------------------------------------------------
 * Define user mesh extrude through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_extrude(cs_mesh_t  *mesh)
{
  if (!cs_gui_file_is_loaded())
    return;

  int n_ext = cs_gui_get_tag_count("/solution_domain/extrusion/extrude_mesh", 1);

  if (n_ext == 0)
    return;

  for (int ext = 0; ext < n_ext; ext++) {
    char *path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "extrusion");
    cs_xpath_add_element_num(&path, "extrude_mesh", ext + 1);
    cs_xpath_add_element(&path, "selector");
    cs_xpath_add_function_text(&path);
    char *value = cs_gui_get_text_value(path);
    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "extrusion");
    cs_xpath_add_element_num(&path, "extrude_mesh", ext + 1);
    cs_xpath_add_element(&path, "layers_number");
    cs_xpath_add_function_text(&path);
    int n_layers;
    cs_gui_get_int(path, &n_layers);
    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "extrusion");
    cs_xpath_add_element_num(&path, "extrude_mesh", ext + 1);
    cs_xpath_add_element(&path, "thickness");
    cs_xpath_add_function_text(&path);
    double thickness;
    cs_gui_get_double(path, &thickness);
    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "extrusion");
    cs_xpath_add_element_num(&path, "extrude_mesh", ext + 1);
    cs_xpath_add_element(&path, "reason");
    cs_xpath_add_function_text(&path);
    double reason;
    cs_gui_get_double(path, &reason);
    BFT_FREE(path);

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;
    BFT_MALLOC(selected_faces, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(value,
                                &n_selected_faces,
                                selected_faces);

    cs_mesh_extrude_constant(mesh,
                             true,
                             n_layers,
                             thickness,
                             reason,
                             n_selected_faces,
                             selected_faces);

#if _XML_DEBUG_
    bft_printf("cs_gui_mesh_extrude==> \n");
    bft_printf("--selector  = %s\n", value);
    bft_printf("--n_layers  = %i\n", n_layers);
    bft_printf("--thickness = %f\n", thickness);
    bft_printf("--reason    = %f\n", reason);
#endif
    BFT_FREE(selected_faces);
    BFT_FREE(value);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
