/*============================================================================
 * Management of the GUI parameters file: mesh related options
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
#include "cs_gui_boundary_conditions.h"
#include "cs_join.h"
#include "cs_join_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_cartesian.h"
#include "cs_mesh_warping.h"
#include "cs_mesh_smoother.h"
#include "cs_mesh_boundary.h"
#include "cs_mesh_extrude.h"
#include "cs_parameters.h"

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
 * Get transformation parameters associated with a translational periodicity
 *
 * parameter:
 *   node   <-- pointer to a periodic faces definition node
 *   trans  --> translation values
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_translation(cs_tree_node_t  *node,
                             double           trans[3])
{
  cs_tree_node_t  *tn = cs_tree_node_get_child(node, "translation");

  if (tn != NULL) {

    const char *names[] = {"translation_x", "translation_y", "translation_z"};

    for (int i = 0; i < 3; i++) {
      const cs_real_t *v = cs_tree_node_get_child_values_real(tn, names[i]);
      if (v != NULL)
        trans[i] = v[0];
    }

  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--translation = [%f %f %f]\n",
             trans[0], trans[1], trans[2]);
#endif
}

/*-----------------------------------------------------------------------------
 * Get transformation parameters associated with a rotational periodicity
 *
 * parameter:
 *   node      <-- pointer to a periodic faces definition node
 *   angle     --> rotation angle
 *   axis      --> rotation axis
 *   invariant --> invariant point
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_rotation(cs_tree_node_t  *node,
                          double          *angle,
                          double           axis[3],
                          double           invariant[3])
{
  cs_tree_node_t  *tn = cs_tree_node_get_child(node, "rotation");

  if (tn != NULL) {

    const cs_real_t *v;

    /* Angle */

    v = cs_tree_node_get_child_values_real(tn, "angle");
    if (v != NULL)
      *angle =  v[0];
    else
      *angle = 0.0;

    /* Axis */

    const char *a_names[] = {"axis_x", "axis_y", "axis_z"};

    for (int i = 0; i < 3; i++) {
      v = cs_tree_node_get_child_values_real(tn, a_names[i]);
      if (v != NULL)
        axis[i] = v[0];
      else
        axis[i] = 0;
    }

    /* Invariant */

    const char *i_names[] = {"invariant_x", "invariant_y", "invariant_z"};

    for (int i = 0; i < 3; i++) {
      v = cs_tree_node_get_child_values_real(tn, i_names[i]);
      if (v != NULL)
        invariant[i] = v[0];
      else
        invariant[i] = 0;
    }

  }
  else
    *angle = 0.0;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--angle = %f\n",
             *angle);
  bft_printf("--axis = [%f %f %f]\n",
             axis[0], axis[1], axis[2]);
  bft_printf("--invariant = [%f %f %f]\n",
             invariant[0], invariant[1], invariant[2]);
#endif
}

/*-----------------------------------------------------------------------------
 * Get transformation parameters associated with a mixed periodicity
 *
 * parameter:
 *   node   <-- pointer to a periodic faces definition node
 *   matrix --> translation values (m11, m12, m13, m14, ..., m33, m34)
 *----------------------------------------------------------------------------*/

static void
_get_periodicity_mixed(cs_tree_node_t  *node,
                       double           matrix[3][4])
{
  cs_tree_node_t  *tn = cs_tree_node_get_child(node, "mixed");

  if (tn != NULL) {

    char c_name[] = "matrix_11";

    size_t coeff_id = strlen("matrix_");
    const char id_str[] = {'1', '2', '3','4'};

    for (int i = 0; i < 3; i++) {
      c_name[coeff_id] = id_str[i];

      for (int j = 0; j < 4; j++) {
        c_name[coeff_id + 1] = id_str[j];

        const cs_real_t *v = cs_tree_node_get_child_values_real(tn, c_name);
        if (v != NULL)
          matrix[i][j] = v[0];
        else {
          if (i != j)
            matrix[i][j] = 0.0;
          else
            matrix[i][j] = 1.0;
        }

      }
    }

  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--matrix = [[%f %f %f %f]\n"
             "            [%f %f %f %f]\n"
             "            [%f %f %f %f]]\n",
             matrix[0][0], matrix[0][1] ,matrix[0][2], matrix[0][3],
             matrix[1][0], matrix[1][1] ,matrix[1][2], matrix[1][3],
             matrix[2][0], matrix[2][1] ,matrix[2][2], matrix[2][3]);
#endif
}

/*----------------------------------------------------------------------------
 * Get cartesian mesh parameters defined with GUI.
 *----------------------------------------------------------------------------*/

static void
_get_cartesian_parameters(int        idim,
                          int        *ip,
                          cs_real_t  *rp)
{
  cs_tree_node_t *tn0
    = cs_tree_get_node(cs_glob_tree,"solution_domain/mesh_cartesian");

  if (tn0 != NULL) {

    cs_tree_node_t *tn = NULL;
    if (idim == 0)
      tn = cs_tree_node_get_child(tn0, "x_direction");
    else if (idim == 1)
      tn = cs_tree_node_get_child(tn0, "y_direction");
    else if (idim == 2)
      tn = cs_tree_node_get_child(tn0, "z_direction");

    const char *law = cs_gui_node_get_tag(tn, "law");
    if (strcmp(law, "constant") == 0)
      ip[0] = 0;
    else if (strcmp(law, "geometric") == 0)
      ip[0] = 1;
    else if (strcmp(law, "parabolic") == 0)
      ip[0] = 2;

    cs_gui_node_get_child_int(tn, "ncells", ip + 1);

    cs_gui_node_get_child_real(tn, "min",  rp);
    cs_gui_node_get_child_real(tn, "max",  rp + 1);
    cs_gui_node_get_child_real(tn, "prog", rp + 2);

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Error: There is no cartesian mesh defined by the XML file.\n"));
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
  const char path0[] = "solution_domain/faces_cutting";

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);

  if (tn == NULL)
    return;

  int cut_warped_faces = false;
  cs_gui_node_get_status_int(tn, &cut_warped_faces);

  if (cut_warped_faces) {

    const cs_real_t *v_r = cs_tree_node_get_child_values_real(tn,
                                                              "warp_angle_max");
    double max_warp_angle = (v_r != NULL) ? v_r[0] : -1;

    /* Apply warp angle options now */

    if (max_warp_angle > 0)
      cs_mesh_warping_set_defaults(max_warp_angle, 0);

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--cut_warped_faces = %d\n"
               "--warp_angle_max   = %f\n",
               cut_warped_faces, max_warp_angle);
#endif

  }
}

/*-----------------------------------------------------------------------------
 * Define joinings using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_joinings(void)
{
  const int *v_i = NULL;
  const cs_real_t *v_r = NULL;

  /* Loop on joining definitions */

  const char path_j[] = "solution_domain/joining/face_joining";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path_j);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char _default_criteria[] = "all[]";

    const char *selector = cs_tree_node_get_child_value_str(tn, "selector");
    if (selector == NULL) selector = _default_criteria;

    v_r = cs_tree_node_get_child_values_real(tn, "fraction");
    double fraction = (v_r != NULL) ? v_r[0] : 0.1;

    v_r = cs_tree_node_get_child_values_real(tn, "plane");
    double plane = (v_r != NULL) ? v_r[0] : 25.0;

    v_i = cs_tree_node_get_child_values_int(tn, "verbosity");
    int verbosity = (v_i != NULL) ? v_i[0] : 1;

    v_i = cs_tree_node_get_child_values_int(tn, "visualization");
    int visualization = (v_i != NULL) ? v_i[0] : 1;

    cs_join_add(selector,
                fraction,
                plane,
                verbosity,
                visualization);

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--selector  = %s\n", selector);
    bft_printf("--fraction  = %s\n", fraction);
    bft_printf("--plane     = %g\n", plane);
    bft_printf("--verbosity = %d\n", verbosity);
    bft_printf("--visualization = %d\n", visualization);
#endif
  }
}

/*-----------------------------------------------------------------------------
 * Define periodicities using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_periodicities(void)
{
  const int *v_i = NULL;
  const cs_real_t *v_r = NULL;

  /* Loop on periodicity definitions */

  int perio_id = 0;

  const char path_p[] = "solution_domain//periodicity/face_periodicity";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path_p);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), perio_id++) {

    double angle, trans[3], axis[3], invariant[3], matrix[3][4];

    /* Get mode associated with each periodicity */

    const char *mode = cs_tree_node_get_tag(tn, "mode");

    if (mode == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("\"%s\" node %d is missing a \"%s\" tag/child."),
                tn->name, perio_id, "mode");

    const char _default_criteria[] = "all[]";

    const char *selector = cs_tree_node_get_child_value_str(tn, "selector");
    if (selector == NULL) selector = _default_criteria;

    v_r = cs_tree_node_get_child_values_real(tn, "fraction");
    double fraction = (v_r != NULL) ? v_r[0] : 0.1;

    v_r = cs_tree_node_get_child_values_real(tn, "plane");
    double plane = (v_r != NULL) ? v_r[0] : 25.0;

    v_i = cs_tree_node_get_child_values_int(tn, "verbosity");
    int verbosity = (v_i != NULL) ? v_i[0] : 1;

    v_i = cs_tree_node_get_child_values_int(tn, "visualization");
    int visualization = (v_i != NULL) ? v_i[0] : 1;

    if (!strcmp(mode, "translation")) {
      _get_periodicity_translation(tn, trans);
      cs_join_perio_add_translation(selector,
                                    fraction,
                                    plane,
                                    verbosity,
                                    visualization,
                                    trans);
    }

    else if (!strcmp(mode, "rotation")) {
      _get_periodicity_rotation(tn, &angle, axis, invariant);
      cs_join_perio_add_rotation(selector,
                                 fraction,
                                 plane,
                                 verbosity,
                                 visualization,
                                 angle,
                                 axis,
                                 invariant);
    }

    else if (!strcmp(mode, "mixed")) {
      _get_periodicity_mixed(tn, matrix);
      cs_join_perio_add_mixed(selector,
                              fraction,
                              plane,
                              verbosity,
                              visualization,
                              matrix);
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Periodicity mode \"%s\" unknown."), mode);

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--selector      = %s\n", selector);
    bft_printf("--fraction      = %g\n", fraction);
    bft_printf("--plane         = %g\n", plane);
    bft_printf("--verbosity     = %d\n", verbosity);
    bft_printf("--visualization = %d\n", visu);
#endif
  }
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
  const char path0[] = "solution_domain/mesh_smoothing";

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);

  if (tn == NULL)
    return;

  int mesh_smoothing = false;
  cs_gui_node_get_status_int(tn, &mesh_smoothing);

  if (mesh_smoothing) {

    const cs_real_t *v_r = cs_tree_node_get_child_values_real(tn,
                                                              "smooth_angle");
    double angle = (v_r != NULL) ? v_r[0] : 25;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--mesh_smoothing = %d\n"
             "--angle          = %f\n",
             mesh_smoothing, angle);
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
}

/*----------------------------------------------------------------------------
 * Define user thin wall through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_boundary(cs_mesh_t  *mesh)
{
  const char path0[] = "/solution_domain/thin_walls/thin_wall";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char _default_criteria[] = "all[]";

    const char *selector = cs_tree_node_get_child_value_str(tn, "selector");
    if (selector == NULL) selector = _default_criteria;

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--selector  = %s\n", value);
#endif

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;
    BFT_MALLOC(selected_faces, mesh->n_i_faces, cs_lnum_t);

    cs_selector_get_i_face_list(selector,
                                &n_selected_faces,
                                selected_faces);

    cs_mesh_boundary_insert(mesh,
                            n_selected_faces,
                            selected_faces);

    BFT_FREE(selected_faces);

  }
}

/*----------------------------------------------------------------------------
 * Define user mesh extrude through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_extrude(cs_mesh_t  *mesh)
{
  const int *v_i = NULL;
  const cs_real_t *v_r = NULL;

  const char path0[] = "solution_domain/extrusion/extrude_mesh";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char _default_criteria[] = "all[]";

    const char *selector = cs_tree_node_get_child_value_str(tn, "selector");
    if (selector == NULL) selector = _default_criteria;

    v_i = cs_tree_node_get_child_values_int(tn, "layers_number");
    int n_layers = (v_i != NULL) ? v_i[0] : 2;

    v_r = cs_tree_node_get_child_values_real(tn, "thickness");
    double thickness = (v_r != NULL) ? v_r[0] : 1;

    v_r = cs_tree_node_get_child_values_real(tn, "reason");
    double reason = (v_r != NULL) ? v_r[0] : 1.5;

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;
    BFT_MALLOC(selected_faces, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(selector,
                                &n_selected_faces,
                                selected_faces);

    cs_mesh_extrude_constant(mesh,
                             true,
                             n_layers,
                             thickness,
                             reason,
                             n_selected_faces,
                             selected_faces);

    BFT_FREE(selected_faces);

#if _XML_DEBUG_
    bft_printf("%s ==>\n", __func__);
    bft_printf("--selector  = %s\n", value);
    bft_printf("--n_layers  = %i\n", n_layers);
    bft_printf("--thickness = %f\n", thickness);
    bft_printf("--reason    = %f\n", reason);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Define mesh save behavior trough the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_save_if_modified(cs_mesh_t  *mesh)
{
  const char path0[] = "solution_domain/save_mesh_if_modified";

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);

  if (tn == NULL)
    return;

  const char *s = cs_tree_node_get_value_str(tn);
  if (s != NULL) {

    if (!strcmp(s, "no"))
      mesh->save_if_modified = 0;
    else if (!strcmp(s, "yes"))
      mesh->save_if_modified = 1;

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--save_mesh_if_modified = %s\n", s);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Define if cartesian mesh is to be built through GUI.
 *----------------------------------------------------------------------------*/

int
cs_gui_mesh_build_cartesian(void)
{
  int retval = 0;
  cs_tree_node_t *tn =
    cs_tree_get_node(cs_glob_tree,"solution_domain/mesh_origin");

  const char *choice = cs_tree_node_get_child_value_str(tn, "choice");
  if (choice != NULL) {
    if (strcmp(choice, "mesh_cartesian") == 0)
      retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Read cartesian mesh parameters defined with GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_cartesian_define(void)
{
  cs_mesh_cartesian_create();
  for (int idim = 0; idim < 3; idim++) {
    int       iparams[2] = {0, 0};
    cs_real_t rparams[3] = {0., 0., 0.};
    _get_cartesian_parameters(idim, iparams, rparams);

    cs_mesh_cartesian_law_t law = CS_MESH_CARTESIAN_CONSTANT_LAW;
    if (iparams[0] == 0)
      law = CS_MESH_CARTESIAN_CONSTANT_LAW;
    else if (iparams[0] == 1)
      law = CS_MESH_CARTESIAN_GEOMETRIC_LAW;
    else if (iparams[0] == 2)
      law = CS_MESH_CARTESIAN_PARABOLIC_LAW;

    cs_mesh_cartesian_define_dir_params(idim,
                                        law,
                                        iparams[1],
                                        rparams[0],
                                        rparams[1],
                                        rparams[2]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
