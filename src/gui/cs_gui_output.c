/*============================================================================
 * Management of the GUI parameters file: output parameters
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "fvm_selector.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_gui.h"
#include "cs_gui_util.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_meg_prototypes.h"
#include "cs_selector.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_function_default.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_output.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * External global variables
 *============================================================================*/

/*============================================================================
 * Private global variables
 *============================================================================*/

/*============================================================================
 * Static local variables
 *============================================================================*/

static char _rij_c_names[6][4] = {"r11", "r22", "r33", "r12", "r23", "r13"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get of node of the given type with a given name tag.
 *
 * parameters:
 *   name    <--  node type (variable, property, ...)
 *   name    <--  associated name tag (sub-node)
 *
 * return:
 *   pointer to node found, or NULL
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_get_node(const char  *node_type,
          const char  *name)
{
  cs_tree_node_t *root = cs_glob_tree;

  for (cs_tree_node_t *tn = cs_tree_find_node_simple(root, node_type);
       tn != NULL;
       tn = cs_tree_find_node_next_simple(root, tn, node_type)) {
    const char *tag = cs_tree_node_get_tag(tn, "name");
    if (tag != NULL) {
      if (strcmp(tag, name) == 0)
        return tn;
    }
  }

  return NULL;
}

/*-----------------------------------------------------------------------------
 * Post-processing options for fields
 *
 * parameters:
 *   f_id <-- field id
 *----------------------------------------------------------------------------*/

static void
_field_post(const char  *field_type,
            int          f_id)
{
  cs_field_t  *f = cs_field_by_id(f_id);

  /* Now check for options */

  int f_post = -999, f_log = -999, f_monitor = -999;
  const int k_log  = cs_field_key_id("log");
  const int k_lbl = cs_field_key_id("label");
  const int k_post = cs_field_key_id("post_vis");

  cs_tree_node_t *tn = _get_node(field_type, f->name);

  if (tn == NULL)
    return;

  /* Listing output */

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "listing_printing"),
                             &f_log);
  if (f_log != -999)
    cs_field_set_key_int(f, k_log, f_log);

  /* As long as boundary output control is not unified, we need
     to handle them in a specific manner (as keywors may have
     been set by the GUI for fields not listed with the commmon
     hierarchy and "postprocessing_recording" and similar tag) */

  bool allow_default_set = true;
  if (   f->location_id == CS_MESH_LOCATION_BOUNDARY_FACES
      && cs_field_is_key_set(f, k_post))
    allow_default_set = false;

  /* Postprocessing outputs */

  cs_gui_node_get_status_int(cs_tree_node_get_child
                               (tn, "postprocessing_recording"),
                             &f_post);
  if (f_post == 1)
    cs_field_set_key_int_bits(f, k_post, CS_POST_ON_LOCATION);
  else if (f_post == 0)
    cs_field_clear_key_int_bits(f, k_post, CS_POST_ON_LOCATION);
  else if (allow_default_set) { /* status unspecified here but property
                                   referenced in tree, could be improved by
                                   depending on field type or flags */
    cs_field_set_key_int_bits(f, k_post, CS_POST_ON_LOCATION);
  }

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "probes_recording"),
                             &f_monitor);
  if (f_monitor == 1)
    cs_field_set_key_int_bits(f, k_post, CS_POST_MONITOR);
  else if (f_monitor == 0)
    cs_field_clear_key_int_bits(f, k_post, CS_POST_MONITOR);
  else if (allow_default_set) { /* status unspecified here but property
                                   referenced in tree, could be improved by
                                   depending on field type or flags */
    if (f->location_id == CS_MESH_LOCATION_CELLS)
      cs_field_set_key_int_bits(f, k_post, CS_POST_MONITOR);
    else
      cs_field_clear_key_int_bits(f, k_post, CS_POST_MONITOR);
  }

  /* Take into account labels */

  const char *label = cs_tree_node_get_tag(tn, "label");
  if (label != NULL)
    cs_field_set_key_str(f, k_lbl, label);
}

/*-----------------------------------------------------------------------------
 * Post-processing options for user defind calculator functions (GUI)
 *
 * parameters:
 *   f_id <-- field id
 *----------------------------------------------------------------------------*/

static void
_function_post(const int id)
{
  cs_function_t *f = cs_function_by_id(id);

  const char path_s[] = "user_functions/calculator/function";
  cs_tree_node_t *tn_s = cs_tree_get_node(cs_glob_tree, path_s);

  cs_tree_node_t *_tn = NULL;
  for (cs_tree_node_t *tn = tn_s;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    const char *name = cs_gui_node_get_tag(tn, "name");
    if (strcmp(name, f->name) == 0) {
      _tn = tn;
      break;
    }
  }

  if (_tn == NULL)
    return;

  /* By default, if function exists in GUI, use parameters defined in it */
  int f_post   = 1;
  int f_log    = 1;
  int f_probes = 1;

  cs_gui_node_get_status_int(cs_tree_node_get_child(_tn, "postprocessing_recording"),
                             &f_post);
  f->post_vis |= CS_POST_ON_LOCATION;
  if (f_post == 0)
    f->post_vis -= CS_POST_ON_LOCATION;

  cs_gui_node_get_status_int(cs_tree_node_get_child(_tn, "probes_recording"),
                             &f_probes);
  f->post_vis |= CS_POST_MONITOR;
  if (f_probes == 0)
    f->post_vis -= CS_POST_MONITOR;

  cs_gui_node_get_status_int(cs_tree_node_get_child(_tn, "listing_printing"),
                             &f_log);
  f->log = f_log;
}

/*----------------------------------------------------------------------------
 * Get postprocessing value status for surfacic variables
 *
 * parameters:
 *   name         --> name of the parameter
 *   default_val  --> default value
 *----------------------------------------------------------------------------*/

static bool
_surfacic_variable_post(const char  *name,
                        bool         default_val)
{
  bool active = default_val;

  cs_tree_node_t *tn = _get_node("property", name);

  if (tn != NULL) {

    /* If base node present but not recording status, default to true */
    active = true;

    cs_gui_node_get_status_bool(cs_tree_node_get_child
                                 (tn, "postprocessing_recording"),
                                &active);

  }

  return active;
}

/*----------------------------------------------------------------------------
 * Return a string value associated with a child node for a profile.
 *
 * If the matching child node is not present, an error is produced
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

static const char *
_get_profile_child_str(cs_tree_node_t  *tn,
                       const char      *child_name)
{
  const char *name = cs_tree_node_get_child_value_str(tn, child_name);

  if (name == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Incorrect setup tree definition for the following node:\n"));
    cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
    bft_error(__FILE__, __LINE__, 0,
              _("Missing child node: %s"), child_name);
  }

  return name;
}

/*----------------------------------------------------------------------------
 * Return an array of integers associated with a child node for a profile.
 *
 * If the matching child node is not present, an error is produced
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   pointer to matching child values
 *----------------------------------------------------------------------------*/

static const int *
_get_profile_child_int(cs_tree_node_t  *tn,
                       const char      *child_name)
{
  const int *v = cs_tree_node_get_child_values_int(tn, child_name);

  if (v == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Incorrect setup tree definition for the following node:\n"));
    cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
    bft_error(__FILE__, __LINE__, 0,
              _("Missing child node: %s"), child_name);
  }

  return v;
}

/*----------------------------------------------------------------------------
 * Get output format for 1D profile
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   1 for CSV, 0 for DAT
 *----------------------------------------------------------------------------*/

static int
_get_profile_format(cs_tree_node_t  *tn)
{
  int format = 0;

  const char *format_s
    = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "format"),
                           "name");

  if (format_s != NULL) {
    if (cs_gui_strcmp(format_s, "CSV"))
      format = 1;
    else if (cs_gui_strcmp(format_s, "DAT"))
      format = 0;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid profile format: %s"), format_s);
  }

  return format;
}

/*----------------------------------------------------------------------------
 * Return the component of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   tn <-- tree node associated with profile variable
 *----------------------------------------------------------------------------*/

static int
_get_profile_v_component(cs_tree_node_t  *tn)
{
  int comp_id = -1;

  const char *c_name = cs_tree_node_get_tag(tn, "component");

  if (c_name != NULL) {
    int n = sscanf(c_name, "%d", &comp_id);
    if (n != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error converting profile component tag %s to integer."),
                c_name);
  }

  return comp_id;
}

/*----------------------------------------------------------------------------
 * Return associated field for a node referring to a field.
 *
 * A node referring to a field must have a "name" tag (i.e. child with a
 * string value).
 *
 * In the case of NEPTUNE_CFD, an additional "field_id" tag requiring
 * addition of the "_<field_id>" extension to match the field's name
 * may be present, so this is tested also.
 *
 * If the category and associated name are defined and a field matching
 * a tag is not present, a message is logged. If those are not defined
 * and a matching field is not found, an error is produced.
 *
 * parameters:
 *   tn            <-- tree node associated with profile variable
 *   category      <-- category for which this function is used, or NULL
 *   category_name <-- name of object in category
 *
 * return:
 *   pointer to field if match, NULL otherwise
 *----------------------------------------------------------------------------*/

static const cs_field_t *
_tree_node_get_field(cs_tree_node_t  *tn,
                     const char      *category,
                     const char      *category_name)
{
  const cs_field_t *f = NULL;

  const char *name = cs_gui_node_get_tag(tn, "name");
  const char *id_name = cs_tree_node_get_tag(tn, "field_id");

  /* Handle phase id (field_id tag in xml) for NEPTUNE_CFD */

  if (id_name != NULL) {
    if (strcmp(id_name, "none") != 0) {
      char buffer[128];
      snprintf(buffer, 127, "%s_%s", name, id_name);
      buffer[127] = '\0';
      if (strlen(buffer) >= 127)
        bft_error(__FILE__, __LINE__, 0,
                  "Local buffer too small to assemble field name with:\n"
                  "name: %s\n"
                  "field_id: %s\n", name, id_name);
      f = cs_field_by_name_try(buffer);
    }
  }

  if (f == NULL) {

    /* Fix time step output */
    if (strcmp(name, "local_time_step") == 0)
      f = CS_F_(dt);

    /* General case */
    else
      f = cs_field_by_name_try(name);

    /* Handle segregated Reynolds tensor solver */
    if (f == NULL) {
      if (strcmp(name, "rij") == 0) {
        int idim = _get_profile_v_component(tn);
        f = cs_field_by_name_try(_rij_c_names[idim]);
      }
    }

  }

  if (f == NULL) {
    if (category != NULL && category_name != NULL)
      bft_printf(_("  For %s \"%s\", field with name \"%s\" not found\n"),
                 category, category_name,  name);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field with name \"%s\" not found"), name);
  }

  return f;
}

/*----------------------------------------------------------------------------
 * Profile definitions
 *----------------------------------------------------------------------------*/

static void
_define_profiles(void)
{
  /* Loop on 1D profile definitions */

  int profile_id = 0;

  const char path0[] = "analysis_control/profiles/profile";

  int n_writers = 0;

  int        *w_i_vals = NULL;
  cs_real_t  *w_r_vals = NULL;

  int writer_id_start = cs_post_get_free_writer_id();

  const char *format_name[] = {"dat", "csv"};

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), profile_id++) {

    const char *name = cs_gui_node_get_tag(tn, "label");

    /* for each profile, check the output frequency */

    int writer_id = 0;
    int output_format = _get_profile_format(tn);
    int output_at_end = 0;
    int output_frequency = -1;
    cs_real_t time_output = -1.;

    const char *output_type = _get_profile_child_str(tn, "output_type");
    if (cs_gui_strcmp(output_type, "time_value")) {
      const cs_real_t *v
        = cs_tree_node_get_child_values_real(tn, "output_frequency");
      if (v != NULL)
        time_output = v[0];
    }
    else if (cs_gui_strcmp(output_type, "frequency")) {
      const int *v
        = cs_tree_node_get_child_values_int(tn, "output_frequency");
      if (v != NULL)
        output_frequency = v[0];
      else
        output_frequency = 1;
      output_at_end = true;  /* Debatable, but consistent with v6.0 behavior */
    }
    else if (cs_gui_strcmp(output_type, "end")) {
      output_at_end = true;
    }

    /* Check if the output frequency matches an existing writer id;
       create a writer if necessary */

    for (int i = 0; i < n_writers; i++) {
      if (   w_i_vals[i*3]   == output_format
          && w_i_vals[i*3+1] == output_at_end
          && w_i_vals[i*3+2] == output_frequency
          && fabs(w_r_vals[i] - time_output) < cs_math_epzero) {
        writer_id = writer_id_start - i;
        break;
      }
    }

    if (writer_id == 0) {  /* Add new writer */

      int i = n_writers;

      n_writers += 1;
      BFT_REALLOC(w_i_vals, n_writers*3, int);
      BFT_REALLOC(w_r_vals, n_writers, cs_real_t);
      w_i_vals[i*3] = output_format;
      w_i_vals[i*3+1] = output_at_end;
      w_i_vals[i*3+2] = output_frequency;
      w_r_vals[i] = time_output;
      writer_id = writer_id_start - i;

      bool _output_at_end = (output_at_end) ? true : false;

      char format_options[64];
      strncpy(format_options, format_name[output_format], 63);
      format_options[63] = '\0';
      if (_output_at_end && output_frequency < 0 && time_output < 0) {
        size_t l = strlen(format_options);
        strncat(format_options, ", no_time_step", 63-l);
        format_options[63] = '\0';
      }

      cs_post_define_writer(writer_id,
                            "",
                            "profiles",
                            "plot",
                            format_options,
                            FVM_WRITER_FIXED_MESH,
                            false, /* output_at_start */
                            _output_at_end,
                            output_frequency,
                            time_output);

    }

    int n_coords = 0;
    const int *v_np = _get_profile_child_int(tn, "points");
    if (v_np != NULL)
      n_coords = v_np[0];

    cs_real_3_t *coords;
    BFT_MALLOC(coords, n_coords, cs_real_3_t);

    /* Be debugger-friendly in case cs_meg_post_profiles does not handle
       this well (if generation is wrong or missing) */
    for (cs_lnum_t i = 0; i < n_coords; i++) {
      coords[i][0] = -1.e30;
      coords[i][1] = -1.e30;
      coords[i][2] = -1.e30;
    }

    cs_meg_post_profiles(name, n_coords, coords);

    cs_probe_set_t *pset
      = cs_probe_set_create_from_array(name,
                                       n_coords,
                                       coords,
                                       NULL);

    BFT_FREE(coords);

    cs_probe_set_assign_curvilinear_abscissa(pset, NULL);

    int writer_ids[1] = {writer_id};
    cs_probe_set_associate_writers(pset, 1, writer_ids);

    cs_probe_set_auto_var(pset, false);

    /* Set snap mode. Default is "SNAP_TO_CENTER" */
    const char *snap_mode = cs_tree_node_get_child_value_str(tn, "snap_mode");
    if (cs_gui_strcmp(snap_mode, "snap_to_vertex"))
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_VERTEX);
    else if (cs_gui_strcmp(snap_mode, "none"))
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_NONE);
    else
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_ELT_CENTER);

    /* Activate interpolation if needed. Default is no */
    const char *activate_interpolation
      = cs_tree_node_get_child_value_str(tn, "interpolation");
    if (cs_gui_strcmp(activate_interpolation, "yes"))
      cs_probe_set_option(pset, "interpolation", "1");

    if (cs_glob_mesh->time_dep > CS_MESH_FIXED)
      cs_probe_set_option(pset, "transient_location", "true");

    /* Output coordinates */

    cs_probe_set_auto_curvilinear_coords(pset, true);
    cs_probe_set_auto_cartesian_coords(pset, true);

    /* Associate fields */

    for (cs_tree_node_t *tn_vp = cs_tree_node_get_child(tn, "var_prop");
         tn_vp != NULL;
         tn_vp = cs_tree_node_get_next_of_name(tn_vp)) {

      const cs_field_t *f = _tree_node_get_field(tn_vp, "profile", name);

      if (f == NULL)
        continue;

      int comp_id = _get_profile_v_component(tn_vp);

      cs_probe_set_associate_field(pset,
                                   writer_id,
                                   f->id,
                                   comp_id);

    }

  }

  BFT_FREE(w_i_vals);
  BFT_FREE(w_r_vals);
}

/*----------------------------------------------------------------------------
 * Boundary zone cells selection wrapper function
 *----------------------------------------------------------------------------*/

static void
_selection_func_boundary_cells(void        *input,
                               cs_lnum_t   *n_elts,
                               cs_lnum_t  **elt_list)
{
  const char *criteria = (const char *)input;

  if (criteria == NULL)
    criteria = "all[]";

  /* Pointer is not allocated when given to this function.
   * It will be deallocated afterwards by the calling function. */

  cs_lnum_t *_cell_list = NULL;
  BFT_MALLOC(_cell_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_cells_list(criteria, n_elts, _cell_list);

  *elt_list = _cell_list;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine general output options.
 *----------------------------------------------------------------------------*/

void
cs_gui_output(void)
{
  const int *v_i = NULL;

  const char path_o[] = "analysis_control/output";
  cs_tree_node_t *tn_o = cs_tree_get_node(cs_glob_tree, path_o);

  v_i = cs_tree_node_get_child_values_int(tn_o,
                                          "listing_printing_frequency");
  if (v_i != NULL) cs_glob_log_frequency = v_i[0];

  const int n_fields = cs_field_n_fields();

  /* temporary field -> moment ids */
  int *moment_id = NULL;
  const int n_moments = cs_time_moment_n_moments();
  if (n_moments > 0) {
    BFT_MALLOC(moment_id, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      moment_id[f_id] = -1;
    for (int m_id = 0; m_id < n_moments; m_id++) {
      const cs_field_t *f = cs_time_moment_get_field(m_id);
      if (f != NULL)
        moment_id[f->id] = m_id;
    }
  }

  /* variable output */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE)
      _field_post("variable", f->id);
    else if (   (f->type & CS_FIELD_PROPERTY)
             || (f->type & CS_FIELD_POSTPROCESS)) {
      if (moment_id != NULL) {
        if (moment_id[f_id] > -1) {
          _field_post("time_average", f->id);
          continue;
        }
      }
      _field_post("property", f->id);
    }
  }

  BFT_FREE(moment_id);

  /* User functions */
  for (int f_id = 0; f_id < cs_function_n_functions(); f_id++) {
    _function_post(f_id);
  }

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
  bft_printf("--ntlist = %i\n", cs_glob_log_frequency);
#endif
}

/*----------------------------------------------------------------------------
 * Determine output boundary fields
 *----------------------------------------------------------------------------*/

void
cs_gui_output_boundary(void)
{
  int k_vis = cs_field_key_id("post_vis");

  /* Surfacic variables output */

  bool ignore_stresses = false;
  bool ignore_yplus = false;

  if (   cs_glob_physical_model_flag[CS_GROUNDWATER] != -1
      || cs_glob_physical_model_flag[CS_SOLIDIFICATION] != -1) {
    ignore_stresses = true;
    ignore_yplus = true;
  }

  /* CDO-based velocity-pressure does not currently involve
     yplus and stresses computation. */

  cs_field_t  *f_vel = cs_field_by_name_try("velocity");
  if (f_vel != NULL) {
    if (f_vel->type & CS_FIELD_CDO) {
      ignore_stresses = true;
      ignore_yplus = true;
    }
  }
  else {
    ignore_stresses = true;
    ignore_yplus = true;
  }

  if (ignore_stresses == false) {
    if (_surfacic_variable_post("stress", true))
      cs_function_define_boundary_stress();
    if (_surfacic_variable_post("stress_tangential", false))
      cs_function_define_boundary_stress_tangential();
    if (_surfacic_variable_post("stress_normal", false))
      cs_function_define_boundary_stress_normal();
  }

    /* TODO: move this following field and function definitions earlier
       (with thermal model), and only handle "post_vis" option here,
       to also allow for logging */

  if (ignore_yplus == false) {
    if (_surfacic_variable_post("yplus", true)) {
      cs_field_t *bf = cs_field_by_name_try("yplus");
      if (bf == NULL) {
        bf = cs_field_create("yplus",
                             CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             1,
                             false);
        cs_field_set_key_int(bf, cs_field_key_id("log"), 1);
        cs_field_set_key_str(bf, cs_field_key_id("label"), "Yplus");
      }
      cs_field_set_key_int(bf, k_vis, 1);
    }

    if (_surfacic_variable_post("thermal_flux", true)) {
      cs_function_define_boundary_thermal_flux();
    }

    if (_surfacic_variable_post("boundary_layer_nusselt", false)) {
      cs_function_define_boundary_nusselt();
    }

    /* Boundary T+/H+/E+ */

    if (cs_glob_thermal_model->itherm != CS_THERMAL_MODEL_NONE) {
      int post_vis = _surfacic_variable_post("tplus", true) ? 1 : 0;
      cs_field_t *bf = cs_field_by_name_try("tplus");
      if (bf != NULL)
        cs_field_set_key_int(bf, k_vis, post_vis);
      else if (post_vis) {
        bf = cs_field_create("tplus",
                             CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             1,
                             false);
        cs_field_set_key_int(bf, k_vis, post_vis);
        switch (cs_glob_thermal_model->itherm) {
        case CS_THERMAL_MODEL_ENTHALPY:
          cs_field_set_key_str(bf, cs_field_key_id("label"), "Hplus");
          break;
        case CS_THERMAL_MODEL_TOTAL_ENERGY:
          cs_field_set_key_str(bf, cs_field_key_id("label"), "Eplus");
          break;
        default:
          cs_field_set_key_str(bf, cs_field_key_id("label"), "Tplus");
          break;
        }
      }
    }

    /* Boundary temperature */

    bool post_b_temp = _surfacic_variable_post("boundary_temperature", true);

    /* Activate by default when using GUI; ignore for non-temperature variable
       when properties are not present in the GUI, or the thermal model is not
       set in the GUI, as this implies the GUI was probably not used and we
       cannot determine easily whether enthalpy to temperature conversion is
       available */

    if (cs_glob_thermal_model->itherm != CS_THERMAL_MODEL_TEMPERATURE) {
      if (   cs_tree_find_node_simple(cs_glob_tree, "property") == NULL
          || cs_gui_thermal_model_code() <= 0)
        post_b_temp = false;
    }

    if (post_b_temp) {
      cs_field_t *bf = cs_parameters_add_boundary_temperature();
      if (bf != NULL)
        cs_field_set_key_int(bf, k_vis, 1);
    }

  }
}

/*----------------------------------------------------------------------------
 * Define postprocessing meshes using an XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_meshes(void)
{
  const int *v_i = NULL;
  const cs_real_t *v_r = NULL;

  const char path_o[] = "analysis_control/output";
  cs_tree_node_t *tn_o = cs_tree_get_node(cs_glob_tree, path_o);

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_o, "mesh");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    v_i = cs_tree_node_get_child_values_int(tn, "id");
    const char *label = cs_tree_node_get_tag(tn, "label");
    const char *type = cs_tree_node_get_tag(tn, "type");

    if (v_i == NULL || label == NULL || type == NULL) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("Incorrect setup tree definition for the following node:\n"));
      cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
      bft_error(__FILE__, __LINE__, 0,
                _("One of the following child (tag) nodes is missing: %s"),
                "id, label, type");
    }

    int id = v_i[0];

    const char *location = cs_tree_node_get_child_value_str(tn, "location");
    if (location == NULL)
      location = "all[]";

    bool add_groups = true;
    bool auto_vars = true;

    cs_gui_node_get_status_bool(cs_tree_node_get_child(tn, "all_variables"),
                                &auto_vars);

    int n_writers = cs_tree_get_node_count(tn, "writer");
    int *writer_ids = NULL;
    BFT_MALLOC(writer_ids, n_writers, int);
    n_writers = 0;
    for (cs_tree_node_t *tn_w = cs_tree_get_node(tn, "writer");
         tn_w != NULL;
         tn_w = cs_tree_node_get_next_of_name(tn_w)) {
      v_i = cs_tree_node_get_child_values_int(tn_w, "id");
      if (v_i != NULL) {
        writer_ids[n_writers] = v_i[0];
        n_writers++;
      }
    }

    if (cs_gui_strcmp(type, "cells")) {
      cs_post_define_volume_mesh(id, label, location, add_groups, auto_vars,
                                 n_writers, writer_ids);
    }
    else if(cs_gui_strcmp(type, "interior_faces")) {
      cs_post_define_surface_mesh(id, label, location, NULL,
                                  add_groups, auto_vars,
                                  n_writers, writer_ids);
    }
    else if(cs_gui_strcmp(type, "boundary_faces")) {
      cs_post_define_surface_mesh(id, label, NULL, location,
                                  add_groups, auto_vars,
                                  n_writers, writer_ids);
    }
    else if(cs_gui_strcmp(type, "boundary_cells")) {
      /* Ensure location string remains available, or pass NULL */
      const char *_location = cs_tree_node_get_child_value_str(tn, "location");
      cs_post_define_volume_mesh_by_func(id, label,
                                         _selection_func_boundary_cells,
                                         (void *)_location,   /* input */
                                         true,
                                         add_groups, auto_vars,
                                         n_writers, writer_ids);
    }
    else if (   cs_gui_strcmp(type, "VolumicZone")
             || cs_gui_strcmp(type, "volume_czone")) {
      const cs_zone_t *z = cs_volume_zone_by_name(location);
      cs_post_define_mesh_by_location(id, label, z->location_id,
                                      add_groups, auto_vars,
                                      n_writers, writer_ids);
    }
    else if (   cs_gui_strcmp(type, "BoundaryZone")
             || cs_gui_strcmp(type, "boundary_zone")) {
      const cs_zone_t *z = cs_boundary_zone_by_name(location);
      cs_post_define_mesh_by_location(id, label, z->location_id,
                                      add_groups, auto_vars,
                                      n_writers, writer_ids);
    }
    else if (   cs_gui_strcmp(type, "BoundaryZone_cells")
             || cs_gui_strcmp(type, "boundary_zone_cells")) {
      const cs_zone_t *z = cs_boundary_zone_by_name(location);
      const char *criteria = NULL;
      /* FIXME: filter will not apply if zone is defined by a function */
      if (z->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
        criteria = cs_mesh_location_get_selection_string(z->location_id);
      cs_post_define_volume_mesh_by_func(id, label,
                                         _selection_func_boundary_cells,
                                         (void *)criteria,  /* input */
                                         true,
                                         add_groups, auto_vars,
                                         n_writers, writer_ids);
    }
    else if(cs_gui_strcmp(type, "interior_face_centers")) {
      cs_post_define_surface_mesh(id, label, location, NULL,
                                  add_groups, auto_vars,
                                  n_writers, writer_ids);
      cs_post_mesh_set_element_centers_only(id, true);
    }
    else if(   cs_gui_strcmp(type, "particles")
            || cs_gui_strcmp(type, "trajectories")) {
      bool trajectory = cs_gui_strcmp(type, "trajectories") ? true : false;
      v_r = cs_tree_node_get_child_values_real(tn, "density");
      double density = (v_r != NULL) ? v_r[0] : 1;
      cs_post_define_particles_mesh(id, label, location,
                                    density, trajectory, auto_vars,
                                    n_writers, writer_ids);
    }

    BFT_FREE(writer_ids);
  }

  /* Probe definitions */

  int n_probes = cs_tree_get_node_count(tn_o, "probe");

  if (n_probes > 0) {

    char **probe_labels;
    BFT_MALLOC(probe_labels, n_probes, char *);
    for (int ii = 0; ii < n_probes; ii++)
      BFT_MALLOC(probe_labels[ii], 128, char);

    const char *coord_node_name[] = {"probe_x", "probe_y", "probe_z"};

    cs_real_3_t *p_coords;
    BFT_MALLOC(p_coords, n_probes, cs_real_3_t);

    int i = 0;
    for (cs_tree_node_t *tn = cs_tree_get_node(tn_o, "probe");
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn), i++) {

      /* Probe coordinates */
      for (int j = 0; j < 3; j++) {
        v_r = cs_tree_node_get_child_values_real(tn, coord_node_name[j]);
        p_coords[i][j] = (v_r != NULL) ? v_r[0] : 0;
      }

      /* Probe name */
      const char *pn = cs_tree_node_get_child_value_str(tn, "name");
      strcpy(probe_labels[i], pn);
    }

    cs_probe_set_t *pset =
      cs_probe_set_create_from_array("probes",
                                     n_probes,
                                     (const cs_real_3_t *)p_coords,
                                     (const char **)probe_labels);

    /* Set snap mode. Default is "SNAP_TO_CENTER" */
    const char *snap_mode
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn_o, "probes_snap"),
                             "choice");
    if (cs_gui_strcmp(snap_mode, "snap_to_vertex"))
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_VERTEX);
    else if (cs_gui_strcmp(snap_mode, "none"))
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_NONE);
    else
      cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_ELT_CENTER);

    /* Activate interpolation if needed. Default is no */
    const char *activate_interpolation
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn_o,
                                                    "probes_interpolation"),
                             "choice");
    if (cs_gui_strcmp(activate_interpolation, "yes"))
      cs_probe_set_option(pset, "interpolation", "1");


    BFT_FREE(p_coords);

    for (int ii = 0; ii < n_probes; ii++)
      BFT_FREE(probe_labels[ii]);
    BFT_FREE(probe_labels);

  }

  /* Profile definitions;
     note that this may lead to additional writer definitions, as
     the GUI does not yet present profiles in a consistent
     writer/mesh logic (TODO) */

  _define_profiles();
}

/*----------------------------------------------------------------------------
 * Define postprocessing writers using an XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_writers(void)
{
  const char path_o[] = "analysis_control/output";
  cs_tree_node_t *tn_o = cs_tree_get_node(cs_glob_tree, path_o);

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_o, "writer");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const int *v_i = cs_tree_node_get_child_values_int(tn, "id");
    const char *label = cs_tree_node_get_tag(tn, "label");

    if (v_i == NULL || label == NULL) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("Incorrect setup tree definition for the following node:\n"));
      cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
      bft_error(__FILE__, __LINE__, 0,
                _("One of the following child (tag) nodes is missing: %s"),
                "id, label");
    }

    int id = v_i[0];

    fvm_writer_time_dep_t  time_dep = FVM_WRITER_FIXED_MESH;
    bool output_at_start = false;
    bool output_at_end = true;

    const char *directory
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "directory"),
                             "name");

    const char *frequency_choice
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "frequency"),
                             "period");
    int time_step = -1;
    cs_real_t time_value = -1.;

    if (cs_gui_strcmp(frequency_choice, "none"))
      time_step = -1;
    else if (cs_gui_strcmp(frequency_choice, "time_step")) {
      v_i = cs_tree_node_get_child_values_int(tn, "frequency");
      if (v_i != NULL) time_step = v_i[0];
    }
    else if (cs_gui_strcmp(frequency_choice, "time_value")) {
      const cs_real_t *v_r = cs_tree_node_get_child_values_real(tn, "frequency");
      if (v_r != NULL)
        time_value = v_r[0];
      else {
        v_r = cs_tree_node_get_child_values_real(tn, "frequency_time");
        if (v_r != NULL)
          time_value = v_r[0];
      }
    }
    else if (cs_gui_strcmp(frequency_choice, "formula"))
      time_step = -1;

    cs_gui_node_get_status_bool(cs_tree_node_get_child(tn, "output_at_start"),
                                &output_at_start);
    cs_gui_node_get_status_bool(cs_tree_node_get_child(tn, "output_at_end"),
                                &output_at_end);

    const char *format_name
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "format"),
                             "name");
    const char *format_options
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "format"),
                             "options");
    const char *time_dependency
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "time_dependency"),
                             "choice");

    if (cs_gui_strcmp(time_dependency, "fixed_mesh"))
      time_dep = FVM_WRITER_FIXED_MESH;
    else if(cs_gui_strcmp(time_dependency, "transient_coordinates"))
      time_dep = FVM_WRITER_TRANSIENT_COORDS;
    else if(cs_gui_strcmp(time_dependency, "transient_connectivity"))
      time_dep = FVM_WRITER_TRANSIENT_CONNECT;

    cs_post_define_writer(id,
                          label,
                          directory,
                          format_name,
                          format_options,
                          time_dep,
                          output_at_start,
                          output_at_end,
                          time_step,
                          time_value);
  }

  /* Probes default writer */

  const int *v_i
    = cs_tree_node_get_child_values_int(tn_o, "probe_recording_frequency");
  int frequency_n = (v_i != NULL) ? v_i[0] : 1;

  const cs_real_t *v_r
    = cs_tree_node_get_child_values_real(tn_o, "probe_recording_frequency_time");
  cs_real_t frequency_t = (v_r != NULL) ? v_r[0] : -1.;

  /* Time plot (probe) format string */
  const char *fmt_opts
    = cs_tree_node_get_tag(cs_tree_node_get_child(tn_o, "probe_format"),
                           "choice");

  cs_post_define_writer(CS_POST_WRITER_PROBES,   /* writer_id */
                        "",                      /* case_name */
                        "monitoring",            /* dir_name */
                        "time_plot",
                        fmt_opts,
                        FVM_WRITER_FIXED_MESH,
                        false,                   /* output_at_start */
                        false,                   /* output_at_end */
                        frequency_n,
                        frequency_t);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
