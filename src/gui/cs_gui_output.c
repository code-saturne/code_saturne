/*============================================================================
 * Management of the GUI parameters file: output parameters
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "mei_evaluate.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_gui.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_selector.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_parameters.h"
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
 * Static global variables
 *============================================================================*/

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

  cs_tree_node_t *tn =_get_node(field_type, f->name);

  if (tn == NULL)
    return;

  /* Listing output */

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "listing_printing"),
                             &f_log);
  if (f_log != -999)
    cs_field_set_key_int(f, k_log, f_log);

  /* Postprocessing outputs */

  cs_gui_node_get_status_int(cs_tree_node_get_child
                               (tn, "postprocessing_recording"),
                             &f_post);
  if (f_post == 1)
    cs_field_set_key_int_bits(f, k_post, CS_POST_ON_LOCATION);
  else if (f_post == 0)
    cs_field_clear_key_int_bits(f, k_post, CS_POST_ON_LOCATION);
  else /* status unspecified here but property referenced in tree,
          could be improved by depending on field type or flags */
    cs_field_set_key_int_bits(f, k_post, CS_POST_ON_LOCATION);

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "probes_recording"),
                             &f_monitor);
  if (f_monitor == 1)
    cs_field_set_key_int_bits(f, k_post, CS_POST_MONITOR);
  else if (f_monitor == 0)
    cs_field_clear_key_int_bits(f, k_post, CS_POST_MONITOR);
  else /* status unspecified here but property referenced in tree,
          could be improved by depending on field type or flags */
    if (f->location_id == CS_MESH_LOCATION_CELLS)
      cs_field_set_key_int_bits(f, k_post, CS_POST_MONITOR);
    else
      cs_field_clear_key_int_bits(f, k_post, CS_POST_MONITOR);

  /* Take into account labels */

  const char *label = cs_tree_node_get_tag(tn, "label");
  if (label != NULL)
    cs_field_set_key_str(f, k_lbl, label);
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

/*-----------------------------------------------------------------------------
 * Initialize mei tree and check for symbols existence.
 *
 * parameters:
 *   formula        <--  mei formula
 *   nt_cur         <--  current time step
 *   t_cur          <--  current time value
 *
 * return:
 *   MEI tree object
 *----------------------------------------------------------------------------*/

static mei_tree_t *
_init_mei_tree(const char      *formula,
               const int        nt_cur,
               const cs_real_t  t_cur)
{
  /* return an empty interpreter */

  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "niter", nt_cur );
  mei_tree_insert(tree, "t", t_cur);

  /* add variable from notebook */
  cs_gui_add_notebook_variables(tree);

  /* try to build the interpreter */
  if (mei_tree_builder(tree))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not interpret expression: %s\n"), tree->string);
  /* check for symbols */
  if (mei_tree_find_symbol(tree, "iactive"))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not find the required symbol: %s\n"), "iactive");

  return tree;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine output boundary fields
 *----------------------------------------------------------------------------*/

void CS_PROCF (cspstb, CSPSTB) (cs_int_t        *ipstdv)
{
  if (!cs_gui_file_is_loaded())
    return;

  /* Surfacic variables output */

  for (int i = 0; i < 5; i++)
    ipstdv[i] = 0;

  if (cs_glob_physical_model_flag[CS_GROUNDWATER] == -1) {
    if (_surfacic_variable_post("stress", true))
      ipstdv[0] += 1;
    if (_surfacic_variable_post("stress_tangential", false))
      ipstdv[0] += 2;
    if (_surfacic_variable_post("stress_normal", false))
      ipstdv[0] += 4;

    if (_surfacic_variable_post("yplus", true))
      ipstdv[1] = 1;
    if (_surfacic_variable_post("tplus", false))
      ipstdv[2] = 1;
    if (_surfacic_variable_post("thermal_flux", true))
      ipstdv[3] = 1;
    if (_surfacic_variable_post("boundary_temperature", true)) {
      cs_field_t *bf = cs_parameters_add_boundary_temperature();
      if (bf != NULL) {
        int k_vis = cs_field_key_id("post_vis");
        cs_field_set_key_int(bf, k_vis, 1);
      }
    }

    if (_surfacic_variable_post("boundary_layer_nusselt", false))
      ipstdv[4] = 1;
  }
}

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
    else if (f->type & CS_FIELD_PROPERTY)
      _field_post("property", f->id);
    else if (moment_id != NULL) {
      if (moment_id[f_id] > -1)
        _field_post("time_average", f->id);
    }
  }

  BFT_FREE(moment_id);

#if _XML_DEBUG_
  bft_printf("==>CSENSO\n");
  bft_printf("--iecaux = %i\n", *iecaux);
  bft_printf("--ntlist = %i\n", cs_glob_log_frequency);
#endif
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
    else if(cs_gui_strcmp(type, "VolumicZone")) {
      const cs_zone_t *z = cs_volume_zone_by_name(location);
      const char *criteria =
        cs_mesh_location_get_selection_string(z->location_id);
      cs_post_define_volume_mesh(id, label, criteria, add_groups, auto_vars,
                                 n_writers, writer_ids);
    }
    else if(cs_gui_strcmp(type, "BoundaryZone")) {
      const cs_zone_t *z = cs_boundary_zone_by_name(location);
      const char *criteria =
        cs_mesh_location_get_selection_string(z->location_id);
      cs_post_define_surface_mesh(id, label, NULL, criteria,
                                  add_groups, auto_vars,
                                  n_writers, writer_ids);
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

    const char *coord_node_name[] = {"probe_x", "probe_y", "probe_z"};

    cs_real_3_t *p_coords;
    BFT_MALLOC(p_coords, n_probes, cs_real_3_t);

    int i = 0;
    for (cs_tree_node_t *tn = cs_tree_get_node(tn_o, "probe");
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn), i++) {

      for (int j = 0; j < 3; j++) {
        v_r = cs_tree_node_get_child_values_real(tn, coord_node_name[j]);
        p_coords[i][j] = (v_r != NULL) ? v_r[0] : 0;
      }
    }

    cs_probe_set_create_from_array("probes",
                                   n_probes,
                                   (const cs_real_3_t *)p_coords,
                                   NULL);

    BFT_FREE(p_coords);

    v_i = cs_tree_node_get_child_values_int(tn_o, "probe_recording_frequency");
    int frequency_n = (v_i != NULL) ? v_i[0] : 1;

    v_r = cs_tree_node_get_child_values_real
           (tn_o, "probe_recording_frequency_time");
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
}

/*----------------------------------------------------------------------------
 * Activate writers depending on a formula
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_activate(void)
{
  const char path_o[] = "analysis_control/output";
  cs_tree_node_t *tn_o = cs_tree_get_node(cs_glob_tree, path_o);

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_o, "writer");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const int *v_i = cs_tree_node_get_child_values_int(tn, "id");
    if (v_i == NULL) continue;

    int id = v_i[0];

    const char *frequency_choice
      = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "frequency"),
                             "period");

    if (cs_gui_strcmp(frequency_choice, "formula")) {
      const char *formula = cs_tree_node_get_child_value_str(tn, "frequency");
      assert(formula != NULL);
      const cs_time_step_t *ts = cs_glob_time_step;
      mei_tree_t *ev_formula = _init_mei_tree(formula, ts->nt_cur, ts->t_cur);
      mei_evaluate(ev_formula);
      int iactive =  mei_tree_lookup(ev_formula, "iactive");
      mei_tree_destroy(ev_formula);
      if (iactive == 1)
        cs_post_activate_writer(id, true);
      else
        cs_post_activate_writer(id, false);
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
