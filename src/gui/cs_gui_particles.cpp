/*============================================================================
 * Management of the GUI parameters file: particles tracking
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "gui/cs_gui.h"
#include "gui/cs_gui_util.h"
#include "gui/cs_gui_boundary_conditions.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_post.h"
#include "lagr/cs_lagr_stat.h"
#include "lagr/cs_lagr_tracking.h"

#include "base/cs_parameters.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "gui/cs_gui_particles.h"

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return value of the particles model
 *
 * parameters:
 *   tn <-- "particles_models" tree node
 *
 * return:
 *   particles model based on "model" tag
 *----------------------------------------------------------------------------*/

static int
_get_particles_model(cs_tree_node_t  *tn)
{
  int retval = 0;

  const char *attr = cs_tree_node_get_tag(tn, "model");

  if (attr != nullptr) {
    if (cs_gui_strcmp(attr, "off"))
      retval = 0;
    else if (cs_gui_strcmp(attr, "thermal"))
      retval = 1;
    else if (cs_gui_strcmp(attr, "coal"))
      retval = 2;
  }

  return retval;
}

/*-----------------------------------------------------------------------------
 * Define postprocessing output status for a given particle attribute.
 *
 * parameters:
 *   tn_o    <-- tree node associated to lagrangian output
 *   attr_id <-- associated attribute id
 *   name    <-- name of attribute in XML, or nullptr for automatic name
 *----------------------------------------------------------------------------*/

static void
_attr_post_status(cs_tree_node_t       *tn_o,
                  cs_lagr_attribute_t   attr_id,
                  const char           *name)
{
  bool status = false;

  if (name != nullptr)
    cs_gui_node_get_status_bool(cs_tree_node_get_child(tn_o, name),
                                &status);
  else {
    cs_lagr_particle_attr_in_range(attr_id);
    const char *_name = cs_lagr_attribute_name[attr_id];
    cs_gui_node_get_status_bool(cs_tree_node_get_child(tn_o, _name),
                                &status);
  }

  cs_lagr_post_set_attr(attr_id, status);
}

/*-----------------------------------------------------------------------------
 * Activate statistics
 *
 * parameters:
 *   tn_s     <--  parent tree node ("volume" or "boundary")
 *----------------------------------------------------------------------------*/

static void
_get_stats_post(cs_tree_node_t  *tn_s)
{
  cs_tree_node_t *tn = nullptr;

  for (tn = cs_tree_node_get_child(tn_s, "property");
       tn != nullptr;
       tn = cs_tree_node_get_next_of_name(tn)) {
    const char *_name = cs_tree_node_get_tag(tn, "name");
    if (_name != nullptr) {
      int stat_type = cs_lagr_stat_type_by_name(_name);
      if (stat_type > -1) {
        int status = 1; /* default if active */
        cs_gui_node_get_status_int(tn, &status);
        if (status < 1)
          cs_lagr_stat_deactivate(stat_type);
        else
          cs_lagr_stat_activate(stat_type);
      }
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define lagrangian model options
 *----------------------------------------------------------------------------*/

void
cs_gui_particles_model(void)
{
  cs_tree_node_t *tn_lagr = cs_tree_get_node(cs_glob_tree, "lagrangian");

  const char *model = cs_tree_node_get_tag(tn_lagr, "model");

  cs_glob_lagr_time_scheme->iilagr = CS_LAGR_OFF;
  if (model != nullptr) {
    if (! strcmp(model, "one_way"))
      cs_glob_lagr_time_scheme->iilagr = CS_LAGR_ONEWAY_COUPLING;
    else if (! strcmp(model, "two_way"))
      cs_glob_lagr_time_scheme->iilagr = CS_LAGR_TWOWAY_COUPLING;
    else if (! strcmp(model, "frozen"))
      cs_glob_lagr_time_scheme->iilagr = CS_LAGR_FROZEN_CONTINUOUS_PHASE;
  }

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
  bft_printf("--iilagr = %i\n", cs_glob_lagr_time_scheme->iilagr);
#endif

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_OFF)
    return;

  /* Global settings */

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr, "restart"),
                             &(cs_glob_lagr_time_scheme->isuila));

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr,
                                                    "carrier_field_stationary"),
                             &(cs_glob_lagr_time_scheme->isttio));

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr,
                                                    "deposition_submodel"),
                             &(cs_glob_lagr_model->deposition));

  /* Particles model */

  cs_tree_node_t *tn_pm = cs_tree_get_node(tn_lagr, "particles_models");

  cs_glob_lagr_model->physical_model = _get_particles_model(tn_pm);

  switch (cs_glob_lagr_model->physical_model) {
  case CS_LAGR_PHYS_HEAT:
    {
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "break_up"),
                          &(cs_glob_lagr_specific_physics->solve_diameter));
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "evaporation"),
                                 &(cs_glob_lagr_specific_physics->solve_mass));
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "thermal"),
                          &(cs_glob_lagr_specific_physics->solve_temperature));
    }
    break;
  case CS_LAGR_PHYS_COAL:
    {
      cs_tree_node_t *tn_cf = cs_tree_node_get_child(tn_pm, "coal_fouling");

      cs_gui_node_get_status_int(tn_cf, &(cs_glob_lagr_model->fouling));

      /* Query various real values */

      const char *attr_name[] = {"threshold_temperature",
                                 "critical_viscosity",
                                 "fouling_coefficient_1",
                                 "fouling_coefficient_2"};
      cs_real_t *attr_val[] = {cs_glob_lagr_encrustation->tprenc,
                               cs_glob_lagr_encrustation->visref,
                               cs_glob_lagr_encrustation->enc1,
                               cs_glob_lagr_encrustation->enc2};

      for (int attr_id = 0; attr_id < 4; attr_id++) {
        for (cs_tree_node_t *tn = cs_tree_node_get_child
                                    (tn_cf, attr_name[attr_id]);
             tn != nullptr;
             tn = cs_tree_node_get_next_of_name(tn)) {

          const int *v_i = cs_tree_node_get_child_values_int(tn, "coal");
          if (v_i == nullptr) continue;
          int icoal = v_i[0] - 1;

          const cs_real_t *v_r = cs_tree_node_get_values_real(tn);
          if (v_r != nullptr) attr_val[attr_id][icoal] = v_r[0];

        }
      }

    }
    break;
  }

  /* Two-way coupling */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
    cs_tree_node_t *tn = cs_tree_node_get_child(tn_lagr, "two_way_coupling");

    cs_gui_node_get_child_int(tn, "iteration_start",
                              (&cs_glob_lagr_source_terms->nstits));

    cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "dynamic"),
                               (&cs_glob_lagr_source_terms->ltsdyn));

    cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "mass"),
                               (&cs_glob_lagr_source_terms->ltsmas));

    cs_gui_node_get_status_int(cs_tree_node_get_child(tn, "thermal"),
                               (&cs_glob_lagr_source_terms->ltsthe));
  }

  /* Numerical modeling */

  const char *choice;

  choice = cs_tree_node_get_tag(cs_tree_node_get_child(tn_lagr, "scheme_order"),
                                "choice");
  if (choice != nullptr)
    cs_glob_lagr_time_scheme->t_order = atoi(choice);

  cs_gui_node_get_status_int
    (cs_tree_node_get_child(tn_lagr, "fluid_particles_turbulent_diffusion"),
     &(cs_glob_lagr_model->idiffl));

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr,
                                                    "deposition_submodel"),
                             &(cs_glob_lagr_model->deposition));

  /* Regular or fluid particles */
  cs_gui_node_get_child_int(tn_lagr, "regular_particles",
                            (&cs_glob_lagr_model->modcpl));

  /* Output */

  cs_tree_node_t *tn_o = cs_tree_node_get_child(tn_lagr, "output");
  if (tn_o != nullptr) {
    _attr_post_status(tn_o, CS_LAGR_VELOCITY, "velocity_particles");
    _attr_post_status(tn_o, CS_LAGR_VELOCITY_SEEN, "velocity_fluid_seen");
    _attr_post_status(tn_o, CS_LAGR_RESIDENCE_TIME, "resident_time");
    _attr_post_status(tn_o, CS_LAGR_DIAMETER, "diameter");
    _attr_post_status(tn_o, CS_LAGR_TEMPERATURE, "temperature");
    _attr_post_status(tn_o, CS_LAGR_MASS, "mass");
    _attr_post_status(tn_o, CS_LAGR_AGGLO_CLASS_ID, "parcel_class");
    _attr_post_status(tn_o, CS_LAGR_STAT_WEIGHT, "stat_weight");

    if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL) {
      _attr_post_status(tn_o, CS_LAGR_SHRINKING_DIAMETER,
                        "shrinking_core_diameter");
      _attr_post_status(tn_o, CS_LAGR_WATER_MASS, "moisture_mass_fraction");
      _attr_post_status(tn_o, CS_LAGR_COAL_MASS, "raw_coal_mass_fraction");
      _attr_post_status(tn_o, CS_LAGR_COKE_MASS, "char_mass_fraction");
    }

    cs_gui_node_get_child_int(tn_o, "listing_printing_frequency",
                              &cs_glob_lagr_log_frequency_n);
  }

  /* Statistics */

  bool volume_stats = false;
  bool boundary_stats = false;

  cs_tree_node_t *tn_s = cs_tree_node_get_child(tn_lagr, "statistics");
  if (tn_s != nullptr) {

    cs_gui_node_get_child_int(tn_s, "statistics_groups_of_particles",
                              &cs_glob_lagr_model->n_stat_classes);

    cs_gui_node_get_child_int(tn_s, "iteration_start",
                              &cs_glob_lagr_stat_options->idstnt);

    cs_gui_node_get_child_int(tn_s, "iteration_steady_start",
                              &cs_glob_lagr_stat_options->nstist);

    cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr,
                                                      "restart"),
                               &cs_glob_lagr_stat_options->isuist);

    cs_gui_node_get_child_real(tn_s, "threshold",
                               &cs_glob_lagr_stat_options->threshold);

    cs_tree_node_t *tn_vs = cs_tree_node_get_child(tn_s, "volume");

    cs_gui_node_get_status_bool(tn_vs, &volume_stats);

    if (volume_stats) {

      /* Default statistics */

      cs_lagr_stat_activate(CS_LAGR_STAT_CUMULATIVE_WEIGHT);
      cs_lagr_stat_activate(CS_LAGR_STAT_VOLUME_FRACTION);
      cs_lagr_stat_activate_attr(CS_LAGR_RESIDENCE_TIME);
      cs_lagr_stat_activate_attr(CS_LAGR_DIAMETER);
      cs_lagr_stat_activate_attr(CS_LAGR_MASS);
      cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

      switch(cs_glob_lagr_model->physical_model) {
      case CS_LAGR_PHYS_HEAT:
        cs_lagr_stat_activate_attr(CS_LAGR_TEMPERATURE);
        break;
      case CS_LAGR_PHYS_COAL:
        cs_lagr_stat_activate_attr(CS_LAGR_WATER_MASS);
        cs_lagr_stat_activate_attr(CS_LAGR_TEMPERATURE);
        cs_lagr_stat_activate_attr(CS_LAGR_COAL_MASS);
        cs_lagr_stat_activate_attr(CS_LAGR_COKE_MASS);
        cs_lagr_stat_activate_attr(CS_LAGR_COAL_DENSITY);
      default:
        break;
      }

      /* XML-defined */
      _get_stats_post(tn_vs);
    }

    cs_tree_node_t *tn_bs = cs_tree_node_get_child(tn_s, "boundary");

    cs_gui_node_get_status_bool(tn_bs, &boundary_stats);

    if (boundary_stats) {

      /* Default statistics */

      cs_lagr_stat_activate(CS_LAGR_STAT_E_CUMULATIVE_WEIGHT);
      cs_lagr_stat_activate(CS_LAGR_STAT_MASS_FLUX);
      cs_lagr_stat_activate(CS_LAGR_STAT_IMPACT_ANGLE);
      cs_lagr_stat_activate(CS_LAGR_STAT_IMPACT_VELOCITY);

      if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL) {
        cs_lagr_stat_activate(CS_LAGR_STAT_FOULING_CUMULATIVE_WEIGHT);
        cs_lagr_stat_activate(CS_LAGR_STAT_FOULING_MASS_FLUX);
        cs_lagr_stat_activate(CS_LAGR_STAT_FOULING_DIAMETER);
        cs_lagr_stat_activate(CS_LAGR_STAT_FOULING_COKE_FRACTION);
      }

      _get_stats_post(tn_bs);
    }

  }

  /* When activating the complete turbulent dispersion model,
   * statistics are required, so activate it after the start time of
   * statistics.
   */
  if (cs_glob_lagr_model->modcpl > 0)
    cs_glob_lagr_model->modcpl
      = cs::max(cs_glob_lagr_model->modcpl, cs_glob_lagr_stat_options->idstnt);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iilagr = %i\n", cs_glob_lagr_time_scheme->iilagr);
  bft_printf("--isuila = %i\n", cs_glob_lagr_time_scheme->isuila);
  bft_printf("--isttio = %i\n", cs_glob_lagr_time_scheme->isttio);
  bft_printf("--idepst = %i\n", cs_glob_lagr_model->deposition);
  bft_printf("--iphyla = %i\n", cs_glob_lagr_model->physical_model);
  switch(cs_glob_lagr_model->physical_model) {
  case CS_LAGR_PHYS_OFF:
    break;
  case CS_LAGR_PHYS_HEAT:
    bft_printf("--solve_diameter = %i\n",
        cs_glob_lagr_specific_physics->solve_diameter);
    bft_printf("--solve_mass = %i\n",
        cs_glob_lagr_specific_physics->solve_mass);
    bft_printf("--solve_temperature = %i\n",
        cs_glob_lagr_specific_physics->solve_temperature);
    break;
  case CS_LAGR_PHYS_COAL:
    bft_printf("--iencra = %i\n", cs_glob_lagr_model->fouling);
    const cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
    for (int icoal = 0; icoal < extra->ncharb; icoal++) {
      bft_printf("--tprenc[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->tprenc[icoal]);
      bft_printf("--visref[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->visref[icoal]);
      bft_printf("--enc1[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->enc1[icoal]);
      bft_printf("--enc2[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->enc2[icoal]);
    }
    break;
  }

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
    bft_printf("--nstits = %i\n", cs_glob_lagr_source_terms->nstits);
    bft_printf("--ltsdyn = %i\n", cs_glob_lagr_source_terms->ltsdyn);
    bft_printf("--ltsmas = %i\n", cs_glob_lagr_source_terms->ltsmas);
    bft_printf("--ltsthe = %i\n", cs_glob_lagr_source_terms->ltsthe);
  }

  bft_printf("--nordre = %i\n", cs_glob_lagr_time_scheme->t_order);
  bft_printf("--idistu = %i\n", cs_glob_lagr_model->idistu);
  bft_printf("--idiffl = %i\n", cs_glob_lagr_model->idiffl);
  bft_printf("--modcpl = %i\n", cs_glob_lagr_model->modcpl);

  bft_printf("--isuist = %i\n", cs_glob_lagr_stat_options->isuist);
  bft_printf("--nbclst = %i\n", cs_glob_lagr_model->n_stat_classes);

  bft_printf("--idstnt = %i\n", cs_glob_lagr_stat_options->idstnt);
  bft_printf("--nstist = %i\n", cs_glob_lagr_stat_options->nstist);
  int v_stats = (volume_stats) ? 1 : 0;
  bft_printf("--vol_stats = %i\n", v_stats);

  int b_stats = (boundary_stats) ? 1 : 0;
  bft_printf("--boundary_output = %i\n", b_stats);
  if (b_stats) {
    bft_printf("--has_part_impact_nbr   = %i\n",
               cs_glob_lagr_boundary_interactions->has_part_impact_nbr);
    bft_printf("--iencnbbd = %i\n", cs_glob_lagr_boundary_interactions->iencnbbd);
    bft_printf("--iencckbd = %i\n", cs_glob_lagr_boundary_interactions->iencckbd);
  }

#endif
}

/*----------------------------------------------------------------------------
 * Define Lagrangian model boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_gui_particles_bcs(void)
{
  cs_lnum_t iphyla = cs_glob_lagr_model->physical_model;
  cs_lagr_zone_data_t *bdy_cond = cs_lagr_get_boundary_conditions();

  /* zone 0 for "all", following zones defined by GUI */

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
#endif

  cs_tree_node_t *tn0 = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  /* First iteration only: memory allocation */

  int zone_id = 1; /* 1 to n */

  for (cs_tree_node_t *tn_bc = cs_tree_node_get_child(tn0, "boundary");
       tn_bc != nullptr;
       tn_bc = cs_tree_node_get_next_of_name(tn_bc), zone_id++) {

    const char *label = cs_tree_node_get_tag(tn_bc, "label");
    const char *nature = cs_tree_node_get_tag(tn_bc, "nature");

    /* Find associated BC description with particle sub-node if present */
    /* Loop on all sub-nodes, of which some describe matching one BC info */

    cs_tree_node_t *tn_p = nullptr;

    for (cs_tree_node_t *tn1 = cs_tree_node_get_child(tn0, nature);
         tn1 != nullptr && tn_p == nullptr;
         tn1 = cs_tree_node_get_next_of_name(tn1)) {

      if (! cs_gui_strcmp(label, cs_tree_node_get_tag(tn1, "label")))
        continue;

      tn_p = cs_tree_node_get_child(tn1, "particles");
    }

    if (tn_p == nullptr)
      continue;

    /* Now we have a particles node for the given BC */

    const char *interaction = cs_tree_node_get_tag(tn_p, "choice");
    if (interaction == nullptr)
      continue;

    if (! strcmp(interaction, "inlet"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_INLET;

    else if (! strcmp(interaction, "outlet"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_OUTLET;

    else if (! strcmp(interaction, "bounce"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_REBOUND;

    else if (! strcmp(interaction, "part_symmetry"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_SYM;

    else if (! strcmp(interaction, "deposit1"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO1;

    else if (! strcmp(interaction, "deposit2"))
      bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO2;

    else if (! strcmp(interaction, "fouling") && iphyla == CS_LAGR_PHYS_COAL)
      bdy_cond->zone_type[zone_id] = CS_LAGR_FOULING;

    else if (   ! strcmp(interaction, "fouling")
             && (iphyla == CS_LAGR_PHYS_OFF  || iphyla == CS_LAGR_PHYS_HEAT))
      bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO_DLVO;

#if _XML_DEBUG_
    bft_printf("--zone_type[%i] = %i \n", zone_id,
               bdy_cond->zone_type[zone_id]);

    bft_printf("--        : label    %s \n", label);
    bft_printf("--        : nature   %s \n", nature);
    bft_printf("--        : p_nature %i \n", bdy_cond->zone_type[zone_id]);
#endif

    /* Additional info for inlet */

    if (bdy_cond->zone_type[zone_id] == CS_LAGR_INLET) {

      /* Loop on injection sets (mislabelled as "class") */

      int set_id = 0;
      for (cs_tree_node_t *tn_i = cs_tree_node_get_child(tn_p, "class");
           tn_i != nullptr;
           tn_i = cs_tree_node_get_next_of_name(tn_i), set_id++) {

        cs_lagr_injection_set_t *zis
          = cs_lagr_get_injection_set(bdy_cond, zone_id, set_id);

        cs_lagr_injection_set_default(zis);

        const char *choice = nullptr;
        const int *v_i;
        const cs_real_t *v_r;
        cs_tree_node_t *tn_c; /* child nodes */

        v_i = cs_tree_node_get_child_values_int(tn_i, "number");
        if (v_i != nullptr) zis->n_inject = v_i[0];

        v_i = cs_tree_node_get_child_values_int(tn_i, "frequency");
        if (v_i != nullptr) zis->injection_frequency = v_i[0];

        v_i = cs_tree_node_get_child_values_int(tn_i, "statistical_groups");
        if (v_i != nullptr) zis->cluster = v_i[0];

#if _XML_DEBUG_
        bft_printf("---number = %llu \n", (unsigned long long)(zis->n_inject));
        bft_printf("---frequency = %i \n", zis->injection_frequency);
        bft_printf("---statistical_groups = %i \n", zis->cluster);
#endif

        /* velocity */

        tn_c = cs_tree_node_get_child(tn_i, "velocity");
        if (tn_c != nullptr) {

          choice = cs_tree_node_get_tag(tn_c, "choice");

          if (choice == nullptr || cs_gui_strcmp(choice, "fluid"))
            zis->velocity_profile = CS_LAGR_IN_IMPOSED_FLUID_VALUE;

          else if (cs_gui_strcmp(choice, "norm")) {
            zis->velocity_profile = CS_LAGR_IN_IMPOSED_NORM;

            v_r = cs_tree_node_get_child_values_real(tn_c, "norm");
            if (v_r != nullptr) zis->velocity_magnitude = v_r[0];
          }
          else if (cs_gui_strcmp(choice, "components")) {
            const char *cname[] = {"velocity_x", "velocity_y", "velocity_z"};
            zis->velocity_profile = CS_LAGR_IN_IMPOSED_COMPONENTS;
            for (int i = 0; i < 3; i++) {
              v_r = cs_tree_node_get_child_values_real(tn_c, cname[i]);
              if (v_r != nullptr) zis->velocity[i] = v_r[0];
            }
          }

#if _XML_DEBUG_
          bft_printf("---velocity choice: %i "
                     " (-1: fluid, 0: norm, 1: components, 2: subroutine)\n",
                     zis->velocity_profile);

          if (zis->velocity_profile == CS_LAGR_IN_IMPOSED_NORM)
            bft_printf("----norm = %f \n", zis->velocity_magnitude);

          else if (zis->velocity_profile == CS_LAGR_IN_IMPOSED_COMPONENTS) {
            bft_printf("----u = %f \n", zis->velocity[0]);
            bft_printf("----v = %f \n", zis->velocity[1]);
            bft_printf("----w = %f \n", zis->velocity[2]);
          }
#endif
        }

        /* statistical_weight, mass_flow_rate*/

        tn_c = cs_tree_node_get_child(tn_i, "statistical_weight");
        if (tn_c != nullptr) {

          choice = cs_tree_node_get_tag(tn_c, "choice");

          if (cs_gui_strcmp(choice, "rate")) {
            v_r = cs_tree_node_get_child_values_real(tn_i, "mass_flow_rate");
            zis->stat_weight = 0;
            if (v_r != nullptr) zis->flow_rate = v_r[0];
          }
          else if (cs_gui_strcmp(choice, "prescribed") || choice == nullptr) {
            v_r = cs_tree_node_get_values_real(tn_c);
            if (v_r != nullptr) zis->stat_weight = v_r[0];
            zis->flow_rate = 0;
          }

#if _XML_DEBUG_
          bft_printf("---statistical weight choice: %s "
                     " (1: prescribed, 2: rate)\n", choice);

          bft_printf("----statistical weight = %f \n", zis->stat_weight);
          bft_printf("----mass flow rate = %f \n", zis->flow_rate);
#endif
        }

        /* diameter */

        v_r = cs_tree_node_get_child_values_real(tn_i, "diameter");
        if (v_r != nullptr) zis->diameter = v_r[0];

        v_r = cs_tree_node_get_child_values_real
                (tn_i, "diameter_standard_deviation");
        if (v_r != nullptr) zis->diameter_variance = v_r[0];

#if _XML_DEBUG_
        bft_printf("----diameter = %f \n", zis->diameter);
        bft_printf("----standard deviation = %f \n", zis->diameter_variance);
#endif

        /* density */
        if (iphyla != CS_LAGR_PHYS_COAL) {
          v_r = cs_tree_node_get_child_values_real(tn_i, "density");
          if (v_r != nullptr) zis->density = v_r[0];

#if _XML_DEBUG_
          bft_printf("---density = %f \n", zis->density);
#endif
        }

        /* Fouling index*/
        v_r = cs_tree_node_get_child_values_real(tn_i, "fouling_index");
        if (v_r != nullptr) zis->fouling_index = v_r[0];

#if _XML_DEBUG_
          bft_printf("---fouling_index = %f \n", zis->fouling_index);
#endif

        if (iphyla == CS_LAGR_PHYS_HEAT) {

          /* temperature, specific_heat, emissivity */

          tn_c = cs_tree_node_get_child(tn_i, "temperature");
          if (tn_c != nullptr) {

            choice = cs_tree_node_get_tag(tn_c, "choice");

            if (cs_gui_strcmp(choice, "prescribed")) {
              zis->temperature_profile = 1;
              v_r = cs_tree_node_get_values_real(tn_c);
              if (v_r != nullptr) zis->temperature = v_r[0];
            }
            else if (cs_gui_strcmp(choice, "fluid") || choice == nullptr) {
              zis->temperature_profile = 0;
              zis->temperature = 0;
            }

          }

          v_r = cs_tree_node_get_child_values_real(tn_i, "specific_heat");
          if (v_r != nullptr) zis->cp = v_r[0];

          v_r = cs_tree_node_get_child_values_real(tn_i, "emissivity");
          if (v_r != nullptr) zis->emissivity = v_r[0];

#if _XML_DEBUG_
          bft_printf("---temperature choice = %s "
                     "(0: fluid, 1: prescribed)\n",
                     choice);

          bft_printf("----temperature = %f \n", zis->temperature);

          bft_printf("---specific heat = %f \n", zis->cp);
          bft_printf("---emissivity = %f \n", zis->emissivity);
#endif
        }

        /* coal */
        else if (iphyla == CS_LAGR_PHYS_COAL) {

          /* Read the coal number */
          v_i = cs_tree_node_get_child_values_int(tn_i, "coal_number");
          if (v_i != nullptr) zis->coal_number = v_i[0];

          /* Data are read in pulverized fuel combustion module profile */
          v_r = cs_tree_node_get_child_values_real(tn_i, "coal_temperature");
          if (v_r != nullptr) zis->temperature = v_r[0];

#if _XML_DEBUG_
          bft_printf("---coal number = %i \n", zis->coal_number);
          bft_printf("---coal temperature = %f \n", zis->temperature);
#endif /* _XML_DEBUG_ */

        }

      } /* End of loop on injection sets */

    } /* End of test on inlet */

  } /* End of loop on zones and interaction */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
