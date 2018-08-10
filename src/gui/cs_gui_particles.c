/*============================================================================
 * Management of the GUI parameters file: particles tracking
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
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui.h"
#include "cs_gui_util.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_prototypes.h"

#include "cs_lagr.h"
#include "cs_lagr_post.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_particles.h"

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

  if (attr != NULL) {
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
 *   name    <-- name of attribute in XML, or NULL for automatic name
 *----------------------------------------------------------------------------*/

static void
_attr_post_status(cs_tree_node_t       *tn_o,
                  cs_lagr_attribute_t   attr_id,
                  const char           *name)
{
  bool status = false;

  if (name != NULL)
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
 * Activate statistics and postprocessing for volume or boundary statistics
 *
 * parameters:
 *   tn_s     <--  parent tree node ("volume" or "boundary")
 *   name     <--  name of statistics in tree
 *   status   <->  associated postprocessing status
 *----------------------------------------------------------------------------*/

static void
_get_stats_post(cs_tree_node_t  *tn_s,
                const char      *name,
                int             *status)
{
  cs_tree_node_t *tn = NULL;

  for (tn = cs_tree_node_get_child(tn_s, "property");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    const char *_name = cs_tree_node_get_tag(tn, "name");
    if (_name != NULL) {
      if (! strcmp(_name, name)) {
        *status = 1; /* default if active */
        cs_gui_node_get_status_int
          (cs_tree_node_get_child(tn, "postprocessing_recording"),
           status);
        break;
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

  cs_glob_lagr_time_scheme->iilagr = 0;
  if (model != NULL) {
    if (! strcmp(model, "one_way"))
      cs_glob_lagr_time_scheme->iilagr = 1;
    else if (! strcmp(model, "two_way"))
      cs_glob_lagr_time_scheme->iilagr = 2;
    else if (! strcmp(model, "frozen"))
      cs_glob_lagr_time_scheme->iilagr = 3;
  }

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
  bft_printf("--iilagr = %i\n", cs_glob_lagr_time_scheme->iilagr);
#endif

  if (cs_glob_lagr_time_scheme->iilagr == 0)
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
  case 1:
    {
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "break_up"),
                                 &(cs_glob_lagr_specific_physics->idpvar));
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "evaporation"),
                                 &(cs_glob_lagr_specific_physics->impvar));
      cs_gui_node_get_status_int(cs_tree_node_get_child(tn_pm, "thermal"),
                                 &(cs_glob_lagr_specific_physics->itpvar));
    }
    break;
  case 2:
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
             tn != NULL;
             tn = cs_tree_node_get_next_of_name(tn)) {

          const int *v_i = cs_tree_node_get_child_values_int(tn, "coal");
          if (v_i == NULL) continue;
          int icoal = v_i[0] - 1;

          const cs_real_t *v_r = cs_tree_node_get_values_real(tn);
          if (v_r != NULL) attr_val[attr_id][icoal] = v_r[0];

        }
      }

    }
    break;
  }

  /* Two-way coupling */

  if (cs_glob_lagr_time_scheme->iilagr == 2) {
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
  if (choice != NULL)
    cs_glob_lagr_time_scheme->t_order = atoi(choice);

  cs_gui_node_get_status_int
    (cs_tree_node_get_child(tn_lagr, "turbulent_dispersion"),
     &(cs_glob_lagr_time_scheme->idistu));

  cs_gui_node_get_status_int
    (cs_tree_node_get_child(tn_lagr, "fluid_particles_turbulent_diffusion"),
     &(cs_glob_lagr_time_scheme->idiffl));

  cs_gui_node_get_status_int(cs_tree_node_get_child(tn_lagr,
                                                    "deposition_submodel"),
                             &(cs_glob_lagr_model->deposition));

  cs_gui_node_get_child_int(tn_lagr, "complete_model",
                            (&cs_glob_lagr_time_scheme->modcpl));

  choice = cs_tree_node_get_tag(cs_tree_node_get_child
                                  (tn_lagr, "complete_model_direction"),
                                "choice");
  if (choice != NULL)
    cs_glob_lagr_time_scheme->idirla = atoi(choice);

  /* Output */

  cs_tree_node_t *tn_o = cs_tree_node_get_child(tn_lagr, "output");
  if (tn_o != NULL) {
    _attr_post_status(tn_o, CS_LAGR_VELOCITY, "velocity_particles");
    _attr_post_status(tn_o, CS_LAGR_VELOCITY_SEEN, "velocity_fluid_seen");
    _attr_post_status(tn_o, CS_LAGR_RESIDENCE_TIME, "resident_time");
    _attr_post_status(tn_o, CS_LAGR_DIAMETER, "diameter");
    _attr_post_status(tn_o, CS_LAGR_TEMPERATURE, "temperature");
    _attr_post_status(tn_o, CS_LAGR_MASS, "mass");

    if (cs_glob_lagr_model->physical_model == 2) {
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
  if (tn_s != NULL) {

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

      int flag = 0;
      _get_stats_post(tn_vs, "Part_vol_frac", &flag);
      if (flag)
        cs_lagr_stat_activate(CS_LAGR_STAT_VOLUME_FRACTION);

      flag = 0;
      _get_stats_post(tn_vs, "Part_velocity", &flag);
      if (flag)
        cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

      flag = 0;
      _get_stats_post(tn_vs, "Part_resid_time", &flag);
      if (flag)
        cs_lagr_stat_activate_attr(CS_LAGR_RESIDENCE_TIME);

      flag = 0;
      _get_stats_post(tn_vs, "Part_stat_weight", &flag);
      if (flag)
        cs_lagr_stat_activate(CS_LAGR_STAT_CUMULATIVE_WEIGHT);

    }

    cs_tree_node_t *tn_bs = cs_tree_node_get_child(tn_s, "boundary");

    cs_gui_node_get_status_bool(tn_bs, &boundary_stats);

    if (boundary_stats) {

      _get_stats_post(tn_bs, "Part_impact_number",
                      &cs_glob_lagr_boundary_interactions->inbrbd);
      _get_stats_post(tn_bs, "Part_bndy_mass_flux",
                      &cs_glob_lagr_boundary_interactions->iflmbd);
      _get_stats_post(tn_bs, "Part_impact_angle",
                      &cs_glob_lagr_boundary_interactions->iangbd);
      _get_stats_post(tn_bs, "Part_impact_velocity",
                      &cs_glob_lagr_boundary_interactions->ivitbd);

      /* Coal fouling statistics*/
      _get_stats_post(tn_bs, "Part_fouled_impact_number",
                      &cs_glob_lagr_boundary_interactions->iencnbbd);
      _get_stats_post(tn_bs, "Part_fouled_mass_flux",
                      &cs_glob_lagr_boundary_interactions->iencmabd);
      _get_stats_post(tn_bs, "Part_fouled_diam",
                      &cs_glob_lagr_boundary_interactions->iencdibd);
      _get_stats_post(tn_bs, "Part_fouled_Xck",
                      &cs_glob_lagr_boundary_interactions->iencckbd);

    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iilagr = %i\n", cs_glob_lagr_time_scheme->iilagr);
  bft_printf("--isuila = %i\n", cs_glob_lagr_time_scheme->isuila);
  bft_printf("--isttio = %i\n", cs_glob_lagr_time_scheme->isttio);
  bft_printf("--idepst = %i\n", cs_glob_lagr_model->deposition);
  bft_printf("--iphyla = %i\n", cs_glob_lagr_model->physical_model);
  switch(cs_glob_lagr_model->physical_model) {
  case 0:
    break;
  case 1:
    bft_printf("--idpvar = %i\n", cs_glob_lagr_specific_physics->idpvar);
    bft_printf("--impvar = %i\n", cs_glob_lagr_specific_physics->impvar);
    bft_printf("--itpvar = %i\n", cs_glob_lagr_specific_physics->itpvar);
    break;
  case 2:
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

  if (cs_glob_lagr_time_scheme->iilagr == 2) {
    bft_printf("--nstits = %i\n", cs_glob_lagr_source_terms->nstits);
    bft_printf("--ltsdyn = %i\n", cs_glob_lagr_source_terms->ltsdyn);
    bft_printf("--ltsmas = %i\n", cs_glob_lagr_source_terms->ltsmas);
    bft_printf("--ltsthe = %i\n", cs_glob_lagr_source_terms->ltsthe);
  }

  bft_printf("--nordre = %i\n", cs_glob_lagr_time_scheme->t_order);
  bft_printf("--idistu = %i\n", cs_glob_lagr_time_scheme->idistu);
  bft_printf("--idiffl = %i\n", cs_glob_lagr_time_scheme->idiffl);
  bft_printf("--modcpl = %i\n", cs_glob_lagr_time_scheme->modcpl);
  bft_printf("--idirla = %i\n", cs_glob_lagr_time_scheme->idirla);

  bft_printf("--isuist = %i\n", cs_glob_lagr_stat_options->isuist);
  bft_printf("--nbclst = %i\n", cs_glob_lagr_model->n_stat_classes);

  bft_printf("--idstnt = %i\n", cs_glob_lagr_stat_options->idstnt);
  bft_printf("--nstist = %i\n", cs_glob_lagr_stat_options->nstist);
  int v_stats = (volume_stats) ? 1 : 0;
  bft_printf("--vol_stats = %i\n", v_stats);

  int b_stats = (boundary_stats) ? 1 : 0;
  bft_printf("--boundary_output = %i\n", b_stats);
  if (b_stats) {
    bft_printf("--inbrbd   = %i\n", cs_glob_lagr_boundary_interactions->inbrbd);
    bft_printf("--iflmbd   = %i\n", cs_glob_lagr_boundary_interactions->iflmbd);
    bft_printf("--iangbd   = %i\n", cs_glob_lagr_boundary_interactions->iangbd);
    bft_printf("--ivitbd   = %i\n", cs_glob_lagr_boundary_interactions->ivitbd);
    bft_printf("--iencnbbd = %i\n", cs_glob_lagr_boundary_interactions->iencnbbd);
    bft_printf("--iencmabd = %i\n", cs_glob_lagr_boundary_interactions->iencmabd);
    bft_printf("--iencdibd = %i\n", cs_glob_lagr_boundary_interactions->iencdibd);
    bft_printf("--iencckbd = %i\n", cs_glob_lagr_boundary_interactions->iencckbd);
  }

#endif
}

/*----------------------------------------------------------------------------
 * Define lagrangian model boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_gui_particles_bcs(void)
{
  cs_lnum_t iphyla = cs_glob_lagr_model->physical_model;
  cs_lagr_zone_data_t *bdy_cond = cs_lagr_get_boundary_conditions();
  cs_lagr_get_internal_conditions();

  /* zone 0 for "all", following zones defined by GUI */

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
#endif

  cs_tree_node_t *tn0 = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  /* First iteration only: memory allocation */

  int zone_id = 1; /* 1 to n */

  for (cs_tree_node_t *tn_bc = cs_tree_node_get_child(tn0, "boundary");
       tn_bc != NULL;
       tn_bc = cs_tree_node_get_next_of_name(tn_bc), zone_id++) {

    const char *label = cs_tree_node_get_tag(tn_bc, "label");
    const char *nature = cs_tree_node_get_tag(tn_bc, "nature");

    /* Find associated BC description with particle sub-node if present */
    /* Loop on all sub-nodes, of which some describe matching one BC info */

    cs_tree_node_t *tn_p = NULL;

    for (cs_tree_node_t *tn1 = cs_tree_node_get_child(tn0, nature);
         tn1 != NULL && tn_p == NULL;
         tn1 = cs_tree_node_get_next_of_name(tn1)) {

      if (! cs_gui_strcmp(label, cs_tree_node_get_tag(tn1, "label")))
        continue;

      tn_p = cs_tree_node_get_child(tn1, "particles");
    }

    if (tn_p == NULL)
      continue;

    /* Now we have a particles node for the given BC */

    const char *interaction = cs_tree_node_get_tag(tn_p, "choice");
    if (interaction == NULL)
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

    else if (! strcmp(interaction, "fouling") && iphyla == 2)
      bdy_cond->zone_type[zone_id] = CS_LAGR_FOULING;

    else if (   ! strcmp(interaction, "fouling")
             && (iphyla == 0  || iphyla == 1))
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
           tn_i != NULL;
           tn_i = cs_tree_node_get_next_of_name(tn_i), set_id++) {

        cs_lagr_injection_set_t *zis
          = cs_lagr_get_injection_set(bdy_cond, zone_id, set_id);

        cs_lagr_injection_set_default(zis);

        const char *choice = NULL;
        const int *v_i;
        const cs_real_t *v_r;
        cs_tree_node_t *tn_c; /* child nodes */

        v_i = cs_tree_node_get_child_values_int(tn_i, "number");
        if (v_i != NULL) zis->n_inject = v_i[0];

        v_i = cs_tree_node_get_child_values_int(tn_i, "frequency");
        if (v_i != NULL) zis->injection_frequency = v_i[0];

        v_i = cs_tree_node_get_child_values_int(tn_i, "statistical_groups");
        if (v_i != NULL) zis->cluster = v_i[0];

#if _XML_DEBUG_
        bft_printf("---number = %llu \n", (unsigned long long)(zis->n_inject));
        bft_printf("---frequency = %i \n", zis->injection_frequency);
        bft_printf("---statistical_groups = %i \n", zis->cluster);
#endif

        /* velocity */

        tn_c = cs_tree_node_get_child(tn_i, "velocity");
        if (tn_c != NULL) {

          choice = cs_tree_node_get_tag(tn_c, "choice");

          if (cs_gui_strcmp(choice, "fluid"))
            zis->velocity_profile = -1;

          else if (cs_gui_strcmp(choice, "norm")) {
            zis->velocity_profile = 0;

            v_r = cs_tree_node_get_child_values_real(tn_c, "norm");
            if (v_r != NULL) zis->velocity_magnitude = v_r[0];
          }
          else if (cs_gui_strcmp(choice, "components")) {
            const char *cname[] = {"velocity_x", "velocity_y", "velocity_z"};
            zis->velocity_profile = 1;
            for (int i = 0; i < 3; i++) {
              v_r = cs_tree_node_get_child_values_real(tn_c, cname[i]);
              if (v_r != NULL) zis->velocity[i] = v_r[0];
            }
          }
          else if (cs_gui_strcmp(choice, "subroutine"))
            zis->velocity_profile = 2;

#if _XML_DEBUG_
          bft_printf("---velocity choice: %i "
                     " (-1: fluid, 0: norm, 1: components, 2: subroutine)\n",
                     zis->velocity_profile);

          if (zis->velocity_profile == 0)
            bft_printf("----norm = %f \n", zis->velocity_magnitude);

          else if (zis->velocity_profile == 1) {
            bft_printf("----u = %f \n", zis->velocity[0]);
            bft_printf("----v = %f \n", zis->velocity[1]);
            bft_printf("----w = %f \n", zis->velocity[2]);
          }
#endif
        }

        /* statistical_weight, mass_flow_rate*/

        tn_c = cs_tree_node_get_child(tn_i, "statistical_weight");
        if (tn_c != NULL) {

          choice = cs_tree_node_get_tag(tn_c, "choice");

          if (cs_gui_strcmp(choice, "rate")) {
            v_r = cs_tree_node_get_child_values_real(tn_i, "mass_flow_rate");
            zis->stat_weight = 0;
            if (v_r != NULL) zis->flow_rate = v_r[0];
          }
          else if (cs_gui_strcmp(choice, "prescribed")) {
            v_r = cs_tree_node_get_values_real(tn_c);
            if (v_r != NULL) zis->stat_weight = v_r[0];
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
        if (v_r != NULL) zis->diameter = v_r[0];

        v_r = cs_tree_node_get_child_values_real
                (tn_i, "diameter_standard_deviation");
        if (v_r != NULL) zis->diameter_variance = v_r[0];

#if _XML_DEBUG_
        bft_printf("----diameter = %f \n", zis->diameter);
        bft_printf("----standard deviation = %f \n", zis->diameter_variance);
#endif

        /* density */
        if (iphyla != 2) {
          v_r = cs_tree_node_get_child_values_real(tn_i, "density");
          if (v_r != NULL) zis->density = v_r[0];

#if _XML_DEBUG_
          bft_printf("---density = %f \n", zis->density);
#endif
        }

        /* Fouling index*/
        v_r = cs_tree_node_get_child_values_real(tn_i, "fouling_index");
        if (v_r != NULL) zis->fouling_index = v_r[0];

#if _XML_DEBUG_
          bft_printf("---fouling_index = %f \n", zis->fouling_index);
#endif

        if (iphyla == 1) {

          /* temperature, specific_heat, emissivity */

          tn_c = cs_tree_node_get_child(tn_i, "temperature");
          if (tn_c != NULL) {

            choice = cs_tree_node_get_tag(tn_c, "choice");

            if (cs_gui_strcmp(choice, "prescribed")) {
              zis->temperature_profile = 1;
              v_r = cs_tree_node_get_values_real(tn_c);
              if (v_r != NULL) zis->temperature = v_r[0];
            }
            else if (cs_gui_strcmp(choice, "fluid")) {
              zis->temperature_profile = 0;
              zis->temperature = 0;
            }

          }

          v_r = cs_tree_node_get_child_values_real(tn_i, "specific_heat");
          if (v_r != NULL) zis->cp = v_r[0];

          v_r = cs_tree_node_get_child_values_real(tn_i, "emissivity");
          if (v_r != NULL) zis->emissivity = v_r[0];

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
        else if (iphyla == 2) {

          /* Read the coal number */
          v_i = cs_tree_node_get_child_values_int(tn_i, "coal_number");
          if (v_i != NULL) zis->coal_number = v_i[0];

          /* Data are read in pulverized fuel combustion module profile */
          v_r = cs_tree_node_get_child_values_real(tn_i, "coal_temperature");
          if (v_r != NULL) zis->temperature = v_r[0];

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
