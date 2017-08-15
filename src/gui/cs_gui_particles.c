/*============================================================================
 * Management of the GUI parameters file: particles tracking
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *----------------------------------------------------------------------------*/

static int
_get_particles_model(const char  *model)
{
  int retval = 0;

  char *path;
  char *attr;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "lagrangian", model);
  cs_xpath_add_attribute(&path, "model");
  attr = cs_gui_get_attribute_value(path);

  if (attr != NULL) {
    if (cs_gui_strcmp(attr, "off"))
      retval = 0;
    else if (cs_gui_strcmp(attr, "thermal"))
      retval = 1;
    else if (cs_gui_strcmp(attr, "coal"))
      retval = 2;
    BFT_FREE(attr);
  }
  BFT_FREE(path);

  return retval;
}

/*-----------------------------------------------------------------------------
 * Return value of the parameter of the character type for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_status(int  *keyword,
            int   nbr,
            ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  int i;
  int result;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_attribute(&path, "status");
  if(cs_gui_get_status(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer parameters for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_int(int *const keyword, const int nbr, ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  int value = 0;
  int i;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_int(path, &value))
    *keyword = value;

  BFT_FREE(path);

}


/*-----------------------------------------------------------------------------
 * Return float parameters for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_double(double *const  keyword,
            const int      nbr,
            ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  double value = 0;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for (int i=0; i < nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value))
    *keyword = value;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of the attribute of the character type for lagrangian
 *
 *   parameters:
 *   param     <--   name of the attribute
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static char*
_get_attr(const char *const param,
          const int nbr,
          ...)
{
  va_list list;

  int i;
  char *elt = NULL;
  char *path;
  char *name;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_attribute(&path, param);

  name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name;
}

/*-----------------------------------------------------------------------------
 * Define postprocessing output status for a given particle attribute.
 *
 * parameters:
 *   attr_id <-- associated attribute id
 *   name    <-- name of attribute in XML, or NULL for automatic name
 *----------------------------------------------------------------------------*/

static void
_attr_post_status(cs_lagr_attribute_t   attr_id,
                  const char           *name)
{
  int key = 0;
  if (name != NULL)
    _get_status(&key, 3, "lagrangian", "output", name);
  else {
    cs_lagr_particle_attr_in_range(attr_id);
    const char *_name = cs_lagr_attribute_name[attr_id];
    _get_status(&key, 3, "lagrangian", "output", _name);
  }

  bool status = (key = 0) ? false : true;

  cs_lagr_post_set_attr(attr_id, status);
}

/*-----------------------------------------------------------------------------
 * Return float parameters for coal parameters
 *
 *   parameters:
 *    param         -->   value to modify
 *    name          -->   name of property
 *    icoal         -->   number of coal
 *----------------------------------------------------------------------------*/

static void
_get_coal_double(double *const param, const char *const name, int icoal)
{
  double result = 0;
  char *path = NULL;
  char scoal[2];

  sprintf(scoal, "%i", icoal);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4,
                        "lagrangian", "particles_models", "coal_fouling", name);
  cs_xpath_add_test_attribute(&path, "coal", scoal);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *param = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return status and label of the property for post treatment
 *
 *   parameters:
 *    type          -->   type of property ('volume' or 'boundary')
 *    name          -->   name of property
 *    list_value    <--   status for listing
 *    record_value  <--   status for post processing
 *----------------------------------------------------------------------------*/

static char*
_get_char_post(const char *const type,
               const char *const name,
               int  *record_value)
{
  char *path, *path1, *path2 = NULL;
  char *label = NULL;
  int result;

  *record_value = 1;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "lagrangian", "statistics", type, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  BFT_MALLOC(path1, strlen(path)+1, char);
  strcpy(path1, path);
  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);
  cs_xpath_add_attribute(&path, "label");
  label = cs_gui_get_attribute_value(path);

  if (cs_gui_strcmp(type, "volume")) {

    cs_xpath_add_element(&path1, "postprocessing_recording");
    cs_xpath_add_attribute(&path1, "status");
    if (cs_gui_get_status(path1, &result))
      *record_value = result;
  }

  else if (cs_gui_strcmp(type, "boundary")) {

    cs_xpath_add_element(&path2, "postprocessing_recording");
    cs_xpath_add_attribute(&path2, "status");
    if (cs_gui_get_status(path2, &result))
      *record_value = result;
  }

  BFT_FREE(path);
  BFT_FREE(path1);
  BFT_FREE(path2);

  return label;
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
  int icoal, ncoals = 0;
  int flag = 0;
  char *attr = NULL;
  char *path1 = NULL;

  attr = _get_attr("model", 1, "lagrangian");
  if (attr == NULL || cs_gui_strcmp(attr, "off")) {
    cs_glob_lagr_time_scheme->iilagr = 0;
    BFT_FREE(attr);
#if _XML_DEBUG_
    bft_printf("%s\n", __func__);
    bft_printf("--iilagr = %i\n", cs_glob_lagr_time_scheme->iilagr);
#endif
    return;
  }
  else if (cs_gui_strcmp(attr, "one_way"))
    cs_glob_lagr_time_scheme->iilagr = 1;
  else if (cs_gui_strcmp(attr, "two_way"))
    cs_glob_lagr_time_scheme->iilagr = 2;
  else if (cs_gui_strcmp(attr, "frozen"))
    cs_glob_lagr_time_scheme->iilagr = 3;
  BFT_FREE(attr);

  /* Global settings */

  _get_status(&(cs_glob_lagr_time_scheme->isuila), 2,
              "lagrangian", "restart");
  _get_status(&(cs_glob_lagr_time_scheme->isttio), 2,
              "lagrangian", "carrier_field_stationary");

  _get_status(&(cs_glob_lagr_model->deposition), 2,
              "lagrangian", "deposition_submodel");

  /* Particles model */

  cs_glob_lagr_model->physical_model = _get_particles_model("particles_models");

  switch (cs_glob_lagr_model->physical_model) {
  case 1:
    _get_status(&(cs_glob_lagr_specific_physics->idpvar), 3,
                "lagrangian", "particles_models", "break_up");
    _get_status(&(cs_glob_lagr_specific_physics->impvar), 3,
                "lagrangian", "particles_models", "evaporation");
    _get_status(&(cs_glob_lagr_specific_physics->itpvar), 3,
                "lagrangian", "particles_models", "thermal");
#if 0
    if (cs_glob_lagr_specific_physics->itpvar == 1) {
      _get_double(tpart,  4, "lagrangian", "particles_models", "thermal",
                  "particle_temperature");
      _get_double(cppart, 4, "lagrangian", "particles_models",
                  "thermal", "particle_specific_heat");
    }
#endif
    break;
  case 2:
    _get_status(&cs_glob_lagr_model->fouling, 3, "lagrangian", "particles_models",
                "coal_fouling");
    path1 = cs_xpath_init_path();
    cs_xpath_add_elements(&path1, 4, "lagrangian", "particles_models",
                          "coal_fouling", "threshold_temperature");
    ncoals = cs_gui_get_nb_element(path1);
    BFT_FREE(path1);

    for (icoal =1; icoal <= ncoals; icoal++) {
      _get_coal_double(&cs_glob_lagr_encrustation->tprenc[icoal-1],
                       "threshold_temperature", icoal);

      _get_coal_double(&cs_glob_lagr_encrustation->visref[icoal-1],
                       "critical_viscosity",    icoal);

      _get_coal_double(&cs_glob_lagr_encrustation->enc1[icoal-1],
                       "fouling_coefficient_1", icoal);

      _get_coal_double(&cs_glob_lagr_encrustation->enc2[icoal-1],
                       "fouling_coefficient_2", icoal);
    }
    break;
  }

  /* Two-way coupling */

  if (cs_glob_lagr_time_scheme->iilagr == 2) {
    _get_int(&cs_glob_lagr_source_terms->nstits, 3, "lagrangian",
             "two_way_coupling", "iteration_start");
    _get_status(&cs_glob_lagr_source_terms->ltsdyn, 3, "lagrangian",
                "two_way_coupling", "dynamic");
    _get_status(&cs_glob_lagr_source_terms->ltsmas, 3, "lagrangian",
                "two_way_coupling", "mass");
    _get_status(&cs_glob_lagr_source_terms->ltsthe, 3, "lagrangian",
                "two_way_coupling", "thermal");
  }

  /* Numerical modeling */

  attr = _get_attr("choice", 2, "lagrangian", "scheme_order");
  if (attr) {
    cs_glob_lagr_time_scheme->t_order = atoi(attr);
    BFT_FREE(attr);
  }
  attr = _get_attr("choice", 2, "lagrangian", "complete_model_direction");
  if (attr) {
    cs_glob_lagr_time_scheme->idirla = atoi(attr);
    BFT_FREE(attr);
  }
  _get_status(&cs_glob_lagr_time_scheme->idistu, 2, "lagrangian",
              "turbulent_dispersion");
  _get_status(&cs_glob_lagr_time_scheme->idiffl, 2, "lagrangian",
              "fluid_particles_turbulent_diffusion");
  _get_int(&cs_glob_lagr_time_scheme->modcpl, 2, "lagrangian",
           "complete_model");

  /* Output */

  _attr_post_status(CS_LAGR_VELOCITY, "velocity_particles");
  _attr_post_status(CS_LAGR_VELOCITY_SEEN, "velocity_fluid_seen");
  _attr_post_status(CS_LAGR_RESIDENCE_TIME, "resident_time");
  _attr_post_status(CS_LAGR_DIAMETER, "diameter");
  _attr_post_status(CS_LAGR_TEMPERATURE, "temperature");
  _attr_post_status(CS_LAGR_MASS, "mass");

  if (cs_glob_lagr_model->physical_model == 2) {
    _attr_post_status(CS_LAGR_SHRINKING_DIAMETER, "shrinking_core_diameter");
    _attr_post_status(CS_LAGR_WATER_MASS, "moisture_mass_fraction");
    _attr_post_status(CS_LAGR_COAL_MASS, "raw_coal_mass_fraction");
    _attr_post_status(CS_LAGR_COKE_MASS, "char_mass_fraction");
  }

  _get_int(&cs_glob_lagr_log_frequency_n,
           3, "lagrangian", "output", "listing_printing_frequency");

  /* Statistics */

  _get_int(&cs_glob_lagr_model->n_stat_classes, 3, "lagrangian",
           "statistics", "statistics_groups_of_particles");
  _get_status(&cs_glob_lagr_stat_options->isuist, 3, "lagrangian",
              "statistics", "restart");

  _get_double(&cs_glob_lagr_stat_options->threshold, 3, "lagrangian",
              "statistics", "threshold");

  _get_int(&cs_glob_lagr_stat_options->idstnt, 3, "lagrangian",
           "statistics", "iteration_start");
  _get_int(&cs_glob_lagr_stat_options->nstist, 3, "lagrangian",
           "statistics", "iteration_steady_start");

  int vol_stats = 0;
  _get_status(&vol_stats, 3, "lagrangian", "statistics", "volume");

  if (vol_stats == 1) {

    /* labels */

    flag = 0;
    _get_char_post("volume", "Part_vol_frac", &flag);
    if (flag)
      cs_lagr_stat_activate(CS_LAGR_STAT_VOLUME_FRACTION);

    _get_char_post("volume", "Part_velocity", &flag);
    if (flag)
      cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

    _get_char_post("volume", "Part_resid_time", &flag);
    if (flag)
      cs_lagr_stat_activate_attr(CS_LAGR_RESIDENCE_TIME);

    _get_char_post("volume", "Part_stat_weight", &flag);
    if (flag)
      cs_lagr_stat_activate(CS_LAGR_STAT_CUMULATIVE_WEIGHT);

  }

  int b_key = 0;
  _get_status(&b_key, 3, "lagrangian", "statistics", "boundary");

  if (b_key) {

    _get_char_post("boundary", "Part_impact_number",
                   &cs_glob_lagr_boundary_interactions->inbrbd);
    _get_char_post("boundary", "Part_bndy_mass_flux",
                   &cs_glob_lagr_boundary_interactions->iflmbd);
    _get_char_post("boundary", "Part_impact_angle",
                   &cs_glob_lagr_boundary_interactions->iangbd);
    _get_char_post("boundary", "Part_impact_velocity",
                   &cs_glob_lagr_boundary_interactions->ivitbd);

    /* Coal fouling statistics*/
    _get_char_post("boundary", "Part_fouled_impact_number",
                   &cs_glob_lagr_boundary_interactions->iencnbbd);
    _get_char_post("boundary", "Part_fouled_mass_flux",
                   &cs_glob_lagr_boundary_interactions->iencmabd);
    _get_char_post("boundary", "Part_fouled_diam",
                   &cs_glob_lagr_boundary_interactions->iencdibd);
    _get_char_post("boundary", "Part_fouled_Xck",
                   &cs_glob_lagr_boundary_interactions->iencckbd);

  }

#if _XML_DEBUG_
  bft_printf("==>UILAG1\n");
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
    for (icoal=1; icoal <= ncoals; icoal++)
    {
      bft_printf("--tprenc[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->tprenc[icoal-1]);
      bft_printf("--visref[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->visref[icoal-1]);
      bft_printf("--enc1[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->enc1[icoal-1]);
      bft_printf("--enc2[%i] = %f\n", icoal,
                 cs_glob_lagr_encrustation->enc2[icoal-1]);
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

  bft_printf("--ivisv1 = %i\n", lagr_post_options->ivisv1);
  bft_printf("--ivisv2 = %i\n", lagr_post_options->ivisv2);
  bft_printf("--ivistp = %i\n", lagr_post_options->ivistp);
  bft_printf("--ivisdm = %i\n", lagr_post_options->ivisdm);
  bft_printf("--iviste = %i\n", lagr_post_options->iviste);
  bft_printf("--ivismp = %i\n", lagr_post_options->ivismp);

  if (cs_glob_lagr_model->physical_model == 2) {
    bft_printf("--ivisdk  = %i\n", lagr_post_options->ivisdk);
    bft_printf("--iviswat = %i\n", lagr_post_options->iviswat);
    bft_printf("--ivisch  = %i\n", lagr_post_options->ivisch);
    bft_printf("--ivisck  = %i\n", lagr_post_options->ivisck);
  }

  bft_printf("--isuist = %i\n", cs_glob_lagr_stat_options->isuist);
  bft_printf("--nbclst = %i\n", cs_glob_lagr_model->n_stat_classes);

  bft_printf("--idstnt = %i\n", cs_glob_lagr_stat_options->idstnt);
  bft_printf("--nstist = %i\n", cs_glob_lagr_stat_options->nstist);
  bft_printf("--vol_stats = %i\n", vol_stats);

  bft_printf("--boundary_output = %i\n", b_key);
  if (b_key == 1) {
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
  char *interaction = NULL;
  char sclass[10];
  char *path1, *path2;
  char *choice;

  cs_lnum_t iphyla = cs_glob_lagr_model->physical_model;
  cs_lagr_zone_data_t *bdy_cond = cs_lagr_get_boundary_conditions();
  cs_lagr_get_internal_conditions();

  /* zone 0 for "all", next zones defined by GUI */

  int n_zones = cs_gui_boundary_zones_number();

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
#endif

  /* First iteration only: memory allocation */

  for (cs_lnum_t izone = 0; izone < n_zones; izone++) {

    /* zone_id is incremented because access is made through real
       zone number (from 1 to nzone), and not through array index */

    int zone_id = izone + 1;

    char *label = cs_gui_boundary_zone_label(zone_id);
    char *nature = cs_gui_boundary_zone_nature(zone_id);

    path2 = cs_xpath_init_path();
    cs_xpath_add_elements(&path2, 2, "boundary_conditions",
                          nature);
    cs_xpath_add_test_attribute(&path2, "label", label);
    cs_xpath_add_test_attribute(&path2, "field_id", "none");
    cs_xpath_add_element(&path2, "particles");

    BFT_MALLOC(path1, strlen(path2)+1, char);
    strcpy(path1, path2);
    cs_xpath_add_attribute(&path1, "choice");
    interaction = cs_gui_get_attribute_value(path1);

    if (interaction != NULL) {

      if (cs_gui_strcmp(interaction, "inlet"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_INLET;

      else if(cs_gui_strcmp(interaction, "outlet"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_OUTLET;

      else if(cs_gui_strcmp(interaction, "bounce"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_REBOUND;

      else if(cs_gui_strcmp(interaction, "part_symmetry"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_SYM;

      else if(cs_gui_strcmp(interaction, "deposit1"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO1;

      else if(cs_gui_strcmp(interaction, "deposit2"))
        bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO2;

      else if(cs_gui_strcmp(interaction, "fouling") && iphyla == 2)
        bdy_cond->zone_type[zone_id] = CS_LAGR_FOULING;

      else if(   cs_gui_strcmp(interaction, "fouling")
              && (iphyla == 0  || iphyla == 1))
        bdy_cond->zone_type[zone_id] = CS_LAGR_DEPO_DLVO;

#if _XML_DEBUG_
      bft_printf("--zone_type[%i] = %i has %i set(s) \n", izone,
                 bdy_cond->zone_type[zone_id],
                 bdy_cond->b_zone_sets[zone_id]);

      bft_printf("--        : label    %s \n", label);
      bft_printf("--        : nature   %s \n", nature);
      bft_printf("--        : p_nature %i \n", bdy_cond->zone_type[zone_id]);
#endif

      /* Additional info for inlet */

      if (bdy_cond->zone_type[zone_id] == CS_LAGR_INLET) {

        strcpy(path1, path2);
        cs_xpath_add_element(&path1, "class");
        int n_injection_sets = cs_gui_get_nb_element(path1);
        strcpy(path1, path2);

        for (int set_id = 0; set_id < n_injection_sets; set_id++) {

          cs_lagr_injection_set_t *zis
            = cs_lagr_get_injection_set(bdy_cond, zone_id, set_id);

          cs_lagr_injection_set_default(zis);

          sprintf(sclass, "class[%i]", set_id+1);
          BFT_REALLOC(path2,
                      ( 20+strlen(nature)
                        +10+strlen(label)
                        +13+strlen(sclass)+1),
                      char);
          strcpy(path2, "");
          sprintf(path2,
                  "boundary_conditions/%s[@label='%s']/particles/%s",
                  nature,
                  label,
                  sclass);

          int itmp0, itmp1, itmp2;
          cs_real_t rtmp0 = 0., rtmp1 = 0., rtmp2;
          _get_int(&(itmp0), 2, path2, "number");
          _get_int(&(itmp1), 2, path2, "frequency");
          _get_int(&(itmp2), 2, path2, "statistical_groups");

          zis->n_inject = itmp0;
          zis->injection_frequency = itmp1;
          zis->cluster = itmp2;

#if _XML_DEBUG_
          bft_printf("---number = %i \n", itmp0);
          bft_printf("---frequency = %i \n", itmp1);
          bft_printf("---statistical_groups = %i \n", itmp2);
#endif
          /* velocity */

          choice = _get_attr("choice", 2, path2, "velocity");

          cs_real_t vel[3] = {0, 0, 0};

          if (cs_gui_strcmp(choice, "fluid"))
            zis->velocity_profile = -1;

          else if (cs_gui_strcmp(choice, "norm")) {
            _get_double(&(vel[0]), 3, path2, "velocity", "norm");
            zis->velocity_profile = 0;
            zis->velocity_magnitude = vel[0];
          }
          else if (cs_gui_strcmp(choice, "components")) {
            zis->velocity_profile = 1;
            _get_double(&(vel[0]), 3, path2, "velocity", "velocity_x");
            _get_double(&(vel[1]), 3, path2, "velocity", "velocity_y");
            _get_double(&(vel[2]), 3, path2, "velocity", "velocity_z");
            for (int i = 0; i < 3; i++)
              zis->velocity[i] = vel[i];
          }
          else if (cs_gui_strcmp(choice, "subroutine"))
            zis->velocity_profile = 2;

#if _XML_DEBUG_
          bft_printf("---velocity choice: %i "
                     " (-1: fluid, 0: norm, 1: components, 2: subroutine)\n",
                     itmp0);

          if (itmp0 == 0)

            bft_printf("----norm = %f \n", vel[0]);

          else if (itmp0 == 1) {

            bft_printf("----u = %f \n", vel[0]);
            bft_printf("----v = %f \n", vel[1]);
            bft_printf("----w = %f \n", vel[2]);
          }
#endif
         BFT_FREE(choice);

          /* statistical_weight, mass_flow_rate*/

          choice = _get_attr("choice", 2, path2, "statistical_weight");

          if (cs_gui_strcmp(choice, "rate")) {
            rtmp0 = 0;
            _get_double(&rtmp1, 2, path2, "mass_flow_rate");
            zis->stat_weight = 0;
            zis->flow_rate = rtmp1;
          }
          else { /* if (cs_gui_strcmp(choice, "prescribed")) */
            _get_double(&rtmp0, 2, path2, "statistical_weight");
            zis->stat_weight = rtmp0;
            zis->flow_rate = 0;
          }

#if _XML_DEBUG_
          bft_printf("---statistical weight choice: %i "
                     " (1: prescribed, 2: rate)\n", itmp0);

          if (itmp0 == 1 || itmp0 == 2) {
            bft_printf("----statistical weight = %f \n", rtmp0);
            bft_printf("----mass flow rate = %f \n", rtmp1);
          }
#endif
          BFT_FREE(choice);

          /* diameter */

          choice = _get_attr("choice", 2, path2, "diameter");

          _get_double(&rtmp0, 2, path2, "diameter");
          _get_double(&rtmp1, 2, path2, "diameter_standard_deviation");
          zis->diameter = rtmp0;
          zis->diameter_variance = rtmp1;

#if _XML_DEBUG_
          if (itmp0 == 1) {
            bft_printf("----diameter = %f \n", rtmp0);
            bft_printf("----standard deviation = %f \n", rtmp1);
          }
#endif
          BFT_FREE(choice);

          /* density */
          if (iphyla != 2) {

            _get_double(&rtmp0, 2, path2, "density");

            zis->density = rtmp0;

#if _XML_DEBUG_
            bft_printf("---density = %f \n", rtmp0);
#endif

          }

          /* Fouling index*/

          _get_double(&rtmp0, 2, path2, "fouling_index");
          zis->fouling_index = rtmp0;

          if (iphyla == 1) {

            /* temperature, specific_heat, emissivity */

            choice = _get_attr("choice", 2, path2, "temperature");

            if (cs_gui_strcmp(choice, "prescribed")) {
              zis->temperature_profile = 1;
              _get_double(&rtmp0, 2, path2, "temperature");
              zis->temperature = rtmp0;
            }
            else if (cs_gui_strcmp(choice, "fluid")) {
              zis->temperature_profile = 0;
              zis->temperature = 0;
            }

            _get_double(&rtmp1, 2, path2, "specific_heat");
            _get_double(&rtmp2, 2, path2, "emissivity");
            zis->cp = rtmp1;
            zis->emissivity = rtmp2;

#if _XML_DEBUG_
            bft_printf("---temperature choice = %i "
                       "(0: fluid, 1: prescribed)\n",
                       itmp0);

            if (itmp0 == 1)
              bft_printf("----temperature = %f \n", rtmp0);

            bft_printf("---specific heat = %f \n", rtmp1);
            bft_printf("---emissivity = %f \n", rtmp2);
#endif
            BFT_FREE(choice);

          }

          /* coal */
          else if (iphyla == 2) {

            /* Read the coal number */

            _get_int(&itmp0, 2, path2, "coal_number");
            zis->coal_number = itmp0;

            /* Data are read in pulverized fuel combustion module profile */

            _get_double(&rtmp0, 2, path2, "coal_temperature");
            zis->temperature = rtmp0;

#if _XML_DEBUG_
            bft_printf("---coal number = %i \n", itmp0);
            bft_printf("---coal composition = %i "
                       "(1: raw coal, 2: user defined)\n", itmp1);
#endif /* _XML_DEBUG_ */

          }

        } /* End of loop on set */

      } /* End of test on inlet */

    } /* End of test on interaction */

    BFT_FREE(path1);
    BFT_FREE(path2);
    BFT_FREE(interaction);

    BFT_FREE(label);
    BFT_FREE(nature);

  } /* End of loop on zones */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
