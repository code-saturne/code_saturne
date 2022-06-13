/*============================================================================
 * Base gas mix data.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_log.h"

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gas_mix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_gas_mix.c
        Base gas mix data.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_gas_mix_t

  \brief Gas mix descriptor.

  Members of this structure are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_gas_mix_t::n_species
        number of species in the gas mix
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main gas mix options and associated pointer */

static cs_gas_mix_t _gas_mix = {
  .n_species = 0,
  .n_species_solved = 0,
  .species_to_field_id = NULL
};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default variable compute options */

static cs_gas_mix_species_prop_t _gas_mix_species_prop =
{
  -1.,   /* molar mass */
  -1.,   /* specific heat */
  -1.,   /* volume diffusion */
  -1.,   /* dynamic viscosity a */
  -1.,   /* dynamic viscosity b */
  -1.,   /* thermal conductivity a */
  -1.,   /* thermal conductivity b */
  -1.,   /* reference viscosity (Sutherland) */
  -1.,   /* reference conductivity (Sutherland) */
  -1.,   /* reference temperature for viscosity */
  -1.,   /* reference temperature for conductivity */
  -1.,   /* Sutherland temperature for viscosity */
  -1.,   /* Sutherland temperature for conductivity */
};

/*============================================================================
 * Global variables
 *============================================================================*/

const cs_gas_mix_t  *cs_glob_gas_mix = &_gas_mix;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_gas_mix_get_pointers(int **nscasp);

int
cs_f_gas_mix_species_to_field_id(int sp_id);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log values of the structure */

static void
_log_func_gas_mix_species_prop(const void *t)
{
  const char fmt[] = N_("      %-19s  %-12.3g\n");
  const cs_gas_mix_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt, "mol_mas ", _t->mol_mas);
  cs_log_printf(CS_LOG_SETUP, fmt, "cp      ", _t->cp);
  cs_log_printf(CS_LOG_SETUP, fmt, "vol_diff", _t->vol_dif);
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_a    ", _t->mu_a);
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_b    ", _t->mu_b);
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_a", _t->lambda_a);
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_b", _t->lambda_b);
  cs_log_printf(CS_LOG_SETUP, fmt, "muref   ", _t->muref);
  cs_log_printf(CS_LOG_SETUP, fmt, "lamref  ", _t->lamref);
  cs_log_printf(CS_LOG_SETUP, fmt, "trefmu  ", _t->trefmu);
  cs_log_printf(CS_LOG_SETUP, fmt, "treflam ", _t->treflam);
  cs_log_printf(CS_LOG_SETUP, fmt, "smu     ", _t->smu);
  cs_log_printf(CS_LOG_SETUP, fmt, "slam    ", _t->slam);
}

/* Log default values of the structure */

static void
_log_func_default_gas_mix_species_prop(const void *t)
{
  const char fmt[] = "      %-19s  %-12.3g %s\n";
  const cs_gas_mix_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt, "mol_mas ", _t->mol_mas,
                _("Molar mass"));
  cs_log_printf(CS_LOG_SETUP, fmt, "cp      ", _t->cp,
                _("Specific heat"));
  cs_log_printf(CS_LOG_SETUP, fmt, "vol_diff", _t->vol_dif,
                _("Volume diffusion"));
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_a    ", _t->mu_a,
                _("Dynamic viscosity a"));
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_b    ", _t->mu_b,
                _("Dynamic viscosity b"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_a", _t->lambda_a,
                _("Thermal conductivity a"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_b", _t->lambda_b,
                _("Thermal conductivity b"));
  cs_log_printf(CS_LOG_SETUP, fmt, "muref   ", _t->muref,
                _("Reference thermal viscosity (Sutherland)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lamref  ", _t->lamref,
                _("Reference thermal conductivity (Sutherland)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "trefmu  ", _t->trefmu,
                _("Reference temperature (Sutherland for viscosity)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "treflam ", _t->treflam,
                _("Reference temperature (Sutherland conductivity)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "smu     ", _t->smu,
                _("Sutherland temperature for viscosity"));
  cs_log_printf(CS_LOG_SETUP, fmt, "slam    ", _t->slam,
                _("Sutherland temperature for conductivity"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * The species to field mapping is usually done in the order of addition
 * of species, but non-resolved (deduced) fields are placed last in the
 * list; the order of solved species relative to each other is maintained,
 * as is that of deduced species. In most cases, we have at most one deduced
 * species, which is mapped last.
 *
 * \param[in]   f   pointer to the field
 */
/*----------------------------------------------------------------------------*/

static void
_map_field(const cs_field_t *f)
{
  /* Check if already mapped */

  for (int i = 0; i < _gas_mix.n_species; i++) {
    if (_gas_mix.species_to_field_id[i] == f->id)
      return;
  }

  int is_solved = f->type & CS_FIELD_VARIABLE;
  int insert_id = (is_solved) ?
    _gas_mix.n_species_solved : _gas_mix.n_species;

  int n_species_ini = _gas_mix.n_species;
  _gas_mix.n_species++;

  if (is_solved)
    _gas_mix.n_species_solved += 1;

  BFT_REALLOC(_gas_mix.species_to_field_id, _gas_mix.n_species, int);

  /* If we need to insert a solved variable with non-solved fields
     already mapped, shift the non-solved field map positions */

  for (int i = n_species_ini; i > insert_id; i--)
    _gas_mix.species_to_field_id[i] = _gas_mix.species_to_field_id[i-1];

  _gas_mix.species_to_field_id[insert_id] = f->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set predefined gas properties for a field based on its name.
 *
 * This only applies to fields for which properties have not already been
 * set (i.e. for which \ref cs_gas_mix_add_species_with_properties
 * has not been called), so as to complete the property definitions
 * based on predefined values.
 *
 * \param[in, out]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

static void
_set_predefined_property(cs_field_t  *f)
{
  int k_id = cs_gas_mix_get_field_key();

  cs_gas_mix_species_prop_t gmp;
  cs_field_get_key_struct(f, k_id, &gmp);

  if (gmp.mol_mas > 0)
    return;

  if (strcmp(f->name, "y_o2") == 0) {
    gmp.mol_mas= 0.032;
    gmp.cp = 930.;
    gmp.vol_dif = 19.70;
    gmp.mu_a = 5.086522e-8;
    gmp.mu_b = 5.512391e-6;
    gmp.lambda_a = 6.2e-5;
    gmp.lambda_b = 8.1e-3;
    gmp.muref = 1.919e-5;
    gmp.lamref = 0.0244;
    gmp.trefmu = 273.;
    gmp.treflam = 273.;
    gmp.smu = 139.;
    gmp.slam = 240.;
  }

  else if (strcmp(f->name, "y_n2") == 0) {
    gmp.mol_mas = 0.028;
    gmp.cp = 1042.;
    gmp.vol_dif = 19.70;
    gmp.mu_a = 4.210130e-8;
    gmp.mu_b = 5.494348e-6;
    gmp.lambda_a = 6.784141e-5;
    gmp.lambda_b = 5.564317e-3;
    gmp.muref = 1.663e-5;
    gmp.lamref = 0.0242;
    gmp.trefmu = 273.;
    gmp.treflam = 273.;
    gmp.smu = 107.;
    gmp.slam = 150.;
  }

  else if (strcmp(f->name, "y_he") == 0) {
    gmp.mol_mas = 0.004;
    gmp.cp = 5194.;
    gmp.vol_dif = 2.67;
    gmp.mu_a = 18.5752e-6;
    gmp.mu_b = 0.0;
    gmp.lambda_a = 0.144;
    gmp.lambda_b = 0.0;
    gmp.muref = 1.874e-5;
    gmp.lamref = 0.647;
    gmp.trefmu = 273.;
    gmp.treflam = 273.;
    gmp.smu = 78.;
    gmp.slam = 78.;
  }

  else if (strcmp(f->name, "y_h2") == 0) {
    gmp.mol_mas = 0.002;
    gmp.cp = 14560.;
    gmp.vol_dif = 6.12;
    gmp.mu_a = 1.93e-9;
    gmp.mu_b = 8.40e-6;
    gmp.lambda_a = 4.431e-4;
    gmp.lambda_b = 5.334e-2;
    gmp.muref = 8.411e-6;
    gmp.lamref = 0.0168;
    gmp.trefmu = 273.;
    gmp.treflam = 273.;
    gmp.smu = 97.;
    gmp.slam = 120.;
  }

  else if (strcmp(f->name, "y_h2o_g") == 0) {
    gmp.mol_mas = 0.018;
    gmp.cp = 2060.;
    gmp.vol_dif = 13.1;
    gmp.mu_a = 3.8496e-8;
    gmp.mu_b = 8.2997e-6;
    gmp.lambda_a = 7.6209e-5;
    gmp.lambda_b = 0.016949;
    gmp.muref = 1.12e-5;
    gmp.lamref = 0.0181;
    gmp.trefmu = 350.;
    gmp.treflam = 300.;
    gmp.smu = 1064.;
    gmp.slam = 2200.;
  }

  else {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: no predefined properties for field %s."),
              __func__, f->name);

  }

  cs_field_set_key_struct(f, k_id, &gmp);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global gas mix structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nscasp --> pointer to cs_glob_gas_mix->n_species
 *----------------------------------------------------------------------------*/

void
cs_f_gas_mix_get_pointers(int  **nscasp)
{
  *nscasp = &(_gas_mix.n_species_solved);
}

/*----------------------------------------------------------------------------
 * Return field id matching a given species id.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   sp_id --> matching species id (1-based,
 *             deduced species if sp_id == nscasp + 1)
 *----------------------------------------------------------------------------*/

int
cs_f_gas_mix_species_to_field_id(int sp_id)
{
  int retval = -1;

  int _sp_id = sp_id-1;

  if (_sp_id >= 0 && _sp_id < _gas_mix.n_species)
    retval = _gas_mix.species_to_field_id[_sp_id];

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the field key for gas mix properties.
 *
 * \return  field key id for gas mix properties
 */
/*----------------------------------------------------------------------------*/

int
cs_gas_mix_get_field_key(void)
{
  static int k_id = -1;

  /* Structure containing physical properties relative to
     species scalars used by the gas mixture modelling */

  if (k_id < 0) {

    const char key[] = "gas_mix_species_prop";

    cs_field_define_key_struct(key,
                               &_gas_mix_species_prop,
                               _log_func_gas_mix_species_prop,
                               _log_func_default_gas_mix_species_prop,
                               NULL,
                               sizeof(cs_gas_mix_species_prop_t),
                               0);

    k_id = cs_field_key_id_try(key);

  }

  return k_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species(int  f_id)
{
  if (cs_glob_physical_model_flag[CS_GAS_MIX] == -1)
    bft_error(__FILE__, __LINE__, 0,
              _("No gas species can be added."
                " The gas mix model is not enabled.\n"));

  cs_field_t *f = cs_field_by_id(f_id);

  if (   strcmp(f->name, "y_o2") != 0
      && strcmp(f->name, "y_n2") != 0
      && strcmp(f->name, "y_he") != 0
      && strcmp(f->name, "y_h2") != 0)
                bft_error(__FILE__, __LINE__, 0,
                          _("Only the species having the following field names "
                            "can be added to a gas mix:\n"
                            "y_o2, y_n2, y_he, y_h2\n"));

  _map_field(f);

  _set_predefined_property(f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]  f_id         id of field representing species mixture fraction.
 * \param[in]  mol_mass     molar mass
 * \param[in]  cp           specific heat
 * \param[in]  col_diff     volume diffusion
 * \param[in]  mu_a         dynamic viscosity a
 * \param[in]  mu_b         dynamic viscosity b
 * \param[in]  lambda_a     thermal conductivity a
 * \param[in]  lambda_b     thermal conductivity b
 * \param[in]  mu_ref       reference viscosity (Sutherland)
 * \param[in]  lambda_ref   reference conductivity (Sutherland)
 * \param[in]  tref_mu      reference temperature for viscosity
 * \param[in]  tref_lambda  reference temperature for conductivity
 * \param[in]  s_mu         Sutherland temperature for viscosity
 * \param[in]  s_lambda     Sutherland temperature for conductivity
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species_with_properties(int        f_id,
                                       cs_real_t  mol_mass,
                                       cs_real_t  cp,
                                       cs_real_t  vol_diff,
                                       cs_real_t  mu_a,
                                       cs_real_t  mu_b,
                                       cs_real_t  lambda_a,
                                       cs_real_t  lambda_b,
                                       cs_real_t  mu_ref,
                                       cs_real_t  lambda_ref,
                                       cs_real_t  tref_mu,
                                       cs_real_t  tref_lambda,
                                       cs_real_t  s_mu,
                                       cs_real_t  s_lambda)
{
  cs_field_t *f = cs_field_by_id(f_id);

  _map_field(f);

  cs_gas_mix_species_prop_t gmp
    = {.mol_mas = mol_mass,
       .cp = cp,
       .vol_dif = vol_diff,
       .mu_a = mu_a,
       .mu_b = mu_b,
       .lambda_a = lambda_a,
       .lambda_b = lambda_b,
       .muref = mu_ref,
       .lamref = lambda_ref,
       .trefmu = tref_mu,
       .treflam = tref_lambda,
       .smu = s_mu,
       .slam = s_lambda};

  int k_id = cs_gas_mix_get_field_key();

  cs_field_set_key_struct(f, k_id, &gmp);

  /* Set model flag.
     Reserve lower values (0-5 currently used) for predefined cases;
     TODO: only 0/1 should be needed here, using specific functions
     or an ENUM for pre-definition */

  cs_glob_physical_model_flag[CS_GAS_MIX] = CS_GAS_MIX_USER;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add property fields specific to a gas mix.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_property_fields(void)
{
  cs_field_t *f;

  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  bool add_deduced = true;

  char name[32], label[32];

  switch(cs_glob_physical_model_flag[CS_GAS_MIX]) {

  case CS_GAS_MIX_AIR_HELIUM:
    strncpy(name, "y_he", 32);
    strncpy(label, "Y_he", 32);
    break;

  case CS_GAS_MIX_AIR_HYDROGEN:
    strncpy(name, "y_h2", 32);
    strncpy(label, "Y_H2", 32);
    break;

  case CS_GAS_MIX_AIR_STEAM:
  case CS_GAS_MIX_AIR_HELIUM_STEAM:
  case CS_GAS_MIX_AIR_HYDROGEN_STEAM:
    strncpy(name, "y_h2o_g", 32);
    strncpy(label, "Y_H20_g", 32);
    break;

  case CS_GAS_MIX_HELIUM_AIR:
    strncpy(name, "y_o2", 32);
    strncpy(label, "Y_O2", 32);
    break;

  default:
    add_deduced = false;  /* should be added by user */

  }

  if (add_deduced) {

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    f = cs_field_create(name,
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        true);

    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, label);

    _map_field(f);
    _set_predefined_property(f);

    /* Add binary diffusion coefficient of steam into non-condensables */
    f = cs_field_create("steam_binary_diffusion",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        false);

    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);

    /* Add molecular weight of non-condensable mixture */
    f = cs_field_create("mol_mas_ncond",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        false);

    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);

    /* Add molecular weight of non-condensable mixture */
    f = cs_field_create("tempk",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        false);

    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);

  }

  /* Add molar mass property */

  f = cs_field_create("mix_mol_mas",
                      field_type,
                      CS_MESH_LOCATION_CELLS,
                      1, /* dim */
                      false);

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free array mapping gas mix species ids to field ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_finalize(void)
{
  BFT_FREE(_gas_mix.species_to_field_id);
  _gas_mix.n_species = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
