/*============================================================================
 * Base gas mix data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <math.h>
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

#include "cs_array.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"

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
  .species_to_field_id = nullptr
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
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dynamic viscosity and
 *        conductivity coefficient associated to each gas species.
 *
 * TODO: replace name with enum to avoid test on string, as this
 *       function is called in a loop.
 *
 * \param[in]     name          name of the field associated to the gas species
 * \param[in]     tk            temperature variable in kelvin
 * \param[in]     spro          constants used for the physcial laws
 * \param[out]    mu            dynamic viscosity associated to the gas species
 * \param[out]    lambda        conductivity coefficient of the gas species
 */
/*----------------------------------------------------------------------------*/

static void
_compute_mu_lambda(const char                       *name,
                   const cs_real_t                   tk,
                   const cs_gas_mix_species_prop_t   spro,
                   cs_real_t                        *mu,
                   cs_real_t                        *lambda)
{
  cs_real_t _mu;
  cs_real_t _lambda;
  const cs_real_t tkelvin = cs_physical_constants_celsius_to_kelvin;

  /*  The viscosity law  and conductivity expressionfor each species is
   * defined as:
   *
   * 1/. for Steam species:
   *      mu = mu_a.(tk - tkelvi), with a law in (°C) unit
   *      lambda = lambda_a .(tk-tkelvi) + lambda_b, with a law in (°C) unit
   * 2/. for Helium species:
   *     mu = mu_a .(tk/tkelvi)**c + mu_b, with nodim. law
   *     lambda = lambda_a .(tk/tkelvi)**c + lambda_b, with nodim. law
   * 3/. for Hydrogen species:
   *      mu = mu_a .(tk-tkelvi) + mu_b, with t (°C)
   *      lambda = lambda_a .tk + lambda_b, with tk (°K)
   * 4/. for Oxygen and Nitrogen species :
   *      mu = mu_a .(tk) + mu_b, with t in (°K)
   *      lambda = lambda_a .(tk) + lambda_b, with t (°K) */

  if (strcmp(name, "y_h2o_g") == 0) {
    _mu     = spro.mu_a    *(tk-tkelvin) + spro.mu_b;
    _lambda = spro.lambda_a*(tk-tkelvin) + spro.lambda_b;
  }
  else if (strcmp(name, "y_he") == 0) {
    _mu     = spro.mu_a     * pow(tk/tkelvin, 0.7);
    _lambda = spro.lambda_a * pow(tk/tkelvin, 0.7);
  }
  else if (strcmp(name, "y_h2") == 0) {
    _mu     = spro.mu_a     * (tk-tkelvin) + spro.mu_b;
    _lambda = spro.lambda_a * tk + spro.lambda_b;
  }
  else if (   strcmp(name, "y_o2") == 0
           || strcmp(name, "y_n2") == 0) {
    _mu     = spro.mu_a     * tk + spro.mu_b;
    _lambda = spro.lambda_a * tk + spro.lambda_b;
  }
  else {
    _mu = -1;
    _lambda = -1;
    bft_error(__FILE__, __LINE__, 0,
              _("%s: no predefined properties for field %s."),
              __func__, name);
  }

  *mu = _mu;
  *lambda = _lambda;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sutherland law for viscosity and thermal conductivity
 *
 * The viscosity law for each specie is defined as below:
 *
 * mu = mu_ref*(T/T_ref)**(3/2)*(T_ref+S1)/(T+S1)
 *
 * The conductivity expression for each specie is defined as:
 *
 *  lambda = lambda_ref*(T/T_ref)**(3/2)*(T_ref+S2)/(T+S2)
 *
 * S1 and S2 are respectively Sutherland temperature for conductivity and
 * Sutherland temperature for viscosity of the considered specie
 * T_ref is a reference temperature, equal to 273K for a perfect gas.
 * For steam (H20), Tref has not the same value in the two formulae.
 * Available species : O2, N2, H2, H20 and  He
 * The values for the parameters come from F.M. White's book
 * "Viscous Fluid Flow"
 *
 * \param[in]     tk            temperature variable in kelvin
 * \param[in]     spro          constants used for the physcial laws
 * \param[out]    mu            dynamic viscosity associated to the gas species
 * \param[out]    lambda        conductivity coefficient of the gas species
 */
/*----------------------------------------------------------------------------*/

static void
_compute_mu_lambda_suth(const cs_real_t                  tk,
                        const cs_gas_mix_species_prop_t  spro,
                        cs_real_t                        *mu,
                        cs_real_t                        *lambda)
{
  const cs_real_t smu = spro.smu;
  const cs_real_t slam = spro.slam;
  const cs_real_t muref = spro.muref;
  const cs_real_t lamref = spro.lamref;
  const cs_real_t trefmu = spro.trefmu;
  const cs_real_t treflam = spro.treflam;

  *mu =  muref * pow(tk / trefmu, 1.5) * ((trefmu+smu) / (tk+smu));
  *lambda = lamref * pow(tk / treflam, 1.5) * ((treflam+slam) / (tk+slam));
}

/*----------------------------------------------------------------------------*/
/* Log values of the structure                                                */
/*----------------------------------------------------------------------------*/

static void
_log_func_gas_mix_species_prop(const void *t)
{
  const char fmt[] = N_("      %-19s  %-12.3g\n");
  const cs_gas_mix_species_prop_t *_t =
    static_cast<const cs_gas_mix_species_prop_t *>(t);
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

/*----------------------------------------------------------------------------*/
/* Log default values of the structure                                        */
/*----------------------------------------------------------------------------*/

static void
_log_func_default_gas_mix_species_prop(const void *t)
{
  const char fmt[] = "      %-19s  %-12.3g %s\n";
  const cs_gas_mix_species_prop_t *_t =
    static_cast<const cs_gas_mix_species_prop_t *>(t);
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
                               nullptr,
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

    /* Add temperature in kelvin */
    f = cs_field_create("tempk",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        false);

    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);

    cs_field_pointer_map(CS_ENUMF_(t_kelvin), f);
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
 * \brief Fills physical properties which are variable in time
 *        for the gas mixtures modelling with or without steam
 *        inside the fluid domain. In presence of steam, this one
 *        is deduced from the noncondensable gases transported
 *        as scalars (by means of the mass fraction of each species).
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_physical_properties(void)
{
  if (   (cs_glob_fluid_properties->icp < 0)
      && (cs_glob_fluid_properties->irovar < 1)
      && (cs_glob_fluid_properties->ivivar < 1)  ) {

    if (cs_glob_fluid_properties->icp < 0)
      bft_error(__FILE__, __LINE__, 0,
              _("Inconsistent calculation data icp = %d\n"
                "The calculation will not be run.\n"
                "Modify cs_user_parameters or cs_user_physical_properties.\n"),
                cs_glob_fluid_properties->icp);

    bft_error(__FILE__, __LINE__, 0,
              _("Inconsistent calculation data for irovar = %d or ivivar = %d\n"),
              cs_glob_fluid_properties->irovar,
              cs_glob_fluid_properties->ivivar);
  }

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_field_t *th_f = cs_thermal_model_field();
  const int kivisl = cs_field_key_id("diffusivity_id");
  int ifcvsl = cs_field_get_key_int(th_f, kivisl);

  /*  Initializations
      --------------- */

  cs_real_t *tempk = nullptr;
  cs_real_t *lambda = nullptr;
  cs_real_t *cpro_rho = nullptr;
  const cs_real_t *cvar_enth = nullptr;

  /* Specific heat value */
  cs_real_t *cpro_cp = CS_F_(cp)->val;

  /* Molecular dynamic viscosity value */
  cs_real_t *cpro_viscl = CS_F_(mu)->val;

  /* Lambda/Cp value */
  cs_real_t *cpro_venth = cs_field_by_id(ifcvsl)->val;

  cs_real_t *steam_binary_diffusion
    = cs_field_by_name("steam_binary_diffusion")->val;
  cs_real_t *mix_mol_mas = cs_field_by_name("mix_mol_mas")->val;
  cs_real_t *mol_mas_ncond = cs_field_by_name("mol_mas_ncond")->val;

  /* Deduce mass fraction (y_d) which is
   * y_h2o_g in presence of steam or
   * y_he/y_h2 with noncondensable gases */
  const cs_field_t *f = nullptr;
  cs_gas_mix_species_prop_t s_d;
  const int k_id = cs_gas_mix_get_field_key();

  if (cs_glob_physical_model_flag[CS_GAS_MIX] == CS_GAS_MIX_AIR_HELIUM)
    f = cs_field_by_name("y_he");
  else if (cs_glob_physical_model_flag[CS_GAS_MIX] == CS_GAS_MIX_AIR_HYDROGEN)
    f = cs_field_by_name("y_h2");
  else if (   (cs_glob_physical_model_flag[CS_GAS_MIX] >= CS_GAS_MIX_AIR_STEAM)
           && (cs_glob_physical_model_flag[CS_GAS_MIX] <  CS_GAS_MIX_HELIUM_AIR))
    f = cs_field_by_name("y_h2o_g");
  else
    f = cs_field_by_name("y_o2");

  cs_real_t *y_d = f->val;
  cs_field_get_key_struct(f, k_id, &s_d);
  /* Storage the previous value of the deduced mass fraction ya_d */
  cs_array_real_copy(n_cells_ext, y_d, f->val_pre);

  /* In compressible, the density is updated after the pressure step */
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
    cvar_enth = th_f->val;
    cpro_rho  = CS_F_(rho)->val;
    tempk     = cs_field_by_name("tempk")->val;
    BFT_MALLOC(lambda, n_cells_ext, cs_real_t);
  }
  else {
    tempk = CS_F_(t_kelvin)->val;
    ifcvsl = cs_field_get_key_int(CS_F_(t_kelvin), kivisl);
    lambda = cs_field_by_id(ifcvsl)->val;
  }

  /* Define the physical properties for the gas mixture with:
   *  - the density (rho_m) and specific heat (cp_m) of the gas mixture
   *    function temperature and species scalar (yk),
   *  - the dynamic viscosity (mu_m) and conductivity coefficient (lbd_m) of
   *    the gas mixture function ot the enthalpy and species scalars,
   *  - the diffusivity coefficients of the scalars (Dk, D_enh).
   ---------------------------------------------------------------------------*/

  cs_real_t pressure = cs_glob_fluid_properties->p0;
  if (cs_glob_velocity_pressure_model->idilat == 3)
    pressure = cs_glob_fluid_properties->pther;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    const cs_real_t x_0 = 0.0;
    const cs_real_t x_1 = 1.0;

    // Initialization
    y_d[c_id] = x_1;

    // Thermal conductivity
    lambda[c_id] = x_0;

    // Mixture specific heat
    cpro_cp[c_id] = x_0;

    // Mixture molecular diffusivity
    cpro_viscl[c_id] = x_0;

    // Mixture molar mass
    mix_mol_mas[c_id] = x_0;

    mol_mas_ncond[c_id] = x_0;
    steam_binary_diffusion[c_id] = x_0;

    /* Mass fraction array of the different species */
    for (int spe_id = 0; spe_id < _gas_mix.n_species_solved; spe_id++) {
      const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
      const cs_field_t *f_spe = cs_field_by_id(f_spe_id);
      const cs_real_t *cvar_yk = f_spe->val;

      cs_gas_mix_species_prop_t s_k;
      cs_field_get_key_struct(f_spe, k_id, &s_k);

      y_d[c_id] -= cvar_yk[c_id];
      mix_mol_mas[c_id] += cvar_yk[c_id]/s_k.mol_mas;
      mol_mas_ncond[c_id] += cvar_yk[c_id] / s_k.mol_mas;
    }

    // Clipping
    y_d[c_id] = cs_math_fmax(y_d[c_id], x_0);

    // Finalize the computation of the Mixture molar mass
    mix_mol_mas[c_id] += y_d[c_id]/s_d.mol_mas;
    mix_mol_mas[c_id] = x_1/mix_mol_mas[c_id];
    mol_mas_ncond[c_id] = (x_1 - y_d[c_id])/mol_mas_ncond[c_id];

    for (int spe_id = 0; spe_id < _gas_mix.n_species_solved + 1; spe_id++) {
      /* Mixture specific heat function of species specific heat (cpk)
       * and mass fraction of each gas species (yk), as below:
       *             -----------------------------
       * - noncondensable gases and the mass fraction deduced:
       *             cp_m[c_id] = Sum( yk.cpk)_k[0, n_species_solved]
       *                        + y_d.cp_d
       *             -----------------------------
       * remark:
       *   The mass fraction is deduced depending of the
       *    modelling chosen by the user, with:
       *      - CS_GAS_MIX = CS_GAS_MIX_AIR_HELIUM or
       *        CS_GAS_MIX_AIR_HYDROGEN, a noncondensable gas
       *      - CS_GAS_MIX > CS_GAS_MIX_AIR_STEAM, a condensable gas (steam) */
      const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
      const cs_field_t *f_spe = cs_field_by_id(f_spe_id);
      if (spe_id == _gas_mix.n_species_solved)
        f_spe = f;
      const cs_real_t *cvar_yk = f_spe->val;

      cs_gas_mix_species_prop_t s_k;
      cs_field_get_key_struct(f_spe, k_id, &s_k);
      cpro_cp[c_id] += cvar_yk[c_id]*s_k.cp;
    }

  }

  /* gas mixture density function of the temperature, pressure
   * and the species scalars with taking into account the
   * dilatable effects, as below:
   * ----------------------------
   * - with inlet/outlet conditions:
   *   [idilat=2]: rho= p0/(R. temp(1/sum[y_i/M_i]))
   * - with only wall conditions:
   *   [idilat=3]: rho= pther/(R. temp(1/sum[y_i/M_i]))
   *                          ----------------------
   *         i ={1, .. ,N} : species scalar number
   *         y_i           : mass fraction of each species
   *         M_i           : molar fraction [kg/mole]
   *         R             : ideal gas constant [J/mole/K]
   *         p0            : atmos. pressure (Pa)
   *         pther         : pressure (Pa) integrated on the
   *                         fluid domain */
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {
      // Evaluate the temperature thanks to the enthalpy
      tempk[c_id] = cvar_enth[c_id]/ cpro_cp[c_id];
      cpro_rho[c_id]
        = pressure*mix_mol_mas[c_id]/(cs_physical_constants_r*tempk[c_id]);
    }

  /* Dynamic viscosity and conductivity coefficient
   * the physical properties associated to the gas
   * mixture with or without condensable gas */

  /* Loop on all species */
  for (int spe_id = 0; spe_id < _gas_mix.n_species_solved + 1; spe_id++) {

    const cs_field_t *f_spe = f;
    if (spe_id < _gas_mix.n_species_solved) {
      const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
      f_spe = cs_field_by_id(f_spe_id);
    }

    cs_gas_mix_species_prop_t s_i;
    const cs_real_t *cvar_yi = f_spe->val;
    cs_field_get_key_struct(f_spe, k_id, &s_i);

    cs_real_t  mu_i, lambda_i,  mu_j, lambda_j;

    const int ivsuth = cs_glob_fluid_properties->ivsuth;
    if (ivsuth == 1) {
      if (   (strcmp(f_spe->name, "y_he")    != 0)
          && (strcmp(f_spe->name, "y_h2")    != 0)
          && (strcmp(f_spe->name, "y_o2")    != 0)
          && (strcmp(f_spe->name, "y_n2")    != 0)
          && (strcmp(f_spe->name, "y_h2o_g") != 0)   )
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: no predefined properties for field %s."),
                  __func__, f_spe->name);
    }

    /* TODO: use local arrays to store s_j and field name or identifiers,
     *       so as to avoid querying them for each field for each cell. */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

      if (ivsuth == 0)
        _compute_mu_lambda(f_spe->name, tempk[c_id], s_i,
                           &mu_i, &lambda_i);
      else
        _compute_mu_lambda_suth(tempk[c_id], s_i,
                                &mu_i, &lambda_i);

      cs_real_t xsum_mu = 0.0, xsum_lambda = 0.0;

      for (int s_id = 0; s_id < _gas_mix.n_species_solved + 1; s_id++) {

        const cs_field_t *f_s = f;
        if (s_id < _gas_mix.n_species_solved) {
          const int f_s_id = _gas_mix.species_to_field_id[s_id];
          f_s = cs_field_by_id(f_s_id);
        }

        cs_gas_mix_species_prop_t s_j;
        const cs_real_t *cvar_yj = f_s->val;
        cs_field_get_key_struct(f_s, k_id, &s_j);

        if (ivsuth == 0)
          _compute_mu_lambda(f_s->name, tempk[c_id], s_j,
                             &mu_j, &lambda_j);
        else
          _compute_mu_lambda_suth(tempk[c_id], s_j,
                                  &mu_j, &lambda_j);

        const cs_real_t phi_mu
          =   (1.0/sqrt(8.0))
            * pow(1.0 + s_i.mol_mas/s_j.mol_mas, -0.5)
            * pow(1.0 + pow(mu_i/mu_j, 0.5)
            * pow(s_j.mol_mas/s_i.mol_mas, 0.25), 2);

        const cs_real_t phi_lambda
          =   (1.0/sqrt(8.0))
            * pow(1.0 + s_i.mol_mas/s_j.mol_mas, -0.5)
            * pow(1.0 + pow(lambda_i/lambda_j, 0.5)
            * pow(s_j.mol_mas / s_i.mol_mas, 0.25), 2);

        const cs_real_t x_k = cvar_yj[c_id]*mix_mol_mas[c_id]/s_j.mol_mas;
        xsum_mu += x_k * phi_mu;
        xsum_lambda += x_k * phi_lambda;
      }

      /* Mixture viscosity defined as function of the scalars
         ----------------------------------------------------- */
      const cs_real_t x_k
        = cvar_yi[c_id]*mix_mol_mas[c_id]/s_i.mol_mas;
      cpro_viscl[c_id] = cpro_viscl[c_id] + x_k * mu_i / xsum_mu;

      lambda[c_id] += x_k * lambda_i / xsum_lambda;

    } //loop on cells
  } // end of loop on species

  /* Dynamic viscosity and conductivity coefficient
   * the physical properties filled for the gas mixture */

  /* Same diffusivity for all the scalars except the enthalpy */
  for (int spe_id = 0; spe_id < _gas_mix.n_species_solved; spe_id++) {

    const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
    const cs_field_t *f_spe = cs_field_by_id(f_spe_id);
    ifcvsl = cs_field_get_key_int(f_spe, kivisl);
    cs_real_t *cpro_vyk = cs_field_by_id(ifcvsl)->val;

    cs_array_real_copy(n_cells, cpro_viscl, cpro_vyk);
  }

  const cs_real_t patm = 101320.0;

  /* Steam binary diffusion */
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    cs_real_t x_ncond_tot = 0.0;
    const cs_real_t ratio_tkpr = pow(tempk[c_id], 1.75)/pressure;

    for (int spe_id = 0; spe_id < _gas_mix.n_species_solved; spe_id++) {
      const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
      const cs_field_t *f_spe = cs_field_by_id(f_spe_id);

      const cs_real_t *cvar_yi = f_spe->val;
      cs_gas_mix_species_prop_t s_i;
      cs_field_get_key_struct(f_spe, k_id, &s_i);

      const cs_real_t y_k = cvar_yi[c_id];
      const cs_real_t x_k = y_k * mix_mol_mas[c_id] / s_i.mol_mas;
      const cs_real_t xmab
        = sqrt(2.0/( 1.0 / (s_d.mol_mas*1000.0) +1.0 / (s_i.mol_mas*1000.0)));
      const cs_real_t xvab
        = pow(pow(s_d.vol_dif, 1.0/3.0) + pow(s_i.vol_dif, 1.0/3.0), 2.0);
      const cs_real_t a1 = 1.43e-7 / (xmab * xvab) * patm;
      steam_binary_diffusion[c_id] += x_k / (a1 * ratio_tkpr);
      x_ncond_tot += x_k;
    }

    steam_binary_diffusion[c_id] = x_ncond_tot / steam_binary_diffusion[c_id];

  }

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
      cpro_venth[c_id] = lambda[c_id]/cpro_cp[c_id];

    BFT_FREE(lambda);
  }
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
/*!
 * \brief Initialization of calculation variables for gas mixture modelling
 *        in presence of the steam gas or another gas used as variable deduced
 *        and not solved.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_initialization(void)
{
  if (cs_glob_fluid_properties->icp < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Inconsistent calculation data icp = %d\n"
                "The calculation will not be run.\n"
                "Modify cs_user_parameters or cs_user_physical_properties.\n"),
              cs_glob_fluid_properties->icp);

  /* If this is restarted computation, do not reinitialize values */
  if (cs_glob_time_step->nt_prev > 0)
    return;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* Initializations
     ---------------- */
  cs_real_t *cvar_enth = nullptr;

  /* Specific heat value */
  cs_real_t *cpro_cp = CS_F_(cp)->val;

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
    cvar_enth = cs_thermal_model_field()->val;

  /* Deduced species (h2o_g) with steam gas
   * or Helium or Hydrogen  with noncondensable gases */
  const cs_field_t *f = nullptr;
  cs_gas_mix_species_prop_t s_d;
  const int k_id = cs_gas_mix_get_field_key();

  if (cs_glob_physical_model_flag[CS_GAS_MIX] == CS_GAS_MIX_AIR_HELIUM)
    f = cs_field_by_name("y_he");
  else if (cs_glob_physical_model_flag[CS_GAS_MIX] == CS_GAS_MIX_AIR_HYDROGEN)
    f = cs_field_by_name("y_h2");
  else if (   (cs_glob_physical_model_flag[CS_GAS_MIX] >= CS_GAS_MIX_AIR_STEAM)
           && (cs_glob_physical_model_flag[CS_GAS_MIX] <  CS_GAS_MIX_HELIUM_AIR))
    f = cs_field_by_name("y_h2o_g");
  else
    f = cs_field_by_name("y_o2");

  cs_real_t *y_d = f->val;
  cs_field_get_key_struct(f, k_id, &s_d);

  cs_real_t *mix_mol_mas = cs_field_by_name("mix_mol_mas")->val;

  int iok = 0;
  cs_real_t vol_d = 0.0;
  const cs_real_t t0 = cs_glob_fluid_properties->t0;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {
    y_d[c_id] = 1.;
    cpro_cp[c_id] = 0.;
    mix_mol_mas[c_id] = 0.;
  }

  for (int spe_id = 0; spe_id < _gas_mix.n_species_solved; spe_id++) {
    const int f_spe_id = _gas_mix.species_to_field_id[spe_id];
    const cs_field_t *f_spe = cs_field_by_id(f_spe_id);
    const cs_real_t *cvar_yk = f_spe->val;

    cs_gas_mix_species_prop_t s_k;
    cs_field_get_key_struct(f_spe, k_id, &s_k);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {
      if ((cvar_yk[c_id] > 1.0) || (cvar_yk[c_id] < 0.0))
        iok++;
      y_d[c_id] -= cvar_yk[c_id];
      cpro_cp[c_id] += cvar_yk[c_id]*s_k.cp;
      mix_mol_mas[c_id] += cvar_yk[c_id]/s_k.mol_mas;
    }
  }

  /* Finalization and check */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    if ((y_d[c_id] > 1.0) || (y_d[c_id] < 0.0))
      iok++;

    y_d[c_id] = cs_math_fmin(cs_math_fmax(y_d[c_id], 0.0), 1.0);

    // specific heat (Cp_m0) of the gas mixture
    cpro_cp[c_id] += y_d[c_id]*s_d.cp;

    if (cvar_enth != nullptr)
      cvar_enth[c_id] = cpro_cp[c_id]*t0;

    mix_mol_mas[c_id] += y_d[c_id]/s_d.mol_mas;
    mix_mol_mas[c_id]  = 1.0/mix_mol_mas[c_id];

    /* Gas deduced and Total gas volumes injected */
    vol_d += cell_vol[c_id]*(y_d[c_id]/s_d.mol_mas)*mix_mol_mas[c_id];

  }

  cs_parall_sum(1, CS_REAL_TYPE, &vol_d);
  const cs_real_t volgas = cs_glob_mesh_quantities->tot_vol;

  /* Print to the log to check the variables intialization
     ----------------------------------------------------- */

  bft_printf("----------------------------------------------------------\n"
             "**     Gas mixture : Check variables initialization     **\n"
             "----------------------------------------------------------\n"
             "   Total   gas Volume: %10.17le\n"
             "   Deduced gas Volume: %10.17le\n", volgas , vol_d);

  if (iok > 0)
     bft_error(__FILE__, __LINE__, 0,
              _("Abort in the variables initialization.\n"
                " The variables initialization is incomplete or\n"
                " incoherent with the parameters value of the calculation.\n"
                "The calculation will not be run (%d,' errors). \n"
                "Refer to the previous warnings for further information."
                "Pay attention to the initialization of, \n "
                "                              the time-step\n"
                "                              the turbulence\n"
                "                              the scalars and variances\n"
                "                              the time averages\n\n"
                "Verify the initialization and the restart file."
                "In the case where the values read in the restart file"
                "are incorrect, they may be modified with"
                "cs_user_initialization.c or with the interface."), iok);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
