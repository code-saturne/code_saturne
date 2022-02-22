#ifndef __CS_GAS_MIX_H__
#define __CS_GAS_MIX_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Gas mix types
   ------------- */

typedef enum {

  CS_GAS_MIX_OFF = -1,                   /*!< gas mix model off */
  CS_GAS_MIX_AIR_HELIUM = 0,             /*!< air/helium, helium deduced */
  CS_GAS_MIX_AIR_HYDROGEN = 1,           /*!< air/hydrogen, hydrogen deduced */
  CS_GAS_MIX_AIR_STEAM = 2,              /*!< air/steam, steam deduced */
  CS_GAS_MIX_AIR_HELIUM_STEAM = 3,       /*!< helium/steam, steam deduced */
  CS_GAS_MIX_AIR_HYDROGEN_STEAM = 4,     /*!< hydrogen/steam, steam deduced */
  CS_GAS_MIX_HELIUM_AIR = 5,             /*!< helium/air, O2 from air deduced */
  CS_GAS_MIX_USER                        /*!< user defined */

} cs_gas_mix_type_t;

/* Gas mix descriptor
   ------------------ */

typedef struct {

  int    n_species;             /*!< number of species in the gas mix */
  int    n_species_solved;      /*!< number of species which
                                  are solved variables */
  int   *species_to_field_id;   /*!< species to field mapping
                                  (solved variables first) */

} cs_gas_mix_t;

/*
 * Gas mix modelling physical properties
 * ------------------------------------- */

typedef struct {

  double  mol_mas;   /* molar mass */
  double  cp;        /* specific heat at constant pressure */
  double  vol_dif;   /* volume diffusion */
  double  mu_a;      /* dynamic viscosity a */
  double  mu_b;      /* dynamic viscosity a */
  double  lambda_a;  /* thermal conductivity a */
  double  lambda_b;  /* thermal conductivity b */
  double  muref;     /* ref. viscosity for Sutherland law */
  double  lamref;    /* ref. thermal conductivity for Sutherland law */
  double  trefmu;    /* ref. temperature for viscosity in Sutherland law */
  double  treflam;   /* ref. temperature for conductivity Sutherland law */
  double  smu;       /* Sutherland temperature for viscosity */
  double  slam;      /* Sutherland temperature for conductivity */

} cs_gas_mix_species_prop_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main physical constants structure */

extern const cs_gas_mix_t  *cs_glob_gas_mix;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the field key for gas mix properties.
 *
 * \return  field key id for gas mix properties
 */
/*----------------------------------------------------------------------------*/

int
cs_gas_mix_get_field_key(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species(int f_id);

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
                                       cs_real_t  s_lambda);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add property fields specific to a gas mix.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_property_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free array mapping gas mix species ids to field ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GAS_MIX_H__ */
