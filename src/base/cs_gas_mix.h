#ifndef CS_GAS_MIX_H
#define CS_GAS_MIX_H

/*============================================================================
 * Base gas mix data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

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
  CS_GAS_MIX_CO2_AIR = 6,                /*!< CO2/air, O2 from air deduced */
  CS_GAS_MIX_NO2_AIR = 7,                /*!< NO2/air, O2 from air deduced */
  CS_GAS_MIX_USER                        /*!< user defined */

} cs_gas_mix_type_t;

/* Known gases
   ----------- */
/*! Enum for gas types used in mixtures */
enum class cs_gas_mix_y_type {
  h2o,     /*!< Water vapor (H2O) */
  he,      /*!< Helium (He) */
  h2,      /*!< Hydrogen (H2) */
  co2,     /*!< CO2 */
  no2,     /*!< NO2 */
  o2,      /*!< Oxygen (O2)*/
  n2,      /*!< Nitrogen (N2) */
  n_gases, /*!< Number of predefined gases */
  unknown  /*!< Uknown gas -> Used for error detection */
};

/* Gas mix descriptor
   ------------------ */

/*! Structure containing the data related to a gas mixture */

typedef struct cs_gas_mix_t{

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor
   */
  /*--------------------------------------------------------------------------*/

  cs_gas_mix_t() = default;

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default destructor
   */
  /*--------------------------------------------------------------------------*/

  ~cs_gas_mix_t() = default;

  /* ------- */
  /* Members */
  /* ------- */

  int n_species {0};        /*!< number of species in the gas mix */
  int n_species_solved {0}; /*!< number of species which
                                 are solved variables */

  cs_array<int> species_to_field_id; /*!< species to field mapping
                                          (solved variables first) */
  cs_array<cs_gas_mix_y_type> gas_type; /*!< Type of predefined gas */

  cs_array_2d<cs_real_t> acp;     /*!< Polynomial coefficients used for Cp calculation */

  cs_array<cs_real_t>  mol_mas;   /*!< molar mass */
  cs_array<cs_real_t>  cp;        /*!< specific heat at constant pressure */
  cs_array<cs_real_t>  vol_dif;   /*!< volume diffusion */
  cs_array<cs_real_t>  mu_a;      /*!< dynamic viscosity a */
  cs_array<cs_real_t>  mu_b;      /*!< dynamic viscosity a */
  cs_array<cs_real_t>  lambda_a;  /*!< thermal conductivity a */
  cs_array<cs_real_t>  lambda_b;  /*!< thermal conductivity b */
  cs_array<cs_real_t>  muref;     /*!< ref. viscosity for Sutherland law */
  cs_array<cs_real_t>  lamref;    /*!< ref. thermal conductivity for Sutherland law */
  cs_array<cs_real_t>  trefmu;    /*!< ref. temperature for viscosity in Sutherland law */
  cs_array<cs_real_t>  treflam;   /*!< ref. temperature for conductivity Sutherland law */
  cs_array<cs_real_t>  smu;       /*!< Sutherland temperature for viscosity */
  cs_array<cs_real_t>  slam;      /*!< Sutherland temperature for conductivity */

  /* ------- */
  /* Methods */
  /* ------- */

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Check if a polynomial formula for Cp is used.
   *
   * \return True if used, false otherwise
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  bool use_polynomial_cp
  () const
  {
    return (_cp_poly_d > 1);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get polynomial degree. Currently is either 1 or 6
   *
   * \return value of Cp polynomial degree
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  int
  cp_poly_d
  () const
  {
    return _cp_poly_d;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Activate the Cp polynomial Cp formula
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  polynomial_cp_activate()
  {
    // Currently only degree 6 is handled, could be changed in the future
    _cp_poly_d = 6;
  }

private:
  int _cp_poly_d {1}; /*!< Polynomial degree used for Cp calculation */

} cs_gas_mix_t;

/*
 * Gas mix modelling physical properties
 * ------------------------------------- */

typedef struct {

  cs_real_t  mol_mas;   /* molar mass */
  cs_real_t  cp;        /* specific heat at constant pressure */
  cs_real_t  vol_dif;   /* volume diffusion */
  cs_real_t  mu_a;      /* dynamic viscosity a */
  cs_real_t  mu_b;      /* dynamic viscosity a */
  cs_real_t  lambda_a;  /* thermal conductivity a */
  cs_real_t  lambda_b;  /* thermal conductivity b */
  cs_real_t  muref;     /* ref. viscosity for Sutherland law */
  cs_real_t  lamref;    /* ref. thermal conductivity for Sutherland law */
  cs_real_t  trefmu;    /* ref. temperature for viscosity in Sutherland law */
  cs_real_t  treflam;   /* ref. temperature for conductivity Sutherland law */
  cs_real_t  smu;       /* Sutherland temperature for viscosity */
  cs_real_t  slam;      /* Sutherland temperature for conductivity */

} cs_gas_mix_species_prop_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main physical constants structure */

extern const cs_gas_mix_t  *cs_glob_gas_mix;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*
 * \brief Finalize setup by creating all data structure and doing final checks
 */
/*--------------------------------------------------------------------------*/

void
cs_gas_mix_setup_finalize(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get the field key for gas mix properties.
 *
 * \return  field key id for gas mix properties
 */
/*----------------------------------------------------------------------------*/

int
cs_gas_mix_get_field_key(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species(int f_id);

/*----------------------------------------------------------------------------*/
/*
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
/*
 * \brief Add variable fields specific to a gas mix.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_variable_fields(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Add property fields specific to a gas mix.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_property_fields(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Free array mapping gas mix species ids to field ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_finalize(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Initialization of calculation variables for gas mixture modelling
 *        in presence of the steam gas or another gas used as variable deduced
 *        and not solved.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_initialization(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Fills physical properties which are variable in time
 *        for the gas mixtures modelling with or without steam
 *        inside the fluid domain. In presence of steam, this one
 *        is deduced from the noncondensable gases transported
 *        as scalars (by means of the mass fraction of each species).
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_physical_properties(void);

/*--------------------------------------------------------------------------*/
/*
 * \brief Update the deduced species fraction based on the solved species
 *        fractions, since the sum must be equal to 1, and fractions between
 *        0 and 1.
 */
/*--------------------------------------------------------------------------*/

void
cs_gas_mix_update_deduced_fraction(void);

/*--------------------------------------------------------------------------*/
/*
 * \brief Activate polynomial formula for Cp
 */
/*--------------------------------------------------------------------------*/

void
cs_gas_mix_use_cp_polynomial_formula(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_GAS_MIX_H */
