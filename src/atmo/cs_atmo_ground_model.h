#ifndef CS_ATMO_GROUND_MODEL_H
#define CS_ATMO_GROUND_MODEL_H

/*============================================================================
 * Atmospheric ground module - Build constants and variables to describe ground model
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \file cs_atmo_ground_model.h */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Plant model options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {
  /*! Height of the canopy */
  cs_real_t h_canopy;
  /*! Plant air exchange height */
  cs_real_t h_canopy_exch;
  /*! Leaf length */
  cs_real_t d_leaf;
  /*! Leaf Area Index */
  cs_real_t leaf_area_index;
  /*! Canopy extinction coefficient for the radiation */
  cs_real_t k_ext_coef;
  /*! Albedo of the plant */
  cs_real_t albedo_plant;
  /*! Emissivity of the plant */
  cs_real_t emi_plant;
  /*! CO2 concentration in the air */
  cs_real_t c_co2_air;
  /*! Constant for the CO2 flux in the ground */
  cs_real_t f_co2_0_ground;
  /*! Reference temperature for the CO2 exchanges */
  cs_real_t temp_ref_co2;
  /*! Activation energy for the CO2 */
  cs_real_t ea_ground;
  /*! Stomatal conductance constant */
  cs_real_t a_pho;
  /*! Stomatal conductance minimum value */
  cs_real_t g0_pho;
  /*! Constant factor for the water stress model */
  cs_real_t sf_psiv;
  /*! Water stress reference value */
  cs_real_t psif_psiv;
  /*! Plant process constant 1 */
  cs_real_t oi;
  /*! Plant process constant 2 */
  cs_real_t koref;
  /*! Plant process constant 3 */
  cs_real_t kcref;
  /*! Activation energy linked to constant 2 */
  cs_real_t eao;
  /*! Activation energy linked to constant 3 */
  cs_real_t eac;
  /*! Gamma const 0 */
  cs_real_t gamma_0;
  /*! Gamma const 1 */
  cs_real_t gamma_1;
  /*! Gamma const 2 */
  cs_real_t gamma_2;
  /*! Reference temperature for the photsynthesis */
  cs_real_t tem_ref_pho;
  /*! Photosynhtesis process 1 constant */
  cs_real_t vlmaxref;
  /*! Leaf allocation coef for the nitrogen */
  cs_real_t knitro;
  /*! Max photon flux */
  cs_real_t jlmaxref;
  /*! Activation energy for photosynthesis process 1 */
  cs_real_t ea_vlmax;
  /*! Desactivation energy for photsynthesis process 1 */
  cs_real_t ed_vlmax;
  /*! Entropy constant for photsynthesis process 1 */
  cs_real_t s_vlmax;
  /*! Activation energy for photosynthesis process 2 */
  cs_real_t ea_jlmax;
  /*! Desactivation energy for photsynthesis process 2 */
  cs_real_t ed_jlmax;
  /*! Entropy constant for photsynthesis process 2 */
  cs_real_t s_jlmax;
  /*! Photon flux constant 1 */
  cs_real_t alpha_pho;
  /*! Photon flux constant 2 */
  cs_real_t theta_pho;
  /*! Speed of increase of gco2 */
  cs_real_t v_gco2_incr;
  /*! Speed of decrease of gco2 */
  cs_real_t v_gco2_decr;
  /*! Hydraulic resistance */
  cs_real_t xiv;
  /*! Root radius */
  cs_real_t r_root;
  /*! Root depth */
  cs_real_t z_root;
  /*! Length of roots in a volume of ground */
  cs_real_t ld_root;
  /*! Choice of water potential model */
  int water_potential_model;
  /*! Ref water potential */
  cs_real_t psi_e;
  /*! Air hydraulic power */
  cs_real_t psi_pow;
  /*! Reference ground hydraulic conductivity */
  cs_real_t kref;
  /*! Canopy mixing length */
  cs_real_t canopy_mix_l;
  /*! Turbulent Prandtl in the canopy */
  cs_real_t turb_prandtl;
  /*! Leaf draf coefficient */
  cs_real_t cdrag_leaf;
  /*! Bare ground rugosity */
  cs_real_t zos;
  /*! Reference height of the plant */
  cs_real_t zref_plant;
  /*! Eddy diffusivity at the top of the canopy */
  cs_real_t k_eddy_hc;
  /*! Ground tortuosity factor */
  cs_real_t tort_factor;
  /*! Ground exchange length */
  cs_real_t ground_exch_l;
  /*! Dry ground porosity */
  cs_real_t dry_ground_porosity;
  /*! Ground constant for the porosity */
  cs_real_t dv;
} cs_plant_option_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to atmo options structure */
extern cs_plant_option_t *cs_glob_plant_option;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build constants and variables to describe ground model
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_ground_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ground - atmosphere parameters computed from a "Land use" file
 *
 * \param[in] call_stage  first pass to set default values,
 *                        second pass to perform some checks and log
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_ground_cat(int call_stage);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute ground and interface values using Deardorff force restore method
 */
/*----------------------------------------------------------------------------*/

void
cs_ground_model(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_ATMO_GROUND_MODEL_H */
