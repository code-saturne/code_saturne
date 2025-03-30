/*============================================================================
 * Base wall condensation model.
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_1d_wall_thermal.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_gas_mix.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_restart.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_wall_condensation_1d_thermal.h"
#include "base/cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_wall_condensation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

//> Constants for the correlation of steam saturated pressure
static const cs_real_t pr_c = 221.2e+5;
static const cs_real_t t_c  = 647.3e0;
static const cs_real_t c_k1 = -7.691234564e0;
static const cs_real_t c_k2 = -26.08023696e0;
static const cs_real_t c_k3 = -168.1706546e0;
static const cs_real_t c_k4 = 64.23285504e0;
static const cs_real_t c_k5 = -118.9646225e0;
static const cs_real_t c_k6 = 4.16711732e0;
static const cs_real_t c_k7 = 20.9750676e0;
static const cs_real_t c_k8 = -1.e+9;
static const cs_real_t c_k9 = 6.e0;

//> Caracteristic length for natural convection models
static const cs_real_t lcar_nc = 1.0;
//> Turbulent Prandtl
static const cs_real_t pr_tur = 0.9;
//> Latent heat of condensation (water)
static const cs_real_t lcond = 2278.0e+3;

static cs_wall_condensation_t _wall_cond
  = {.icondb             = -1,
     .icondv             = -1,
     .nztag1d            = -1,
     .natural_conv_model = CS_WALL_COND_MODEL_COPAIN,
     .forced_conv_model  = CS_WALL_COND_MODEL_WALL_LAW, // fixed for now
     .mixed_conv_model   = CS_WALL_COND_MIXED_MAX,      // fixed for now

     /* Surface wall condensation */

     // Mesh related quantities
     // TODO: clean unnecessary quantities
     .nfbpcd                    = 0,
     .ifbpcd                    = nullptr,
     .itypcd                    = nullptr,
     .izzftcd                   = nullptr,
     .spcond                    = nullptr,
     .hpcond                    = nullptr,
     .twall_cond                = nullptr,
     .thermal_condensation_flux = nullptr,
     .convective_htc            = nullptr,
     .condensation_htc          = nullptr,
     .total_htc                 = nullptr,
     .flthr                     = nullptr,
     .dflthr                    = nullptr,

     // Zone related quantities
     .nzones    = -1,
     .izcophc   = nullptr,
     .izcophg   = nullptr,
     .iztag1d   = nullptr,
     .ztpar     = nullptr,
     .zxrefcond = nullptr,
     .zprojcond = nullptr,

     /* Volumetric wall condensation */

     .ncmast                    = 0,
     .nvolumes                  = -1,
     .ltmast                    = nullptr,
     .itypst                    = nullptr,
     .izmast                    = nullptr,
     .svcond                    = nullptr,
     .flxmst                    = nullptr,
     .itagms                    = nullptr};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

const cs_wall_condensation_t *cs_glob_wall_condensation = &_wall_cond;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute saturation pressure at given temperature (in Kelvins)
 *----------------------------------------------------------------------------*/

static cs_real_t
_compute_psat(cs_real_t temperature)
{
  cs_real_t dtheta      = temperature / t_c;
  cs_real_t dtheta_comp = 1.0 - dtheta;
  cs_real_t psat
    = pr_c
      * exp((1.0 / dtheta)
              * (c_k1 * dtheta_comp + c_k2 * cs_math_pow2(dtheta_comp)
                 + c_k3 * cs_math_pow3(dtheta_comp)
                 + c_k4 * cs_math_pow4(dtheta_comp)
                 + c_k5 * pow(dtheta_comp, 5.0))
              / (1.0 + c_k6 * dtheta_comp + c_k7 * cs_math_pow2(dtheta_comp))
            - dtheta_comp / (c_k8 * cs_math_pow2(dtheta_comp) + c_k9));
  return psat;
}

/*---------------------------------------------------------------------------
 * Compute mole fraction from mass fraction
 *----------------------------------------------------------------------------*/

static cs_real_t
_compute_mole_fraction(cs_real_t  mass_fraction,
                       cs_real_t  mix_mol_mas,
                       cs_real_t  mol_mas)
{
  return mass_fraction * mix_mol_mas / mol_mas;
}

/*---------------------------------------------------------------------------
 * Get wall temperature at a given face
 *----------------------------------------------------------------------------*/

static cs_real_t
_get_wall_temperature(cs_lnum_t  ieltcd)
{
  cs_real_t temperature;
  int       izone = _wall_cond.izzftcd[ieltcd];
  if (_wall_cond.iztag1d[izone] == 1) {
    if ((cs_glob_time_step->nt_cur == 1) && !cs_restart_present()) {
      temperature = cs_glob_wall_cond_1d_thermal->ztpar0[izone];
    }
    else {
      temperature = cs_glob_wall_cond_1d_thermal->ztmur[ieltcd];
    }
    temperature += cs_physical_constants_celsius_to_kelvin;
  }
  else if (_wall_cond.iztag1d[izone] == 2) {
    cs_lnum_t ifac = _wall_cond.ifbpcd[ieltcd];
    temperature    = cs_field_by_name("boundary_temperature")->val[ifac];
  }
  else {
    temperature = _wall_cond.ztpar[izone];
    temperature += cs_physical_constants_celsius_to_kelvin;
  }

  _wall_cond.twall_cond[ieltcd] = temperature;

  return temperature;
}

/*----------------------------------------------------------------------------
 * Compute Mac Adams natural convection correlation for mass or
 * heat transfer exchange coefficient
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_compute_mac_adams(cs_real_t theta,
                   cs_real_t grashof,
                   cs_real_t schmidt_or_prandtl)
{
  constexpr cs_real_t c_1ov3 = 1./3.;

  return theta * 0.13 * pow(grashof * schmidt_or_prandtl, c_1ov3);
}

/*----------------------------------------------------------------------------
 * Compute Schlichting number
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_compute_schlichting(cs_real_t theta,
                     cs_real_t reynolds,
                     cs_real_t schmidt_or_prandtl)
{
  constexpr cs_real_t c_1ov3 = 1./3.;

  return   theta * 0.0296 * pow(reynolds, 0.8)
         * pow(schmidt_or_prandtl, c_1ov3);
}

/*----------------------------------------------------------------------------
 * Compute Grashof number
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_compute_grashof(cs_real_t  gravity,
                 cs_real_t  drho,
                 cs_real_t  length,
                 cs_real_t  kin_viscosity)
{
  return   gravity * fabs(drho) * cs_math_pow3(length)
         / cs_math_pow2(kin_viscosity);
}

/*----------------------------------------------------------------------------
 * Compute characteristic length for schlichting model
 *----------------------------------------------------------------------------*/

static cs_real_t
_compute_characteristic_length(const cs_real_t  point[3],
                               const cs_real_t  ref_point[3],
                               const cs_real_t  ref_direction[3])
{
  cs_real_t lcar = 0.0;
  for (int idir = 0; idir < 3; idir++) {
    lcar += (point[idir] - ref_point[idir]) * ref_direction[idir];
  }
  return lcar;
}

/*----------------------------------------------------------------------------*/
/* Compute tangential velocity (for Schlichting model)                        */
/*----------------------------------------------------------------------------*/

static cs_real_t
_compute_tangential_velocity(const cs_real_t  velocity[3],
                             const cs_real_t  normal_vector[3],
                             const cs_real_t  coeff)
{
  cs_real_t u_square = 0.0;
  cs_real_t u_normal = 0.0;
  for (int idir = 0; idir < 3; idir++) {
    u_normal += velocity[idir] * normal_vector[idir] * coeff;
    u_square += velocity[idir] * velocity[idir];
  }
  return sqrt(u_square - u_normal * u_normal);
}

/*----------------------------------------------------------------------------
 * Compute forced convection exchange
 *
 * TODO: this function uses costly "cs_field_by_name" calls at each
 *       element. The loop on faces should be inside this function,
 *       not outside.
 *----------------------------------------------------------------------------*/

static void
_compute_exchange_forced_convection(cs_lnum_t   ieltcd,
                                    cs_real_t  *hconv,
                                    cs_real_t  *hcond)
{
  // Mesh related indices
  const cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;
  const cs_lnum_t  ifac   = _wall_cond.ifbpcd[ieltcd];
  const cs_lnum_t  iel    = ifabor[ifac];
  const int        iz     = _wall_cond.izzftcd[ieltcd];

  // Prandtl number
  const cs_real_t dyn_visc = cs_field_by_name("molecular_viscosity")->val[iel];
  const cs_real_t lambda_over_cp
    = cs_field_by_name("thermal_diffusivity")->val[iel];
  const cs_real_t pr_lam = dyn_visc / lambda_over_cp;

  // Schmidt number
  const cs_real_t mass_diffusion
    = cs_field_by_name("steam_binary_diffusion")->val[iel];
  const cs_real_t rho      = cs_field_by_name("density")->val[iel];
  const cs_real_t kin_visc = dyn_visc / rho;
  const cs_real_t Sch      = kin_visc / mass_diffusion;

  // Mass and mole fractions of steam, in bulk and at interface
  const cs_real_t t_wall   = _get_wall_temperature(ieltcd);
  cs_real_t       psat     = _compute_psat(t_wall);
  const cs_real_t pressure = (cs_glob_velocity_pressure_model->idilat == 3)
                               ? cs_glob_fluid_properties->pther
                               : cs_glob_fluid_properties->p0;
  cs_real_t x_vap_int = psat / pressure;

  cs_field_t *              f_vap = cs_field_by_name("y_h2o_g");
  cs_gas_mix_species_prop_t s_vap;
  const int                 k_id = cs_gas_mix_get_field_key();
  cs_field_get_key_struct(f_vap, k_id, &s_vap);

  cs_real_t y_vap_core  = f_vap->val[iel];
  cs_real_t mol_mas_vap = s_vap.mol_mas;
  cs_real_t mix_mol_mas = cs_field_by_name("mix_mol_mas")->val[iel];
  cs_real_t x_vap_core
    = _compute_mole_fraction(y_vap_core, mix_mol_mas, mol_mas_vap);
  cs_real_t kcond;

  switch (_wall_cond.forced_conv_model) {

  case CS_WALL_COND_MODEL_SCHLICHTING:
    {
      cs_real_3_t *cdgfbo = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
      const cs_real_3_t *b_face_normal
        = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
      const cs_real_t *  surfbn = cs_glob_mesh_quantities->b_face_surf;
      const cs_real_3_t *velocity
        = (cs_real_3_t *)cs_field_by_name("velocity")->val;
      const cs_real_3_t *n_ref = (cs_real_3_t *)_wall_cond.zprojcond;
      const cs_real_3_t *x_ref = (cs_real_3_t *)_wall_cond.zxrefcond;
      cs_real_3_t        n_ref_norm;
      cs_math_3_normalize(n_ref[iz], n_ref_norm);
      const cs_real_t lcar
        = _compute_characteristic_length(cdgfbo[ifac], x_ref[iz], n_ref_norm);
      const cs_real_t u_ref
        = _compute_tangential_velocity(velocity[iel], b_face_normal[ifac],
                                       1. / surfbn[ifac]);
      // Reynolds number
      const cs_real_t re    = rho * u_ref * lcar / dyn_visc;
      cs_real_t       theta = 1.0;
      if (x_vap_int < x_vap_core) {
        // here a choice between the different suction coefficients must be made
        // => we choose to use the updated value of Benteboula and Dabbene
        theta = 0.8254 + 0.616 * (x_vap_core - x_vap_int) / (1.0 - x_vap_int);
      }
      const cs_real_t nu  = _compute_schlichting(theta, re, pr_lam);
      *hconv              = lambda_over_cp * nu / lcar;
      const cs_real_t She = _compute_schlichting(theta, re, Sch);
      kcond               = mass_diffusion * She / lcar;
    }
    break;

  case CS_WALL_COND_MODEL_WALL_LAW:
  default:
    {
      const cs_real_t   rough_t = 0.0;
      cs_real_t         uk      = 0.0;
      const cs_field_t *f_ustar = cs_field_by_name_try("boundary_ustar");
      if (f_ustar != nullptr) {
        uk = f_ustar->val[ifac];
      }
      const cs_real_t          dplus  = 0.0;
      const cs_real_t          yplus  = cs_field_by_name("yplus")->val[ifac];
      const cs_wall_f_s_type_t iwalfs = cs_glob_wall_functions->iwalfs;
      cs_real_t                hpflui, ypth;
      cs_wall_functions_scalar(iwalfs,
                               kin_visc,
                               pr_lam,
                               pr_tur,
                               rough_t,
                               uk,
                               yplus,
                               dplus,
                               &hpflui,
                               &ypth);

      const cs_real_t distbf = cs_glob_mesh_quantities->b_dist[ifac];
      *hconv                 = lambda_over_cp * hpflui / distbf;
      cs_wall_functions_scalar(iwalfs, kin_visc, Sch, pr_tur, rough_t,
                               uk, yplus, dplus, &hpflui, &ypth);
      kcond = mass_diffusion * hpflui / distbf;
    }
    break;
  }

  // Spalding number (only for condensation)
  cs_real_t mol_mas_ncond = cs_field_by_name("mol_mas_ncond")->val[iel];
  cs_real_t mol_mas_int
    = x_vap_int * mol_mas_vap + (1.0 - x_vap_int) * mol_mas_ncond;
  cs_real_t y_vap_int = x_vap_int * mol_mas_vap / mol_mas_int;
  cs_real_t spalding  = 0.0;
  if (y_vap_int < y_vap_core)
    spalding = (y_vap_core - y_vap_int) / (1.0 - y_vap_int);

  cs_real_t t_inf      = cs_field_by_name("tempk")->val[iel];
  cs_real_t delta_temp = t_inf - t_wall;

  *hcond = 0.0;
  if (fabs(delta_temp) > cs_math_epzero)
    *hcond = kcond * rho * spalding * lcond / delta_temp;

  cs_real_t cp = cs_field_by_name("specific_heat")->val[iel];
  *hconv *= cp;
}

/*----------------------------------------------------------------------------
 * Compute natural convection exchange
 *
 * TODO: this function uses costly "cs_field_by_name" calls at each
 *       element. The loop on faces should be inside this function,
 *       not outside.
 *----------------------------------------------------------------------------*/

static void
_compute_exchange_natural_convection(cs_lnum_t  ieltcd,
                                     cs_real_t *hconv,
                                     cs_real_t *hcond)
{
  constexpr cs_real_t c_2ov3 = 2./3.;

  cs_lnum_t *     ifabor   = cs_glob_mesh->b_face_cells;
  const cs_real_t pressure =   (cs_glob_velocity_pressure_model->idilat == 3)
                             ? cs_glob_fluid_properties->pther
                             : cs_glob_fluid_properties->p0;

  // Mesh quantities
  cs_lnum_t ifac = _wall_cond.ifbpcd[ieltcd];
  cs_lnum_t iel  = ifabor[ifac];

  // Laminar Prandtl
  cs_real_t lambda_over_cp = cs_field_by_name("thermal_diffusivity")->val[iel];
  cs_real_t rho            = cs_field_by_name("density")->val[iel];
  cs_real_t mu             = cs_field_by_name("molecular_viscosity")->val[iel];
  cs_real_t Pr             = mu / lambda_over_cp;

  // Schmidt number
  cs_real_t diff = cs_field_by_name("steam_binary_diffusion")->val[iel];
  cs_real_t Sc   = mu / (rho * diff);

  cs_field_t *              f_vap = cs_field_by_name("y_h2o_g");
  cs_gas_mix_species_prop_t s_vap;
  const int                 k_id = cs_gas_mix_get_field_key();
  cs_field_get_key_struct(f_vap, k_id, &s_vap);

  // Spalding number
  cs_real_t t_wall        = _get_wall_temperature(ieltcd);
  cs_real_t psat          = _compute_psat(t_wall);
  cs_real_t x_vap_int     = cs_math_fmin(1.0, psat / pressure);
  cs_real_t mol_mas_ncond = cs_field_by_name("mol_mas_ncond")->val[iel];
  cs_real_t mol_mas_vap   = s_vap.mol_mas;
  cs_real_t mol_mas_int
    = x_vap_int * mol_mas_vap + (1.0 - x_vap_int) * mol_mas_ncond;
  cs_real_t y_vap_int  = x_vap_int * mol_mas_vap / mol_mas_int;
  cs_real_t y_vap_core = f_vap->val[iel];
  cs_real_t spalding   = 0.0;
  if (y_vap_int < y_vap_core)
    spalding = (y_vap_core - y_vap_int) / (1.0 - y_vap_int);

  // h_cond / h_conv ratio (Chilton-Colburn analogy)
  cs_real_t cp           = cs_field_by_name("specific_heat")->val[iel];
  cs_real_t t_inf        = cs_field_by_name("tempk")->val[iel];
  cs_real_t delta_temp   = (t_inf - t_wall);
  cs_real_t hcd_over_hcv = 0.0; // default : no condensation
  if (fabs(delta_temp) > cs_math_epzero)
    hcd_over_hcv
      = 1.0 / cp * pow(Pr / Sc, c_2ov3) * lcond / delta_temp * spalding;

  // h_conv computation
  cs_real_t       theta   = 1.0;
  const cs_real_t gravity = cs_math_3_norm(cs_glob_physical_constants->gravity);
  switch (_wall_cond.natural_conv_model) {

  case CS_WALL_COND_MODEL_COPAIN:
    {
      cs_real_t drho = fabs((t_inf - t_wall) / t_inf);
      if (y_vap_int < y_vap_core) { // condensation
        cs_real_t mix_mol_mas = cs_field_by_name("mix_mol_mas")->val[iel];
        cs_real_t x_vap_core
          = _compute_mole_fraction(y_vap_core, mix_mol_mas, mol_mas_vap);
        theta = 1.0 + 0.625 * (x_vap_core - x_vap_int) / (1.0 - x_vap_int);
      }
      cs_real_t gr = _compute_grashof(gravity, drho, lcar_nc, mu / rho);
      cs_real_t nu = _compute_mac_adams(theta, gr, Pr);
      *hconv       = cp * lambda_over_cp * nu / lcar_nc;
    }
    break;

  case CS_WALL_COND_MODEL_COPAIN_BD:
    {
      cs_real_t drho
        = fabs(1. - t_wall / t_inf
               + (y_vap_core - y_vap_int)
               / (  mol_mas_ncond / (mol_mas_ncond - mol_mas_vap)
                  - 1.0 + y_vap_int));
      if (y_vap_int < y_vap_core) { // condensation
        cs_real_t mix_mol_mas = cs_field_by_name("mix_mol_mas")->val[iel];
        cs_real_t x_vap_core
          = _compute_mole_fraction(y_vap_core, mix_mol_mas, mol_mas_vap);
        theta = 0.8254 + 0.616 * (x_vap_core - x_vap_int) / (1.0 - x_vap_int);
      }
      cs_real_t gr = _compute_grashof(gravity, drho, lcar_nc, mu / rho);
      cs_real_t nu = _compute_mac_adams(theta, gr, Pr);
      *hconv       = cp * lambda_over_cp * nu / lcar_nc;
    }
    break;

  case CS_WALL_COND_MODEL_UCHIDA:
    {
      cs_real_t h_uchida = 0.0;
      if ((y_vap_core > 0)
          && (y_vap_core < 1)) // Uchida correlation is developed for steam
                               // fractions between 0.10 and 0.80
        h_uchida = 380.0 * pow(y_vap_core / (1.0 - y_vap_core), 0.7);
      _wall_cond.total_htc[ieltcd] = h_uchida;
      *hconv                       = h_uchida / (1.0 + hcd_over_hcv);
    }
    break;

  case CS_WALL_COND_MODEL_DEHBI:
    {
      cs_real_t rho_wall
        = pressure * mol_mas_int / (cs_physical_constants_r * t_wall);
      cs_real_t drho = fabs(rho_wall - rho) / rho; // different from Fortran
      if (y_vap_int < y_vap_core) {
        cs_real_t B = (y_vap_int - y_vap_core) / (1.0 - y_vap_int);
        theta = 1.33 * log(1.0 + B) / B;
      }
      cs_real_t gr = _compute_grashof(gravity, drho, lcar_nc, mu / rho);
      cs_real_t nu = _compute_mac_adams(theta, gr, Pr);
      *hconv       = cp * lambda_over_cp * nu / lcar_nc;
    }
    break;
  }

  *hcond = *hconv * hcd_over_hcv;
}

/*----------------------------------------------------------------------------
 * Compute mixed convection exchange
 *
 * TODO: this function uses costly "cs_field_by_name" calls at each
 *       element. The loop on faces should be inside this function,
 *       not outside.
 *----------------------------------------------------------------------------*/

static void
_compute_exchange_mixed_convection(cs_lnum_t  ieltcd,
                                   cs_real_t *hconv,
                                   cs_real_t *hcond)
{
  cs_real_t hconv_fc, hcond_fc, hconv_nc, hcond_nc;
  _compute_exchange_forced_convection(ieltcd, &hconv_fc, &hcond_fc);
  _compute_exchange_natural_convection(ieltcd, &hconv_nc, &hcond_nc);

  switch (_wall_cond.mixed_conv_model) {
  case CS_WALL_COND_MIXED_MAX:
    *hcond = cs::max(hcond_fc, hcond_nc);
    *hconv = cs::max(hconv_fc, hconv_nc);
    break;
  case CS_WALL_COND_MIXED_INCROPERA:;
    const cs_lnum_t *  ifabor = cs_glob_mesh->b_face_cells;
    const cs_lnum_t    ifac   = _wall_cond.ifbpcd[ieltcd];
    const cs_lnum_t    iel    = ifabor[ifac];
    const cs_real_3_t *velocity
      = (cs_real_3_t *)cs_field_by_name("velocity")->val;
    const cs_real_t g_dot_u = cs_math_3_dot_product(
      cs_glob_physical_constants->gravity, velocity[iel]);
    if (g_dot_u > 0) { // buoyancy-aided flow
      *hcond = fabs(hcond_fc - hcond_nc);
      *hconv = fabs(hconv_fc - hconv_nc);
    }
    else { // buoyance-opposed flow
      *hcond = fabs(hcond_fc + hcond_nc);
      *hconv = fabs(hconv_fc + hconv_nc);
    }
    break;
  }
}

/*----------------------------------------------------------------------------
 * Compute exchange coefficients
 *
 * TODO: this function uses costly "cs_field_by_name" calls at each
 *       element. The loop on faces should be inside this function,
 *       not outside.
 *----------------------------------------------------------------------------*/

static void
_compute_exchange_coefficients(cs_lnum_t  ieltcd,
                               cs_real_t *hconv,
                               cs_real_t *hcond)
{
  const int iz = _wall_cond.izzftcd[ieltcd];
  // For now, izcophc and izcophg are kept for backward compatibility
  if (_wall_cond.izcophc[iz] != _wall_cond.izcophg[iz]) {
    bft_printf("  -- Warning: izcophc and izcophg have different values for "
               "zone: %d\n",
               iz);
    bft_printf("              izcophg is set to the value of izcophg: %d\n",
               _wall_cond.izcophc[iz]);
    _wall_cond.izcophg[iz] = _wall_cond.izcophc[iz];
  }
  switch (_wall_cond.izcophc[iz]) {
  case 1:
    _compute_exchange_forced_convection(ieltcd, hconv, hcond);
    break;
  case 2:
    _compute_exchange_natural_convection(ieltcd, hconv, hcond);
    break;
  case 3:
    _compute_exchange_mixed_convection(ieltcd, hconv, hcond);
    break;
  }
}

/*--------------------------------------------------------------------------- */
/*!
 * \brief  Compute natural convection exchange for volume wall
 *
 * \param[out] hcond  array of condensation exchange coefficients
 * \param[out] hconv  array of convective exchange coefficients
 */
/*----------------------------------------------------------------------------*/

static void
_compute_exchange_natural_convection_volume_structure(cs_real_t   *hcond,
                                                      cs_real_t   *hconv)
{
  constexpr cs_real_t c_2ov3 = 2./3.;

  const cs_real_t pressure =   (cs_glob_velocity_pressure_model->idilat == 3)
                             ? cs_glob_fluid_properties->pther
                             : cs_glob_fluid_properties->p0;

  cs_real_t *cpro_venth = cs_field_by_name("thermal_diffusivity")->val;
  cs_real_t *cpro_rho   = cs_field_by_name("density")->val;
  cs_real_t *cpro_viscl = cs_field_by_name("molecular_viscosity")->val;
  cs_real_t *cpro_stmbd = cs_field_by_name("steam_binary_diffusion")->val;
  cs_real_t *cpro_mol_mass_ncond = cs_field_by_name("mol_mas_ncond")->val;
  cs_real_t *cpro_cp = cs_field_by_name("specific_heat")->val;
  cs_real_t *cpro_tempk = cs_field_by_name("tempk")->val;
  cs_real_t *cpro_mix_mol_mass = cs_field_by_name("mix_mol_mas")->val;

  cs_field_t *f_vap = cs_field_by_name("y_h2o_g");
  cs_gas_mix_species_prop_t s_vap;
  const int  k_id = cs_gas_mix_get_field_key();
  cs_field_get_key_struct(f_vap, k_id, &s_vap);

  for (cs_lnum_t ii = 0; ii < _wall_cond.ncmast; ii++) {

    /* Mesh quantities  */
    cs_lnum_t iel  = _wall_cond.ltmast[ii];
    cs_lnum_t vol_id = _wall_cond.izmast[ii];

    /* Laminar Prandtl */
    cs_real_t lambda_over_cp = cpro_venth[iel];
    cs_real_t rho            = cpro_rho[iel];
    cs_real_t mu             = cpro_viscl[iel];
    cs_real_t Pr             = mu / lambda_over_cp;

    /* Schmidt number */
    cs_real_t diff = cpro_stmbd[iel];
    cs_real_t Sc   = mu / (rho * diff);

    /* Spalding number */
    cs_real_t t_wall = (_wall_cond.itagms[vol_id] == 1) ?
                       cs_glob_wall_cond_0d_thermal->volume_t[ii][0]:
                       cs_glob_wall_cond_0d_thermal->volume_t0[vol_id];

    t_wall += cs_physical_constants_celsius_to_kelvin;

    cs_real_t psat          = _compute_psat(t_wall);
    cs_real_t x_vap_int     = cs::min(1.0, psat / pressure);
    cs_real_t mol_mas_ncond = cpro_mol_mass_ncond[iel];
    cs_real_t mol_mas_vap   = s_vap.mol_mas;
    cs_real_t mol_mas_int
      = x_vap_int * mol_mas_vap + (1.0 - x_vap_int) * mol_mas_ncond;
    cs_real_t y_vap_int  = x_vap_int * mol_mas_vap / mol_mas_int;
    cs_real_t y_vap_core = f_vap->val[iel];
    cs_real_t spalding   = 0.0;
    if (y_vap_int < y_vap_core)
      spalding = (y_vap_core - y_vap_int) / (1.0 - y_vap_int);

    /* h_cond / h_conv ratio (Chilton-Colburn analogy) */
    cs_real_t cp           = cpro_cp[iel];
    cs_real_t t_inf        = cpro_tempk[iel];
    cs_real_t delta_temp   = (t_inf - t_wall);
    cs_real_t hcd_over_hcv = 0.0; // default : no condensation
    if (fabs(delta_temp) > cs_math_epzero)
      hcd_over_hcv
        = 1.0 / cp * pow(Pr / Sc, c_2ov3) * lcond / delta_temp * spalding;

    /* h_conv computation */
    cs_real_t theta = 1.0;
    cs_real_t gravity = cs_math_3_norm(cs_glob_physical_constants->gravity);

    cs_real_t drho = fabs((t_inf - t_wall) / t_inf);
    if (y_vap_int < y_vap_core) { // condensation
      cs_real_t mix_mol_mas = cpro_mix_mol_mass[iel];
      cs_real_t x_vap_core
        = _compute_mole_fraction(y_vap_core, mix_mol_mas, mol_mas_vap);
      theta = 1.0 + 0.625 * (x_vap_core - x_vap_int) / (1.0 - x_vap_int);
    }
    cs_real_t gr = _compute_grashof(gravity, drho, lcar_nc, mu / rho);
    cs_real_t nu = _compute_mac_adams(theta, gr, Pr);
    hconv[ii] = cp * lambda_over_cp * nu / lcar_nc;
    hcond[ii] = hconv[ii] * hcd_over_hcv;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_condensation_t *
cs_get_glob_wall_condensation(void)
{
  return &_wall_cond;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the wall condensation model
 *
 * \param[in] model    integer corresponding to the desired model
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_model(cs_wall_cond_natural_conv_model_t  model)
{
  _wall_cond.natural_conv_model = model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the onoff state of wall condensation modeling
 *
 * \param[in] icondb  integer corresponding to the onoff state (-1: off, 0: on)
 * \param[in] icondv  integer corresponding to the onoff state with
 *                    metal structures (-1: off, 0: on)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_onoff_state(int  icondb,
                                     int  icondv)
{
  _wall_cond.icondb = icondb;
  _wall_cond.icondv = icondv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_create(void)
{
  int n_var = 0;
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    n_var += f->dim;
  }

  if (_wall_cond.nzones < 1)
    _wall_cond.nzones = 1;
  if (_wall_cond.nvolumes < 1)
    _wall_cond.nvolumes = 1;

  // Mesh related quantities
  CS_MALLOC_HD(_wall_cond.ifbpcd, _wall_cond.nfbpcd, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.itypcd, _wall_cond.nfbpcd * n_var, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.izzftcd, _wall_cond.nfbpcd, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.spcond, _wall_cond.nfbpcd * n_var, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.hpcond, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.twall_cond, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.thermal_condensation_flux, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.convective_htc, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.condensation_htc, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.total_htc, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.flthr, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.dflthr, _wall_cond.nfbpcd, cs_real_t, cs_alloc_mode);

  // Zone related quantities
  CS_MALLOC(_wall_cond.izcophc, _wall_cond.nzones, cs_lnum_t);
  CS_MALLOC(_wall_cond.izcophg, _wall_cond.nzones, cs_lnum_t);
  CS_MALLOC(_wall_cond.iztag1d, _wall_cond.nzones, cs_lnum_t);
  CS_MALLOC(_wall_cond.ztpar, _wall_cond.nzones, cs_real_t);
  CS_MALLOC(_wall_cond.zxrefcond, 3 * _wall_cond.nzones, cs_real_t);
  CS_MALLOC(_wall_cond.zprojcond, 3 * _wall_cond.nzones, cs_real_t);

  for (cs_lnum_t i = 0; i < _wall_cond.nfbpcd; i++) {
    _wall_cond.ifbpcd[i]                    = 0;
    _wall_cond.hpcond[i]                    = 0.0;
    _wall_cond.twall_cond[i]                = 0.0;
    _wall_cond.thermal_condensation_flux[i] = 0.0;
    _wall_cond.convective_htc[i]            = 0.0;
    _wall_cond.condensation_htc[i]          = 0.0;
    _wall_cond.total_htc[i]                 = 0.0;
    _wall_cond.flthr[i]                     = 0.0;
    _wall_cond.dflthr[i]                    = 0.0;
    if (_wall_cond.nzones <= 1) {
      _wall_cond.izzftcd[i] = 0;
    }
    else {
      _wall_cond.izzftcd[i] = -1;
    }
    for (cs_lnum_t j = 0; j < n_var; j++) {
      _wall_cond.itypcd[n_var * i + j] = 0;
      _wall_cond.spcond[n_var * i + j] = 0.0;
    }
  }

  for (cs_lnum_t i = 0; i < _wall_cond.nzones; i++) {
    _wall_cond.izcophc[i]           = 0;
    _wall_cond.izcophg[i]           = 0;
    _wall_cond.iztag1d[i]           = 0;
    _wall_cond.ztpar[i]             = -1.0;
    _wall_cond.zxrefcond[3 * i]     = 0.0;
    _wall_cond.zxrefcond[3 * i + 1] = 0.0;
    _wall_cond.zxrefcond[3 * i + 2] = 0.0;
    _wall_cond.zprojcond[3 * i]     = 0.0;
    _wall_cond.zprojcond[3 * i + 1] = 0.0;
    _wall_cond.zprojcond[3 * i + 2] = 0.0;
  }

  CS_MALLOC_HD(_wall_cond.ltmast, _wall_cond.ncmast, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.itypst, _wall_cond.ncmast * n_var, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.svcond, _wall_cond.ncmast * n_var, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.izmast, _wall_cond.ncmast, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(_wall_cond.flxmst, _wall_cond.ncmast, cs_real_t, cs_alloc_mode);
  CS_MALLOC(_wall_cond.itagms, _wall_cond.nvolumes, cs_lnum_t);

  cs_array_lnum_fill_zero(_wall_cond.ncmast, _wall_cond.ltmast);

  if (_wall_cond.nvolumes <= 1)
    cs_array_lnum_set_value(_wall_cond.ncmast, 0, _wall_cond.izmast);
  else
    cs_array_lnum_set_value(_wall_cond.ncmast, -1, _wall_cond.izmast);

  cs_array_lnum_fill_zero(_wall_cond.ncmast * n_var, _wall_cond.itypst);
  cs_array_real_fill_zero(_wall_cond.ncmast * n_var, _wall_cond.svcond);
  cs_array_real_fill_zero(_wall_cond.ncmast, _wall_cond.flxmst);
  cs_array_lnum_set_value(_wall_cond.nvolumes, 1, _wall_cond.itagms);

  cs_base_at_finalize(cs_wall_condensation_free);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_free(void)
{
  CS_FREE(_wall_cond.ifbpcd);
  CS_FREE(_wall_cond.itypcd);
  CS_FREE(_wall_cond.izzftcd);
  CS_FREE(_wall_cond.spcond);
  CS_FREE(_wall_cond.svcond);
  CS_FREE(_wall_cond.hpcond);
  CS_FREE(_wall_cond.twall_cond);
  CS_FREE(_wall_cond.thermal_condensation_flux);
  CS_FREE(_wall_cond.convective_htc);
  CS_FREE(_wall_cond.condensation_htc);
  CS_FREE(_wall_cond.total_htc);
  CS_FREE(_wall_cond.flthr);
  CS_FREE(_wall_cond.dflthr);

  CS_FREE(_wall_cond.izcophc);
  CS_FREE(_wall_cond.izcophg);
  CS_FREE(_wall_cond.iztag1d);
  CS_FREE(_wall_cond.ztpar);
  CS_FREE(_wall_cond.zxrefcond);
  CS_FREE(_wall_cond.zprojcond);

  CS_FREE(_wall_cond.ltmast);
  CS_FREE(_wall_cond.itypst);
  CS_FREE(_wall_cond.izmast);
  CS_FREE(_wall_cond.svcond);
  CS_FREE(_wall_cond.flxmst);
  CS_FREE(_wall_cond.itagms);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize wall condensation models.
 *
 * This includes building the associated meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_initialize(void)
{
  cs_wall_condensation_t *wall_cond = cs_get_glob_wall_condensation();
  cs_wall_cond_1d_thermal_t *wall_thermal = cs_get_glob_wall_cond_1d_thermal();
  cs_wall_cond_0d_thermal_t *wall_0d_thermal = cs_get_glob_wall_cond_0d_thermal();

  if (wall_thermal->nzones < 1)
    wall_thermal->nzones = 1;
  if (wall_0d_thermal->nvolumes < 1)
    wall_0d_thermal->nvolumes = 1;

  cs_wall_condensation_1d_thermal_create(wall_cond->nzones);

  cs_wall_condensation_0d_thermal_create(wall_0d_thermal->nvolumes,
                                         wall_cond->ncmast);

  /* We call cs_user_wall_condensation when there are cells with condensation
     source terms on at least one processor. We only fill the cell indirection
     array. However, cs_user_condensation is called on all processors in case
     the user has implemented global operations. */

  cs_user_wall_condensation(2);

  cs_wall_condensation_set_model(wall_cond->natural_conv_model);

  /* Verification if wall temperature is computed with
     a 1-D thermal model with implicit numerical scheme */
  wall_cond->nztag1d = 0;
  for (cs_lnum_t iz = 0; iz < wall_thermal->nzones; iz++) {
    if (wall_cond->iztag1d[iz] == 1)
      wall_cond->nztag1d
        = cs::max(wall_cond->iztag1d[iz], wall_cond->nztag1d);
  }

  cs_parall_max(1, CS_INT_TYPE, &wall_cond->nztag1d);

  if (wall_cond->nztag1d == 1) {
    // Compute maximal number of discretized points
    wall_thermal->znmurx = 0;
    for (cs_lnum_t iz = 0; iz < wall_cond->nzones; iz++) {
      wall_thermal->znmurx
        = cs::max(wall_thermal->znmur[iz], wall_thermal->znmurx);
    }

    cs_parall_max(1, CS_INT_TYPE, &wall_thermal->znmurx);

    cs_wall_condensation_1d_thermal_mesh_create(wall_thermal->znmurx,
                                                wall_cond->nfbpcd,
                                                wall_thermal->nzones);

    //1-D mesh generated and temperature initialization
    cs_wall_condensation_1d_thermal_mesh_initialize();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return condensing volume structures surface at each cell.
 *
 * \param[out]  surf  array of volume structure surfaces at each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_volume_exchange_surf_at_cells(cs_real_t  *surf)
{
  const cs_wall_cond_0d_thermal_t *wall_cond_0d_thermal
    = cs_glob_wall_cond_0d_thermal;

  const cs_real_t *vol_surf = wall_cond_0d_thermal->volume_surf;
  const cs_real_t *vol_measure = wall_cond_0d_thermal->volume_measure;

  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_lnum_t ncmast = _wall_cond.ncmast;
  const cs_lnum_t *ltmast = _wall_cond.ltmast;
  const cs_lnum_t *izmast = _wall_cond.izmast;

  for (cs_lnum_t ii = 0; ii < ncmast; ii++) {
    const cs_lnum_t c_id = ltmast[ii];
    const cs_lnum_t vol_id = izmast[ii];
    surf[ii] = vol_surf[vol_id]*cell_vol[c_id]/vol_measure[vol_id];
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Reset values
 *
 * \param[in]  wall_cond  wall_condensation strucutre
 * \param[in]  n_var number of variable
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_reset(cs_wall_condensation_t *wall_cond, const int n_var)
{
  for (int ii = 0; ii < wall_cond->nfbpcd; ii++) {
    wall_cond->hpcond[ii] = 0.0;
    for (int ivar = 0; ivar < n_var; ivar++) {
      wall_cond->itypcd[ivar * wall_cond->nfbpcd + ii] = 0;
      wall_cond->spcond[ivar * wall_cond->nfbpcd + ii] = 0.0;
    }
  }

  for (cs_lnum_t ii = 0; ii < wall_cond->ncmast; ii++) {
    wall_cond->flxmst[ii] = 0.0;
    for (int ivar = 0; ivar < n_var; ivar++) {
      wall_cond->itypst[ivar * wall_cond->ncmast + ii] = 0;
      wall_cond->svcond[ivar * wall_cond->ncmast + ii] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the wall condensation source terms.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_compute(cs_real_t  total_htc[])
{
  cs_field_t *f          = cs_field_by_name("pressure");
  const int   var_id_key = cs_field_key_id("variable_id");
  const int   ipr        = cs_field_get_key_int(f, var_id_key) - 1;

  cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;

  const cs_real_t *cpro_cp = cs_field_by_name("specific_heat")->val;
  const cs_real_t *cpro_tempk = cs_field_by_name("tempk")->val;

  for (int ii = 0; ii < _wall_cond.nfbpcd; ii++) {

    cs_lnum_t ifac = _wall_cond.ifbpcd[ii];
    cs_lnum_t iel  = ifabor[ifac];

    /* Exchange coefficients */

    cs_real_t hconv, hcond;
    _compute_exchange_coefficients(ii, &hconv, &hcond);
    _wall_cond.convective_htc[ii]   = hconv;
    _wall_cond.condensation_htc[ii] = hcond;
    _wall_cond.total_htc[ii]        = hconv + hcond;
    cs_real_t cp          = cpro_cp[iel];
    _wall_cond.hpcond[ii] = hconv / cp;
    total_htc[ii]         = _wall_cond.total_htc[ii];

    /* Thermal flux */

    const cs_real_t t_wall = _get_wall_temperature(ii);
    const cs_real_t t_inf  = cpro_tempk[iel];
    cs_real_t       flux   = _wall_cond.total_htc[ii] * (t_inf - t_wall);
    _wall_cond.thermal_condensation_flux[ii] = flux;

    _wall_cond.spcond[ipr * _wall_cond.nfbpcd + ii]
      -= hcond * (t_inf - t_wall) / lcond;
    int iz = _wall_cond.izzftcd[ii];
    if (_wall_cond.iztag1d[iz] == 1) {
      _wall_cond.flthr[ii]  = flux;
      _wall_cond.dflthr[ii] = 0.0;
    }
  }

  /* Compute condensation on volume structures */

  cs_wall_cond_0d_thermal_t *_0d_thermal = cs_get_glob_wall_cond_0d_thermal();

  if ((cs_glob_time_step->nt_cur == 1) && !cs_restart_present()) {
    for (cs_lnum_t ii = 0; ii < _wall_cond.ncmast; ii++) {

      cs_lnum_t vol_id = _wall_cond.izmast[ii];
      if (_wall_cond.itagms[vol_id] == 1) {
        _0d_thermal->volume_t[ii][0]=_0d_thermal->volume_t0[vol_id];
        _0d_thermal->volume_t[ii][1]=_0d_thermal->volume_t0[vol_id];
      }
    }

  }

  cs_real_t *vol_hcond = nullptr, *vol_hconv = nullptr;
  CS_MALLOC(vol_hcond, _wall_cond.ncmast, cs_real_t);
  CS_MALLOC(vol_hconv, _wall_cond.ncmast, cs_real_t);

  _compute_exchange_natural_convection_volume_structure(vol_hcond, vol_hconv);

  for (cs_lnum_t ii = 0; ii < _wall_cond.ncmast; ii++) {

    cs_lnum_t vol_id = _wall_cond.izmast[ii];
    cs_lnum_t c_id = _wall_cond.ltmast[ii];

    cs_real_t t_wall = (_wall_cond.itagms[vol_id] == 1) ?
                       _0d_thermal->volume_t[ii][0]:
                       _0d_thermal->volume_t0[vol_id];

    t_wall += cs_physical_constants_celsius_to_kelvin;

    cs_real_t t_inf = cpro_tempk[c_id];

    _wall_cond.svcond[ipr*_wall_cond.ncmast + ii] -= vol_hcond[ii]
                                                   * (t_inf - t_wall)/lcond;

    if (_wall_cond.itagms[vol_id] == 1) {
      _wall_cond.flxmst[ii] = vol_hconv[ii]*(t_inf - t_wall);
    }
  }

  CS_FREE(vol_hcond);
  CS_FREE(vol_hconv);

  if (cs_log_default_is_active())
    cs_wall_condensation_log();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output statistics about wall condensation source terms (for user log)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_log(void)
{
  const cs_field_t *f          = cs_field_by_name("pressure");
  const int         var_id_key = cs_field_key_id("variable_id");
  const int         ipr        = cs_field_get_key_int(f, var_id_key) - 1;

  if (_wall_cond.icondb == 0) {

    cs_real_t gamma_cond = 0.0;
    cs_real_t h_conv_min = cs_math_infinite_r;
    cs_real_t h_cond_min = cs_math_infinite_r;
    cs_real_t flux_min   = cs_math_infinite_r;
    cs_real_t h_conv_max = -cs_math_infinite_r;
    cs_real_t h_cond_max = -cs_math_infinite_r;
    cs_real_t flux_max   = -cs_math_infinite_r;

    for (cs_lnum_t ii = 0; ii < _wall_cond.nfbpcd; ii++) {
      gamma_cond += _wall_cond.spcond[ipr * _wall_cond.nfbpcd + ii];
      h_conv_min = cs::min(h_conv_min, _wall_cond.convective_htc[ii]);
      h_conv_max = cs::max(h_conv_max, _wall_cond.convective_htc[ii]);
      h_cond_min = cs::min(h_cond_min, _wall_cond.condensation_htc[ii]);
      h_cond_max = cs::max(h_cond_max, _wall_cond.condensation_htc[ii]);
      flux_min = cs::min(flux_min,
                         _wall_cond.thermal_condensation_flux[ii]);
      flux_max = cs::max(flux_max,
                         _wall_cond.thermal_condensation_flux[ii]);
    }

    if (cs_glob_n_ranks > 1) {

      cs_real_t r_buffer[6] = {
        -h_conv_min, h_conv_max,
        -h_cond_min, h_cond_max,
        -flux_min, flux_max
      };

      cs_parall_max(6, CS_REAL_TYPE, r_buffer);

      h_conv_min = - r_buffer[0];
      h_conv_max =   r_buffer[1];
      h_cond_min = - r_buffer[2];
      h_cond_max =   r_buffer[3];
      flux_min   = - r_buffer[4];
      flux_max   =   r_buffer[5];

      cs_parall_sum(1, CS_REAL_TYPE, &gamma_cond);

    }

    bft_printf
      (" Minmax values of convective HTC [W/m2/K]   : %15.12e    %15.12e\n",
       h_conv_min, h_conv_max);
    bft_printf
      (" Minmax values of condensation HTC [W/m2/K] : %15.12e    %15.12e\n",
       h_cond_min, h_cond_max);
    bft_printf
      (" Minmax values of thermal flux [W/m2]       : %15.12e    %15.12e\n",
       flux_min, flux_max);
    bft_printf(" Total mass sink term [kg/s/m2]             : %15.12e\n",
               gamma_cond);

  }

  if (_wall_cond.icondv == 0) {

    cs_real_t gamma_cond = 0.0;
    cs_real_t flux_min   = cs_math_infinite_r;
    cs_real_t flux_max   = -cs_math_infinite_r;

    for (cs_lnum_t ii = 0; ii < _wall_cond.ncmast; ii++) {
      gamma_cond += _wall_cond.svcond[ipr * _wall_cond.ncmast + ii];
      flux_min = cs::min(flux_min, _wall_cond.flxmst[ii]);
      flux_max = cs::max(flux_max, _wall_cond.flxmst[ii]);
    }

    if (cs_glob_n_ranks > 1) {

      cs_real_t r_buffer[3] = {-flux_min, flux_max};

      cs_parall_max(2, CS_REAL_TYPE, r_buffer);

      flux_min = - r_buffer[0];
      flux_max =   r_buffer[1];

      cs_parall_sum(1, CS_REAL_TYPE, &gamma_cond);

    }

    bft_printf(" Minmax values of thermal flux [W/m2] on volume structure:"
               " %15.12e    %15.12e\n", flux_min, flux_max);
    bft_printf(" Total mass sink term [kg/s/m2] on volume structure: %15.12e\n",
               gamma_cond);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Explicit and implicit sources terms from sources
 *        condensation computation.
 *
 * \param[in]      f         pointer to field structure
 * \param[in]      xcpp      array of specific heat (Cp)
 * \param[in]      pvara     variable value at time step beginning
 * \param[in,out]  st_exp    explicit source term part linear in the variable
 * \param[in,out]  st_imp    associated value with \c tsexp
 *                           to be stored in the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_source_terms(const cs_field_t  *f,
                                  const cs_real_t    xcpp[],
                                  const cs_real_t    pvara[],
                                  cs_real_t          st_exp[],
                                  cs_real_t          st_imp[])
{
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  const cs_real_t *restrict b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const cs_wall_condensation_t *wall_cond = cs_glob_wall_condensation;
  const cs_lnum_t nfbpcd = wall_cond->nfbpcd;
  const cs_lnum_t *ifbpcd = wall_cond->ifbpcd;
  const cs_lnum_t *itypcd = wall_cond->itypcd;
  const cs_real_t *spcond = wall_cond->spcond;

  const cs_lnum_t ncmast = wall_cond->ncmast;
  const cs_lnum_t *ltmast = wall_cond->ltmast;
  const cs_lnum_t *itypst = wall_cond->itypst;
  const cs_real_t *svcond = wall_cond->svcond;
  const cs_real_t *flxmst = wall_cond->flxmst;

  const int var_id_key = cs_field_key_id("variable_id");
  const int ivar = cs_field_get_key_int(f, var_id_key) - 1;

  const int ipr = cs_field_get_key_int(CS_F_(p), var_id_key) - 1;

  if (wall_cond->icondb == 0) {

    const cs_real_t *gamma = spcond + ipr*nfbpcd;
    const cs_real_t *spcondp = spcond + ivar*nfbpcd;

    /*  Compute condensation sourrece terms associated to surface zones */
    for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
      const cs_lnum_t f_id = ifbpcd[ii];
      const cs_lnum_t c_id = b_face_cells[f_id];

      /* Compute the explicit and implicit term of the condensation model */
      if (gamma[ii] <= 0)
        st_exp[c_id] += -b_face_surf[f_id]*gamma[ii]*pvara[c_id]*xcpp[c_id];
      else
        st_imp[c_id] += b_face_surf[f_id]*gamma[ii]*xcpp[c_id];

      if (itypcd[ivar*nfbpcd + ii] == 1)
        st_exp[c_id] += b_face_surf[f_id]*gamma[ii]*spcondp[ii]*xcpp[c_id];
    }
  }

  if (wall_cond->icondv == 0) {

    const cs_real_t *gamma = svcond + ipr*ncmast;
    const cs_real_t *svcondp = svcond + ivar*ncmast;

    /* Compute condensation source terms associated to volume zones
       with the metal mass structures modelling */

    cs_real_t *surfbm = nullptr;
    CS_MALLOC(surfbm, ncmast, cs_real_t);

    cs_wall_condensation_volume_exchange_surf_at_cells(surfbm);

    for (cs_lnum_t ii = 0; ii < ncmast; ii++) {
      const cs_lnum_t c_id = ltmast[ii];

      /* Compute the explicit and implicit term of the condensation model */
      if (gamma[ii] <= 0)
        st_exp[c_id] -= surfbm[ii]*gamma[ii]*pvara[c_id]*xcpp[c_id];
      else
        st_imp[c_id] += surfbm[ii]*gamma[ii]*xcpp[c_id];

      if (itypst[ivar*ncmast + ii] == 1) {
        if (f == CS_F_(h))
          st_exp[c_id] += surfbm[ii]*(gamma[ii]*svcondp[ii]
                        - flxmst[ii]);
        else
          st_exp[c_id] += surfbm[ii]*gamma[ii]*svcondp[ii]*xcpp[c_id];
      }
    }

    CS_FREE(surfbm);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
