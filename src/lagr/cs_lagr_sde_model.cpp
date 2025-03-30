/*============================================================================
 * Methods for lagrangian equations
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "bft/bft_mem.h"

#include "base/cs_base.h"
#include "comb/cs_coal.h"
#include "comb/cs_coal_ht_convert.h"
#include "base/cs_math.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"

#include "mesh/cs_mesh.h"
#include "base/cs_thermal_model.h"
#include "turb/cs_turbulence_model.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_tracking.h"
#include "lagr/cs_lagr_prototypes.h"

#include "lagr/cs_lagr_sde.h"

#include "base/cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_sde_model.h"
#include "atmo/cs_air_props.h"
#include "ctwr/cs_ctwr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_sde_model.cpp
        Integration of Lagrangian stochastic diferential equations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Constants */

static const cs_real_t  _c_stephan = 5.6703e-8;
static const cs_real_t  _tkelvi = 273.15;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square of a real value.
 *
 * \param[in]  x  value
 *
 * \return the square of the given value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_pow2(cs_real_t  x)
{
  return x*x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cube of a real value.
 *
 * \param[in]  x  value
 *
 * \return the square of the given value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_pow3(cs_real_t  x)
{
  return x*x*x;
}

/*----------------------------------------------------------------------------
 * Evolution of a coal or biomass particle
 *
 * - compute thermal fluxes (radiation and conduction)
 * - solve the heat diffusion equation inside a coal or biomass particle
 *   accounting for volume fluxes (chemical reaction):
 *   rho cp dT/dt = div( lambda grad(T) ) + phi
 *
 * parameters:
 *   npt         <--  particle id
 *   dt_part     <--  time step associated to the particle
 *   nor         <--  current step id (for 2nd order scheme)
 *   layer_vol   <--  volume occuppied by one layer
 *   tempct      <--  characteristic thermal time
 *   radius      <--  radius of each layer
 *   mlayer      <--  mass of each layer
 *   phith       <--  thermal source term due to chemical reaction (1/layer)
 *   temp        -->  particle temperature after evolution
 *----------------------------------------------------------------------------*/

static void
_lagtmp(cs_lnum_t         npt,
        const cs_real_t   dt_part,
        int               nor,
        cs_real_t         layer_vol,
        const cs_real_2_t tempct,
        const cs_real_t   radius[],
        const cs_real_t   mlayer[],
        const cs_real_t   phith[],
        cs_real_t         temp[])
{
  int l_id;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  int nlayer = cs_glob_lagr_const_dim->nlayer;
  int l_id_max = nlayer - 1;

  cs_real_t delray[nlayer], radiusd[nlayer], rho[nlayer];

  cs_real_t a[nlayer], b[nlayer], c[nlayer], d[nlayer];
  cs_real_t w1[nlayer], w2[nlayer];

  /* Particles management */

  cs_lagr_particle_set_t   *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  /* Initialization */

  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

  cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_DIAMETER);
  cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_MASS);
  cs_real_t p_init_diam = cs_lagr_particle_get_real(particle, p_am,
                                                    CS_LAGR_INITIAL_DIAMETER);
  cs_real_t p_shrink_diam = cs_lagr_particle_get_real(particle, p_am,
                                                      CS_LAGR_SHRINKING_DIAMETER);
  cs_real_t part_cp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP);

  const cs_coal_model_t  *coal_model = cs_glob_coal_model;

  const cs_real_t *part_temp
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_TEMPERATURE);
  const cs_real_t *prev_part_temp
    = cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                 CS_LAGR_TEMPERATURE);

  cs_real_t dd2 = cs_math_sq(p_diam);

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);
  cs_lnum_t co_id = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_COAL_ID);

  /* Multiple-layer resolution
     ------------------------- */

  if (nlayer > 1) {

    /* Compute radii and radii deltas */

    radiusd[0]  = radius[0] / 2.0;
    delray[0]   = radius[1] / 2.0;
    for (l_id = 1; l_id < l_id_max; l_id++) {
      radiusd[l_id]  = (radius[l_id - 1] + radius[l_id]) / 2.0;
      delray[l_id]  = (radius[l_id + 1] - radius[l_id - 1]) / 2.0;
    }
    radiusd[l_id_max]  = (radius[l_id_max - 1] + radius[l_id_max]) / 2.0;

    /* Compute density of layers */

    for (l_id = 0; l_id < nlayer; l_id++) {

      rho[l_id] = mlayer[l_id] / layer_vol;

      if (rho[l_id] <= 0.)
        bft_error(__FILE__, __LINE__, 0,
                  _("Particle layer %d with nonpositive (%g) density detected."),
                  l_id, rho[l_id]);

    }

    /* Conduction inside particle  */

    cs_real_t lambda = coal_model->thcdch[co_id];

    cs_real_t diamp2
      =          coal_model->xashch[co_id]  * cs_math_sq(p_init_diam)
        + (1.0 - coal_model->xashch[co_id]) * cs_math_sq(p_shrink_diam);

    cs_real_t tpscara = tempct[0] * diamp2 / dd2;

    cs_real_t coefh
      =   cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_MASS)
        * cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_CP)
        / (tpscara * cs_math_pi * diamp2);

    /* Equivalent radiative temperature */

    cs_real_t temprayo = pow(extra->rad_energy->val[cell_id]
                             / (4.0 * _c_stephan), 0.25);

    /* Build system (phith given by _lagich is in W) */

    /* layer 0 */

    cs_real_t prev_part_cp
      = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_CP);

    a[0]  = 0; /* unused */

    b[0]  =  1.0 + 4.0 * (lambda * dt_part) / (rho[0] * prev_part_cp)
                       * (  1.0 + 1.0 / (radius[1] * radius[0])
                         + 2.0 / (radius[1] * (radius[0] + radius[1])));

    c[0]  =  - 4.0 * (lambda * dt_part) / (rho[0] * prev_part_cp)
                   * (  1.0 + 1.0 / (radius[1] * radius[0])
                      + 2.0 / (radius[1] * (radius[0] + radius[1])));

    d[0]  = part_temp[0] + (phith[0] * dt_part) / (mlayer[0] * prev_part_cp);

    /* interior layers */

    for (l_id = 1; l_id < nlayer-1; l_id++) {

      cs_real_t f = (lambda * dt_part) / (  rho[l_id] * prev_part_cp
                                      * delray[l_id - 1] * delray[l_id]);

      a[l_id]  =  -f * (  2.0 * delray[l_id]
                        / (   delray[l_id - 1] + delray[l_id])
                           - (delray[l_id] / radiusd[l_id]));

      b[l_id]  = 1.0 + f * (2.0 - (  (delray[l_id] - delray[l_id - 1])
                                   / radiusd[l_id]));

      c[l_id]  =  -f * (  2.0 * delray[l_id - 1]
                              / (delray[l_id - 1] + delray[l_id])
                        + (delray[l_id - 1] / radiusd[l_id]));

      d[l_id]  =   part_temp[l_id]
                 + (phith[l_id] * dt_part) / (mlayer[l_id] * prev_part_cp);

    }

    /* last layer */

    l_id = nlayer - 1;

    cs_real_t f =   _c_stephan
                  * (cs_math_sq(temprayo) + cs_math_sq(part_temp[l_id]))
                  * (temprayo + part_temp[l_id]);

    cs_real_t  t_fluid_l
      =   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_TEMPERATURE_SEEN)
        + _tkelvi;

    a[l_id] = - (lambda * dt_part)
              / (rho[l_id] * prev_part_cp * delray[l_id-1])
              * (1.0 / delray[l_id-1] - 1.0 / radiusd[l_id]);

    b[l_id]  =  1.0 + (lambda * dt_part)
               / (rho[l_id] * prev_part_cp * delray[l_id-1])
               * (1.0 / delray[l_id-1] - 1.0 / radiusd[l_id])
               + (  dt_part * (coefh+f) / (rho[l_id]*prev_part_cp)
                  * (1.0 / delray[l_id-1] + 1.0 / radiusd[l_id]));

    d[l_id]  =   part_temp[l_id]
               + dt_part / (mlayer[l_id] * prev_part_cp)
               * (  phith[l_id] + (coefh*t_fluid_l + f*temprayo)
                  * layer_vol * (1.0 / delray[l_id-1] + 1.0/radiusd[l_id]));

    /* Solve system; we apply the Thomas algorithm:
       a_i T_i-1 + b_i T_i + c_i T_i+1 = d_i
       T_i = w2_i - w1_i T_i+1 */

    /* Relation between T_1 and T_2 is known */

    w1[0] = c[0] / b[0];
    w2[0] = d[0] / b[0];

    /* Compute w1_i and w2_i by recurrence */

    for (l_id = 1; l_id < nlayer; l_id++) {
      w1[l_id] = c[l_id] / (b[l_id] - w1[l_id - 1] * a[l_id]);
      w2[l_id] =   (d[l_id] - w2[l_id - 1] * a[l_id])
                 / (b[l_id] - w1[l_id - 1] * a[l_id]);
    }

    l_id = nlayer-1;

    temp[l_id] = w2[l_id];

    for (l_id = nlayer-2; l_id >= 0; l_id--)
      temp[l_id] = w2[l_id] - w1[l_id] * temp[l_id + 1];

  }

  /* Single-layer resolution
     ----------------------- */

  else if (nlayer == 1) {

    cs_real_t diamp2
      =          coal_model->xashch[co_id]  * cs_math_sq(p_init_diam)
        + (1.0 - coal_model->xashch[co_id]) * cs_math_sq(p_shrink_diam);

    cs_real_t tpscara = tempct[0] * diamp2 / dd2;

    /* Radiation */

    cs_real_t phirayo   =    extra->rad_energy->val[cell_id] / 4.0
                          - _c_stephan * pow(part_temp[0], 4);

    cs_real_t aux1      =  cs_lagr_particle_get_real(particle, p_am,
                                                     CS_LAGR_TEMPERATURE_SEEN)
                         + _tkelvi
                         + tpscara * (phirayo * cs_math_pi * diamp2 + phith[0])
                         / (p_mass * part_cp);
    cs_real_t aux2      = exp (-dt_part / tpscara);

    /* Source term */

    cs_real_t *part_ts_temp = nullptr;
    if (p_set->p_am->source_term_displ != nullptr) {
      if (p_set->p_am->source_term_displ[CS_LAGR_TEMPERATURE] >= 0)
        part_ts_temp = cs_lagr_particles_source_terms(p_set, npt,
                                                      CS_LAGR_TEMPERATURE);
    }

    if (nor == 1) {

      if (part_ts_temp != nullptr)
        part_ts_temp[0] =  0.5 * prev_part_temp[0] * aux2
                         + ( -aux2 + (1.0 - aux2) * tpscara / dt_part) * aux1;

      temp[0] = prev_part_temp[0] * aux2 + (1.0 - aux2) * aux1;

    }
    else if (nor == 2) {

      temp[0] =  part_ts_temp[0]
               + 0.5 * prev_part_temp[0] * aux2
               + (1.0 - (1.0 - aux2) * tpscara / dt_part) * aux1;

    }

  }

}

/*----------------------------------------------------------------------------
 * Compute particle evaporation using a pressure equilibrium model.
 *
 * - compute ath saturating partial pressure at particle temperature
 * - compute vapor flux leaving the particle
 *
 * Possibly, the flux is limited (particle has a fixed behavior over time:
 * either if vaporizes, either it condenses).
 *
 * parameters:
 *   npt         <--  particle id
 *   dt_part     <--  time step associated to the particle
 *   nor         <--  current step id (for 2nd order scheme)
 *   diftl0      <--  reference diffusivity for Enthalpy
 *   layer_vol   <--  volume occuppied by one layer
 *   tempct      <--  characteristic thermal time
 *   radius      <--  radius of each layer
 *   mwat_max    <--  maximum water mass of each layer
 *   sherw       <--  particle's Sherwod number
 *   mlayer      <--  mass of each layer
 *   fwat        -->  drying flux (kg/s) for each layer
 *----------------------------------------------------------------------------*/

static void
_lagsec(cs_lnum_t         npt,
        const cs_real_t   dt_part,
        int               nor,
        cs_real_t         diftl0,
        cs_real_t         layer_vol,
        cs_real_t         mwat_max,
        cs_real_t         sherw,
        cs_real_t         radius[],
        const cs_real_2_t tempct,
        cs_real_t         mlayer[],
        cs_real_t         mwater[],
        cs_real_t         fwat[])
{
  /* Particles management */

  cs_lagr_particle_set_t        *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t *p_am  = p_set->p_am;

  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

  cs_real_t prev_p_diam
    = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_DIAMETER);

  cs_real_t prev_p_cp
    = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_CP);

  cs_real_t *ptsvar = nullptr;
  if (p_set->p_am->source_term_displ != nullptr) {
    if (p_set->p_am->source_term_displ[CS_LAGR_TEMPERATURE] >= 0)
      ptsvar = cs_lagr_particles_source_terms(p_set, npt, CS_LAGR_TEMPERATURE);
  }

  const cs_coal_model_t *coal_model = cs_glob_coal_model;

  const int ih2o = coal_model->ih2o -1;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  cs_lnum_t nlayer = cs_glob_lagr_const_dim->nlayer;

  cs_real_t aux1, aux2, aux3;
  cs_real_t fwatsat[nlayer];
  cs_real_t phith[nlayer], temp[nlayer], tssauv[nlayer];

  cs_real_t precis = 1e-15;
  cs_real_t lv     = 2263000.0;
  cs_real_t tebl   = 100.0 + _tkelvi;
  cs_real_t tlimit = 302.24;
  cs_real_t tmini  = tlimit * (   1. - tlimit * cs_physical_constants_r
                               / (lv * coal_model->wmole[ih2o]));

  /* Compute water flux for layer l_id
     --------------------------------- */

  cs_real_t fwattot = 0.0;

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {
    fwat[l_id] = 0.0;
    fwatsat[l_id] = 0.0;
  }

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);

  /* find layer */

  cs_lnum_t l_id_wat = 0;

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {
    if (mwater[l_id] > 0.)
      l_id_wat = l_id;
  }

  cs_real_t *part_temp
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_TEMPERATURE);
  cs_real_t tpk = part_temp[l_id_wat];

  /* Compute mass fraction of saturating water */

  if (tpk >= tmini) {

    if (tpk >= tlimit) {

      aux1 = coal_model->wmole[ih2o] / extra->x_m->val[cell_id];
      aux2 = aux1 * exp(  lv * coal_model->wmole[ih2o]
                        * (1.0 / tebl - 1.0 / tpk)
                        / cs_physical_constants_r);

    }
    else {

      /* Linearize mass fraction of saturating water between tmini and Tlimit;
       * At Tlimit, the saturating water mass fraction is zero */

      aux1 = coal_model->wmole[ih2o] / extra->x_m->val[cell_id];
      aux2 =  aux1
            * exp(  lv * coal_model->wmole[ih2o]
                  * (1.0 / tebl - 1.0 / tlimit)
                  / cs_physical_constants_r)
            * lv * coal_model->wmole[ih2o]
            / (cs_physical_constants_r * _pow2(tlimit))
            * (tpk - tmini);

    }

    /* Compute diffused water source term */

    aux3      = cs::max(1.0 - aux2, precis);
    fwattot   =   cs_math_pi * prev_p_diam * diftl0 * sherw
                * log((1.0 - extra->x_eau->val[cell_id]) / aux3);

  }
  else
    fwattot = 0.0;  /* fwattot */

  /* Distribute this flux over neighboring cells */

  cs_real_t fwat_remain = fwattot;

  if (fwattot >= 0.0) { /* Particle dries towards its core */

    for (cs_lnum_t l_id = l_id_wat; l_id >= 0; l_id--) {

      /* cannot dry more than water present */
      fwat[l_id] = cs::min(mwater[l_id] / dt_part , fwat_remain);

      /* Update flux remaining to evaporate */
      fwat_remain = cs::max(0.0, fwat_remain - fwat[l_id]);

    }

  }
  else { /* Particle condenses towards exterior */

    for (cs_lnum_t l_id = l_id_wat; l_id < nlayer; l_id++) {

      if (l_id == nlayer - 1)
        /* With nlayer, condense all remaining flux */
        fwat[l_id] = fwat_remain;

      else
        /* Cannot condense more than water on 1 layer */
        fwat[l_id] = cs::max(-(mwat_max - mwater[l_id]) / dt_part,
                             fwat_remain);

      /* Flux remain a condenser  */
      fwat_remain = cs::min(0.0, fwat_remain - fwat[l_id]);

    }

  }

  /* Compute saturation fluxes
     ------------------------- */

  /* Limit flux relative to saturation temperature:

     we limit the drying flux so that at the end of the time step, the particle's
     enthalpy is sufficiently high so that its water saturation pressure is
     greater than the partial water pressure in the surrounding air. */

  /* Compute de tsat, saturating temperature at the air partial fraction */

  cs_real_t tsat;

  if (extra->x_eau->val[cell_id] > precis) {

    aux1 = coal_model->wmole[ih2o] / extra->x_m->val[cell_id];
    tsat = 1 / (  1 / tebl
                - cs_physical_constants_r
                  * log (extra->x_eau->val[cell_id] / aux1)
                  / (lv * coal_model->wmole[ih2o]));

    if (tsat < tlimit)
      tsat =  tmini
            + extra->x_eau->val[cell_id]
              / (aux1 * exp (  lv * coal_model->wmole[ih2o]
                             * (1.0 / tebl - 1.0 / tlimit)
                             / cs_physical_constants_r)
              * (lv * coal_model->wmole[ih2o])
              / (cs_physical_constants_r * _pow2(tlimit)));

  }
  else
    tsat = tmini;

  /* Compute temperature at the end of the time step with no chemistry */

  /* ignore all volume thermal source terms for this computation */
  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    phith[l_id]  = 0.0;

  /* Save the correction array for second order */
  if (ptsvar != nullptr) {
    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
      tssauv[l_id] = ptsvar[l_id];
  }

  _lagtmp(npt, dt_part, nor, layer_vol, tempct, radius, mlayer, phith, temp);

  /* Use array for second order correction */
  if (ptsvar != nullptr) {
    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
      ptsvar[l_id] = tssauv[l_id];
  }

  /* Compute evaporation/condensation flus so that  T_i = Tsat */

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    fwatsat[l_id] = mlayer[l_id] * prev_p_cp * (temp[l_id] - tsat)
                  / (lv * dt_part);

  /* Vapor limitation if needed
     -------------------------- */

  bool limiter = false;

  if (fwattot >= 0.0) {

    /* Particle dries towards its core */

    for (cs_lnum_t l_id = nlayer - 1; l_id >= 0; l_id--) {

      if (limiter == false) {

        /* Check that layer does not have an opposite behavior */

        if (fwatsat[l_id] < 0.0)
          limiter = true; /* block all following layers */

        /* limit flux */
        if (fwat[l_id] > fwatsat[l_id])
          fwat[l_id]  = cs::max(0.0, fwatsat[l_id]);

      }
      else
        fwat[l_id] = 0.0; /* limiter blocks interior layers */

    }

  }
  else if (fwattot < 0.0) {

    /* Check that exterior layers do not have an opposite behavior */

    for (cs_lnum_t l_id = nlayer - 1; l_id <= l_id_wat; l_id++) {

      if (fwatsat[l_id] > 0.0)
        limiter = true;  /* block all following layers */

    }

    /* Particle condenses towards exterior */

    for (cs_lnum_t l_id = l_id_wat; l_id < nlayer; l_id++) {

      if (!limiter) { /* limit flux */

        if (fwatsat[l_id] > fwat[l_id])
          fwat[l_id]  = cs::min(0.0, fwatsat[l_id]);

      }
      else
        fwat[l_id] = 0.0; /* limiter blocks */

    }

  }

}

/*----------------------------------------------------------------------------
 * Integrate stochastic differential equations by particles mass.
 *----------------------------------------------------------------------------*/

static void
_lagimp(const cs_lnum_t       npt,
        const cs_real_t       dt_part,
        int                   nor)
{
  /* Particles management */

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  cs_real_t tcarac;
  cs_real_t pip;

  /* FIXME check this: not called for fixed particle
   * yet a fixed particle's mass might still evolve */

  tcarac = 1.0;
  pip    = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_MASS);
  cs_lagr_sde_attr(CS_LAGR_MASS, npt, nor, dt_part, tcarac, pip);
}

/*----------------------------------------------------------------------------
 * Integrate SDE's for particle diameters
 *----------------------------------------------------------------------------*/

static void
_lagidp(const cs_lnum_t       npt,
        const cs_real_t       dt_part,
        int                   nor)
{
  /* Particles management */

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  cs_real_t tcarac;
  cs_real_t pip;

  /* FIXME check this: not called for fixed particle
   * yet a fixed particle's mass might still evolve */

  tcarac = 1.0;
  pip    = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DIAMETER);
  cs_lagr_sde_attr(CS_LAGR_DIAMETER, npt, nor, dt_part, tcarac, pip);
}

/*----------------------------------------------------------------------------
 * Integration of SDE's for particle temperature.
 *----------------------------------------------------------------------------*/

static void
_lagitp(const cs_lnum_t       npt,
        const cs_real_t       dt_part,
        int                   nor,
        const cs_real_2_t     tempct)
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

  cs_real_t tcarac;
  cs_real_t pip;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  /* ==========================================================================
   * REMPLISSAGE DU TEMPS CARACTERISTIQUE ET DU "PSEUDO SECOND MEMBRE"
   *=========================================================================== */

  /* FIXME check this: not called for fixed particle
   * yet a fixed particle's mass might still evolve */

  tcarac = tempct[0];

  pip = cs_lagr_particles_get_real_n(p_set, npt, 2 - nor,
                                         CS_LAGR_TEMPERATURE_SEEN);

  /* ======================================================================
  * Account for radiation of present
  * ======================================================================= */
  if (extra->radiative_model > 0) {
    /* use previous step for t_order == 1 or prediction step
     * and current one for correction step */
    cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                    CS_LAGR_CELL_ID);

    cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
    cs_real_t p_cp   = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP);
    cs_real_t p_eps  = cs_lagr_particle_get_real(particle, p_am,
                                                 CS_LAGR_EMISSIVITY);

    if (nor == 1) {

      cs_real_t prev_p_diam = cs_lagr_particle_get_real_n(particle, p_am,
                                                          1, CS_LAGR_DIAMETER);
      cs_real_t prev_p_temp = cs_lagr_particle_get_real_n(particle, p_am,
                                                          1, CS_LAGR_TEMPERATURE);

      cs_real_t srad =    cs_math_pi * pow(prev_p_diam, 2.0) * p_eps
                        * (extra->rad_energy->val[cell_id]
                      - 4.0 * _c_stephan * pow (prev_p_temp,4));
      pip =   cs_lagr_particle_get_real_n(particle, p_am, 1,
                                          CS_LAGR_TEMPERATURE_SEEN)
           + tcarac * srad / p_cp / p_mass;

    }
    else {

      cs_real_t p_diam = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                     CS_LAGR_DIAMETER);
      cs_real_t p_temp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                     CS_LAGR_TEMPERATURE);

      cs_real_t srad =    cs_math_pi * pow(p_diam, 2.0) * p_eps
                        * (extra->rad_energy->val[cell_id]
                      - 4.0 * _c_stephan *  pow(p_temp , 4));
      pip =  cs_lagr_particle_get_real(particle, p_am,
                                       CS_LAGR_TEMPERATURE_SEEN)
            + tcarac * srad / p_cp /p_mass;

    }

  }

  cs_lagr_sde_attr(CS_LAGR_TEMPERATURE, npt, nor, dt_part, tcarac, pip);
}

/*----------------------------------------------------------------------------
 * Integration of SDEs for temperature seen by particles.
 *----------------------------------------------------------------------------*/

static void
_lagitf(const cs_lnum_t       npt,
        const cs_real_t       dt_part,
        int                   nor,
        cs_lagr_attribute_t  *iattr)
{
  const cs_real_3_t *cell_cen = (cs_real_3_t*)cs_glob_mesh_quantities->cell_cen;
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;
  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);

  /* Initialize variables to avoid compiler warnings */

  int ltsvar = 0;

  if (p_set->p_am->source_term_displ != nullptr) {
    if (p_set->p_am->source_term_displ[*iattr] >= 0)
      ltsvar   = 1;
  }

  /* Mean temperature seen in degrees C
   * =================================== */
  cs_real_t t_shift = 0.;

  if (    extra->temperature != nullptr
       && cs_glob_thermal_model->temperature_scale
              == CS_TEMPERATURE_SCALE_KELVIN)
    t_shift = cs_physical_constants_celsius_to_kelvin;

  else if (extra->temperature == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("%s (Lagrangian module):\n\n"
                "Temperature field is not defined or mapped."), __func__);

  /* FIXME take value at previous time step if nor == 1 and current otherwise*/
  cs_real_t loc_tempf = extra->temperature->val[cell_id] - t_shift;
  if (cs_glob_lagr_time_scheme->interpol_field > 0
      && extra->grad_tempf != nullptr) {
    cs_real_t *old_part_coords =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                 CS_LAGR_COORDS);
    /* linear interpolation */
    for (int i = 0; i < 3; i++)
      loc_tempf += extra->grad_tempf[cell_id][i]
        * (old_part_coords[i] - cell_cen[cell_id][i]);
  }

  cs_real_t daux1 =   (0.5 + 0.75 * cs_turb_crij_c0) / cs_turb_crij_ct
                    * extra->lagr_time->val[cell_id] / dt_part;
  cs_real_t aux2 = 0.;
  if (daux1 > cs_math_epzero)
    aux2 = exp(-1. / daux1);

  /* Integration of the SDE over particles
   * ===================================== */

  if (nor == 1) {

    cs_real_t ter1 = aux2 *
      cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_TEMPERATURE_SEEN);
    cs_real_t ter2 = loc_tempf * (1.0 - aux2);

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE_SEEN,
                              ter1 + ter2);

    /* Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2  */
    if (ltsvar) {
      cs_real_t *part_ts_fluid_t
        = cs_lagr_particles_source_terms(p_set, npt,
                                         CS_LAGR_TEMPERATURE_SEEN);
      *part_ts_fluid_t = 0.5 * ter1 - loc_tempf
                             * (aux2 + (aux2 - 1.0) * daux1);
    }

  }
  else if (nor == 2) {

  /* FIXME check this: not called for fixed particle
   * yet a fixed particle's mass might still evolve */

    if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_REBOUND_ID) != 0 ) {

      cs_real_t ter1
        = 0.5 * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                            CS_LAGR_TEMPERATURE_SEEN) * aux2;
      cs_real_t ter2   = loc_tempf * (1.0 + (aux2 - 1.0) * daux1);
      cs_real_t *part_ts_fluid_t
        = cs_lagr_particles_source_terms(p_set, npt, CS_LAGR_TEMPERATURE_SEEN);

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE_SEEN,
                                *part_ts_fluid_t + ter1 + ter2);
    }
  }
}

/*----------------------------------------------------------------------------
 * Integrate SDE's for coal and compute shrinking coal diameter
 *
 * parameters:
 *   npt      <--  particle id
 *   dt_part <--  time step associated to the particle
 *   nor     <--  current step id (for 2nd order scheme)
 *   tempct  <-- thermal characteristic times
 *   cpgd1   <-> devolatilization term 1
 *   cpgd2   <-> devolatilization term 1
 *   cpght   <-> term for heterogeneous combustion (coal with thermal
 *               return coupling)
 *----------------------------------------------------------------------------*/

static void
_lagich(const cs_lnum_t       npt,
        const cs_real_t       dt_part,
        int                   nor,
        const cs_real_2_t     tempct,
        cs_real_t            *cpgd1,
        cs_real_t            *cpgd2,
        cs_real_t            *cpght)
{
  /* Particles management */

  cs_lagr_particle_set_t        *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t *p_am  = p_set->p_am;
  const cs_coal_model_t *coal_model = cs_glob_coal_model;

  const int ico = coal_model->ico - 1;
  const int io2 = coal_model->io2 - 1;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  const cs_real_t xcal2j = 4.1855e0;

  /* Small value (for numeric precision) */
  cs_real_t precis = 1e-15;

  /* Latent heat in J/kg */
  cs_real_t lv = 2263000.0;

  /* Check for models activation */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("In Lagrangian module:\n\n"
                "Lagrangian transport of coal particles is activated\n"
                "but the appropriate coal combustion model is not."));

  /* Numerical variables */
  cs_real_t d6spi  = 6.0 / cs_math_pi;
  cs_real_t dpis6  = cs_math_pi / 6.0;
  cs_real_t d1s3   = 1.0 / 3.0;
  cs_real_t d2s3   = 2.0 / 3.0;

  /* File variables */
  cs_real_t coef  = 0.0;

  /* If thermal return coupling  */

  if (cs_glob_lagr_source_terms->ltsthe == 1) {

    coef = 1.0 / cs_glob_lagr_time_scheme->t_order;

    if (nor == 1) {
        *cpgd1 = 0.0;
        *cpgd2 = 0.0;
        *cpght = 0.0;
    }

  }

  if (nor == 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian coal/shrinking diameter model is incompatible\n"
                "with second-order model."));

  cs_real_t diftl0 = -1;
  diftl0 = cs_field_get_key_double(cs_field_by_name("enthalpy"),
                                   cs_field_key_id("diffusivity_ref"));

  /* Main loop on particles
   * ====================== */

  cs_lnum_t nlayer = cs_glob_lagr_const_dim->nlayer;
  cs_real_t f_nlayer = nlayer;

  /* FIXME check this: not called for fixed particle
   * yet a fixed particle's mass might still evolve */

  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);

  /* local variables*/
  cs_real_t aux1, aux2, aux3, aux4, aux5;

  /* Variables generiques */
  cs_real_t diam           = cs_lagr_particle_get_real(particle, p_am,
                                                       CS_LAGR_DIAMETER);
  cs_real_t init_diam      = cs_lagr_particle_get_real(particle, p_am,
                                                       CS_LAGR_INITIAL_DIAMETER);
  cs_real_t shrink_diam    = cs_lagr_particle_get_real(particle, p_am,
                                                       CS_LAGR_SHRINKING_DIAMETER);

  cs_real_t *part_vel_seen =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY_SEEN);
  cs_real_t *part_vel      =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY);

  cs_real_t *part_temp     =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_TEMPERATURE);

  cs_real_t *part_coke_mass      =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                               0, CS_LAGR_COKE_MASS);
  cs_real_t *prev_part_coke_mass =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                               1, CS_LAGR_COKE_MASS);

  cs_real_t *part_coal_mass      =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                               0, CS_LAGR_COAL_MASS);
  cs_real_t *prev_part_coal_mass =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                               1, CS_LAGR_COAL_MASS);

  cs_real_t *part_coal_density   =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_COAL_DENSITY);

  cs_lnum_t co_id = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_COAL_ID);

  cs_real_t layer_vol  = dpis6 * _pow3(init_diam) / f_nlayer;

  /* Reynolds number */
  aux1 = cs_math_3_distance(part_vel_seen, part_vel);

  cs_real_t rom  = extra->cromf->val[cell_id];
  cs_real_t xnul = extra->viscl->val[cell_id] / rom;
  cs_real_t rep  = aux1 * diam / xnul;

  /* Prandtl and Sherwood numbers */
  cs_real_t xrkl;
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 2)
    xrkl = diftl0 / rom;
  else if (extra->cpro_viscls != nullptr )
    xrkl = extra->cpro_viscls->val[cell_id] / rom;
  else
    xrkl = extra->visls0 / rom;

  cs_real_t prt   = xnul / xrkl;
  cs_real_t sherw = 2 + 0.55 * pow (rep, 0.5) * pow (prt, (d1s3));

  /* Compute discretization radii */
  cs_real_t radius[nlayer];
  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {
    cs_real_t f_l = l_id+1;
    radius[l_id] = (init_diam/2.0) * pow(f_l/f_nlayer, d1s3);
  }

  cs_real_t mp0  = dpis6 * _pow3(init_diam) * coal_model->rho0ch[co_id];
  cs_real_t mwat_max  = coal_model->xwatch[co_id] * mp0 / nlayer;

  /* Compute water quantity on each layer */
  aux1 = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_WATER_MASS);

  cs_real_t mwater[nlayer];

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {

    if (l_id == nlayer - 1)
      mwater[l_id] = cs::max(0.0, aux1);

    else
      mwater[l_id] = cs::max(0.0, cs::min(aux1, mwat_max));

    aux1 -= mwater[l_id];

  }

  /* Mass on each layer */
  cs_real_t mlayer[nlayer];
  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    mlayer[l_id] =   coal_model->xashch[co_id] * mp0 / nlayer
                   + mwater[l_id]
                   + part_coal_mass[l_id]
                   + part_coke_mass[l_id];

  /* Compute avaporating water mass; we assume for the computation
     of active coal density that drying occurs at constant volume */

  /* Vapor flux for particle */

  cs_real_t fwat[nlayer];
  _lagsec(npt, dt_part, nor, diftl0, layer_vol, mwat_max, sherw, radius, tempct,
          mlayer, mwater, fwat);

  /* Compute velocity constants SPK1 of SPK2 of the mass transfer by
     devolatilization with the Arrhenius laws */

  cs_real_t skp1[nlayer], skp2[nlayer];

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {

    aux1  = 1.0 / (cs_physical_constants_r * part_temp[l_id]);

    skp1[l_id] =   coal_model->a1ch[co_id]
                 * exp(-coal_model->e1ch[co_id] * aux1);
    skp2[l_id] =   coal_model->a2ch[co_id]
                 * exp(-coal_model->e2ch[co_id] * aux1);

    aux1  = skp1[l_id] * coal_model->y1ch[co_id] * part_coal_mass[l_id];
    aux2  = skp2[l_id] * coal_model->y2ch[co_id] * part_coal_mass[l_id];

    /* Thermal return coupling */

    if (cs_glob_lagr_source_terms->ltsthe == 1) {
      *cpgd1 = *cpgd1 + coef * aux1;
      *cpgd2 = *cpgd2 + coef * aux2;
    }

  }

  /* Compute global heterogeneous diffusion combustion constant
     --------------------------_------------------------------- */

  /* Locate layer where heterogeneous combustion occurs */

  cs_lnum_t l_id_het = 0;

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {
    if (prev_part_coal_mass[l_id] > 0.0)
      l_id_het = l_id;
  }

  /* Check if ck remains on a more external layer */
  for (cs_lnum_t l_id = l_id_het; l_id < nlayer; l_id++) {
    if (prev_part_coke_mass[l_id] > 0.0)
      l_id_het = l_id;
  }

  /* Chemical cinetics coeffcient for CO formation, in (kg.m-2.s-1.atm(-n))
     conversion (kcal/mol -> J/mol) */

  aux1 = coal_model->ehetch[co_id] * 1000.0 * xcal2j;
  aux2 = coal_model->ahetch[co_id] * exp (-aux1 / (  cs_physical_constants_r
                                                   * part_temp[l_id_het]));

  /* Diffusion coefficient in (Kg/m2/s/atm) and global reaction constant */

  cs_real_t skglob;
  if (cs_physical_constants_r * shrink_diam > precis) {

    /* Constant 2.53e-7 is explained in tome 5 of report on code_saturne
       specific physics (HI-81/04/003/A) equation 80 */
    aux3 = sherw * 2.53e-07 * (pow (extra->temperature->val[cell_id], 0.75))
                            / shrink_diam;
    skglob = (aux2 * aux3) / (aux2 + aux3);

  }
  else
    skglob = aux2;

  /* Compute GAMMAhet
     ---------------- */

  /* Compute O2 partial pressure
   *
   *     PO2 = RHO1*cs_physical_constants_r*T*YO2/MO2
   *                                                      */
  aux1 =   extra->cromf->val[cell_id] * cs_physical_constants_r
         * extra->temperature->val[cell_id]
         * extra->x_oxyd->val[cell_id] / coal_model->wmole[coal_model->io2]
           / cs_coal_prefth;

  /* Compute working surface: SE */
  aux2 =  cs_math_pi * (1.0 - coal_model->xashch[co_id]) * pow(shrink_diam, 2);

  /* No heterogeneous combustion if Mch/Mp >= 1.e-3 */
  cs_real_t gamhet;
  if (prev_part_coal_mass[0] <= (0.001 * mlayer[0]))
    gamhet = 0.0;
  else
    /* Compute GamHET */
    gamhet = aux1 * aux2 * skglob;

  /* Thermal return coupling */
  if (cs_glob_lagr_source_terms->ltsthe == 1)
    *cpght = *cpght + coef * gamhet;

  /* Compute 0.5(MO2/MC)*(HO2(Tp)-HO2(TF))
   * ------------------------------------- */

  /* Compute Hc(Tp)-Mco/Mc Hco2(Tp)+0.5Mo2/Mc Ho2(Tf) */

  cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

  /* Compute Hcoke(TP) */
  aux1  =    coal_model->h02ch[co_id]
          +   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP)
            * (part_temp[l_id_het] - cs_coal_trefth);

  /* Compute MCO/MC HCO(TP)  */
  for (cs_lnum_t iii = 0;
       iii < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
       iii++)
    coefe[iii] = 0.0;

  coefe[ico] =   coal_model->wmole[ico]
               / coal_model->wmolat[cs_coal_atom_id_c];

  aux2 = cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(part_temp[l_id_het],
                                                         coefe);

  /* Compute MO2/MC/2. HO2(TF)    */
  for (cs_lnum_t iii = 0;
       iii < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
       iii++)
    coefe[iii]   = 0.0;

  coefe[io2] =   coal_model->wmole[io2]
               / coal_model->wmolat[cs_coal_atom_id_c] / 2.0;

  aux3 = cs_lagr_particle_get_real(particle, p_am,
                                   CS_LAGR_TEMPERATURE_SEEN) + _tkelvi;

  aux4 = cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(aux3, coefe);

  cs_real_t deltah = aux2 - aux4 - aux1;

  /* Integaration of water mass
     -------------------------- */

  if (nor == 1) {

    aux1   = 0.0;

    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
      aux1 += fwat[l_id] * dt_part;

    cs_real_t mwat = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_WATER_MASS) - aux1;

    /* Clipping */
    if (mwat < precis)
      mwat = 0.0;

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS, mwat);

  }
  else if (nor == 2) {

    aux1   = 0.0;
    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
      aux1 += fwat[l_id] * dt_part;

    cs_real_t mwat
      = 0.5 * (  cs_lagr_particle_get_real_n(particle, p_am,
                                             0, CS_LAGR_WATER_MASS)
               + cs_lagr_particle_get_real_n(particle, p_am,
                                             1, CS_LAGR_WATER_MASS)
               - aux1);

    /* Clipping */
    if (mwat < precis)
      mwat = 0.0;

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS, mwat);

  }

  /* Integration of reactive coal mass
     --------------------------------- */

  if (nor == 1) {

    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {

      aux1 = exp(-(skp1[l_id] + skp2[l_id]) * dt_part);

      part_coal_mass[l_id] = prev_part_coal_mass[l_id] * aux1;

      /* Clipping */
      if (part_coal_mass[l_id] < precis)
        part_coal_mass[l_id] = 0.0;

    }

  }
  else if (nor == 2) {

    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {

      aux1 = exp ( -(skp1[l_id] + skp2[l_id]) * dt_part);

      part_coal_mass[l_id] = 0.5 * (  part_coal_mass[l_id]
                                    + prev_part_coal_mass[l_id] * aux1);

      /* Clipping   */
      if (part_coal_mass[l_id] < precis)
        part_coal_mass[l_id] = 0.0;

    }

  }

  /* Integration of coke mass
     ------------------------ */

  cs_real_t fcoke[nlayer];

  if (nor == 1) {

    /* Initialize efective flux of heterogeneous combustion */
    for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
      fcoke[l_id]  = 0.0;

    /* Loop on all cells which have reactive coke or coal */
    for (cs_lnum_t l_id = 0; l_id < l_id_het; l_id++) {

      aux1 =  (skp1[l_id] * (1.0 - coal_model->y1ch[co_id]) + skp2[l_id]
                          * (1.0 - coal_model->y2ch[co_id]))
            / (skp1[l_id] + skp2[l_id]);
      aux2 = exp(-(skp1[l_id] + skp2[l_id]) * dt_part);
      aux3 = aux1 * prev_part_coal_mass[l_id] * (1.0 - aux2) / dt_part;

      if (l_id == l_id_het) {

        /* Compute equivalent coke mass */
        aux4 =  dpis6 * (1.00 - coal_model->xashch[co_id])
              * _pow3(shrink_diam)
              * part_coal_density[l_id];

        if (aux4 > precis) {

          /* Account for heterogeneous combustion */
          aux5 = dt_part * aux4 * (-gamhet + aux3)
               / (d2s3 * gamhet * dt_part + aux4);
          fcoke[l_id] = gamhet;

        }
        else
          /* Ignore heterogeneous combustion */
          aux5 = dt_part * aux3;

      }
      else
        /* Ignore heterogeneous combustion */
        aux5 = dt_part * aux3;

      part_coke_mass[l_id] = prev_part_coke_mass[l_id] + aux5;

    }

    /* If gamhet is too large, spread over several layers */
    for (cs_lnum_t l_id = l_id_het-1; l_id >= 0; l_id--) {

      if (part_coke_mass[l_id] < 0) {

        /* Limit heterogeneous combustion */
        fcoke[l_id] += part_coke_mass[l_id];

        /* Possibly start combustion of next layer */
        if (l_id > 1 )
          part_coke_mass[l_id - 1]
            = part_coke_mass[l_id - 1] + part_coke_mass[l_id];

        /* Limit coke mass */
        part_coke_mass[l_id] = 0.0;

      }

    }

  }
  else if (nor == 2) {
    /* No second order for now */
    assert(0);
  }

  /* Integrate temperature of coal grains */

  cs_real_t  phith[nlayer];
  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    /* Thermal source terms layer by layer:
       thermal exchanges with the exterior are computed directly by _lagtmp */
    phith[l_id] = (-fcoke[l_id] * deltah) - fwat[l_id] * lv;

  cs_real_t temp[nlayer];
  _lagtmp(npt, dt_part, nor,  layer_vol, tempct, radius, mlayer, phith, temp);

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    part_temp[l_id]  = temp[l_id];

  /* Update coke density */

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {

    if (prev_part_coal_mass[l_id] >= 0.001 * mlayer[l_id]) {

      /* mv represents mlayer which left the grain (drying + pyrolysis) */
      cs_real_t mv =  mp0 * (1 - coal_model->xashch[co_id]) / nlayer
                    - part_coal_mass[l_id] - part_coke_mass[l_id]
                    - (mwater[l_id] - fwat[l_id] * dt_part);

      /* density of coke only */
      part_coal_density[l_id]
        =   coal_model->rho0ch[co_id] - mv
          / (layer_vol * (1.0 - coal_model->xashch[co_id]));

    }

  }

  /* Update diameter of shrinking core
     --------------------------------- */

  /* Locate most external layer with ch */

  l_id_het = 0;
  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++) {
    if (part_coal_mass[l_id] > 0.0)
      l_id_het = l_id;
  }

  /* Check for remaining coke on a more external layer */

  for (cs_lnum_t l_id = l_id_het; l_id < nlayer; l_id++) {
    if (part_coke_mass[l_id] > 0.0)
      l_id_het = l_id;
  }

  if (part_coal_mass[l_id_het] >= 0.001 * mlayer[l_id_het]) {
    /* Pyrolysis is not finished, char has initial diameter */
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_SHRINKING_DIAMETER,
                              2.0 * radius[l_id_het]);
  }
  else {

    /* Distribute char in uniform manner */
    if (l_id_het == 0) {

      aux5 = pow(  d6spi
                 / (1.0 - coal_model->xashch[co_id])
                 * (  part_coal_mass[l_id_het] / coal_model->rho0ch[co_id]
                    + part_coke_mass[l_id_het] / part_coal_density[l_id_het]),
                   d1s3);

      /* Clipping   */
      if (aux5 > 2.0 * radius[l_id_het])
        aux5 = 2.0 * radius[l_id_het];

      else if (aux5 < 0.0)
        aux5 = 0.0;

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_SHRINKING_DIAMETER,
                                aux5);

    }
    else {

      cs_real_t f1 = part_coal_mass[l_id_het] / coal_model->rho0ch[co_id];
      cs_real_t f2 = 0;
      if (part_coal_density[l_id_het] > 0.0)
        f2 = part_coke_mass[l_id_het] / part_coal_density[l_id_het];
      aux5 = pow(  cs_math_pow3(2.0 * radius[l_id_het - 1])
                 + (  d6spi / (1.0 - coal_model->xashch[co_id])
                    * (f1 + f2)), d1s3);

      /* Clipping   */
      if (aux5 > 2.0 * radius[l_id_het])
        aux5 = 2.0 * radius[l_id_het];

      else if (aux5 < 2.0 * radius[l_id_het - 1])
        aux5 = 2.0 * radius[l_id_het - 1];

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_SHRINKING_DIAMETER,
                                aux5);

    }

  }

  shrink_diam = cs_lagr_particle_get_real(particle, p_am,
                                          CS_LAGR_SHRINKING_DIAMETER);

  /* Compute diameter of coal grains
   * ------------------------------- */

  aux5 = sqrt(         coal_model->xashch[co_id]  * cs_math_pow2(init_diam)
              + (1.0 - coal_model->xashch[co_id]) * cs_math_pow2(shrink_diam));

  /* Compute mass of coal grains
   * --------------------------- */

  aux1 = 0.0;

  for (cs_lnum_t l_id = 0; l_id < nlayer; l_id++)
    aux1 += part_coal_mass[l_id] + part_coke_mass[l_id];

  cs_real_t mwat = cs_lagr_particle_get_real(particle, p_am,
                                             CS_LAGR_WATER_MASS);

  aux1 += mwat + coal_model->xashch[co_id] * mp0;

  cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, aux1);

}

/*----------------------------------------------------------------------------
 * Integrate SDE's for cooling tower model
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Compute the Lewis factor used for the evaluation of the heat transfer
 * phase change source terms
 *
 * parameters:
 *   evap_model  <-- Evaporation model: CS_CTWR_POPPE or CS_CTWR_MERKEL
 *   molmassrat  <-- Dry air to water vapor molecular mass ratio
 *   x           <-- Humidity
 *   x_s_tl      <-- Saturation humidity at the temperature of the liquid
 *
 * returns:
 *   xlew        --> Lewis factor
 *----------------------------------------------------------------------------*/

static cs_real_t
_lewis_factor(const int        evap_model,
              const cs_real_t  molmassrat,
              const cs_real_t  x,
              const cs_real_t  x_s_tl)
{
  /* Merkel Model
     Hypothesis of unity Lewis factor */
  cs_real_t xlew = 1.;

  if (evap_model == CS_CTWR_POPPE) {
    /* Poppe evaporation model
       Compute Lewis factor using Bosnjakovic hypothesis
       NB: clippings ensuring xi > 1 and xlew > 0 */
    cs_real_t xi = (molmassrat + x_s_tl)/(molmassrat + cs::min(x, x_s_tl));
    if ((xi - 1.) < 1.e-15)
      xlew = pow(0.866,(2./3.));
    else
      xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);
  }

  return xlew;
}

/*----------------------------------------------------------------------------
 * Compute the heat and mass tranfer due to evaporation of a single droplet.
 *
 * The mathematical models used to describe these
 * phenomena are similar to those used for the rain zone of the cooling tower
 * model (see the cs_ctwr.c subroutine).
 *
 * parameters:
 *   npt         <--  particle id
 *   dt_part     <--  time step associated to the particle
 *   nor         <--  current step id (for 2nd order scheme)
 *
 *----------------------------------------------------------------------------*/

static void
_sde_i_ct(const cs_lnum_t       npt,
          const cs_real_t       dt_part,
          int                   nor)
{
  /* Adressing structures of the Lagrangian module */
  cs_lagr_particle_set_t        *p_set  = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t *p_am   = p_set->p_am;
  cs_lagr_extra_module_t        *extra  = cs_glob_lagr_extra_module;

  /* Adressing structures of the cooling tower module */
  cs_ctwr_option_t              *ct_opt = cs_get_glob_ctwr_option();
  int                            evap_model= ct_opt->evap_model;

  /* Adressing field values of the Eulerian cooling tower model */
  cs_real_t *x   = cs_field_by_name("humidity")->val;
  cs_real_t *x_s = cs_field_by_name("x_s")->val;

  /* User defined air properties of the cooling tower model*/
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;
  cs_real_t lambda_h             = air_prop->lambda_h;
  cs_real_t hv0                  = air_prop->hv0;
  cs_real_t rho_l                = air_prop->rho_l;

  /* General fluid properties*/
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  cs_real_t p0                       = fluid_props->p0;

  /* Constant local numerical parameters*/
  cs_real_t molmassrat = air_prop->molmass_rat;
  cs_real_t precis = 1e-15; /* clipping purposes */
  // cs_real_t r = 461.5; /* J /kg K */

  unsigned char *particle = p_set->p_buffer + p_am->extents * npt;
  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);
  /* Calculating the current particle surface */
  cs_real_t dia = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_DIAMETER);
  //cs_real_t surf_p = cs_math_pi*cs_math_pow2(dia) * 0.25;

  /* Drop diameter based Reynolds number */
  const cs_real_t *vel_p
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY);
  const cs_real_t *vel_f
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY_SEEN);
  cs_real_t rho_h = extra->cromf->val[cell_id];
  cs_real_t visc  = extra->viscl->val[cell_id];

  /* local Reynolds number */
  cs_real_t rey   = rho_h * cs_math_3_distance(vel_p, vel_f) * dia / visc;

  /* Prandtl number */
  cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);
  cs_real_t pr   = cp_h * visc / lambda_h;

  cs_real_t t_l_p
    = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_TEMPERATURE);
  cs_real_t temp_h
    = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_TEMPERATURE_SEEN);
  cs_real_t x_s_tl = cs_air_x_sat(t_l_p, p0);
  cs_real_t x_s_th = cs_air_x_sat(temp_h, p0);
  cs_real_t xlew   = _lewis_factor(evap_model, molmassrat,
                                 x[cell_id], x_s_tl);

  /* Nusselt number correlations */
  /* Ranz-Marshall or Hughmark when rey <= 776.06 && pr <= 250. */
  cs_real_t nusselt = 2.+0.6*sqrt(rey)*pow(pr,(1./3.));
  /* Hughmark when rey > 776.06 && pr <= 250. */
  if (rey > 776.06 && pr <= 250.) {
    nusselt = 2. + 0.27*pow(rey, 0.62)*pow(pr,(1./3.));
  }

  cs_real_t a_c = (nusselt * lambda_h) / dia;
  cs_real_t beta_x = a_c / (xlew * cp_h);
  cs_real_t beta_x_ai = beta_x * cs_math_pi*cs_math_pow2(dia);
  cs_real_t mass_source = 0;
  if (x[cell_id] <= x_s_th) {
    mass_source = beta_x_ai * (x_s_tl - x[cell_id]);
  }
  else {
    mass_source = beta_x_ai * (x_s_tl - x_s_th);
  }

  /* Integration of droplet mass
     -------------------------- */

  cs_real_t aux1 = 0.0;
  aux1 +=  mass_source * dt_part;

  if (nor == 1) {

    cs_real_t mass_p = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                   CS_LAGR_MASS) - aux1;

    /* Clipping */
    if (mass_p < precis)
      mass_p = 0.0;

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mass_p);

  }
  else if (nor == 2) {

    cs_real_t mass_p = 0.5 * (  cs_lagr_particle_get_real_n(particle, p_am,
                                                          0, CS_LAGR_MASS)
                            + cs_lagr_particle_get_real_n(particle, p_am,
                                                            1, CS_LAGR_MASS)
                            - aux1);

    /* Clipping */
    if (mass_p < precis)
      mass_p = 0.0;

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mass_p) ;

  }

  /* Droplet diameter */

  cs_real_t mass_p = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_MASS) - aux1;
  cs_real_t droplet_dia = pow(6*mass_p/(cs_math_pi*rho_l),1./3);
  cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, droplet_dia);

  /* Integration of droplet enthalpy
     -------------------------- */
  /* Energy loss due to evaporation */
  aux1 = 0.0;
  aux1 += mass_source * hv0 * dt_part; /* */

  /* Energy loss due to heat conduction */
  cs_real_t aux2 =   (lambda_h * nusselt / dia)
                    * cs_math_pi * cs_math_pow2(dia) *dt_part;
  //aux2 += a_c * surf_p * dtp;

  if (nor == 1) {

    /* Evaporation term */
    cs_real_t temp_p = (cs_lagr_particle_get_real_n(particle, p_am, 1,
                                         CS_LAGR_TEMPERATURE) + _tkelvi) *
                        cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                    CS_LAGR_CP)          *
                        cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                    CS_LAGR_MASS) - aux1;
    /* Heat conduction term */
    /* temperature difference in C is the same as for K */
    temp_p -= aux2 * ( cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                   CS_LAGR_TEMPERATURE)  -
                       cs_lagr_particle_get_real_n(particle, p_am, 1,
                                            CS_LAGR_TEMPERATURE_SEEN) );

    /* FIXME: Shall we apply a sort of clipping to the enthalpy before
       converting it into a temperature value ? */
    temp_p /= ( cs_lagr_particle_get_real_n(particle, p_am,1,CS_LAGR_CP) *
        cs_lagr_particle_get_real_n(particle, p_am,1,CS_LAGR_MASS) );

    temp_p -= _tkelvi;
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE, temp_p);

  }
  else if (nor == 2) {

    /* Evaporation term */
    cs_real_t temp_p
      = 0.5 * (   (  (cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                  CS_LAGR_TEMPERATURE) + _tkelvi)
                   * cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                 CS_LAGR_CP)
                   * cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                 CS_LAGR_MASS))
               +
                  (  (cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                  CS_LAGR_TEMPERATURE) + _tkelvi)
                   * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_CP)
                   * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_MASS))) - aux1;

    /* Heat conduction term */
    /* temperature difference in C is the same as for K */
    temp_p
      += aux2 * (  0.5 * (  cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                        CS_LAGR_TEMPERATURE)
                          + cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                        CS_LAGR_TEMPERATURE))
                 - 0.5 * (  cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                      CS_LAGR_TEMPERATURE_SEEN)
                          + cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                        CS_LAGR_TEMPERATURE_SEEN)));

    /* FIXME: Shall we apply a sort of clipping to the enthalpy before
       converting it into a temperature value ? */
    temp_p
      /=  0.5 * (  (   cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                CS_LAGR_CP)
                     * cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                CS_LAGR_MASS) )
                 +  (  cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_CP)
                     * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                 CS_LAGR_MASS)));
    temp_p -= _tkelvi;
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, temp_p) ;

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle stochastic differential equations
 *        for specific physical models.
 *
 * - temperature seen by particles,
 * - particle temperature,
 * - particle diameter
 * - particle mass
 * - variables related to coal grains (Temp, MCH, MCK)
 * - additional user parameters
 *
 * \param[in]  npt           particle id
 * \param[in]  dt_part      time step associated to the particle
 * \param[in]  nor          current step id (for 2nd order scheme)
 * \param[in]  tempct       thermal characteristic time
 * \param[out] cpgd1,cpgd2  devolatilisation terms 1 and 2
 * \param[out] cpght        heterogeneos combusion terms (coal with thermal
 *                          return coupling)
 */
/*------------------------------------------------------------------------- */

void
cs_lagr_sde_model(const cs_lnum_t    npt,
                  const cs_real_t    dt_part,
                  int                nor,
                  const cs_real_2_t  tempct,
                  cs_real_t         *cpgd1,
                  cs_real_t         *cpgd2,
                  cs_real_t         *cpght)
{
  cs_lagr_attribute_t temp_seen = CS_LAGR_TEMPERATURE_SEEN;

  /* Integration of temperature seen by particles */

  if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
      || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR
      || (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
          && cs_glob_lagr_specific_physics->solve_temperature_seen == 1))
    _lagitf(npt, dt_part, nor, &temp_seen);

  /* Integration of particles temperature */

  if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
      && cs_glob_lagr_specific_physics->solve_temperature == 1)
    _lagitp(npt, dt_part, nor, tempct);

  /* Integration of particles diameter */

  if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
      && cs_glob_lagr_specific_physics->solve_diameter == 1)
    _lagidp(npt, dt_part, nor);

  /* Integration of particles mass */

  if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
      && cs_glob_lagr_specific_physics->solve_mass == 1)
    _lagimp(npt, dt_part, nor);

  /* Integration of coal equations: hp, mch, mck */

  if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL)
    _lagich(npt , dt_part, nor, tempct, cpgd1, cpgd2, cpght);

  /* Integration of cooling tower model equations*/

  if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR)
    _sde_i_ct(npt, dt_part, nor);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
