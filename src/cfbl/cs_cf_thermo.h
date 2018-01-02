#ifndef __CS_CF_THERMO_H__
#define __CS_CF_THERMO_H__

/*============================================================================
 * Thermodynamic laws for the compressible module
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute \f$\gamma\f$.
 *
 * \param[in]     cp      array of isobaric specific heat values
 * \param[in]     cv      array of isochoric specific heat values
 * \param[out]    gamma   array of values of ratio of specific heat
 * \param[in]     l_size  size of the array
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_cf_thermo_gamma(cs_real_t *cp,
                   cs_real_t *cv,
                   cs_real_t *gamma,
                   cs_lnum_t l_size)
{
  /*  Local variables */
  int ieos = cs_glob_fluid_properties->ieos;

  /*  Gamma is supposed to be superior or equal to 1.
      It is computed at each call, even if this may seem costly,
      to be coherent with the "constant gamma" case for which this
      constant is not saved. */

  /* single ideal gas - constant gamma
     or ideal gas mix - gamma for the mixture */
  if (ieos == 1 || ieos == 3) {
    for (cs_lnum_t ii = 0; ii < l_size; ii++) {
      gamma[ii] = cp[ii]/cv[ii];
      if (gamma[ii] < 1.)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error in thermodynamics computations for "
                    "compressible flows:\n"
                    "Value of gamma smaller to 1. encountered.\n"
                    "Gamma (specific heat ratio) must be a real number "
                    "greater or equal to 1.\n"));
    }
  }
  /* stiffened gas - constant gamma (parameter of the law) */
  else if (ieos == 2) {
    for (cs_lnum_t ii = 0; ii < l_size; ii++)
      gamma[ii] = cs_glob_fluid_properties->gammasg;
  }
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set variability of isobaric specific heat and isochoric specific heat
 * according to the chosen thermodynamic law.
 *----------------------------------------------------------------------------*/

void
cs_cf_set_thermo_options(void);

/*----------------------------------------------------------------------------
 * Initialize density, total energy and isochoric specific heat
 * according to the chosen thermodynamic law using the default parameters.
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_default_init(void);

// TODO: the check function should be generalized (pass the name as argument).

/*----------------------------------------------------------------------------
 * Check the positivity of the pressure.
 *
 * parameters:
 *   pres   <-- array of pressure values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_check_pressure(cs_real_t *pres,
                     cs_lnum_t l_size);

/*----------------------------------------------------------------------------
 * Check the positivity of the internal energy.
 *
 * parameters:
 *   ener   <-- array of total energy values
 *   l_size <-- l_size of the array
 *   vel    <-- array of velocity values
 *----------------------------------------------------------------------------*/

void
cs_cf_check_internal_energy(cs_real_t   *ener,
                            cs_lnum_t    l_size,
                            cs_real_3_t *vel);

/*----------------------------------------------------------------------------
 * Check the positivity of the density given by the user.
 *
 * parameters:
 *   dens   <-- array of density values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_check_density(cs_real_t *dens,
                    cs_lnum_t l_size);

/*----------------------------------------------------------------------------
 * Check strict positivity of temperature (Celsius) given by the user.
 *
 * parameters:
 *   temp   <-- array of temperature values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_check_temperature(cs_real_t *temp,
                        cs_lnum_t l_size);

/*----------------------------------------------------------------------------
 * Compute temperature and total energy from density and pressure.
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   pres   <-- array of pressure values
 *   dens   <-- array of density values
 *   temp   --> array of temperature values
 *   ener   --> array of total energy values
 *   vel    <-- array of velocity component values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_te_from_dp(cs_real_t   *cp,
                        cs_real_t   *cv,
                        cs_real_t   *pres,
                        cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size);

/*----------------------------------------------------------------------------
 * Compute density and total energy from pressure and temperature
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   pres   <-- array of pressure values
 *   temp   <-- array of temperature values
 *   dens   --> array of density values
 *   ener   --> array of total energy values
 *   vel    <-- array of velocity component values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_de_from_pt(cs_real_t   *cp,
                        cs_real_t   *cv,
                        cs_real_t   *pres,
                        cs_real_t   *temp,
                        cs_real_t   *dens,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size);

/*----------------------------------------------------------------------------
 * Compute density and temperature from pressure and total energy.
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   pres   <-- array of pressure values
 *   ener   <-- array of total energy values
 *   dens   --> array of density values
 *   temp   --> array of temperature values
 *   vel    <-- array of velocity component values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_dt_from_pe(cs_real_t   *cp,
                        cs_real_t   *cv,
                        cs_real_t   *pres,
                        cs_real_t   *ener,
                        cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size);

/*----------------------------------------------------------------------------
 * Compute pressure and total energy from density and temperature
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   dens   <-- array of density values
 *   temp   <-- array of temperature values
 *   pres   --> array of pressure values
 *   ener   --> array of total energy values
 *   vel    <-- array of velocity component values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_pe_from_dt(cs_real_t   *cp,
                        cs_real_t   *cv,
                        cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_t   *pres,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size);

/*----------------------------------------------------------------------------
 * Compute pressure and temperature from density and total energy.
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   dens   <-- array of density values
 *   ener   <-- array of total energy values
 *   pres   --> array of pressure values
 *   temp   --> array of temperature values
 *   vel    <-- array of velocity component values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_pt_from_de(cs_real_t   *cp,
                        cs_real_t   *cv,
                        cs_real_t   *dens,
                        cs_real_t   *ener,
                        cs_real_t   *pres,
                        cs_real_t   *temp,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size);

/*----------------------------------------------------------------------------
 * Compute square of sound velocity for perfect gas.
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   pres   <-- array of pressure values
 *   dens   <-- array of density values
 *   c2     --> array of the values of the square of sound velocity
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_c_square(cs_real_t *cp,
                      cs_real_t *cv,
                      cs_real_t *pres,
                      cs_real_t *dens,
                      cs_real_t *c2,
                      cs_lnum_t  l_size);

/*----------------------------------------------------------------------------
 * Compute the thermal expansion coefficient for a perfect gas.
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   dens   <-- array of density values
 *   beta   --> array of beta values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_beta(cs_real_t *cp,
                  cs_real_t *cv,
                  cs_real_t *dens,
                  cs_real_t *beta,
                  cs_lnum_t  l_size);

/*----------------------------------------------------------------------------
 * Compute the isochoric specific heat:
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   xmasml <-- array of molar mass values
 *   cv     --> array of isochoric specific heat values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_cv(cs_real_t *cp,
                cs_real_t *xmasml,
                cs_real_t *cv,
                cs_lnum_t  l_size);

/*----------------------------------------------------------------------------
 * Compute entropy from pressure and density:
 *
 * parameters:
 *   cp     <-- array of isobaric specific heat values
 *   cv     <-- array of isochoric specific heat values
 *   dens   <-- array of density values
 *   pres   <-- array of pressure values
 *   entr   --> array of total energy values
 *   l_size <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_s_from_dp(cs_real_t *cp,
                       cs_real_t *cv,
                       cs_real_t *dens,
                       cs_real_t *pres,
                       cs_real_t *entr,
                       cs_lnum_t  l_size);

/*----------------------------------------------------------------------------
 * Compute wall boundary condition values.
 *
 * parameters:
 *   wbfa    --> output work array
 *   wbfb    --> output work array
 *   face_id <-- boundary face index
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_wall_bc(cs_real_t *wbfa,
                     cs_real_t *wbfb,
                     cs_lnum_t  face_id);

/*----------------------------------------------------------------------------
 * Compute subsonic outlet boundary conditions values.
 *
 * parameters:
 *   bc_en   <--> total energy values at boundary faces
 *   bc_pr   <--> pressure values at boundary faces
 *   bc_vel  <--> velocity values at boundary faces
 *   face_id <--  boundary face index
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_subsonic_outlet_bc(cs_real_t   *bc_en,
                                cs_real_t   *bc_pr,
                                cs_real_3_t *bc_vel,
                                cs_lnum_t    face_id);

/*----------------------------------------------------------------------------
 * Compute inlet boundary condition with total pressure and total
 * enthalpy imposed.
 *
 * parameters:
 *   bc_en   <--> total energy values at boundary faces
 *   bc_pr   <--> pressure values at boundary faces
 *   bc_vel  <--> velocity values at boundary faces
 *   face_id <--  boundary face number
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_ph_inlet_bc(cs_real_t   *bc_en,
                         cs_real_t   *bc_pr,
                         cs_real_3_t *bc_vel,
                         cs_lnum_t    face_id);

/*----------------------------------------------------------------------------
 * Compute epsilon sup for perfect gas.
 *
 * parameters:
 *   dens    <-- array of density values
 *   eps_sup --> epsilon sup array
 *   l_size  <-- l_size of the array
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo_eps_sup(cs_real_t *dens,
                     cs_real_t *eps_sup,
                     cs_lnum_t  l_size);

/*----------------------------------------------------------------------------
 * This subroutine is a driver allowing to call the appropriate
 * thermodynamical functions depending on the quantities provided by the user.
 * Hence it is only used during the initialization step and at the boundaries
 * of type supersonic inlet. It is described in the following how to select the
 * quantity to be returned.
 *
 * When calling the user subroutine, the integer 'iccfth' specifies which
 * calculation has to be performed (and which quantity has to be returned).
 * The values for 'iccfth' for each case are provided below.
 *
 *   The variables are referred to using a different index i:
 *
 *     - pressure: 2
 *     - density: 3
 *     - temperature: 5
 *     - internal energy: 7
 *     - entropy: 13
 *
 *   iccfth is as follows, depending on which quantity needs to be computed:
 *     - variables at cell centers from variable i and variable j (i<j):
 *           iccfth = i*j*10000
 *     - variables at boundary faces from variable i and variable j (i<j):
 *           iccfth = i*j*10000+900
 *
 * Detailed values of iccfth and corresponding computations:
 *
 *   Values at the cell centers:
 *
 *     - temperature and energy from pressure and density: iccfth =  60000
 *     - density and energy from pressure and temperature: iccfth =  100000
 *     - density and temperature from pressure and energy: iccfth =  140000
 *     - pressure and energy from density and temperature: iccfth =  150000
 *     - pressure and temperature from density and energy: iccfth =  210000
 *
 *   Values at the faces for boundary conditions:
 *     - temperature and energy from pressure and density: iccfth = 60900
 *     - density and energy from pressure and temperature: iccfth = 100900
 *     - density and temperature from pressure and energy: iccfth = 140900
 *     - pressure and energy from density and temperature: iccfth = 150900
 *     - pressure and temperature from density and energy: iccfth = 210900
 *
 * parameters:
 *   iccfth  --> id of computation
 *   face_id --> face index if the computation is for a B.C.
 *   bc_en   <-- total energy values at boundary faces
 *   bc_pr   <-- pressure values at boundary faces
 *   bc_tk   <-- temperature values at boundary faces
 *   bc_vel  <-- velocity values at boundary faces
 *----------------------------------------------------------------------------*/

void
cs_cf_thermo(const int    iccfth,
             cs_lnum_t    face_id,
             cs_real_t   *bc_en,
             cs_real_t   *bc_pr,
             cs_real_t   *bc_tk,
             cs_real_3_t *bc_vel);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CF_THERMO_H__ */
