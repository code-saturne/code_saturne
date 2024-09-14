#ifndef __CS_THERMAL_MODEL_H__
#define __CS_THERMAL_MODEL_H__

/*============================================================================
 * Base thermal model data.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Thermal model type
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_THERMAL_MODEL_INIT = -999,
  CS_THERMAL_MODEL_NONE = 0,
  CS_THERMAL_MODEL_TEMPERATURE,
  CS_THERMAL_MODEL_ENTHALPY,
  CS_THERMAL_MODEL_TOTAL_ENERGY,
  CS_THERMAL_MODEL_INTERNAL_ENERGY,
  CS_THERMAL_MODEL_N_TYPES

} cs_thermal_model_variable_t;

typedef enum {

  CS_TEMPERATURE_SCALE_NONE = 0,
  CS_TEMPERATURE_SCALE_KELVIN = 1,
  CS_TEMPERATURE_SCALE_CELSIUS = 2

} cs_temperature_scale_t;

/* thermal model descriptor */
/*--------------------------*/

typedef struct {

  union {
    cs_thermal_model_variable_t  thermal_variable;   /* Thermal variable */
    int                          itherm;
  };

  union {
    cs_temperature_scale_t       temperature_scale;  /* Temperature scale */
    int                          itpscl;
  };

  /* Has kinetic source terme correction */
  int           has_kinetic_st;
  int           cflt;           /* compute the thermal cfl condition */
  int           cflp;           /* compute the pressure cfl condition */
  bool          has_pdivu;
  bool          has_dissipation;

} cs_thermal_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to thermal model structure */

extern const cs_thermal_model_t  *cs_glob_thermal_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_thermal_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_thermal_model_t *
cs_get_glob_thermal_model(void);

/*----------------------------------------------------------------------------
 * Return thermal field (temperature, enthalpy, total energy according to
 * thermal model).
 *
 * returns:
 *   pointer to thermal field
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_thermal_model_field(void);

/*----------------------------------------------------------------------------
 * Print the thermal model structure to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_thermal_model_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize thermal variables if needed
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of the square of sound velocity multiplied
 *        by gamma.
 *
 * \param[in]      cp      array of isobaric specific heat values for dry air
 * \param[in]      temp    array of temperature values
 * \param[in]      pres    array of pressure values
 * \param[in,out]  fracv   array of volume fraction values
 * \param[in,out]  fracm   array of mass fraction values
 * \param[in,out]  frace   array of energy fraction values
 * \param[out]     dc2     array of the values of the square of sound velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_c_square(const cs_real_t  cp[],
                          const cs_real_t  temp[],
                          const cs_real_t  pres[],
                          const cs_real_t  fracv[],
                          const cs_real_t  fracm[],
                          const cs_real_t  frace[],
                          cs_real_t        dc2[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the
 *        temperature at constant pressure.
 *
 * \param[in]  pres   array of pressure values
 * \param[in]  temp   array of temperature values (in Kelvin)
 * \param[in]  yw     array of the total water mass fraction
 * \param[in]  rvsra  ratio gas constant h2o / dry air
 * \param[in]  cva    difference between heat capacity of the dry air
 *                    and r_pg_const
 * \param[in]  cvv    difference beteween heat capacity of the water
 *                    in its gaseous phase and r_v_cnst
 * \param[in]  cpl    heat capacity of the water in its liquid phase
 * \param[in]  l00    water latent heat
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE cs_real_t
cs_thermal_model_demdt(const cs_real_t  pres,
                       const cs_real_t  temp,
                       const cs_real_t  yw,
                       const cs_real_t  rvsra,
                       const cs_real_t  cva,
                       const cs_real_t  cvv,
                       const cs_real_t  cpl,
                       const cs_real_t  l00);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the
 *        temperature at constant internal energy.
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     temp    array of temperature values (in Kelvin)
 * \param[in]     yw      array of the total water mass fraction
 * \param[in]     rvsra   ratio gas constant h2o / dry air
 * \param[in]     cva     difference between heat capacity of the dry air
 *                        and r_pg_const
 * \param[in]     cvv     difference beteween heat capacity of the water
 *                        in its gaseous phase and r_v_cnst
 * \param[in]     cpl     heat capacity of the water in its liquid phase
 * \param[in]     l00     water latent heat
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE cs_real_t
cs_thermal_model_demdt_ecsnt(const cs_real_t  pres,
                             const cs_real_t  temp,
                             const cs_real_t  yw,
                             const cs_real_t  rvsra,
                             const cs_real_t  cva,
                             const cs_real_t  cvv,
                             const cs_real_t  cpl,
                             const cs_real_t  l00);

/*----------------------------------------------------------------------------*/
/*!
 * \brief First pass to compute the contribution of the kinetic energy based
 *        source term from the prediction step
 *
 * \param[in]       imasfl    inner mass flux used in the momentum equation
 * \param[in]       bmasfl    boundary mass flux used in the momentum equation
 * \param[in]       vela      velocity at previous time step
 * \param[in]       vel       velocity at iteration k
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_kinetic_st_prepare(const cs_real_t  imasfl[],
                                    const cs_real_t  bmasfl[],
                                    const cs_real_t  vela[][3],
                                    const cs_real_t  vel[][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the computation of the kinetic energy based source term
 *
 * \param[in]       cromk1    density values at time n+1/2,k-1
 * \param[in]       cromk     density values at time n+1/2,k
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_kinetic_st_finalize(const cs_real_t  cromk1[],
                                     const cs_real_t  cromk[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the kinetic source term if needed
 *
 * \param[in, out]  smbrs  RHS of the thermal equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_add_kst(cs_real_t  smbrs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the pressure equation.
 *
 * \param[in]       croma     density values at the last time iteration
 * \param[in]       trav2     predicted velocity
 * \param[in]       cvara_pr  pressure values at the last time iteration
 * \param[in]       imasfl    face mass fluxes
 * \param[in, out]  cflp      CFL condition related to the pressure equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cflp(const cs_real_t  croma[],
                      const cs_real_t  trav2[][3],
                      const cs_real_t  cvara_pr[],
                      const cs_real_t  imasfl[],
                      cs_real_t        cflp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the thermal equation
 *
 * \param[in]     croma       array of density values at the last time iteration
 * \param[in]     tempk       array of the temperature
 * \param[in]     tempka      array of the temperature at the previous time step
 * \param[in]     xcvv        array of the isochoric heat capacity
 * \param[in]     vel         array of the velocity
 * \param[in]     i_massflux  array of the inner faces mass fluxes
 * \param[in]     b_massflux  array of the boundary faces mass fluxes
 * \param[in]     cflt        CFL condition related to thermal equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cflt(const cs_real_t  croma[],
                      const cs_real_t  tempk[],
                      const cs_real_t  tempka[],
                      const cs_real_t  xcvv[],
                      const cs_real_t  vel[][3],
                      const cs_real_t  i_massflux[],
                      const cs_real_t  b_massflux[],
                      cs_real_t        cflt[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the isochoric heat capacity
 *
 * \param[in]     xcvv      isobaric heat capacity
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cv(cs_real_t  *xcvv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and add the dissipation term of the thermal equation to
 *        its right hand side.
 *
 * \param[in]      vistot  array for the total viscosity
 * \param[in]      gradv   tensor for the velocity gradient
 * \param[in,out]  smbrs   array of equation right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_dissipation(const cs_real_t  vistot[],
                             const cs_real_t  gradv[][3][3],
                             cs_real_t        smbrs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the Newton method to compute the temperature from the
 *        internal energy
 *
 * \param[in]       method    method used to compute the temperature
 * \param[in]       th_scal   internal energy values
 * \param[in]       pk1       pressure values at the last inner iteration
 * \param[in]       cvar_pr   pressure values
 * \param[in]       cvara_pr  pressure values at the last time iteration
 * \param[in]       yw        total water mass fraction
 * \param[in, out]  yv        vapor of water mass fraction
 * \param[in, out]  temp      temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_newton_t(int               method,
                          const cs_real_t  *pk1,
                          const cs_real_t   th_scal[],
                          const cs_real_t   cvar_pr[],
                          const cs_real_t   cvara_pr[],
                          const cs_real_t   yw[],
                          cs_real_t         yv[],
                          cs_real_t         temp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the term pdivu to the thermal equation rhs.
 *
 * \param[in, out]  smbrs     array of the right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_pdivu(cs_real_t  smbrs[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_THERMAL_MODEL_H__ */
