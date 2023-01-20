#ifndef __CS_THERMAL_MODEL_H__
#define __CS_THERMAL_MODEL_H__

/*============================================================================
 * Base thermal model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

  CS_THERMAL_MODEL_NONE,
  CS_THERMAL_MODEL_TEMPERATURE,
  CS_THERMAL_MODEL_ENTHALPY,
  CS_THERMAL_MODEL_TOTAL_ENERGY,
  CS_THERMAL_MODEL_INTERNAL_ENERGY,
  CS_THERMAL_MODEL_N_TYPES

} cs_thermal_model_variable_t;

typedef enum {

  CS_TEMPERATURE_SCALE_NONE,
  CS_TEMPERATURE_SCALE_KELVIN,
  CS_TEMPERATURE_SCALE_CELSIUS

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
  int		        has_kinetic_st;
  int           cflt;         /* compute the thermal cfl condition */
  int           cflp;         /* compute the pressure cfl condition */
  bool          has_pdivu;
  bool          has_dissipation;
  int           unstd_multiplicator;
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
cs_thermal_model_ini(void);



/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of the square of sound velocity multiplied by gamma:
 *
 *
 * \param[in]     cp      array of isobaric specific heat values for dry air
 * \param[in]     cpv     array of isobaric specific heat values for moist air
 * \param[in]     l00     latent heat
 * \param[in]     temp    array of temperature values
 * \param[in]     pres    array of pressure values
 * \param[in,out] fracv   array of volume fraction values
 * \param[in,out] fracm   array of mass fraction values
 * \param[in,out] frace   array of energy fraction values
 * \param[out]    dc2      array of the values of the square of sound velocity
 * \param[in]     l_size  l_size of the array
 */

/*----------------------------------------------------------------------------*/

void
cs_thermal_model_c_square(cs_real_t *cp,
                          cs_real_t cpv,
                          cs_real_t cpl,
                          cs_real_t l00,
                          cs_real_t *temp,
                          cs_real_t *pres,
                          cs_real_t *fracv,
                          cs_real_t *fracm,
                          cs_real_t *frace,
                          cs_real_t *c2,
                          cs_lnum_t  l_size);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the temperature
 *        at constant pressure
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     temp    array of temperature values (in Kelvin)
 * \param[in]     yw      array of the total water mass fraction
 * \param[in]     cpa     heat capacity of the dry air
 * \param[in]     cpv     heat capacity of the water in its gaseous phase
 * \param[in]     cpl     heat capacity of the water in its liquid phase
 * \param[in]     l00     water latent heat
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_demdt(cs_real_t pres,
                       cs_real_t temp,
		                   cs_real_t yw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the temperature
 *        at constant internal energy
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     temp    array of temperature values
 * \param[in,out] demdt   array of the partial derivative of the internal energy
 * 		          related to the temperature
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_demdt_ecsnt(cs_real_t pres,
                             cs_real_t temp,
		                         cs_real_t yw,
		                         cs_real_t cpa,
		                         cs_real_t cpv,
		                         cs_real_t cpl,
		                         cs_real_t l00);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the kinetic energy based source term
 *
 * \param[in]     croma     array of density values at the last time iteration
 * \param[in]     cromaa    array of density values at the n-2 time iteration
 * \param[in]     crom_eos  density value
 * \param[in]     vel       array of velocity
 * \param[in]     vela      array of ancient velocity
 * \param[in]     sk        kinetic source term
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_compute_kinetic_st(cs_real_t      *croma,
		                                cs_real_t      *cromaa,
                                    cs_real_t      *crom_eos,
                                    cs_real_3_t    *vel,
                                    cs_real_3_t    *vela,
                                    cs_real_t      *sk);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the kinetic energy related source term in the thermal equation
 *
 * \param[in]     croma     array of density values at the last time iteration
 */
 /* ------------------------------------------------------------------------- */
cs_real_t
cs_thermal_model_add_kst(cs_real_t  *smbrs);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the pressure equation
 *
 * \param[in]     croma     array of density values at the last time iteration
 * \param[in]     trav2     array of the predicted velocity
 * \param[in]     cvara_pr  array of pressure values at the last time iteration
 * \param[in]     imasfl    array of the faces mass fluxes
 * \param[in]     cflp      CFL condition related to the pressure equation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_cflp(cs_real_t     *croma,
		                  cs_real_3_t   *trav2,
                      cs_real_t     *cvara_pr,
                      cs_real_t     *imasfl,
                      cs_real_t     *cflp);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the Newton method to compute the temperature from the
 * internal energy
 *
 * \param[in]     yw        array of total water mass fraction
 * \param[in]     yv        array of vapor of water mass fraction
 * \param[in]     temp      array of temperature values
 * \param[in]     scalt     array of internal energy values
 * \param[in]     pk1       array of pressure values at the last
 *                          inner iteration
 * \param[in]     cvar_pr   array of pressure values
 * \param[in]     cvara_pr  array of pressure values at the last time iteration
 * \param[in]     method    method used to compute the temperature
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_newton_t (cs_real_t      *yw,
                           cs_real_t      *yv,
                           cs_real_t      *temp,
                           cs_real_t      *scalt,
                           cs_real_t      *pk1,
                           cs_real_t      *cvar_pr,
                           cs_real_t      *cvara_pr,
                           int            method);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of the square
 * of sound velocity multiplied by gamma:
 *
 *
 * \param[in]       temp_     array of temperature
 * \param[in]       tempa_    array of temperature at the previous time step
 * \param[in]       cvar_var  array of the internal energy
 * \param[in]       cvara_var array of the internal energy at the previous
 *                            time step
 * \param[in]       tempa_    array of temperature at the previous time step
 * \param[in]       thetv     theta parameter
 * \param[in]       vel       array of the velocity
 * \param[in]       xcvv      array of the isobaric heat capacity
 * \param[in]       cpro_yw   array of the total water mass fraction
 * \param[in]       cpro_ywa  array of the total water mass fraction at the
 *                            previous time step
 * \param[in]       cpro_yv   array of the vapor of water mass fraction
 * \param[in]       cpro_yva  array of the vapor of water mass fraction at the
 *                            previous time step
 * \param[in]       gradp     array of the pressure gradient
 * \param[in]       gradphi   array of the pressure increment gradient
 * \param[in, out]  smbrs     array of the right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_pdivu(cs_real_t   *temp_,
                       cs_real_t   *tempa_,
                       cs_real_t   *cvar_var,
                       cs_real_t   *cvara_var,
                       cs_real_t   thetv,
                       cs_real_3_t *vel,
                       cs_real_t   *xcvv,
                       cs_real_t   *cpro_yw,
                       cs_real_t   *cpro_ywa,
                       cs_real_t   *cpro_yv,
                       cs_real_t   *cpro_yva,
                       cs_real_3_t *gradp,
                       cs_real_3_t *gradphi,
                       cs_real_t   *smbrs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and add the dissipation term of the thermal equation to
 *        its right hand side.
 *
 *
 * \param[in]     vistot  array for the total viscosity
 * \param[in]     gradv   tensor for the velocity gradient
 * \param[in,out] smbrs   array of equation right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_dissipation(cs_real_t    *vistot,
                             cs_real_33_t *gradv,
                             cs_real_t    *smbrs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the thermal equation
 *
 * \param[in]     croma     array of density values at the last time iteration
 * \param[in]     tempk     array of the temperature
 * \param[in]     tempka    array of the temperature at the previous time step
 * \param[in]     xcvv      array of the isochoric heat capacity
 * \param[in]     vel       array of the velocity
 * \param[in]     imasfl    array of the faces mass fluxes
 * \param[in]     cflt      CFL condition related to thermal equation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_cflt (cs_real_t      *croma,
                       cs_real_t      *tempk,
                       cs_real_t      *tempka,
                       cs_real_t      *xcvv,
		                   cs_real_3_t    *vel,
                       cs_real_t      *imasfl,
                       cs_real_t      *cflt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the isobaric heat capacity
 *
 * \param[in]     xcvv      isobaric heat capacity
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_model_cv(cs_real_t *xcvv);

END_C_DECLS
#endif /* __CS_THERMAL_MODEL_H__ */
