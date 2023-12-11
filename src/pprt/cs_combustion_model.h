#ifndef __CS_COMBUSTION_MODEL_H__
#define __CS_COMBUSTION_MODEL_H__

/*============================================================================
 * Combustion model parameters.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_coal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*! Maximum number of globals species */
#define  CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES  25

/*! Maximum number of elementary gas components */
#define  CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS  20

/*! Maximum number of tabulation points */
#define  CS_COMBUSTION_MAX_TABULATION_POINTS  500

/*! Maximum number of global reactions in gas phase*/
#define CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS 1

/*! Maximum number of oxydants */
#define CS_COMBUSTION_MAX_OXYDANTS 3

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Gas combustion model parameters structure */
/*--------------------------------------------*/

typedef struct {

  int     iic;                     /*!< rank of C in gas composition (base 1) */

  double  hinfue;                  /*!< input mass enthalpy for fuel */
  double  tinfue;                  /*!< input temperature for fuel en K */

  double  xsoot;                   /*!< soot fraction production (isoot = 0) */
  double  rosoot;                  /*!< soot density */

  /*! molar mass of global species */
  double  wmolg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];

  /*! mass fraction conversion coefficients from global species to
      elementary species */
  double  coefeg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES]
                [CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  /*! mole fraction conversion coefficients from global species to
      elementary species */
  double  compog[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES]
                [CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  /*! Mixing rate at the stoichiometry */
  double fs[CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS];

  /*! cpgazg[j][i] is the massic calorific capacity
      (J/kg/K) of the i-th global species at temperature */
  double cpgazg[CS_COMBUSTION_MAX_TABULATION_POINTS]
               [CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];

} cs_combustion_gas_model_t;

/*! Combustion model parameters structure */
/*----------------------------------------*/

typedef struct {

  cs_combustion_gas_model_t  *gas;   /*!< gas combustion model parameters */
  cs_coal_model_t            *coal;  /*!< coal combustion model parameters */

  int     n_gas_el_comp;             /*!< number of elementary gas components */
  int     n_gas_species;             /*!< number of global species */
  int     n_atomic_species;          /*!< number of atomic species */

  int     n_reactions;               /*!< number of global reactions
                                      *   in gas phase */

  int     n_tab_points;              /*!< number of tabulation points */

  int     idrift;                    /*!< drift (0: off, 1: on) */

  int     ieqco2;                    /*!< kinetic model for CO <=> CO2
                                       - 0  unused (maximal conversion
                                            in turbulent model)
                                       - 1  transport of CO2 mass fraction
                                       - 2  transport of CO mass fraction  */

  int     ieqnox;                    /*!< NOx model (0: off; 1: on) */

  int     isoot;                     /*!< soot production modeling flag */

  int     io2;                       /*!< index of o2 in wmole */
  int     in2;                       /*!< index of n2 in wmole */
  int     ico2;                      /*!< index of co2 in wmole */
  int     ih2o;                      /*!< index of h2o in wmole */

  double  ckabs0;                    /*!< absorption coefficient of gas mix */
  double  diftl0;                    /*!< molecular diffusivity for the enthalpy
                                       for gas or coal combustion
                                       (\ref diffusivity_ref is automatically set
                                       to \ref diftl0 for the enthalpy). */
  double  pcigas;                    /*!< combustible reaction enthalpy
                                       (Lower Calorific Value)*/
  double  xco2;                      /*!< molar coefficient of CO2 */
  double  xh2o;                      /*!< molar coefficient of H2O */
  double  hinoxy;                    /*!< input mass enthalpy for the oxidant */
  double  tinoxy;                    /*! input temperature for oxydant */

  /*! molar mass of an elementary gas component */
  double  wmole[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  /*! composition of oxidants in O2 */
   double oxyo2[CS_COMBUSTION_MAX_OXYDANTS];

  /*! composition of N2 oxidants */
  double oxyn2[CS_COMBUSTION_MAX_OXYDANTS];

  /*! composition of H2O oxidants */
  double oxyh2o[CS_COMBUSTION_MAX_OXYDANTS];

  /*! composition of CO2 oxidants */
  double oxyco2[CS_COMBUSTION_MAX_OXYDANTS];

  /*! temperature in K */
  double th[CS_COMBUSTION_MAX_TABULATION_POINTS];

} cs_combustion_model_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Combustion model parameters structure */

extern cs_combustion_model_t  *cs_glob_combustion_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize combustion model based on active physical models
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize combustion model based on active physical models
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_MODEL_H__ */
