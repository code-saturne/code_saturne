#ifndef __CS_COMBUSTION_MODEL_H__
#define __CS_COMBUSTION_MODEL_H__

/*============================================================================
 * Combustion model parameters.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

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

/*! Maximum number of globals species */
#define  CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES  25

/*! Maximum number of elementary gas components */
#define  CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS  20

/*! Maximum number of coals */
#define  CS_COMBUSTION_MAX_COALS  5

/*! Maximum number of coal classes */
#define  CS_COMBUSTION_MAX_CLASSES_PER_COAL  20

/*! Maximum number of coal classes */
#define  CS_COMBUSTION_MAX_COAL_CLASSES    CS_COMBUSTION_MAX_COALS \
                                         * CS_COMBUSTION_MAX_CLASSES_PER_COAL

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

} cs_combustion_gas_model_t;

/*! Coal combustion model parameters structure */
/*---------------------------------------------*/

typedef struct {

  int     n_coals;                  /*< number of coal types */
  int     nclacp;                   /*< number of coal classes */

  /*! number of classes per coal */
  int     n_classes_per_coal[CS_COMBUSTION_MAX_COALS];

  /*! coal id if considered class belongs to coal ich[1, 2, ...] */
  int     ichcor[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! ashes concentration (kg/kg) */
  double  xashch[CS_COMBUSTION_MAX_COALS];

  /*! initial diameter (m) */
  double  diam20[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! minimum diameter (m) */
  double  dia2mn[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! initial density (kg/m^3) */
  double  rho20[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! minimal density (kg/m^3) */
  double  rho2mn[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! initial particle mass (kg) */
  double  xmp0[CS_COMBUSTION_MAX_COAL_CLASSES];

  /*! particle char mass (kg) */
  double  xmasch[CS_COMBUSTION_MAX_COAL_CLASSES];

} cs_combustion_coal_model_t;

/*! Coal combustion model parameters structure */
/*---------------------------------------------*/

typedef struct {

  int     nclafu;                   /*! number of fuel classes */
  double  hinfue;                   /*! input mass enthalpy for fuel */
  double  h02fol;                   /*! H0 of liquid fuel oil */
  double  cp2fol;                   /*! fuel oil liquid CP */

} cs_combustion_fuel_model_t;

/*! Combustion model parameters structure */
/*----------------------------------------*/

typedef struct {

  cs_combustion_gas_model_t   gas;   /*!< gas combustion model parameters */
  cs_combustion_coal_model_t  coal;  /*!< coal combustion model parameters */
  cs_combustion_fuel_model_t  fuel;  /*!< fuel combustion model parameters */

  int     n_gas_el_comp;             /*!< number of elementary gas components */
  int     n_gas_species;             /*!< number of global species */
  int     n_atomic_species;          /*!< number of atomic species */

  int     n_reactions;               /*!< number of global reactions
                                      *   in gas phase */

  int     isoot;                     /*!< soot production modeling flag */

  int     io2;                       /*!< index of o2 in wmole */
  int     in2;                       /*!< index of n2 in wmole */
  int     ico2;                      /*!< index of co2 in wmole */
  int     ih2o;                      /*!< index of h2o in wmole */

  double  ckabs0;                    /*!< absorption coefficient of gas mix */

  double  xco2;                      /*!< molar coefficient of CO2 */
  double  xh2o;                      /*!< molar coefficient of H2O */
  double  hinoxy;                    /*!< input mass enthalpy for the oxidant */

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
 * \brief Print the combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_MODEL_H__ */
