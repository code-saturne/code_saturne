#ifndef __CS_COMBUSTION_GAS_H__
#define __CS_COMBUSTION_GAS_H__

/*============================================================================
 * Gas combustion model.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*! Maximum number of global species */
#define  CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES  25

/*! Maximum number of atomic species */
#define  CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES  5

/*! Maximum number of oxydants */
#define CS_COMBUSTION_GAS_MAX_OXYDANTS 3

/*! Maximum number of elementary gas components */
#define  CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS  20

/*! Maximum number of global reactions in gas phase*/
#define CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS 1

/*! Maximum number of tabulation points */
#define  CS_COMBUSTION_GAS_MAX_TABULATION_POINTS  50

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Gas combustion model type */
/*------------------------------*/

typedef enum {

  CS_COMBUSTION_GAS_NONE = -1,

  /*! Infinitely fast 3-point combustion model, adiabatic */
  CS_COMBUSTION_3PT_ADIABATIC = 100,

  /*! Infinitely fast 3-point combustion model, permeatic */
  CS_COMBUSTION_3PT_PERMEATIC = 101,

  /*! Steady laminar flamelet model, adiabatic conditions */
  CS_COMBUSTION_SLFM_STEADY_ADIABATIC = 200,

  /*! Steady laminar flamelet model, enthalpy transport with heat loss */
  CS_COMBUSTION_SLFM_STEADY_PERMEATIC = 201,

  /*! Flamelet/progress variable model, adiabatic conditions */
  CS_COMBUSTION_SLFM_PROGRESS_ADIABATIC = 202,

  /*! Flamelet/progress variable model, enthalpy transport with heat loss */
  CS_COMBUSTION_SLFM_PROGRESS_PERMEATIC = 203,

  /*! Eddy Break Up pre-mixed flame combustion model,
    adiabatic conditions at constant richness */
  CS_COMBUSTION_EBU_CONSTANT_ADIABATIC = 300,

  /*! Eddy Break Up pre-mixed flame combustion model,
    adiabatic conditions at constant richness */
  CS_COMBUSTION_EBU_CONSTANT_PERMEATIC = 301,

  /*! Eddy Break Up pre-mixed flame combustion model,
    adiabatic conditions at variable richness */
  CS_COMBUSTION_EBU_VARIABLE_ADIABATIC = 302,

  /*! Eddy Break Up pre-mixed flame combustion model,
    adiabatic conditions at variable richness */
  CS_COMBUSTION_EBU_VARIABLE_PERMEATIC = 303,

  /*! Libby-Williams pre-mixed flame combustion,
    two peak model with adiabatic conditions */
  CS_COMBUSTION_LW_2PEAK_ADIABATIC = 400,

  /*! Libby-Williams pre-mixed flame combustion,
    two peak model with permeatic conditions */
  CS_COMBUSTION_LW_2PEAK_PERMEATIC = 401,

  /*! Libby-Williams pre-mixed flame combustion,
    three peak model with adiabatic conditions */
  CS_COMBUSTION_LW_3PEAK_ADIABATIC = 402,

  /*! Libby-Williams pre-mixed flame combustion,
    three peak model with permeatic conditions */
  CS_COMBUSTION_LW_3PEAK_PERMEATIC = 403,

  /*! Libby-Williams pre-mixed flame combustion,
    four peak model with adiabatic conditions */
  CS_COMBUSTION_LW_4PEAK_ADIABATIC = 404,

  /*! Libby-Williams pre-mixed flame combustion,
    four peak model with permeatic conditions */
  CS_COMBUSTION_LW_4PEAK_PERMEATIC = 405

} cs_combustion_gas_model_type_t;

/*! Gas combustion model parameters structure */
/*--------------------------------------------*/

typedef struct {

  /* Generic members
     ---------------
     (keep aligned with gas combustion, so that we can
     move to an object inheritance model in the future) */

  int     n_gas_el_comp;             /*!< number of elementary gas components */
  int     n_gas_species;             /*!< number of global species */
  int     n_atomic_species;          /*!< number of atomic species */

  int     n_reactions;               /*!< number of global reactions
                                      *   in gas phase */

  int     n_tab_points;              /*!< number of tabulation points */

  double  pcigas;                    /*!< combustible reaction enthalpy
                                       (Lower Calorific Value)*/
  double  xco2;                      /*!< molar coefficient of CO2 */
  double  xh2o;                      /*!< molar coefficient of H2O */

  /*! molar mass of an elementary gas component */
  double  wmole[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  /*! molar mass of atomic species */
  double  wmolat[CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES];

  /*! composition of N2 oxidants */
  double oxyn2[CS_COMBUSTION_GAS_MAX_OXYDANTS];

  /*! composition of H2O oxidants */
  double oxyh2o[CS_COMBUSTION_GAS_MAX_OXYDANTS];

  /*! composition of CO2 oxidants */
  double oxyco2[CS_COMBUSTION_GAS_MAX_OXYDANTS];

  /*! temperature in K */
  double th[CS_COMBUSTION_GAS_MAX_TABULATION_POINTS];

  /* Members specific to the gas combustion model
     -------------------------------------------- */

  cs_combustion_gas_model_type_t  type;  /*!< combustion model type */

  char   *data_file_name;          /* Name of thermochemistry data file */

  bool    use_janaf;               /*! use NIST-JANAF tables for
                                     enthalpy/temperature conversion */

  int     iic;                     /*!< rank of C in gas composition (base 1) */
  int     isoot;                   /*!< soot production modeling flag */

  double  hinfue;                  /*!< input mass enthalpy for fuel */
  double  tinfue;                  /*!< input temperature for fuel en K */
  double  tinoxy;                  /*!< input temperature for oxydant */
  double  hinoxy;                  /*!< input mass enthalpy for the oxidant */

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
  double cpgazg[CS_COMBUSTION_GAS_MAX_TABULATION_POINTS]
               [CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];

  /* Numerical parameters
     -------------------- */

  double srrom;  /*!< sub-relaxation coefficient for the density:
                   \f$\rho^{n+1}$\,=\,srrom\,$\rho^n$+(1-srrom)\,$\rho^{n+1}\f$
                   (hence, with a zero value, there is no sub-relaxation) */

  /*! Libby Williams parameters */
  double fmin_lwc;
  double fmax_lwc;
  double hmin_lwc;
  double hmax_lwc;

} cs_combustion_gas_model_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Combustion model parameters structure */

extern cs_combustion_gas_model_t  *cs_glob_combustion_gas_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Activate gas combustion model.
 *
 * \return  pointer to gas combustion model structure.
 *
 * \param[in]  type  gas combustion model type
 */
/*----------------------------------------------------------------------------*/

cs_combustion_gas_model_t  *
cs_combustion_gas_set_model(cs_combustion_gas_model_type_t  type);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set the thermochemical data file name.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_set_thermochemical_data_file(const char  *file_name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute molar and mass fractions of
 *        elementary species Ye, Xe (fuel, O2, CO2, H2O, N2) from
 *        global species Yg (fuel, oxidant, products).
 *
 * \param[in]     yg            global mass fraction
 * \param[out]    ye            elementary mass fraction
 * \param[out]    xe            elementary molar fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_yg2xye(const cs_real_t  yg[],
                         cs_real_t        ye[],
                         cs_real_t        xe[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Print the gas combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_GAS_H__ */
