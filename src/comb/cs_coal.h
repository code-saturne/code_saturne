#ifndef __CS_COAL_H__
#define __CS_COAL_H__

/*============================================================================
 * Coal combustion model.
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

/*! Maximum number of coals */
#define  CS_COMBUSTION_MAX_COALS  5

/*! Maximum number of coal classes pae coal */
#define  CS_COMBUSTION_MAX_CLASSES_PER_COAL  20

/*! Maximum total number of coal classes */
#define  CS_COMBUSTION_COAL_MAX_CLASSES    CS_COMBUSTION_MAX_COALS \
                                         * CS_COMBUSTION_MAX_CLASSES_PER_COAL

/*! Maximum number of tabulation points */
#define  CS_COMBUSTION_COAL_MAX_TABULATION_POINTS  500

/*! Maximum number of solid constituants */
#define  CS_COMBUSTION_COAL_MAX_SOLIDS  CS_COMBUSTION_MAX_COALS * 4

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*! Coal combustion model parameters structure */
/*---------------------------------------------*/

typedef struct {

  int     n_coals;     /*!< number of coal types */
  int     nclacp;      /*!< number of coal classes */

  int     nsolim;      /*!< number of solid constituents
                         (reactive coal, coke, ash) */

  /*! number of classes per coal */
  int     n_classes_per_coal[CS_COMBUSTION_MAX_COALS];

  /* Enthalpy of reactive coal, coke, and ash */

  /*! position in ehsoli array for reactive coal */
  int     ich[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for coke */
  int     ick[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for ash */
  int     iash[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for humidity */
  int     iwat[CS_COMBUSTION_MAX_COALS];

  /*! mass enthalpy (J/kg) at temperature T of solid constituent S */
  double  ehsoli[CS_COMBUSTION_COAL_MAX_TABULATION_POINTS]
                [CS_COMBUSTION_COAL_MAX_SOLIDS];

  /*! molar mass of solid constituents */
  double  wmols[CS_COMBUSTION_COAL_MAX_SOLIDS];

  /*! formation enthalpy (J/kg) of solid constituents */
  double  eh0sol[CS_COMBUSTION_COAL_MAX_SOLIDS];

  /*! coal id if considered class belongs to coal ich[1, 2, ...] */
  int     ichcor[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! coal specific heat */
  double  cp2ch[CS_COMBUSTION_MAX_COALS];

  /*! ashes concentration (kg/kg) */
  double  xashch[CS_COMBUSTION_MAX_COALS];

  /*! humidity (kg/kg) */
  double  xwatch[CS_COMBUSTION_MAX_COALS];

  /*! initial diameter (m) */
  double  diam20[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! minimum diameter (m) */
  double  dia2mn[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! initial density (kg/m^3) */
  double  rho20[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! minimal density (kg/m^3) */
  double  rho2mn[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! initial particle mass (kg) */
  double  xmp0[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! particle char mass (kg) */
  double  xmasch[CS_COMBUSTION_COAL_MAX_CLASSES];

} cs_coal_model_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create coal model object.
 *
 * \return  pointer to coal model structure.
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy coal model object.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_model_destroy(cs_coal_model_t  **cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take in account the radiative source terms in the particle equation
 *        of a given class for pulverized coal flame.
 *
 * \param[in]      f       pointer to scalar field
 * \param[in, out] smbrs   right and side (explicit ST contribution)
 * \param[in, out] rovsdt  system diagonal (implicit ST contribution)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_rad_transfer_st(const cs_field_t  *f,
                        cs_real_t         *smbrs,
                        cs_real_t         *rovsdt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute gas enthalpy
 *        Function with gas temperature and concentrations
 *
 * \param[in] xesp      mass fraction of species
 * \param[in] f1mc      average f1
 * \param[in] f2mc      average f2
 * \param[in] tp        gas temperature (in kelvin)
 *
 * \return   gas enthalpy (in \f$ j . kg^{-1}) \f$ of mixed gas
 */
/*----------------------------------------------------------------------------*/

double
cs_coal_thconvers1(const double  xesp[],
                   const double  f1mc[],
                   const double  f2mc[],
                   cs_real_t     tp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute particles enthalpy
 *         Function with temperature and concentrations

 * \param[in]     class_id      class id
 * \param[in]     xsolid        mass fraction of components
 * \param[in,out] temper        temperature (in kelvin)

 * \return   mass enthalpy (in \f$ j . kg^{-1}) \f$
 */
/*----------------------------------------------------------------------------*/

double
cs_coal_thconvers2(int           class_id,
                   const double  xsolid[],
                   cs_real_t     temper);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COAL__ */
