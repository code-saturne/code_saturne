#ifndef __CS_LAGR_OPTIONS_H__
#define __CS_LAGR_OPTIONS_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian module initialization
 *   1) Initialisation par defaut du parametrage du module
 *   lagrangien
 *   2) Lecture du parametrage utilisateur
 *   3) Verifications du parametrage utilisateur et
 *   controles de coherence
 *   4) Initialisation des variables en COMMON et des pointeurs
 *   sur les tableaux lies aux particules, aux statistiques,
 *   aux conditions aux limites, aux variables parietales,
 *   aux donnees pour le couplage retour.
 *
 * \param[in]  isuite
 * \param[in]  iccvfg
 * \param[in]  iscalt
 * \param[in]  dtref
 */
/*----------------------------------------------------------------------------*/

void
CS_PROCF (lagopt, LAGOPT) (cs_int_t   *isuite,
                           cs_int_t   *iccvfg,
                           cs_int_t   *iscalt,
                           cs_real_t *dtref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian module initialization
 *   1) Initialisation par defaut du parametrage du module
 *   lagrangien
 *   2) Lecture du parametrage utilisateur
 *   3) Verifications du parametrage utilisateur et
 *   controles de coherence
 *   4) Initialisation des variables en COMMON et des pointeurs
 *   sur les tableaux lies aux particules, aux statistiques,
 *   aux conditions aux limites, aux variables parietales,
 *   aux donnees pour le couplage retour.
 *
 * \param[in]  isuite
 * \param[in]  iccvfg
 * \param[in]  iscalt
 * \param[in]  dtref
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_option_definition(cs_int_t   *isuite,
                          cs_int_t   *iccvfg,
                          cs_int_t   *iscalt,
                          cs_real_t *dtref);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_OPTIONS_H__ */
