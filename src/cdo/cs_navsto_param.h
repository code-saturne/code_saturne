#ifndef __CS_NAVSTO_PARAM_H__
#define __CS_NAVSTO_PARAM_H__

/*============================================================================
 * Routines to handle cs_navsto_param_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_NAVSTO_PARAM_MODEL_STOKES,
  CS_NAVSTO_PARAM_MODEL_NAVIER_STOKES, /* (U, p) standard formulation */
  CS_NAVSTO_PARAM_N_MODELS

} cs_navsto_param_model_t;

typedef enum {

  /* Algebraic param solved in one block */
  CS_NAVSTO_PARAM_COUPLING_FULL,

  /* (U,p) coupling using an artificial compressibility algorithm */
  CS_NAVSTO_PARAM_COUPLING_ARTIFICAL_COMPRESSIBILITY,

  /* (U,p) coupling using an incremental projection algorithm */
  CS_NAVSTO_PARAM_COUPLING_PROJECTION,

  CS_NAVSTO_PARAM_N_COUPLINGS

} cs_navsto_param_coupling_t;

typedef struct {

  cs_navsto_param_model_t       model;
  cs_navsto_param_coupling_t    coupling;


} cs_navsto_param_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes system
 *
 * \param[in]  model     model related to the system to solve
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(cs_navsto_param_model_t        model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  param    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *param);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_PARAM_H__ */
