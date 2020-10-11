#ifndef __CS_CDO_TURBULENCE_H__
#define __CS_CDO_TURBULENCE_H__

/*============================================================================
 * Routines to handle the resolution of the turbulence modelling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_cdo_turbulence_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         tubulence modelling. Several members are structures defined in
 *         cs_turbulence_model.h
 *
 */

typedef struct {

  /*!
   * @name Turbulence modelling
   * Set of parameters to handle turbulence modelling.
   * @{
   */

  /*! \var model_param
   * Main set of parameters to handle turbulence modelling. This
   * structure is shared with the legacy part.
   */

  cs_turb_model_t            *model_param;

  /*! \var rans_modelling
   * Main set of parameters to handle RANS modelling. This
   * structure is shared with the legacy part.
   * RANS means Reynolds Average Navier-Stokes
   */

  cs_turb_rans_model_t       *rans_modelling;

  /*! \var les_modelling
   * Main set of parameters to handle LES modelling. This
   * structure is shared with the legacy part.
   * LES means Large Eddy Simulation
   */

  cs_turb_les_model_t        *les_modelling;

  /*! \var reference_values
   *  Set of reference values associated to the turbulence modelling
   */

  cs_turb_ref_values_t       *reference_values;

} cs_cdo_turbulence_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure managing the turbulence modelling
 *
 * \return a pointer to a new allocated cs_cdo_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_turbulence_t *
cs_cdo_turbulence_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure managing the turbulence modelling
 *
 * \param[in, out]  p_turb_struct   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_turbulence_free(cs_cdo_turbulence_t   **p_turb_struct);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the structure managing the turbulence modelling
 *
 * \param[in, out]  turb_struct   pointer to the structure to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_turbulence_init(cs_cdo_turbulence_t   *turb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_TURBULENCE_H__ */
