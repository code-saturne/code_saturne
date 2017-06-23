#ifndef __CS_DOMAIN_POST_H__
#define __CS_DOMAIN_POST_H__

/*============================================================================
 * Manage specific post-processing related to a computational domain
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_property.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the generic post-processing related to a domain
 *
 * \param[in]  dt        reference time step value
 * \param[in]  quant     pointer to a cs_cdo_quantities_t
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_init(double                dt,
                    cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the hidden view of the domain dedicated for post-processing
 *
 * \param[in]    dt      current value of the time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_update(double    dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate writers and output meshes if needed
 *
 * \param[in]  time_step    pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_activate(cs_time_step_t    *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the hidden view of the domain dedicated for post-processing
 *
 * \param[in]  time_step    pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post(cs_time_step_t    *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize post-processing related to the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_POST_H__ */
