#ifndef __CS_LAGR_POST_H__
#define __CS_LAGR_POST_H__

/*============================================================================
 * Lagrangian module postprocessing
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "assert.h"
#include "cs_base.h"
#include "cs_field.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
  Global variables
  ============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Lagrangian postprocessing.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_post_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactive postprocessing for a given particle attribute.
 *
 * \param[in]  attr_id  associated attribute id
 *
 * \return     true if output of given attribute is active, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_lagr_post_get_attr(cs_lagr_attribute_t  attr_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactive postprocessing for a given particle attribute.
 *
 * \param[in]  attr_id  associated attribute id
 * \param[in]  active   true if postprocessing is required, false otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_post_set_attr(cs_lagr_attribute_t  attr_id,
                      bool                 active);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_POST_H__ */
