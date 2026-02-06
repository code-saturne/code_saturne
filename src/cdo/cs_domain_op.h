#ifndef __CS_DOMAIN_OP_H__
#define __CS_DOMAIN_OP_H__

/*============================================================================
 * Manage specific post-processing related to a computational domain and
 * restart files
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "cdo/cs_domain.h"
#include "base/cs_time_step.h"

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
 * \brief Add a time control dedicated to the extra-operations.
 *        The lifecycle of the time control is now managed by the domain
 *        structure
 *
 * \param[in,out] domain  pointer to a cs_domain_t structure
 * \param[in]     tc      pointer to a time control structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_extra_op_time_control(cs_domain_t        *domain,
                                    cs_time_control_t  *tc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the generic post-processing related to a domain
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_init(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extra-operations of the computational domain after a resolution.
 *        Mostly used inside CDO schemes
 *        User-defined extra-operations are also called inside this function.
 *
 * \param[in, out] domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_extra_op(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a restart file for the CDO/HHO module
 *
 * \param[in, out]  domain     pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_read_restart(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a restart file for the CDO/HHO module
 *
 * \param[in]  domain     pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_write_restart(const cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_OP_H__ */
