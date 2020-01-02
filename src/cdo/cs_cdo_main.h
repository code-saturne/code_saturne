#ifndef __CS_CDO_MAIN_H__
#define __CS_CDO_MAIN_H__

/*============================================================================
 * Routines for solving equations with CDO discretizations
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_domain.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the computational domain when CDO/HHO schemes are
 *         activated and cs_user_model() has been called
 *         At this stage of the settings, mesh quantities and adjacencies are
 *         not defined. Only the major moddeling options are set. The related
 *         equations and main properties have been added.
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_setup(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the structures related to the computational domain when
 *         CDO/HHO schemes are activated
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  m        pointer to a cs_mesh_t struct.
 * \param[in]       mq       pointer to a cs_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_structures(cs_domain_t           *domain,
                             cs_mesh_t             *m,
                             cs_mesh_quantities_t  *mq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free all structures allocated during the resolution of CDO/HHO
 *          schemes
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_finalize(cs_domain_t    *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main program for running a simulation with the CDO kernel
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_main(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_MAIN_H__ */
