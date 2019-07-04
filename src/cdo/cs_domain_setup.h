#ifndef __CS_DOMAIN_SETUP_H__
#define __CS_DOMAIN_SETUP_H__

/*============================================================================
 * Manage the definition/setting of a computation
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_advection_field.h"
#include "cs_domain.h"
#include "cs_equation.h"
#include "cs_gwf.h"
#include "cs_param.h"
#include "cs_property.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set to true the automatic update of all advection fields
 * \brief  Set auxiliary parameters related to the way output is done
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in]       nt_interval  frequency for the restart process
 * \param[in]       nt_list      output frequency into the log
 * \param[in]       verbosity    level of information displayed
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_output_param(cs_domain_t       *domain,
                           int                nt_interval,
                           int                nt_list,
                           int                verbosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters for unsteady computations: the max number of time
 *         steps or the final physical time of the simulation
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 * \param[in]       nt_max    max. number of time step iterations
 * \param[in]       t_max     final physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_time_param(cs_domain_t       *domain,
                         int                nt_max,
                         double             t_max);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set time step parameters for unsteady computations when this is not
 *         already done. This situation should occur when the GUI is used to
 *         set a constant time step.
 *
 * \param[in, out]  domain    pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_automatic_time_step_settings(cs_domain_t       *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the time step thanks to a predefined function
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      func        pointer to a cs_time_func_t function
 * \param[in]      func_input  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_function(cs_domain_t        *domain,
                                    cs_time_func_t     *func,
                                    void               *func_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the time step.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        dt        value of the constant time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_value(cs_domain_t   *domain,
                                 double         dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  First setup stage of the cs_domain_t structure
 *         Define extra domain boundaries
 *         Setup predefined equations
 *         Create fields
 *         Define cs_sles_t structures for variable fields
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_initialize_setup(cs_domain_t                 *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  After having read the mesh and the first setup stage build the
 *         connectivities and mesh quantities related to CDO/HHO schemes
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_init_cdo_structures(cs_domain_t              *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last setup stage of the cs_domain_t structure
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_finalize_setup(cs_domain_t         *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize systems of equations and their related field values
 *         according to the user settings
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_initialize_systems(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main domain settings
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_log(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_SETUP_H__ */
