#ifndef __CS_DOMAIN_H__
#define __CS_DOMAIN_H__

/*============================================================================
 * Manage a computational domain
 *  - equations, settings, fields, connectivities and geometrical quantities
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

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_domain_boundary.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Flag related to the activation (or not) of the CDO schemes */
#define CS_DOMAIN_CDO_MODE_OFF     -1  /* CDO schemes are not used */
#define CS_DOMAIN_CDO_MODE_WITH_FV  1  /* CDO and legacy FV schemes are used */
#define CS_DOMAIN_CDO_MODE_ONLY     2  /* CDO schemes are exclusively used */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  /* Mode for CDO: activated, switched off... */
  int                       mode;

  /* Metadata for CDO */
  bool                      force_advfield_update;

  /* Flag to know if scalar or vector equations are requested and which kind
     of numerical schemes is requested to solve these equations */
  cs_flag_t                 fb_scheme_flag;
  cs_flag_t                 vb_scheme_flag;
  cs_flag_t                 vcb_scheme_flag;
  cs_flag_t                 hho_scheme_flag;

} cs_domain_cdo_context_t;

/* Structure storing the main feature */
typedef struct {

  /* Code_Saturne mesh and mesh quantities structures already computed */
  const  cs_mesh_t              *mesh;
  const  cs_mesh_quantities_t   *mesh_quantities;

  /* CDO structures:
     - cs_cdo_connect_t contains additional information about connectivity
     - cs_cdo_quantities_t contains additional information on mesh quantities
  */
  cs_cdo_connect_t              *connect;
  cs_cdo_quantities_t           *cdo_quantities;

  /* Time step management */
  bool                      only_steady;
  bool                      is_last_iter;     // true or false
  double                    dt_cur;           // current time step
  cs_xdef_t                *time_step_def;    // Definition of the time_step
  cs_time_step_t           *time_step;        // time step descriptor
  cs_time_step_options_t    time_options;     // time step options

  /* Output options */
  int        output_nt;   /* Log information every nt iterations */
  int        restart_nt;  /* Write a restart every nt iterations */
  int        verbosity;   /* Level of details given in log */
  bool       profiling;   /* Activate a set of timer statistics (details
                             differ according to the verbosity level) */

  /* Specific context structure related to the numerical schemes */
  cs_domain_cdo_context_t   *cdo_context;

  /* Monitoring */
  cs_timer_counter_t    tcp; /* Cumulated elapsed time for extra-operations
                                and post-processing */
  cs_timer_counter_t    tcs; /* Cumulated elapsed time for setup operations */

} cs_domain_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_domain_t *cs_glob_domain; /* Pointer to main computational domain
                                       used in CDO/HHO schmes */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize by default a cs_domain_t structure
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_domain_t structure
 *
 * \param[in, out]   p_domain    pointer of pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_free(cs_domain_t   **p_domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the global variable storing the mode of activation to apply
 *          to CDO/HHO schemes
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        mode      type of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_cdo_mode(cs_domain_t    *domain,
                       int             mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the mode of activation for the CDO/HHO schemes
 *
 * \param[in]   domain       pointer to a cs_domain_t structure
 *
 * \return the mode of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

int
cs_domain_get_cdo_mode(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set auxiliary parameters related to the way output is done
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in]       nt_interval  frequency for the restart process
 * \param[in]       nt_list      output frequency into the listing
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
 * \brief  Activate profiling output
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_activate_profiling(cs_domain_t       *domain);

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
 * \brief  Define the value of the time step thanks to a predefined function
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      func        pointer to a cs_timestep_func_t function
 * \param[in]      func_input  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_function(cs_domain_t          *domain,
                                    cs_timestep_func_t   *func,
                                    void                 *func_input);

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
 * \brief  Check if one needs to continue iterations in time
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_iteration(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if an ouput is requested according to the domain setting
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_log(const cs_domain_t      *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time_step(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a welcome message indicating which mode of CDO is activated
 *
 * \param[in]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_cdo_log(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_H__ */
