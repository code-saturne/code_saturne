#ifndef __CS_DOMAIN_H__
#define __CS_DOMAIN_H__

/*============================================================================
 * Manage a computational domain
 *  - equations, settings, fields, connectivities and geometrical quantities
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 *  Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_boundary.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
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

/*! \enum cs_domain_stage_t

 *  \brief Indicator describing at which stage is the computation.
 *
 * This indicator is useful to know at which stage the computation is when a
 * user-defined function is called (for instance \ref
 * cs_user_physical_properties). Up to now, this information is only used in
 * the CDO module.
 *
 * \var CS_DOMAIN_STAGE_BEFORE_STEADY_COMPUTATION
 * The computation run is at a stage before computing the steady-state equations
 *
 * \var CS_DOMAIN_STAGE_BEFORE_TIME_LOOP
 * The computation run is at a stage before starting the time loop and after
 * solving all steady-state equations
 *
 * \var CS_DOMAIN_STAGE_TIME_STEP_BEGIN
 * The computation run is inside the time loop and at the beginning of a
 * time step
 *
 * \var CS_DOMAIN_STAGE_TIME_STEP_SUB_ITERATION
 * The computation run is inside the time loop and inside a sub-iteration
 * process
 *
 * \var CS_DOMAIN_STAGE_TIME_STEP_END
 * The computation run is inside the time loop and at the end of an iteration
 * but before the computation of user-defined equation (not triggered by a
 * user)
 *
 * \var CS_DOMAIN_STAGE_AFTER_TIME_LOOP
 * The computation run is at a stage after the time loop (finalize stage).
 * This stage is useful to free memory allocated during the computation.
 */

typedef enum {

  CS_DOMAIN_STAGE_BEFORE_STEADY_COMPUTATION,
  CS_DOMAIN_STAGE_BEFORE_TIME_LOOP,
  CS_DOMAIN_STAGE_TIME_STEP_BEGIN,
  CS_DOMAIN_STAGE_TIME_STEP_SUB_ITERATION,
  CS_DOMAIN_STAGE_TIME_STEP_END,
  CS_DOMAIN_STAGE_AFTER_TIME_LOOP,

  CS_DOMAIN_N_STAGES

} cs_domain_stage_t;


/*! \struct cs_domain_cdo_context_t
 *  \brief  High-level metadata for handling CDO/HHO schemes
 */

typedef struct {

  /* Mode for CDO: activated, switched off... */

  int                       mode;

  /* Flag to know if scalar or vector equations are requested and which kind
     of numerical schemes is requested to solve these equations */

  cs_flag_t                 eb_scheme_flag;
  cs_flag_t                 fb_scheme_flag;
  cs_flag_t                 vb_scheme_flag;
  cs_flag_t                 vcb_scheme_flag;
  cs_flag_t                 hho_scheme_flag;

} cs_domain_cdo_context_t;

/*! \struct cs_domain_t
 *  \brief  Structure storing the main features of the computational domain
 *  and pointers to the main geometrical structures
 */

typedef struct {

  /* code_saturne mesh and mesh quantities structures already computed */

  cs_mesh_t                *mesh;
  cs_mesh_quantities_t     *mesh_quantities;

  /* CDO structures:
   * - cs_cdo_connect_t contains additional information about connectivity
   * - cs_cdo_quantities_t contains additional information on mesh quantities
   */

  cs_cdo_connect_t         *connect;
  cs_cdo_quantities_t      *cdo_quantities;

  /* Boundary of the computational domain */

  cs_boundary_t            *boundaries;
  cs_boundary_t            *ale_boundaries;

  /* Time step management */

  bool                      only_steady;
  bool                      is_last_iter;     /* true or false */

  cs_time_step_t           *time_step;        /* time step descriptor */
  cs_time_step_options_t    time_options;     /* time step options */

  cs_domain_stage_t         stage;   /* store the stage of the computation */

  /* Output options */

  int                       output_nt;   /* Logging done every nt iterations */
  int                       restart_nt;  /* Restart done every nt iterations */
  int                       verbosity;   /* Level of details given in log */

  /* Specific context structure related to the numerical schemes */

  cs_domain_cdo_context_t   *cdo_context;

  /* Monitoring */

  cs_timer_counter_t    tcp; /* Cumulated elapsed time for extra-operations
                                and post-processing */
  cs_timer_counter_t    tca; /* Cumulated elapsed time for all operations */

} cs_domain_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_domain_t  *cs_glob_domain; /* Pointer to main computational domain */

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_domain_increment_time_step(cs_domain_t  *domain)
{
  cs_time_step_t  *ts = domain->time_step;

  /* Increment time iteration */
  ts->nt_cur++;
}

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
 * \brief  Set the computation stage in the domain structure
 *
 * \param[in, out] domain    pointer to a cs_domain_t structure
 * \param[in]      stage     stage in the computation run
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_stage(cs_domain_t         *domain,
                    cs_domain_stage_t    stage);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the computation stage from the domain structure
 *
 * \param[in] domain    pointer to a cs_domain_t structure
 *
 * \return the current stage in the computation run
 */
/*----------------------------------------------------------------------------*/

cs_domain_stage_t
cs_domain_get_stage(const cs_domain_t    *domain);

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
 * \brief  Check if an output is requested according to the domain setting
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 * \param[in]   oneplus   add or not plus one to the current time step
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_log(const cs_domain_t      *domain,
                    bool                    oneplus);

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
