#ifndef __CS_NAVSTO_SYSTEM_H__
#define __CS_NAVSTO_SYSTEM_H__

/*============================================================================
 * Routines to handle cs_navsto_system_t structure
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

#include "cs_advection_field.h"
#include "cs_field.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_mesh.h"
#include "cs_navsto_param.h"
#include "cs_time_step.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_navsto_system.h

  \brief  Routines to handle the cs_navsto_system_t structure
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the context structure related to a given
 *         discretization scheme for the resolution of the Navier-Stokes system
 *
 * \param[in]      nsp        pointer to a cs_navsto_param_t structure
 * \param[in, out] nscc       Navier-Stokes coupling context: pointer to a
 *                            structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_navsto_init_scheme_context_t)(const cs_navsto_param_t    *nsp,
                                  const void                 *nscc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the context structure related to a given
 *         discretization scheme for the resolution of the Navier-Stokes system
 *
 * \param[in]      nsp        pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_navsto_free_scheme_context_t)(const cs_navsto_param_t      *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the
 *         Navier-Stokes system. This means that equations are built and then
 *         solved.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      nsp        pointer to a cs_navsto_param_t structure
 * \param[in, out] nscc       Navier-Stokes coupling context: pointer to a
 *                            structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_navsto_compute_t)(const cs_mesh_t              *mesh,
                      double                        dt_cur,
                      const cs_navsto_param_t      *nsp,
                      void                         *nscc);

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*! \struct cs_navsto_system_t
 *  \brief  Structure managing the Navier-Stokes system
 *
 */

typedef struct {

  /*! \var param
   *  Set of parameters to handle the Navier-Stokes system
   */
  cs_navsto_param_t  *param;

  /*!
   * @name Fields
   * Set of fields (resolved variables): fields are created according to the
   * choice of model for Navier-Stokes
   * @{
   */

  cs_adv_field_t     *adv_field;
  cs_field_t         *velocity;
  cs_field_t         *pressure;
  cs_field_t         *temperature;

  /*!
   * @}
   * @name Physical properties
   * Set of properties: properties and their related fields are allocated
   * according to the choice of model for Navier-Stokes
   * @{
   */

  cs_property_t      *density;
  cs_property_t      *lami_viscosity;

  /*! \var context
   * Additional structure allocated if needed i.e. according to the choice
   * of model for the Navier-Stokes system and the way the equations are
   * solved
   */
  void               *context;

  /*!
   * @}
   * @name Pointer to functions handling specific tasks
   * @{
   */

  /*! \var init
   *  Pointer of functions related to the initialization of the context
   *  structure related to a given discretization scheme for the resolution
   *  of the Navier-Stokes system
   */
  cs_navsto_init_scheme_context_t   *init;

  /*! \var free
   *  Pointer of functions related to the destruction of the context
   *  structure related to a given discretization scheme for the resolution
   *  of the Navier-Stokes system
   */
  cs_navsto_free_scheme_context_t   *free;

  /*! \var compute
   *  Pointer of functions related to resolution of the Navier-Stokes system
   *  Handle the build of the system and its resolution
   */
  cs_navsto_compute_t               *compute;

  /*!
   * @}
   * @name Initial conditions
   *
   * Set of parameters used to take into account the initial condition on the
   * pressure and/or the velocity.
   * CAUTION: so far, there is no check if the different IC are compatible
   * @{
   */

  /*! \var n_velocity_ic_defs
   *  Number of initial conditions associated to the velocity
   */
   int  n_velocity_ic_defs;

  /*! \var n_pressure_ic_defs
   *  Number of initial conditions associated to the pressure
   */
   int  n_pressure_ic_defs;

  /*! \var velocity_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  velocity
   */
   cs_xdef_t  **velocity_ic_defs;

  /*! \var pressure_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  pressure
   */
   cs_xdef_t  **pressure_ic_defs;

  /*! @} */

} cs_navsto_system_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the Navier-Stokes system has been
 *        activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_navsto_system_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the Navier-Stokes (NS) system
 *
 * \param[in] model          type of model related to the NS system
 * \param[in] time_state     state of the time for the NS equations
 * \param[in] algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(cs_navsto_param_model_t        model,
                          cs_navsto_param_time_state_t   time_state,
                          cs_navsto_param_coupling_t     algo_coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Recover the structure storing the parameters for the Navier--Stokes
 *         system
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_system_get_param(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes system
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the Navier-Stokes system
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_finalize_setup(const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns..
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_velocity_ic_by_value(const char    *z_name,
                                   cs_real_t     *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns..
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_pressure_ic_by_value(const char    *z_name,
                                   cs_real_t     *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_velocity_ic_by_analytic(const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_pressure_ic_by_analytic(const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context structure used to build the algebraic system
 *         This is done after the setup step.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_initialize(const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system in case of a
 *         steady-state approach
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute_steady_state(const cs_mesh_t      *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      dt_cur     current value of the time step
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute(const cs_mesh_t              *mesh,
                         double                        dt_cur);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the Navier-Stokes system.
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_extra_post(void                      *input,
                            int                        mesh_id,
                            int                        cat_id,
                            int                        ent_flag[5],
                            cs_lnum_t                  n_cells,
                            cs_lnum_t                  n_i_faces,
                            cs_lnum_t                  n_b_faces,
                            const cs_lnum_t            cell_ids[],
                            const cs_lnum_t            i_face_ids[],
                            const cs_lnum_t            b_face_ids[],
                            const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_SYSTEM_H__ */
