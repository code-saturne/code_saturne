#ifndef __CS_NAVSTO_SYSTEM_H__
#define __CS_NAVSTO_SYSTEM_H__

/*============================================================================
 * Routines to handle cs_navsto_system_t structure
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
 *         with a specified coupling algorithm
 *
 * \param[in]      nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in]      fb_type    type of boundary for each boundary face
 * \param[in, out] nscc       Navier-Stokes coupling context: pointer to a
 *                            structure cast on-the-fly
 *
 * \return a pointer to a new allocated \ref cs_cdofb_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_navsto_init_scheme_context_t)(const cs_navsto_param_t    *nsp,
                                  cs_boundary_type_t         *fb_type,
                                  void                       *nscc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a given discretization scheme
 *         for the resolution of the Navier-Stokes system with a specified
 *         coupling algorithm
 *
 * \param[in, out]  scheme_context   pointer to a structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_navsto_free_scheme_context_t)(void     *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  According to the model, coupling algorithm and the space
 *         discretization, initialize the field values which are not associated
 *         to a \ref cs_equation_t structure (which manages the initialization)
 *
 * \param[in]       nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in]       quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]       ts       pointer to a \ref cs_time_step_t structure
 * \param[in, out]  field    pointer to a \ref cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_navsto_init_values_t)(const cs_navsto_param_t     *nsp,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *ts,
                          cs_field_t                  *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the
 *         Navier-Stokes system. This means that equations are built and then
 *         solved.
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_navsto_compute_t)(const cs_mesh_t              *mesh,
                      const cs_navsto_param_t      *nsp,
                      void                         *scheme_context);

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
  cs_navsto_param_t     *param;

  /*! \var boundary_type
   * Array storing the type of boundary for each boundary face
   */
  cs_boundary_type_t    *bf_type;

  /*!
   * @name Fields
   * Set of fields (resolved variables): fields are created according to the
   * choice of model for Navier-Stokes
   * @{
   */

  /*! \var adv_field
   *  Advection field, pointer to \ref cs_adv_field_t
   */

  cs_adv_field_t       *adv_field;

  /*! \var velocity
   *  Velocity, vector-valued, pointer to \ref cs_field_t
   */

  cs_field_t           *velocity;

  /*! \var velocity_divergence
   *  Divergence of the velocity fied, scalar-valued
   *  pointer to \ref cs_field_t
   */

  cs_field_t           *velocity_divergence;

  /*! \var pressure
   *  Pressure, scalar-valued, pointer to \ref cs_field_t
   */

  cs_field_t           *pressure;

  /*! \var temperature
   *  Temperature, scalar-var, pointer to \ref cs_field_t. NULL if no
   *  energy-like equation is defined
   */

  cs_field_t           *temperature;

  /*! \var coupling_context
   * Additional structure storing information according to the way equations
   * of model for the Navier-Stokes system are coupled and thus solved
   */
  void                 *coupling_context;

  /*! \var scheme_context
   * Additional structure storing information according to the space
   * discretization scheme used for solving the model for the Navier-Stokes
   * system
   */
  void                 *scheme_context;

  /*!
   * @}
   * @name Pointer to functions handling specific tasks
   * @{
   */

  /*! \var init_scheme_context
   *  Pointer of functions related to the initialization of the context
   *  structure related to a given discretization scheme for the resolution
   *  of the Navier-Stokes system
   */
  cs_navsto_init_scheme_context_t   *init_scheme_context;

  /*! \var free_scheme_context
   *  Pointer of functions related to the destruction of the context
   *  structure related to a given discretization scheme for the resolution
   *  of the Navier-Stokes system
   */
  cs_navsto_free_scheme_context_t   *free_scheme_context;

  /*! \var init_velocity
   *  Pointer of functions related to the initialization of variable values
   *  Case of the velocity
   */
  cs_navsto_init_values_t           *init_velocity;

  /*! \var init_pressure
   *  Pointer of functions related to the initialization of variable values
   *  Case of the pressure
   */
  cs_navsto_init_values_t           *init_pressure;

  /*! \var compute_steady
   *  Pointer of functions related to resolution of the Navier-Stokes steady
   *  system. Handle the build of the system and its resolution
   */
  cs_navsto_compute_t               *compute_steady;

  /*! \var compute
   *  Pointer of functions related to resolution of the Navier-Stokes unsteady
   *  system. Handle the build of the system and its resolution
   */
  cs_navsto_compute_t               *compute;

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
 * \param[in] boundaries     pointer to the domain boundaries
 * \param[in] model          type of model related to the NS system
 * \param[in] time_state     state of the time for the NS equations
 * \param[in] algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(const cs_boundary_t           *boundaries,
                          cs_navsto_param_model_t        model,
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
 * \brief  Retrieve the structure storing the parameters for the Navier--Stokes
 *         system
 *
 * \return NULL or the pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_system_get_param(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the equation related to the momentum equation
 *
 * \return NULL or the pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_navsto_system_get_momentum_eq(void);

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
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_finalize_setup(const cs_mesh_t            *mesh,
                                const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant,
                                const cs_time_step_t       *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the settings for SLES related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_set_sles(void);

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
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute_steady_state(const cs_mesh_t        *mesh,
                                      const cs_time_step_t   *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute(const cs_mesh_t         *mesh,
                         const cs_time_step_t    *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Navier-Stokes system
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_extra_op(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq);

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
