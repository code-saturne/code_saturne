#ifndef __CS_NAVSTO_COUPLING_H__
#define __CS_NAVSTO_COUPLING_H__

/*============================================================================
 * Routines to handle structures used as a context when solving the
 * Navier-Stokes equations. Structures are cast on-the-fly according to the
 * type of coupling.
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

#include "cs_equation.h"
#include "cs_field.h"
#include "cs_navsto_param.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_navsto_coupling.h

  \brief  Routines to handle structures used as a context when solving the
          Navier-Stokes equations.
          Structures are cast on-the-fly according to the type of coupling.
          Routines to handle the settings of coupling algorithms
          - Artificial Compressibility algorithm
          - Its variant VVP (Vector Projection Penalty) algorithm
          - Monolithic algorithm
          - Projection algorithm
          - Uzawa-Augmented Lagrangian algorithm
 */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Predefined context structures depending on the settings */
/* ======================================================= */

/*! \struct cs_navsto_ac_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         the "artificial compressibility" algorithm
 *
 *  All equations are not always created. It depends on the choice of the model
 */

typedef struct {

  cs_equation_t  *momentum; /*!< Momentum balance equation (vector-valued) */

  cs_property_t  *zeta;     /*!< Coefficient for the artificial compressibility
                                 attached to the grad-div stabilization term */

} cs_navsto_ac_t;

/*! \struct cs_navsto_ac_vpp_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         the "artificial compressibility" solved by the VPP_eps algorithm
 *
 *  All equations are not always created. It depends on the choice of the model
 */

typedef struct {

  cs_equation_t  *momentum; /*!< Momentum balance equation (vector-valued) */
  cs_equation_t  *graddiv;  /*!< Second equation of the VPP_eps method, that is
                                 where the grad-div operator is used
                                 (vector-valued) */

  cs_property_t  *zeta;    /*!< Parameter (Artificial Compressibility) VPP
                                algorithm attached to the grad-div stabilization
                                term */

} cs_navsto_ac_vpp_t;

/*! \struct cs_navsto_monolithic_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         a fully coupled monolithic algorithm
 */

typedef struct {

  cs_equation_t  *momentum; /*!< Momentum equation (vector-valued) */

} cs_navsto_monolithic_t;

/*! \struct cs_navsto_projection_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         an incremental projection algorithm
 *
 *  All equations are not always created. It depends on the choice of the model
 */

typedef struct {

  cs_equation_t  *prediction; /*!< Velocity prediction step related to the
                                   momentum balance equation (vector-valued) */
  cs_equation_t  *correction; /*!< Pressure correction step related to the mass
                                   balance equation (scalar-valued) */

  /*! \var div_st
   * Source term on the correction step stemming from the divergence of the
   * predicted velocity */
  cs_real_t      *div_st;

  /*! \var bdy_pressure_incr
   * Pressure increment at the boundary. Used as an array to set the boundary
   * condition arising from a Dirichlet on the pressure. */
  cs_real_t      *bdy_pressure_incr;

  /*! \var predicted_velocity
   * Predicted velocity field (value of the velocity at cells after the
   * prediction step). This values may not be divergence-free.
   */
  cs_field_t     *predicted_velocity;

} cs_navsto_projection_t;

/*! \struct cs_navsto_uzawa_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         a fully coupled algorithm using a Uzawa algorithm and an Augmented
 *         Lagrangian approach inside each sub-iteration.
 *
 *  All equations are not always created. It depends on the choice of the model.
 */

typedef struct {

  cs_equation_t  *momentum; /*!< Momentum balance equation (vector-valued) */
  cs_equation_t  *energy;   /*!< Energy balance equation (scalar-valued) */

  cs_property_t  *zeta;     /*!< Coefficient for the augmented Lagrangian
                                 attached to the grad-div stabilzation term */

} cs_navsto_uzawa_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Uzawa-Augmented Lagrangian approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_uzawa_create_context(cs_navsto_param_t    *nsp,
                               cs_param_bc_type_t    bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Uzawa-Augmented Lagrangian
 *         approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_uzawa_free_context(const cs_navsto_param_t    *nsp,
                             void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_uzawa_init_setup(const cs_navsto_param_t    *nsp,
                           void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_uzawa_last_setup(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_navsto_param_t     *nsp,
                           void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Artificial Compressibility approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_create_context(cs_navsto_param_t    *nsp,
                            cs_param_bc_type_t    bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_free_context(const cs_navsto_param_t    *nsp,
                          void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_init_setup(const cs_navsto_param_t    *nsp,
                        void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_last_setup(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_navsto_param_t     *nsp,
                        void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Artificial Compressibility - VPP approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_vpp_create_context(cs_navsto_param_t    *nsp,
                                cs_param_bc_type_t    bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         with the VPP approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_vpp_free_context(const cs_navsto_param_t    *nsp,
                              void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an Artificial
 *         Compressibility with VPP algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_vpp_init_setup(const cs_navsto_param_t    *nsp,
                            void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_vpp_last_setup(const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            const cs_navsto_param_t     *nsp,
                            void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using a monolithic approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_monolithic_create_context(cs_navsto_param_t    *nsp,
                                    cs_param_bc_type_t    bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a monolithic approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_monolithic_free_context(const cs_navsto_param_t    *nsp,
                                  void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_monolithic_init_setup(const cs_navsto_param_t    *nsp,
                                void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a monolithic
 *         algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_monolithic_last_setup(const cs_cdo_connect_t      *connect,
                                const cs_cdo_quantities_t   *quant,
                                const cs_navsto_param_t     *nsp,
                                void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an incremental Projection approach in the
 *         the rotational form (see Minev & Guermond, 2006, JCP)
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_create_context(cs_navsto_param_t    *nsp,
                                    cs_param_bc_type_t    bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a Projection approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_free_context(const cs_navsto_param_t    *nsp,
                                  void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a projection
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in]      loc_id        id related to a mesh location
 * \param[in]      has_previous  values at different time steps (true/false)
 * \param[in, out] context       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_init_setup(const cs_navsto_param_t    *nsp,
                                int                         loc_id,
                                _Bool                       has_previous,
                                void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a
 *         projection algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_last_setup(const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant,
                                const cs_navsto_param_t    *nsp,
                                void                       *context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_COUPLING_H__ */
