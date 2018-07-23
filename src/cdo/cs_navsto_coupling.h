#ifndef __CS_NAVSTO_COUPLING_H__
#define __CS_NAVSTO_COUPLING_H__

/*============================================================================
 * Routines to handle structures used as a context when solving the
 * Navier-Stokes equations. Structures are cast on-the-fly according to the
 * type of coupling.
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

#include "cs_equation.h"
#include "cs_field.h"
#include "cs_param.h"
#include "cs_mesh.h"
#include "cs_navsto_param.h"
#include "cs_time_step.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 * \file cs_navsto_coupling.h
 *
 * \brief  Routines to handle structures used as a context when solving the
 * Navier-Stokes equations. Structures are cast on-the-fly according to the
 * type of coupling.
 */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Predefined context structures depending on the settings */
/* ======================================================= */

/*! \struct cs_navsto_coupling_uzawa_t
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
  cs_real_t       relax;    /*!< Coefficient for the Uzawa algorithm attacehd
                                 to the update of the multiplier */

} cs_navsto_coupling_uzawa_t;

/*! \struct cs_navsto_coupling_ac_t
 *  \brief Set of parameters specific for solving the Navier-Stokes system with
 *         the "artificial compressibility" algorithm
 *
 *  All equations are not always created. It depends on the choice of the model
 */

typedef struct {

  cs_equation_t  *momentum; /*!< Momentum balance equation (vector-valued) */

  cs_property_t  *zeta;     /*!< Coefficient for the artificial compressibility
                                 attached to the grad-div stabilization term */

} cs_navsto_coupling_ac_t;

/*! \struct cs_navsto_coupling_ac_vpp_t
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

} cs_navsto_coupling_ac_vpp_t;

/*! \struct cs_navsto_coupling_projection_t
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

} cs_navsto_coupling_projection_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_COUPLING_H__ */
