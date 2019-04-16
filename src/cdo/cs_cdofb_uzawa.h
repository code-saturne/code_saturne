#ifndef __CS_CDOFB_UZAWA_H__
#define __CS_CDOFB_UZAWA_H__

/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an Augmented Lagrangian-Uzawa algorithm
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_mesh.h"
#include "cs_navsto_param.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values of the velocity on the faces
 *
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 *
 * \return a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t *
cs_cdofb_uzawa_get_face_velocity(void    *scheme_context)
{
  CS_UNUSED(scheme_context);

  return cs_equation_get_face_values(cs_equation_by_name("momentum"));
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_init_common(const cs_cdo_quantities_t     *quant,
                           const cs_cdo_connect_t        *connect,
                           const cs_time_step_t          *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_cdofb_uzawa_t structure
 *
 * \param[in] nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in] fb_type    type of boundary for each boundary face
 * \param[in] nsc_input  pointer to a \ref cs_navsto_uzawa_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_uzawa_init_scheme_context(const cs_navsto_param_t    *nsp,
                                   cs_boundary_type_t         *fb_type,
                                   void                       *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_uzawa_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_uzawa_free_scheme_context(void   *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an ALU algorithm
 *         is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_set_sles(const cs_navsto_param_t    *nsp,
                        void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach.
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_steady(const cs_mesh_t              *mesh,
                              const cs_navsto_param_t      *nsp,
                              void                         *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach and an Euler time scheme
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_implicit(const cs_mesh_t              *mesh,
                                const cs_navsto_param_t      *nsp,
                                void                         *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach and a theta time scheme
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_theta(const cs_mesh_t              *mesh,
                             const cs_navsto_param_t      *nsp,
                             void                         *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach. It builds the matrix
 *         at each iteration
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_steady_rebuild(const cs_mesh_t         *mesh,
                                      const cs_navsto_param_t *nsp,
                                      void                    *scheme_context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_UZAWA_H__ */
