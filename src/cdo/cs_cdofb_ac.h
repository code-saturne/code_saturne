#ifndef __CS_CDOFB_AC_H__
#define __CS_CDOFB_AC_H__

/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an artificial compressibility algorithm
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
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_source_term.h"
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
 * Inline static function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values of the velocity on the faces
 *
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of cs_real_t (size = 3*n_faces)
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t *
cs_cdofb_ac_get_face_velocity(bool   previous)
{
  return cs_equation_get_face_values(cs_equation_by_name("momentum"),
                                     previous);
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
cs_cdofb_ac_init_common(const cs_cdo_quantities_t     *quant,
                        const cs_cdo_connect_t        *connect,
                        const cs_time_step_t          *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_ac_t structure
 *
 * \param[in] nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in] adv_field   pointer to \ref cs_adv_field_t structure
 * \param[in] mflux       current values of the mass flux across primal faces
 * \param[in] mflux_pre   current values of the mass flux across primal faces
 * \param[in] fb_type     type of boundary for each boundary face
 * \param[in] nsc_input   pointer to a \ref cs_navsto_ac_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_ac_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_init_scheme_context(const cs_navsto_param_t   *nsp,
                                cs_adv_field_t            *adv_field,
                                cs_real_t                 *mflux,
                                cs_real_t                 *mflux_pre,
                                cs_boundary_type_t        *fb_type,
                                void                      *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_ac_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_free_scheme_context(void   *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an AC algorithm
 *         is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_set_sles(const cs_navsto_param_t    *nsp,
                     void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and an implicit Euler
 *         time scheme
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute_implicit(const cs_mesh_t              *mesh,
                             const cs_navsto_param_t      *nsp,
                             void                         *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and an implicit Euler
 *         time scheme
 *         For Picard - Navier--Stokes problems
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute_implicit_nl(const cs_mesh_t              *mesh,
                                const cs_navsto_param_t      *nsp,
                                void                         *scheme_context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_AC_H__ */
