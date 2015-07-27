#ifndef __CS_CDOFB_CODITS_H__
#define __CS_CDOFB_CODITS_H__

/*============================================================================
 * Build an algebraic CDO face-based system for convection/diffusion equation
 * with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include "cs_param_eq.h"
#include "cs_cdo_bc.h"
#include "cs_hodge.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO face-based discretization */
typedef struct _cdofb_codits_t cs_cdofb_codits_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the required number of scalar equations based on a face
 *         based discretization
 *
 * \param[in]    n_scal_systems   number of scalar equations
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_create_all(int  n_scal_systems);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_codits_t
 *
 * \param[in] eq       pointer to a structure storing parameters of an eq.
 * \param[in] m        pointer to a mesh structure
 * \param[in] eq_id    id related to the equation to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_init(const cs_param_eq_t  *eq,
                     const cs_mesh_t      *m,
                     int                   eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_cdofb_codits_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a scalar convection/diffusion equation with a CDO
 *         face-based scheme.
 *
 * \param[in] m        pointer to a cs_mesh_t structure
 * \param[in] connect  pointer to a cs_cdo_connect_t structure
 * \param[in] quant    pointer to a cs_cdo_quantities_t structure
 * \param[in] tcur     current physical time of the simulation
 * \param[in] sys_id   pointer to a cs_cdofb_codits_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_solve(const cs_mesh_t            *m,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      double                      tcur,
                      int                         sys_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO face-based scheme
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]  sys_id    id of the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_post(const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     int                         sys_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at each face
 *
 * \param[in]  eq_id    id related to a cs_param_eq_t struct.
 *
 * \return  a pointer to an array of double (size n_faces)
 */
/*----------------------------------------------------------------------------*/

const double *
cs_cdofb_codits_get_face_values(int     eq_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_CODITS_H__ */
