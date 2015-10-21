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
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO face-based discretization */
typedef struct _cs_cdofb_codits_t cs_cdofb_codits_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_codits_t structure
 *
 * \param[in]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_codits_init(const cs_equation_param_t  *eqp,
                     const cs_mesh_t            *mesh,
                     const cs_cdo_connect_t     *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_codits_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_codits_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_codits_free(void   *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      m           pointer to a cs_mesh_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in, out] builder     pointer to a cs_cdofb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_compute_source(const cs_mesh_t            *m,
                               const cs_cdo_connect_t     *connect,
                               const cs_cdo_quantities_t  *quant,
                               const cs_time_step_t       *time_step,
                               void                       *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO face-based scheme.
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdofb_codits_t structure
 * \param[in, out] rhs        pointer to a right-hand side array pointer
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_build_system(const cs_mesh_t             *m,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_real_t             *field_val,
                             const cs_time_step_t        *time_step,
                             double                       dt_cur,
                             void                        *builder,
                             cs_real_t                  **rhs,
                             cs_sla_matrix_t            **sla_mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO face-based scheme
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      solu       solution array
 * \param[in, out] builder    pointer to cs_cdofb_codits_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_update_field(const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             const cs_time_step_t       *time_step,
                             const cs_real_t            *solu,
                             void                       *builder,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at each face
 *
 * \param[in]  builder    pointer to a cs_cdofb_codits_t structure
 * \param[in]  field      pointer to a cs_field_t structure
 *
 * \return  a pointer to an array of double (size n_faces)
 */
/*----------------------------------------------------------------------------*/

const double *
cs_cdofb_codits_get_face_values(const void         *builder,
                                const cs_field_t   *field);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_CODITS_H__ */
