#ifndef __CS_CDOEB_SCALEQ_H__
#define __CS_CDOEB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO edge-based system. Degrees of freedom are defined as
 * a circulation. Degrees of freedom are scalar-valued but the equation to
 * solve is vector-valued
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_cdoeb_priv.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_restart.h"
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

/* Algebraic system for CDO edge-based discretization */
typedef struct _cs_cdoeb_t cs_cdoeb_vecteq_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Eb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdoeb_vecteq_is_initialized(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffers and general structures related to CDO
 *           edge-based schemes. Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_init_common(const cs_cdo_quantities_t    *quant,
                            const cs_cdo_connect_t       *connect,
                            const cs_time_step_t         *time_step,
                            const cs_matrix_structure_t  *ms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *         in case of scalar-valued edge-based scheme
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO edge-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_finalize_common(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdoeb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdoeb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdoeb_vecteq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdoeb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdoeb_vecteq_free_context(void   *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Eb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Eb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh edges (the DoFs)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context     pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_vecteq_get_edge_values(void      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh cells from a reconstruction of edge values.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context     pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_vecteq_get_cell_values(void      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_read_restart(cs_restart_t    *restart,
                             const char      *eqname,
                             void            *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOEB_VECTEQ_H__ */
