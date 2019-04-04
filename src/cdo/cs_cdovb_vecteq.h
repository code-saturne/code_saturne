#ifndef __CS_CDOVB_VECTEQ_H__
#define __CS_CDOVB_VECTEQ_H__

/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of vector-valued equations with source terms
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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Vb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovb_vecteq_is_initialized(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vector-valued vertex-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_init_common(const cs_cdo_quantities_t    *quant,
                            const cs_cdo_connect_t       *connect,
                            const cs_time_step_t         *time_step,
                            const cs_matrix_structure_t  *ms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_finalize_common(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_vecteq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_vecteq_free_context(void   *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of vector-valued CDO-Vb schemes.
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
cs_cdovb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to cs_cdovb_vecteq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_initialize_system(const cs_equation_param_t  *eqp,
                                  cs_equation_builder_t      *eqb,
                                  void                       *data,
                                  cs_matrix_t               **system_matrix,
                                  cs_real_t                 **system_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a vector CDO vertex-based scheme.
 *
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out] context     pointer to the scheme context (cast on-the-fly)
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_set_dir_bc(cs_real_t                     t_eval,
                           const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           void                         *context,
                           cs_real_t                     field_val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovb_vecteq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_update_field(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_vertex_values(void      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an array of values at mesh cells by interpolating the
 *         variable field associated to the given context located at mesh
 *         vertices
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_cell_values(void      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_VECTEQ_H__ */
