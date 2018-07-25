#ifndef __CS_HHO_VECTEQ_H__
#define __CS_HHO_VECTEQ_H__

/*============================================================================
 * Build an algebraic system for vector conv./diff. eq. with Hybrid High Order
 * space discretization
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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_hho_builder.h"
#include "cs_matrix.h"
#include "cs_matrix_assembler.h"
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

/* Algebraic system for HHO discretization */
typedef struct _cs_hho_vecteq_t cs_hho_vecteq_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to HHO schemes
 *         Set shared pointers
 *
 * \param[in]  scheme_flag  flag to identify which kind of numerical scheme is
 *                          requested to solve the computational domain
 * \param[in]  quant        additional mesh quantities struct.
 * \param[in]  connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  ms0          pointer to a cs_matrix_structure_t structure (P0)
 * \param[in]  ms1          pointer to a cs_matrix_structure_t structure (P1)
 * \param[in]  ms2          pointer to a cs_matrix_structure_t structure (P2)
*/
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_init_common(cs_flag_t                      scheme_flag,
                         const cs_cdo_quantities_t     *quant,
                         const cs_cdo_connect_t        *connect,
                         const cs_time_step_t          *time_step,
                         const cs_matrix_structure_t   *ms0,
                         const cs_matrix_structure_t   *ms1,
                         const cs_matrix_structure_t   *ms2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys    pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb      pointer to a pointer on a cs_cell_builder_t structure
 * \param[out]  hhob    pointer to a pointer on a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_get(cs_cell_sys_t       **csys,
                  cs_cell_builder_t   **cb,
                  cs_hho_builder_t    **hhob);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers and generic structures related to HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_finalize_common(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_hho_vecteq_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux__id  id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated \ref cs_hho_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_vecteq_init_context(const cs_equation_param_t   *eqp,
                           int                          var_id,
                           int                          bflux_id,
                           cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_vecteq_t structure
 *
 * \param[in, out]  data    pointer to a cs_hho_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_vecteq_free_context(void   *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a cs_equation_builder_t structure
 * \param[in, out] data     pointer to a cs_hho_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_compute_source(const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         data structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to generic data structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of cs_real_t
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_initialize_system(const cs_equation_param_t  *eqp,
                                cs_equation_builder_t      *eqb,
                                void                       *data,
                                cs_matrix_t               **system_matrix,
                                cs_real_t                 **system_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a HHO scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_vecteq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_build_system(const cs_mesh_t            *mesh,
                           const cs_real_t            *field_val,
                           double                      dt_cur,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *rhs,
                           cs_matrix_t                *matrix);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values required for hybrid discretization
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_vecteq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_update_field(const cs_real_t            *solu,
                           const cs_real_t            *rhs,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at faces (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in, out]  data    pointer to a data structure cast-on-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_hho_vecteq_get_face_values(void          *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at cells (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in, out]  data    pointer to a data structure cast-on-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_hho_vecteq_get_cell_values(void          *data);

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
cs_hho_vecteq_read_restart(cs_restart_t    *restart,
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
cs_hho_vecteq_write_restart(cs_restart_t    *restart,
                            const char      *eqname,
                            void            *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_hho_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_vecteq_extra_op(const char                 *eqname,
                       const cs_field_t           *field,
                       const cs_equation_param_t  *eqp,
                       cs_equation_builder_t      *eqb,
                       void                       *data);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HHO_VECTEQ_H__ */
