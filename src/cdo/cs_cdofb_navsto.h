#ifndef __CS_CDOFB_NAVSTO_H__
#define __CS_CDOFB_NAVSTO_H__

/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members for CDO face-based
 *         schemes
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  sma         pointer to a cs_matrix_assembler_t structure (scalar)
 * \param[in]  sms         pointer to a cs_matrix_structure_t structure (scalar)
 * \param[in]  vma         pointer to a cs_matrix_assembler_t structure (vector)
 * \param[in]  vms         pointer to a cs_matrix_structure_t structure (vector)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_common(const cs_cdo_quantities_t     *quant,
                            const cs_cdo_connect_t        *connect,
                            const cs_time_step_t          *time_step,
                            const cs_matrix_assembler_t   *sma,
                            const cs_matrix_structure_t   *sms,
                            const cs_matrix_assembler_t   *vma,
                            const cs_matrix_structure_t   *vms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of a
 *         Uzawa-Augmented Lagrangian approach
 *
 * \param[in] nsp        pointer to a cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_uzawa_context(const cs_navsto_param_t     *nsp,
                                   const void                  *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of an
 *         Artificial Compressibility approach
 *
 * \param[in] nsp    pointer to a cs_navsto_param_t structure
 * \param[in] nsc    pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_ac_context(const cs_navsto_param_t   *nsp,
                                const void                *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of an
 *         incremental Projection approach
 *
 * \param[in] nsp        pointer to a cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_proj_context(const cs_navsto_param_t    *nsp,
                                  const void                 *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_navsto_t structure
 *
 * \param[in]      nsp        pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_free_context(const cs_navsto_param_t      *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         a Uzawa-Lagrangian Augmented approach
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_uzawa_compute(const cs_mesh_t              *mesh,
                              double                        dt_cur,
                              const cs_cdo_connect_t       *connect,
                              const cs_cdo_quantities_t    *quant,
                              const cs_navsto_param_t      *nsp,
                              void                         *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         an Artificial Compressibility approach.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_ac_compute(const cs_mesh_t              *mesh,
                           double                        dt_cur,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                           const cs_navsto_param_t      *nsp,
                           void                         *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         an incremental correction-projection approach
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_proj_compute(const cs_mesh_t              *mesh,
                             double                        dt_cur,
                             const cs_cdo_connect_t       *connect,
                             const cs_cdo_quantities_t    *quant,
                             const cs_navsto_param_t      *nsp,
                             void                         *nsc_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdofb_navsto_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_update_field(const cs_real_t              *solu,
                             const cs_real_t              *rhs,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *data,
                             cs_real_t                    *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_navsto_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_NAVSTO_H__ */
