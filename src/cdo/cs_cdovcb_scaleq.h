#ifndef __CS_CDOVCB_SCALEQ_H__
#define __CS_CDOVCB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO vertex+cell-based system for unsteady convection
 * diffusion reaction of scalar-valued equations with source terms
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
#include "cs_time_step.h"
#include "cs_mesh.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
#include "cs_restart.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO vertex-based discretization */
typedef struct _cs_cdovcb_scaleq_t cs_cdovcb_scaleq_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-vertex+cell
 *           based scheme are allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovcb_scaleq_is_initialized(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vertex+cell-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_init_common(const cs_cdo_quantities_t    *quant,
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
cs_cdovcb_scaleq_get(cs_cell_sys_t       **csys,
                     cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers and generic structures related to CDO vertex+cell-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_finalize_common(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovcb_scaleq_t structure storing data useful
 *         for building and  managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovcb_scaleq_init_context(const cs_equation_param_t   *eqp,
                              int                          var_id,
                              int                          bflux_id,
                              cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovcb_scaleq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdovcb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovcb_scaleq_free_context(void   *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         data structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] context        pointer to cs_cdovcb_scaleq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_initialize_system(const cs_equation_param_t    *eqp,
                                   cs_equation_builder_t        *eqb,
                                   void                         *context,
                                   cs_matrix_t                **system_matrix,
                                   cs_real_t                  **system_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a scalar CDO vertex+cell-based scheme.
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
cs_cdovcb_scaleq_set_dir_bc(cs_real_t                     t_eval,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context,
                            cs_real_t                     field_val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_steady_state(double                      dt_cur,
                                    const cs_mesh_t            *mesh,
                                    const int                   field_id,
                                    const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         Time scheme is an implicit Euler
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_implicit(double                      dt_cur,
                                const cs_mesh_t            *mesh,
                                const int                   field_id,
                                const cs_equation_param_t  *eqp,
                                cs_equation_builder_t      *eqb,
                                void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         Time scheme is a theta scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_theta(double                      dt_cur,
                             const cs_mesh_t            *mesh,
                             const int                   field_id,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex+cell-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the vertex field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovcb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_build_system(const cs_mesh_t            *mesh,
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
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to data structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_update_field(const cs_real_t            *solu,
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
cs_cdovcb_scaleq_get_vertex_values(void      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at mesh cells from the inverse operation
 *         w.r.t. the static condensation (DoF used in the linear system are
 *         located at primal vertices and field related to the structure
 *         equation is also attached to primal vertices)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovcb_scaleq_get_cell_values(void     *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_compute_flux_across_plane(const cs_real_t             normal[],
                                           const cs_real_t            *pdi,
                                           int                         ml_id,
                                           const cs_equation_param_t  *eqp,
                                           cs_equation_builder_t      *eqb,
                                           void                       *data,
                                           double                     *d_flux,
                                           double                     *c_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_cellwise_diff_flux(const cs_real_t             *values,
                                    const cs_equation_param_t   *eqp,
                                    cs_real_t                    t_eval,
                                    cs_equation_builder_t       *eqb,
                                    void                        *context,
                                    cs_flag_t                    location,
                                    cs_real_t                   *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]       v_values    discrete values for the potential at vertices
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  v_gradient  gradient at vertices
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_vtx_gradient(const cs_real_t         *v_values,
                              cs_equation_builder_t   *eqb,
                              void                    *context,
                              cs_real_t               *v_gradient);

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
cs_cdovcb_scaleq_read_restart(cs_restart_t    *restart,
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
cs_cdovcb_scaleq_write_restart(cs_restart_t    *restart,
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
 * \param[in, out]  context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_extra_op(const char                 *eqname,
                          const cs_field_t           *field,
                          const cs_equation_param_t  *eqp,
                          cs_equation_builder_t      *eqb,
                          void                       *context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVCB_SCALEQ_H__ */
