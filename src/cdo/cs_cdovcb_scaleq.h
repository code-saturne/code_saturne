#ifndef __CS_CDOVCB_SCALEQ_H__
#define __CS_CDOVCB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO vertex+cell-based system for unsteady convection
 * diffusion reaction of scalar-valued equations with source terms
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_time_step.h"
#include "cs_mesh.h"
#include "cs_equation_builder.h"
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
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_init_sharing(const cs_cdo_quantities_t    *quant,
                              const cs_cdo_connect_t       *connect,
                              const cs_time_step_t         *time_step);

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
cs_cdovcb_scaleq_finalize_sharing(void);

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
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-VCb schemes.
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
cs_cdovcb_scaleq_init_values(cs_real_t                     t_eval,
                             const int                     field_id,
                             const cs_mesh_t              *mesh,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         equation with a CDO-VCb scheme. Use for interpolation purpose from
 *         cell values to vertex values.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      cell_values  array of cell values
 * \param[in]      field_id     id of the variable field
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in, out] eqb          pointer to a cs_equation_builder_t structure
 * \param[in, out] context      pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_interpolate(const cs_mesh_t            *mesh,
                             const cs_real_t            *cell_values,
                             const int                   field_id,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_steady_state(bool                        cur2prev,
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
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_implicit(bool                        cur2prev,
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
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_theta(bool                        cur2prev,
                             const cs_mesh_t            *mesh,
                             const int                   field_id,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return a pointer to an array of cs_real_t (size: n_vertices)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovcb_scaleq_get_vertex_values(void      *context,
                                   bool       previous);

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
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovcb_scaleq_get_cell_values(void     *context,
                                 bool      previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex of a boundary face, the portion of diffusive
 *         flux across the boundary face. The surface attached to each vertex
 *         corresponds to the intersection of its dual cell (associated to
 *         a vertex of the face) with the face.
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       t_eval     time at which one performs the evaluation
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       pot_v      pointer to an array of field values at vertices
 * \param[in]       pot_c      pointer to an array of field values at cells
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to a scheme builder structure
 * \param[in, out]  vf_flux    pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                    const cs_equation_param_t   *eqp,
                                    const cs_real_t             *pot_v,
                                    const cs_real_t             *pot_c,
                                    cs_equation_builder_t       *eqb,
                                    void                        *context,
                                    cs_real_t                   *vf_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *         Case of scalar-valued CDO-Vcb schemes
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_flux_across_plane(const cs_real_t             normal[],
                                   const cs_real_t            *pdi,
                                   const cs_equation_param_t  *eqp,
                                   int                         ml_id,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context,
                                   double                     *d_flux,
                                   double                     *c_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux in each cells.
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_diff_flux_in_cells(const cs_real_t             *values,
                                    const cs_equation_param_t   *eqp,
                                    cs_real_t                    t_eval,
                                    cs_equation_builder_t       *eqb,
                                    void                        *context,
                                    cs_real_t                   *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across dual faces
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_diff_flux_dfaces(const cs_real_t             *values,
                                  const cs_equation_param_t   *eqp,
                                  cs_real_t                    t_eval,
                                  cs_equation_builder_t       *eqb,
                                  void                        *context,
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
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                     cs_equation_builder_t      *eqb,
                                     void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVCB_SCALEQ_H__ */
