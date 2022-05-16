#ifndef __CS_CDOVB_SCALEQ_H__
#define __CS_CDOVB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of scalar-valued equations with source terms
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
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_toolbox.h"
#include "cs_equation_builder.h"
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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the generic structures for building a CDO-Vb scheme are
 *         allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovb_scaleq_is_initialized(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vertex-based schemes. Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_init_sharing(const cs_cdo_quantities_t    *quant,
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
cs_cdovb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_finalize_sharing(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_scaleq_t structure
 *
 * \param[in, out]  scheme_context   pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scaleq_free_context(void   *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag     pointer to an array of BC flag for each vtx
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_setup(cs_real_t                      t_eval,
                      const cs_mesh_t               *mesh,
                      const cs_equation_param_t     *eqp,
                      cs_equation_builder_t         *eqb,
                      cs_flag_t                      vtx_bc_flag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main properties before the main loop on cells.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      t_id      thread id if > 0
 * \param[in]      t_eval    time at which one evaluates BCs
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_init_properties(int                           t_id,
                                cs_real_t                     t_eval,
                                const cs_equation_param_t    *eqp,
                                cs_equation_builder_t        *eqb,
                                void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the cell system for the given cell id when the build occurs
 *         in a coupled system -- block in (row_id, col_id)
 *         Case of scalar-valued CDO-Vb schemes.
 *
 *         Warning: Treatment of BCs differs from the "standard" case.
 *         Up to now, one assumes a Dirichlet or a Neumann for all equations
 *         (i.e. all blocks) and only an algebraic treatment is performed.
 *
 * \param[in]      t_id      thread id if openMP is used
 * \param[in]      c_id      cell id
 * \param[in]      f_val     current field values
 * \param[in]      row_id    id related to the row block
 * \param[in]      col_id    id related to the col block
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme context structure
 * \param[in, out] cb        cell builder structure
 * \param[in, out] csys      cell system structure
 *
 * \return the value of the rhs_norm for the cellwise system
 */
/*----------------------------------------------------------------------------*/

double
cs_cdovb_scaleq_build_block_implicit(int                           t_id,
                                     cs_lnum_t                     c_id,
                                     const cs_real_t               f_val[],
                                     int                           row_id,
                                     int                           col_id,
                                     const cs_equation_param_t    *eqp,
                                     cs_equation_builder_t        *eqb,
                                     void                         *context,
                                     cs_cell_builder_t            *cb,
                                     cs_cell_sys_t                *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Vb schemes.
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
cs_cdovb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         One works cellwise and then process to the assembly.
 *
 *         Variant with an incremental approach. The system is modified to fit
 *         the incremental form (unknows are the increments and rhs corresponds
 *         to a residual). This is useful when the resolution is embedded into
 *         a non-linear process.
 *
 * \param[in]      cur2prev   Not used. Should be done before if needed.
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_steady_state_incr(bool                        cur2prev,
                                        const cs_mesh_t            *mesh,
                                        const int                   field_id,
                                        const cs_equation_param_t  *eqp,
                                        cs_equation_builder_t      *eqb,
                                        void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Implicit time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_implicit(bool                        cur2prev,
                               const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Implicit time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 *         Variant with an incremental approach. The system is modified to fit
 *         the incremental form (unknows are the increments and rhs corresponds
 *         to a residual). This is useful when the resolution is embedded into
 *         a non-linear process.
 *
 * \param[in]      cur2prev   Not used. Should be done before if needed.
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_implicit_incr(bool                        cur2prev,
                                    const cs_mesh_t            *mesh,
                                    const int                   field_id,
                                    const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Theta time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_theta(bool                        cur2prev,
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
 * \return  a pointer to an array of cs_real_t (size: n_vertices)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_vertex_values(void        *context,
                                  bool         previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an array of values at mesh cells by interpolating the
 *         variable field associated to the given context located at mesh
 *         vertices
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
cs_cdovb_scaleq_get_cell_values(void      *context,
                                bool       previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *         Case of scalar-valued CDO vertex-based scheme
 *
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context  pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdovb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cellwise stiffness matrix associated to the property
 *         given as a parameter and apply it to the pot array to define
 *         the resulting array associated to entities defined at loc_res
 *         Case of scalar-valued CDO vertex-based scheme
 *
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context  pointer to a scheme builder structure
 * \param[in]      property pointer to the property related to the stiffness op.
 * \param[in]      pot      array to multiply with the stiffness matrix
 * \param[in]      loc_res  location of entities in the resulting array
 * \param[in, out] res      resulting array
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_apply_stiffness(const cs_equation_param_t     *eqp,
                                cs_equation_builder_t         *eqb,
                                void                          *context,
                                const cs_property_t           *property,
                                const cs_real_t               *pot,
                                cs_flag_t                      loc_res,
                                cs_real_t                     *res);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex of a boundary face, the portion of diffusive
 *         flux across the boundary face. The surface attached to each vertex
 *         corresponds to the intersection of its dual cell (associated to
 *         a vertex of the face) with the face.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       t_eval     time at which one performs the evaluation
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       pdi        pointer to an array of field values
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to a scheme builder structure
 * \param[in, out]  vf_flux    pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                   const cs_equation_param_t   *eqp,
                                   const cs_real_t             *pdi,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context,
                                   cs_real_t                   *vf_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *         Case of scalar-valued CDO-Vb schemes
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
cs_cdovb_scaleq_flux_across_plane(const cs_real_t             normal[],
                                  const cs_real_t            *pdi,
                                  const cs_equation_param_t  *eqp,
                                  int                         ml_id,
                                  cs_equation_builder_t      *eqb,
                                  void                       *context,
                                  double                     *d_flux,
                                  double                     *c_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of an approximation of a constant diffusive
 *         flux (a vector) in each cell.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_diff_flux_in_cells(const cs_real_t             *values,
                                   const cs_equation_param_t   *eqp,
                                   cs_real_t                    t_eval,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context,
                                   cs_real_t                   *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux accross dual faces
 *         (a scalar) in each cell.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_diff_flux_dfaces(const cs_real_t             *values,
                                 const cs_equation_param_t   *eqp,
                                 cs_real_t                    t_eval,
                                 cs_equation_builder_t       *eqb,
                                 void                        *context,
                                 cs_real_t                   *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_SCALEQ_H__ */
