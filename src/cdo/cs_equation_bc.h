#ifndef __CS_EQUATION_BC_H__
#define __CS_EQUATION_BC_H__

/*============================================================================
 * Functions to handle the evaluation of boundary conditions when building the
 * algebraic system in CDO/HHO/MAC schemes
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "cdo/cs_cdo_bc.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_local.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_equation_param.h"
#include "cdo/cs_macfb_builder.h"
#include "base/cs_time_step.h"
#include "cdo/cs_xdef_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a boundary condition for a given face (inlet, outlet, wall,
 *        sliding wall, symmetry...)
 *
 * \param[in]       f       face id in the cell mesh numbering
 * \param[in]       eqp     pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm      pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty     pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb      pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys    structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_apply_boundary_t)(short int                     f,
                          const cs_equation_param_t    *eqp,
                          const cs_cell_mesh_t         *cm,
                          const cs_property_data_t     *pty,
                          cs_cell_builder_t            *cb,
                          cs_cell_sys_t                *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enforcement of a boundary condition (Dirichlet, Robin, sliding...)
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_enforce_bc_t)(const cs_equation_param_t      *eqp,
                      const cs_cell_mesh_t           *cm,
                      cs_face_mesh_t                 *fm,
                      cs_hodge_t                     *hodge,
                      cs_cell_builder_t              *cb,
                      cs_cell_sys_t                  *csys);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the values for the normal boundary flux stemming from the Neumann
 *        boundary conditions (zero is left where a Dirichlet is set. This can
 *        be updated later on)
 *
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] values   pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_init_boundary_flux(cs_real_t                     t_eval,
                                  const cs_cdo_quantities_t    *cdoq,
                                  const cs_equation_param_t    *eqp,
                                  cs_real_t                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the BC into a cellwise view of the current system.
 *        Case of vertex-based schemes.
 *
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      face_bc      pointer to a cs_cdo_bc_face_t structure
 * \param[in]      vtx_bc_flag  BC flags associated to vertices
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_cw_vb(const cs_cell_mesh_t         *cm,
                         const cs_equation_param_t    *eqp,
                         const cs_cdo_bc_face_t       *face_bc,
                         const cs_flag_t               vtx_bc_flag[],
                         const cs_real_t               dir_values[],
                         cs_real_t                     t_eval,
                         cs_cell_sys_t                *csys,
                         cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of edge-based schemes
 *
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      face_bc      pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_cw_eb(const cs_cell_mesh_t         *cm,
                         const cs_equation_param_t    *eqp,
                         const cs_cdo_bc_face_t       *face_bc,
                         const cs_real_t               dir_values[],
                         cs_cell_sys_t                *csys,
                         cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of CDO Face-based schemes
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_cw_fb(const cs_cell_mesh_t         *cm,
                         const cs_equation_param_t    *eqp,
                         const cs_cdo_bc_face_t       *face_bc,
                         const cs_real_t               dir_values[],
                         cs_cell_sys_t                *csys,
                         cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of MAC Face-based schemes
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in, out] macb        pointer to a cs_macfb_builder_t strucuture
 * \param[in, out] csys        pointer to a cellwise view of the system
 */
/*----------------------------------------------------------------------------*/

void cs_equation_bc_set_cw_macfb(const cs_cell_mesh_t      *cm,
                                 const cs_equation_param_t *eqp,
                                 const cs_cdo_bc_face_t    *face_bc,
                                 const cs_real_t            dir_values[],
                                 cs_macfb_builder_t        *macb,
                                 cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of Face-based schemes
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_cw_cb(const cs_cell_mesh_t         *cm,
                         const cs_equation_param_t    *eqp,
                         const cs_cdo_bc_face_t       *face_bc,
                         const cs_real_t               dir_values[],
                         cs_cell_sys_t                *csys,
                         cs_cell_builder_t            *cb);


/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define an array of flags for each vertex collecting the flags
 *          of associated boundary faces
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      face_bc   pointer to a structure collecting boundary
 *                           conditions applied to faces
 * \param[in, out] vflag     BC flag on vertices to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_vertex_flag(const cs_cdo_connect_t     *connect,
                               const cs_cdo_bc_face_t     *face_bc,
                               cs_flag_t                  *vflag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define an array of flags for each edge collecting the flags
 *          of associated boundary faces
 *
 * \param[in]      connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      face_bc     pointer to a structure collecting boundary
 *                             conditions applied to faces
 * \param[in, out] edge_flag   BC flag on edges to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_set_edge_flag(const cs_cdo_connect_t     *connect,
                             const cs_cdo_bc_face_t     *face_bc,
                             cs_flag_t                  *edge_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the Dirichlet BCs when DoFs are attached to
 *        vertices
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in, out] bcflag      pointer to an array storing type of BC
 * \param[in, out] values      pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_dirichlet_at_vertices(cs_real_t                   t_eval,
                                     const cs_mesh_t            *mesh,
                                     const cs_cdo_quantities_t  *quant,
                                     const cs_cdo_connect_t     *connect,
                                     const cs_equation_param_t  *eqp,
                                     const cs_cdo_bc_face_t     *face_bc,
                                     cs_flag_t                  *bcflag,
                                     cs_real_t                  *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the Dirichlet BCs at boundary faces.
 *        This can be applied to CDO face-based schemes (DoFs are attached to
 *        primal faces), to CDO cell-based schemes or even to FV schemes.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      face_bc    pointer to a cs_cdo_bc_face_t structure
 * \param[in]      t_eval     time at which one evaluates the boundary cond.
 * \param[in, out] values     pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_dirichlet_at_faces(const cs_mesh_t            *mesh,
                                  const cs_cdo_quantities_t  *quant,
                                  const cs_cdo_connect_t     *connect,
                                  const cs_equation_param_t  *eqp,
                                  const cs_cdo_bc_face_t     *face_bc,
                                  cs_real_t                   t_eval,
                                  cs_real_t                  *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of the Neumann BCs when DoFs are scalar-valued
 *         and attached to a vertex-based schemes (Vb or VCb)
 *         Case of the Neumann BCs i.e. Neumann is defined by a scalar
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing the Neumann values for all DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_svb(cs_real_t                   t_eval,
                                short int                   def_id,
                                short int                   f,
                                const cs_equation_param_t  *eqp,
                                const cs_cell_mesh_t       *cm,
                                double                     *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of the Neumann BCs when DoFs are scalar-valued
 *         and attached to a vertex-based schemes (Vb or VCb)
 *         Case of the full Neumann BCs i.e. Neumann is defined by a vector
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing the Neumann values for all DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_full_neumann_svb(cs_real_t                   t_eval,
                                     short int                   def_id,
                                     short int                   f,
                                     const cs_equation_param_t  *eqp,
                                     const cs_cell_mesh_t       *cm,
                                     double                     *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          the face f.
 *          Case of scalar-valued equation (not full Neumann BCs)
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values for all DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sfb(cs_real_t                    t_eval,
                                short int                    def_id,
                                short int                    f,
                                const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                double                      *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          the face f.
 *          Case of scalar-valued equation with a full Neumann BC definition.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values for all DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_full_neumann_sfb(cs_real_t                    t_eval,
                                     short int                    def_id,
                                     short int                    f,
                                     const cs_equation_param_t   *eqp,
                                     const cs_cell_mesh_t        *cm,
                                     double                      *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of the Neumann BCs at the face f when DoFs are
 *         attached to faces.
 *         Case of vector-valued equation (not the full Neumann)
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_vfb(cs_real_t                    t_eval,
                                short int                    def_id,
                                short int                    f,
                                const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                double                      *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the Robin BCs for a face (cell-wise compute
 *        relying on the cs_cell_mesh_t structure)
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] rob_values  array storing Robin values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_cw_robin(cs_real_t                    t_eval,
                        short int                    def_id,
                        short int                    f,
                        const cs_equation_param_t   *eqp,
                        const cs_cell_mesh_t        *cm,
                        double                      *rob_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the Robin BCs for a face (cell-wise compute
 *        relying on the cs_cell_mesh_t structure)
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      nu          laminar kinematic viscosity
 * \param[in]      k           turbulent kinetic energy
 * \param[in]      hfc         distance from cell center to the wall
 * \param[in, out] rob_values  array storing Robin values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_cw_turb_smooth_wall(cs_real_t                  t_eval,
                                   short int                  def_id,
                                   short int                  f,
                                   const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const double               nu,
                                   const double               k,
                                   const double               hfc,
                                   double                    *rob_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the circulation along primal edges lying on the
 *        domain boundary (the integral of the tangential component of
 *        vector-valued field). This is used for CDO edge-based schemes where
 *        DoFs are attached to (primal) edge-based schemes.
 *
 * \param[in]      t_eval     time at which one evaluates the boundary cond.
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in, out] values     pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_circulation_at_edges(cs_real_t                    t_eval,
                                    const cs_mesh_t             *mesh,
                                    const cs_cdo_quantities_t   *quant,
                                    const cs_cdo_connect_t      *connect,
                                    const cs_equation_param_t   *eqp,
                                    cs_real_t                   *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the boundary conditions to fullfill the constraint when
 *         an incremental solve is set
 *
 * \param[in, out] csys     pointer to the cell system structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_bc_update_for_increment(cs_cell_sys_t  *csys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_BC_H__ */
