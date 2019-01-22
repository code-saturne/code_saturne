#ifndef __CS_EQUATION_BC_H__
#define __CS_EQUATION_BC_H__

/*============================================================================
 * Routines to handle the evaluation of boundary conditions when building the
 * algebraic system in CDO/HHO schemes
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

#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_xdef_eval.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_local.h"
#include "cs_equation_param.h"
#include "cs_time_step.h"

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
 * \brief  Apply a boundary condition for a given face (inlet, outlet, wall,
 *         sliding wall, symmetry...)
 *
 * \param[in]       f       face id in the cell mesh numbering
 * \param[in]       eqp     pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm      pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb      pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys    structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_apply_boundary_t)(short int                     f,
                          const cs_equation_param_t    *eqp,
                          const cs_cell_mesh_t         *cm,
                          cs_cell_builder_t            *cb,
                          cs_cell_sys_t                *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Enforcement of a boundary condition (Dirichlet, Robin, sliding...)
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_enforce_bc_t)(const cs_equation_param_t      *eqp,
                      const cs_cell_mesh_t           *cm,
                      cs_face_mesh_t                 *fm,
                      cs_cell_builder_t              *cb,
                      cs_cell_sys_t                  *csys);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the values for the normal boundary flux stemming from the
 *         Neumann boundary conditions (zero is left where a Dirichlet is
 *         set. This can be updated later one)
 *
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]       eqp      pointer to a cs_equation_param_t structure
 * \param[in, out]  values   pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_boundary_flux_from_bc(cs_real_t                   t_eval,
                                       const cs_cdo_quantities_t  *cdoq,
                                       const cs_equation_param_t  *eqp,
                                       cs_real_t                  *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of vertex-based schemes
 *
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
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
cs_equation_vb_set_cell_bc(const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
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
 *          Case of Face-based schemes
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_fb_set_cell_bc(const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                           const cs_equation_param_t    *eqp,
                           const cs_cdo_bc_face_t       *face_bc,
                           const cs_real_t               dir_values[],
                           cs_real_t                     t_eval,
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          vertices
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] bcflag      pointer to an array storing type of BC
 * \param[in, out] values      pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_dirichlet_vb(cs_real_t                   t_eval,
                                 const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_face_t     *face_bc,
                                 cs_cell_builder_t          *cb,
                                 cs_flag_t                  *bcflag,
                                 cs_real_t                  *bcvals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          CDO face-based schemes
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      face_bc    pointer to a cs_cdo_bc_face_t structure
 * \param[in]      t_eval     time at which one evaluates the boundary cond.
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_dirichlet_fb(const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_face_t     *face_bc,
                                 cs_real_t                   t_eval,
                                 cs_cell_builder_t          *cb,
                                 cs_real_t                  *values);

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
cs_equation_set_vertex_bc_flag(const cs_cdo_connect_t     *connect,
                               const cs_cdo_bc_face_t     *face_bc,
                               cs_flag_t                  *vflag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to vertices.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing the Neumann values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sv(cs_real_t                   t_eval,
                               short int                   def_id,
                               short int                   f,
                               const cs_cdo_quantities_t  *quant,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               double                     *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          faces.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_fb(cs_real_t                    t_eval,
                               short int                    def_id,
                               short int                    f,
                               const cs_cdo_quantities_t   *quant,
                               const cs_equation_param_t   *eqp,
                               const cs_cell_mesh_t        *cm,
                               double                      *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Robin BCs
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] rob_values  array storing Robin values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_robin(cs_real_t                    t_eval,
                          short int                    def_id,
                          short int                    f,
                          const cs_cdo_quantities_t   *quant,
                          const cs_equation_param_t   *eqp,
                          const cs_cell_mesh_t        *cm,
                          double                      *rob_values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_BC_H__ */
