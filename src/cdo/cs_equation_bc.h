#ifndef __CS_EQUATION_BC_H__
#define __CS_EQUATION_BC_H__

/*============================================================================
 * Routines to handle the evaluation of boundary conditions when building the
 * algebraic system in CDO/HHO schemes
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of Face-based schemes
 *
 * \param[in]      bf_id       id of the border face
 * \param[in]      f           id of the current face in a cellwise numbering
 * \param[in]      face_flag   metadata about the current face
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    Definition id related to each Neumann face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_vb_set_cell_bc(cs_lnum_t                     bf_id,
                           short int                     f,
                           cs_flag_t                     face_flag,
                           const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                            const cs_time_step_t        *time_step,
                           const cs_equation_param_t    *eqp,
                           const cs_real_t               dir_values[],
                           const short int               neu_tags[],
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of Face-based schemes
 *
 * \param[in]      bf_id       id of the border face
 * \param[in]      f           id of the current face in a cellwise numbering
 * \param[in]      face_flag   metadata about the current face
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    Definition id related to each Neumann face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_fb_set_cell_bc(cs_lnum_t                     bf_id,
                           short int                     f,
                           cs_flag_t                     face_flag,
                           const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                            const cs_time_step_t        *time_step,
                           const cs_equation_param_t    *eqp,
                           const cs_real_t               dir_values[],
                           const short int               neu_tags[],
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are scalar-valued
 *          and attached to vertices
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step   pointer to a time step structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      dir         pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 *
 * \return a pointer to a new allocated array storing the dirichlet values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_compute_dirichlet_sv(const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_time_step_t       *time_step,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_list_t     *dir,
                                 cs_cell_builder_t          *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          CDO face-based schemes
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      dir        pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 *
 * \return a pointer to a new allocated array storing the dirichlet values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_compute_dirichlet_fb(const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_time_step_t       *time_step,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_list_t     *dir,
                                 cs_cell_builder_t          *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Tag each face related to a Neumann BC with its definition id.
 *          Default tag is -1 (not a Neumann face)
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t

 * \return an array with prescribed tags
 */
/*----------------------------------------------------------------------------*/

short int *
cs_equation_tag_neumann_face(const cs_cdo_quantities_t    *quant,
                             const cs_equation_param_t    *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to vertices.
 *
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing the Neumann values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sv(short int                   def_id,
                               short int                   f,
                               const cs_cdo_quantities_t  *quant,
                               const cs_time_step_t       *time_step,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               double                     *neu_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          faces.
 *
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_fb(short int                    def_id,
                               short int                    f,
                               const cs_cdo_quantities_t   *quant,
                               const cs_time_step_t        *time_step,
                               const cs_equation_param_t   *eqp,
                               const cs_cell_mesh_t        *cm,
                               double                      *neu_values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_BC_H__ */
