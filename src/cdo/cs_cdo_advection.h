#ifndef __CS_CDO_ADVECTION_H__
#define __CS_CDO_ADVECTION_H__

/*============================================================================
 * Build discrete convection operators for CDO schemes
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
#include "cs_property.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* ---------------------------------------------------------------------------
 * Function pointers for CDO face-based schemes
 * -------------------------------------------------------------------------- */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform preprocessing such as the computation of the advection flux
 *         at the expected location in order to be able to build the advection
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      csys     pointer to a cs_cell_sys_t structure
 * \param[in, out] input    NULL or pointer to a structure cast on-the-fly
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_adv_open_hook_t)(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           const cs_cell_sys_t         *csys,
                           void                        *input,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator in CDO-Fb schemes. Case of an
 *          operator that should be assemble into a matrix. Boundary conditions
 *          are also handled in this function.
 *
 * \param[in]      dim     dimension of the variable (1 or 3)
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      csys    pointer to a cs_cell_sys_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_adv_scheme_t)(int                        dim,
                        const cs_cell_mesh_t      *cm,
                        const cs_cell_sys_t       *csys,
                        cs_cell_builder_t         *cb,
                        cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the cellwise advection operator for CDO-Fb schemes
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          This pointer of function manages how is build the advection term
 *          in a cell and it relies on the lower-level function
 *          \ref cs_cdofb_adv_scheme_t
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      csys         pointer to a cs_cell_sys_t structure
 * \param[in]      scheme_func  function pointer to the scheme definition
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_adv_build_t)(const cs_equation_param_t   *eqp,
                       const cs_cell_mesh_t        *cm,
                       const cs_cell_sys_t         *csys,
                       cs_cdofb_adv_scheme_t       *scheme_func,
                       cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Operation done after the matrix related to the advection term has
 *          been defined.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to the local advection matrix
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_adv_close_hook_t)(const cs_equation_param_t   *eqp,
                            const cs_cell_mesh_t        *cm,
                            cs_cell_sys_t               *csys,
                            cs_cell_builder_t           *cb,
                            cs_sdm_t                    *adv);

/* ---------------------------------------------------------------------------
 * Function pointers for CDO vertex-based schemes
 * -------------------------------------------------------------------------- */

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_advection_t)(const cs_equation_param_t   *eqp,
                       const cs_cell_mesh_t        *cm,
                       const cs_property_data_t    *diff_pty,
                       cs_face_mesh_t              *fm,
                       cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the advection operator in CDO
 *          vertex-based (or vertex+cell-based) schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_advection_bc_t)(const cs_cell_mesh_t       *cm,
                          const cs_equation_param_t  *eqp,
                          cs_real_t                   t_eval,
                          cs_face_mesh_t             *fm,
                          cs_cell_builder_t          *cb,
                          cs_cell_sys_t              *csys);

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the stabilization coefficient used in CIP scheme
 *
 * \param[in]  new_value     value of the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_set_cip_coef(double     new_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the value of the stabilization coefficient used in CIP scheme
 *
 * \return   the value the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_advection_get_cip_coef(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform preprocessing such as the computation of the advection flux
 *         at the expected location in order to be able to build the advection
 *         matrix. Follow the prototype given by cs_cdofb_adv_open_hook_t
 *         Default case.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      csys     pointer to a cs_cell_sys_t structure
 * \param[in, out] input    NULL or pointer to a structure cast on-the-fly
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_open_default(const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                const cs_cell_sys_t         *csys,
                                void                        *input,
                                cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_cdofb_adv_close_hook_t
 *         Default scalar-valued case.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to the local advection matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_close_default_scal(const cs_equation_param_t   *eqp,
                                      const cs_cell_mesh_t        *cm,
                                      cs_cell_sys_t               *csys,
                                      cs_cell_builder_t           *cb,
                                      cs_sdm_t                    *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_cdofb_adv_close_hook_t
 *         Default vector-valued case.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to the local advection matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_close_default_vect(const cs_equation_param_t   *eqp,
                                      const cs_cell_mesh_t        *cm,
                                      cs_cell_sys_t               *csys,
                                      cs_cell_builder_t           *cb,
                                      cs_sdm_t                    *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_cdofb_adv_close_hook_t
 *         Explicit treatment without extrapolation for scalar-valued DoFs.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to the local advection matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_close_exp_none_scal(const cs_equation_param_t   *eqp,
                                       const cs_cell_mesh_t        *cm,
                                       cs_cell_sys_t               *csys,
                                       cs_cell_builder_t           *cb,
                                       cs_sdm_t                    *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_cdofb_adv_close_hook_t
 *         Explicit treatment without extrapolation for vector-valued DoFs.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to the local advection matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_close_exp_none_vect(const cs_equation_param_t   *eqp,
                                       const cs_cell_mesh_t        *cm,
                                       cs_cell_sys_t               *csys,
                                       cs_cell_builder_t           *cb,
                                       cs_sdm_t                    *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main function to build the cellwise advection operator for CDO-Fb
 *         schemes The local matrix related to this operator is stored in
 *         cb->loc
 *
 *         Case of an advection term without a diffusion operator. In this
 *         situation, a numerical issue may arise if an internal or a border
 *         face is such that there is no advective flux. A specil treatment
 *         is performed to tackle this issue.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      csys         pointer to a cs_cell_sys_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_no_diffusion(const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                const cs_cell_sys_t         *csys,
                                cs_cdofb_adv_scheme_t       *scheme_func,
                                cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main function to build the cellwise advection operator for CDO
 *          face-based schemes.
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          One assumes that a diffusion term is present so that there is no
 *          need to perform additional checkings on the well-posedness of the
 *          operator.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      csys         pointer to a cs_cell_sys_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection(const cs_equation_param_t   *eqp,
                   const cs_cell_mesh_t        *cm,
                   const cs_cell_sys_t         *csys,
                   cs_cdofb_adv_scheme_t       *scheme_func,
                   cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme
 *         - non-conservative formulation beta \cdot \grad
 *         - upwind scheme
 *         Rely on the work performed during R. Milani's PhD
 *
 *         A scalar-valued version is built. Only the enforcement of the
 *         boundary condition depends on the variable dimension.
 *         Remark: Usually the local matrix called hereafter adv is stored
 *         in cb->loc
 *
 * \param[in]      dim     dimension of the variable (1 or 3)
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      csys    pointer to a cs_cell_sys_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_upwnoc(int                        dim,
                          const cs_cell_mesh_t      *cm,
                          const cs_cell_sys_t       *csys,
                          cs_cell_builder_t         *cb,
                          cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme
 *         - conservative formulation div(\beta )
 *         - upwind scheme
 *         Rely on the work performed during R. Milani's PhD
 *
 *         A scalar-valued version is built. Only the enforcement of the
 *         boundary condition depends on the variable dimension.
 *         Remark: Usually the local matrix called hereafter adv is stored
 *         in cb->loc
 *
 * \param[in]      dim     dimension of the variable (1 or 3)
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      csys    pointer to a cs_cell_sys_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_upwcsv(int                        dim,
                          const cs_cell_mesh_t      *cm,
                          const cs_cell_sys_t       *csys,
                          cs_cell_builder_t         *cb,
                          cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme
 *         - non-conservative formulation beta.grad
 *         - centered scheme
 *         Rely on the work performed during R. Milani's PhD
 *
 *         A scalar-valued version is built. Only the enforcement of the
 *         boundary condition depends on the variable dimension.
 *         Remark: Usually the local matrix called hereafter adv is stored
 *         in cb->loc
 *
 * \param[in]      dim     dimension of the variable (1 or 3)
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      csys    pointer to a cs_cell_sys_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_cennoc(int                        dim,
                          const cs_cell_mesh_t      *cm,
                          const cs_cell_sys_t       *csys,
                          cs_cell_builder_t         *cb,
                          cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme
 *         - conservative formulation div(beta )
 *         - centered scheme
 *         Rely on the work performed during R. Milani's PhD
 *
 *         A scalar-valued version is built. Only the enforcement of the
 *         boundary condition depends on the variable dimension.
 *         Remark: Usually the local matrix called hereafter adv is stored
 *         in cb->loc
 *
 * \param[in]      dim     dimension of the variable (1 or 3)
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      csys    pointer to a cs_cell_sys_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_cencsv(int                        dim,
                          const cs_cell_mesh_t      *cm,
                          const cs_cell_sys_t       *csys,
                          cs_cell_builder_t         *cb,
                          cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme with an upwind scheme and a conservative
 *          formulation. The portion of upwinding relies on an evaluation
 *          of the weigth associated to the given property
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv_wpty(const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                const cs_property_data_t    *diff_pty,
                                cs_face_mesh_t              *fm,
                                cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion and an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           const cs_property_data_t    *diff_pty,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cencsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           const cs_property_data_t    *diff_pty,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a mixed centered/upwind scheme with
 *          a conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_mcucsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           const cs_property_data_t    *diff_pty,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme with an upwind scheme and a conservative
 *          formulation. The portion of upwinding relies on an evaluation
 *          of the weigth associated to the given property.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc_wpty(const cs_equation_param_t   *eqp,
                                const cs_cell_mesh_t        *cm,
                                const cs_property_data_t    *diff_pty,
                                cs_face_mesh_t              *fm,
                                cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion when an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           const cs_property_data_t    *diff_pty,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a non-conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *          Predefined prototype to match the function pointer
 *          cs_cdovb_advection_t
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to a \ref cs_property_data_t structure
 * \param[in, out] fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cennoc(const cs_equation_param_t    *eqp,
                           const cs_cell_mesh_t         *cm,
                           const cs_property_data_t     *diff_pty,
                           cs_face_mesh_t               *fm,
                           cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme when the advection field is cellwise
 *          constant
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to the property associated to diffusion
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_cw_cst(const cs_equation_param_t   *eqp,
                            const cs_cell_mesh_t        *cm,
                            const cs_property_data_t    *diff_pty,
                            cs_face_mesh_t              *fm,
                            cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_pty  pointer to the property associated to diffusion
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb(const cs_equation_param_t   *eqp,
                     const cs_cell_mesh_t        *cm,
                     const cs_property_data_t    *diff_pty,
                     cs_face_mesh_t              *fm,
                     cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_bc(const cs_cell_mesh_t       *cm,
                       const cs_equation_param_t  *eqp,
                       cs_real_t                   t_eval,
                       cs_face_mesh_t             *fm,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator with CDO
 *          V+C schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_bc(const cs_cell_mesh_t        *cm,
                        const cs_equation_param_t   *eqp,
                        cs_real_t                    t_eval,
                        cs_face_mesh_t              *fm,
                        cs_cell_builder_t           *cb,
                        cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value of the upwinding coefficient in each cell knowing
 *          the related Peclet number
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in]      scheme    type of scheme used for the advection term
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_cell_upwind_coef(const cs_cdo_quantities_t    *cdoq,
                                  cs_param_advection_scheme_t   scheme,
                                  cs_real_t                     coefval[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_ADVECTION_H__ */
