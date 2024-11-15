#ifndef __CS_MACFB_ADVECTION_H__
#define __CS_MACFB_ADVECTION_H__

/*============================================================================
 * Build discrete convection operators for MAC schemes
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include "cs_macfb_builder.h"
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
 * Function pointers for MAC face-based schemes
 * -------------------------------------------------------------------------- */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform preprocessing such as the computation of the advection flux
 *         at the expected location in order to be able to build the advection
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      csys     pointer to a cs_cell_sys_t structure
 * \param[in, out] input    null or pointer to a structure cast on-the-fly
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void(cs_macfb_adv_open_hook_t)(const cs_equation_param_t *eqp,
                                       const cs_cell_mesh_t      *cm,
                                       const cs_macfb_builder_t  *macb,
                                       const cs_cell_sys_t       *csys,
                                       void                      *input,
                                       cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator in MAC-Fb schemes. Case of an
 *          operator that should be assemble into a matrix. Boundary conditions
 *          are also handled in this function.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

typedef void(cs_macfb_adv_scheme_t)(const cs_cell_mesh_t     *cm,
                                    const cs_macfb_builder_t *macb,
                                    cs_cell_builder_t        *cb,
                                    cs_sdm_t                 *adv,
                                    cs_real_t                *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the cellwise advection operator for MAC-Fb schemes
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          This pointer of function manages how is build the advection term
 *          in a cell and it relies on the lower-level function
 *          \ref cs_macfb_adv_scheme_t
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      scheme_func  function pointer to the scheme definition
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

typedef void(cs_macfb_adv_build_t)(const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const cs_macfb_builder_t  *macb,
                                   cs_macfb_adv_scheme_t     *scheme_func,
                                   cs_cell_sys_t             *csys,
                                   cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Operation done after the matrix related to the advection term has
 *          been defined.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void(cs_macfb_adv_close_hook_t)(const cs_cell_mesh_t     *cm,
                                        const cs_macfb_builder_t *macb,
                                        const cs_cell_builder_t  *cb,
                                        cs_cell_sys_t            *csys);

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform preprocessing such as the computation of the advection flux
 *         at the expected location in order to be able to build the advection
 *         matrix. Follow the prototype given by cs_macfb_adv_open_hook_t
 *         Default case.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      macb     pointer to a cs_macfb_builder_t structure
 * \param[in]      csys     pointer to a cs_cell_sys_t structure
 * \param[in, out] input    null or pointer to a structure cast on-the-fly
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_open_default(const cs_equation_param_t *eqp,
                                     const cs_cell_mesh_t      *cm,
                                     const cs_macfb_builder_t  *macb,
                                     const cs_cell_sys_t       *csys,
                                     void                      *input,
                                     cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_macfb_adv_close_hook_t
 *         Default vector-valued case.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_close_default_vect(const cs_cell_mesh_t     *cm,
                                           const cs_macfb_builder_t *macb,
                                           const cs_cell_builder_t  *cb,
                                           cs_cell_sys_t            *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_macfb_adv_close_hook_t
 *         Explicit treatment without extrapolation for vector-valued DoFs.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_close_exp_none_vect(const cs_cell_mesh_t     *cm,
                                            const cs_macfb_builder_t *macb,
                                            const cs_cell_builder_t  *cb,
                                            cs_cell_sys_t            *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main function to build the cellwise advection operator for MAC-Fb
 *         schemes The local matrix related to this operator is stored in
 *         cb->loc
 *
 *         Case of an advection term without a diffusion operator. In this
 *         situation, a numerical issue may arise if an internal or a border
 *         face is such that there is no advective flux. A special treatment
 *         is performed to tackle this issue.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_no_diffusion(const cs_equation_param_t *eqp,
                                     const cs_cell_mesh_t      *cm,
                                     const cs_macfb_builder_t  *macb,
                                     cs_macfb_adv_scheme_t     *scheme_func,
                                     cs_cell_sys_t             *csys,
                                     cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main function to build the cellwise advection operator for MAC
 *          face-based schemes.
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          One assumes that a diffusion term is present so that there is no
 *          need to perform additional checkings on the well-posedness of the
 *          operator.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection(const cs_equation_param_t *eqp,
                        const cs_cell_mesh_t      *cm,
                        const cs_macfb_builder_t  *macb,
                        cs_macfb_adv_scheme_t     *scheme_func,
                        cs_cell_sys_t             *csys,
                        cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - non-conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - upwind scheme

 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_upwnoc(const cs_cell_mesh_t     *cm,
                               const cs_macfb_builder_t *macb,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *adv,
                               cs_real_t                *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - conservative formulation \f$ \nabla\cdot(\beta \otimes u) \f$
 *         - upwind scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_upwcsv(const cs_cell_mesh_t     *cm,
                               const cs_macfb_builder_t *macb,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *adv,
                               cs_real_t                *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - non-conservative formulation beta.grad
 *         - centered scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_cennoc(const cs_cell_mesh_t     *cm,
                               const cs_macfb_builder_t *macb,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *adv,
                               cs_real_t                *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - centered scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_advection_cencsv(const cs_cell_mesh_t     *cm,
                               const cs_macfb_builder_t *macb,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *adv,
                               cs_real_t                *rhs);

END_C_DECLS

#endif /* __CS_MACFB_ADVECTION_H__ */
