#ifndef __CS_CDO_CONVECTION_H__
#define __CS_CDO_CONVECTION_H__

/*============================================================================
 * Build discrete convection operators for CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_mesh.h"

#include "cs_cdo.h"
#include "cs_sla.h"
#include "cs_param.h"
#include "cs_cdo_toolbox.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _convection_builder_t  cs_convection_builder_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection flux accross the dual face df(e) lying
 *          inside the cell c and associated to the edge e.
 *          This function is associated to vertex-based discretization.
 *
 * \param[in]  quant    pointer to the cdo quantities structure
 * \param[in]  a_info   set of options for the advection term
 * \param[in]  tcur     value of the current time
 * \param[in]  xc       center of the cell c
 * \param[in]  qe       quantities related to edge e in E_c
 * \param[in]  qdf      quantities to the dual face df(e)
 *
 * \return the value of the convective flux accross the triangle
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_convection_vbflux_compute(const cs_cdo_quantities_t   *quant,
                             const cs_param_advection_t   a_info,
                             double                       tcur,
                             const cs_real_3_t            xc,
                             const cs_quant_t             qe,
                             const cs_dface_t             qdf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure for the convection operator
 *
 * \param[in]  connect       pointer to the connectivity structure
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  a_info        set of options for the advection term
 * \param[in]  do_diffusion  true is diffusion is activated
 * \param[in]  d_info        set of options for the diffusion term
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

cs_convection_builder_t *
cs_convection_builder_init(const cs_cdo_connect_t      *connect,
                           const cs_time_step_t        *time_step,
                           const cs_param_advection_t   a_info,
                           bool                         do_diffusion,
                           const cs_param_hodge_t       d_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a builder structure for the convection operator
 *
 * \param[in, out] b   pointer to a cs_convection_builder_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_convection_builder_t *
cs_convection_builder_free(cs_convection_builder_t  *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator for pure convection
 *
 * \param[in]      c_id       cell id
 * \param[in]      connect    pointer to the connectivity structure
 * \param[in]      quant      pointer to the cdo quantities structure
 * \param[in]      loc_ids    store the local entity ids for this cell
 * \param[in, out] builder    pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_convection_build_local_op(cs_lnum_t                    c_id,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_lnum_t             *loc_ids,
                             cs_convection_builder_t     *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator for pure convection
 *
 * \param[in]      mesh         pointer to a cs_msesh_t structure
 * \param[in]      connect      pointer to the connectivity structure
 * \param[in]      quant        pointer to the cdo quantities structure
 * \param[in]      dir_vals     values of the Dirichlet boundary condition
 * \param[in, out] builder      pointer to a convection builder structure
 * \param[in, out] rhs_contrib  array storing the contribution for the rhs
 * \param[in, out] diag_contrib array storing the contribution for the diagonal
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_get_bc_contrib(const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_real_t             *dir_vals,
                             cs_convection_builder_t     *builder,
                             cs_real_t                    rhs_contrib[],
                             cs_real_t                    diag_contrib[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell in a given direction
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in]      a_info   set of options for the advection term
 * \param[in]      d_info   set of options for the diffusion term
 * \param[in]      dir_vect  direction in which we estimate the Peclet number
 * \param[in]      tcur      value of the current time
 * \param[in, out] peclet    pointer to the pointer of real numbers to fill
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_get_peclet_cell(const cs_cdo_quantities_t   *cdoq,
                              const cs_param_advection_t   a_info,
                              const cs_param_hodge_t       d_info,
                              const cs_real_3_t            dir_vect,
                              cs_real_t                    tcur,
                              cs_real_t                   *p_peclet[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value in each cell of the upwinding coefficient given
 *          a related Peclet number
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_get_upwind_coef_cell(const cs_cdo_quantities_t   *cdoq,
                                   const cs_param_advection_t   a_info,
                                   cs_real_t                    coefval[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_CONVECTION_H__ */
