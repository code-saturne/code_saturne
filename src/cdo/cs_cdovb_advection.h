#ifndef __CS_CDOVB_ADVECTION_H__
#define __CS_CDOVB_ADVECTION_H__

/*============================================================================
 * Build discrete convection operators for CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "cs_cdo.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_cdovb_adv_t  cs_cdovb_adv_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure for the convection operator
 *
 * \param[in]  connect       pointer to the connectivity structure
 * \param[in]  adv_field     pointer to a cs_adv_field_t structure
 * \param[in]  a_info        set of options for the advection term
 * \param[in]  do_diffusion  true is diffusion is activated
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_adv_t *
cs_cdovb_advection_builder_init(const cs_cdo_connect_t      *connect,
                                const cs_adv_field_t        *adv,
                                const cs_param_advection_t   a_info,
                                bool                         do_diffusion);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a builder structure for the convection operator
 *
 * \param[in, out] b   pointer to a cs_cdovb_adv_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_adv_t *
cs_cdovb_advection_builder_free(cs_cdovb_adv_t  *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      diffmat   tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_advection_build_local(const cs_cell_mesh_t      *cm,
                               const cs_real_33_t         diffmat,
                               cs_cdovb_adv_t            *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      quant         pointer to the cdo quantities structure
 * \param[in]      cm            pointer to a cs_cell_mesh_t struct.
 * \param[in, out] advb          pointer to a convection builder structure
 * \param[in, out] ls            cell-wise structure sotring the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_add_bc(const cs_cdo_quantities_t   *quant,
                          cs_cell_mesh_t              *cm,
                          cs_cdovb_adv_t              *advb,
                          cs_cdo_locsys_t             *ls);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell in a given direction
 *
 * \param[in]      cdoq           pointer to the cdo quantities structure
 * \param[in]      adv            pointer to the advection field struct.
 * \param[in]      diff_property  pointer to the diffusion property struct.
 * \param[in]      dir_vect       direction for estimating the Peclet number
 * \param[in, out] peclet         pointer to the pointer of real numbers to fill
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_get_peclet_cell(const cs_cdo_quantities_t   *cdoq,
                                   const cs_adv_field_t        *adv,
                                   const cs_property_t         *diff_property,
                                   const cs_real_3_t            dir_vect,
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
cs_cdovb_advection_get_upwind_coef_cell(const cs_cdo_quantities_t   *cdoq,
                                        const cs_param_advection_t   a_info,
                                        cs_real_t                    coefval[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_ADVECTION_H__ */
