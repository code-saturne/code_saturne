#ifndef __CS_GWF_SSPF_H__
#define __CS_GWF_SSPF_H__

/*============================================================================
 * Main functions to handle single-phase flows in a saturated porous media
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "gwf/cs_gwf_hydraulic_model.h"

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
 * \brief Allocate and initialize the model context structure in case of
 *        single-phase flows in a saturated porous media
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_sspf_t *
cs_gwf_sspf_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the model context structure in case of a saturated single-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_free(cs_gwf_sspf_t   **p_mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of saturated single-phase
 *        flows in porous media
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_log_setup(cs_gwf_sspf_t   *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a saturated single-phase flows model in porous media
 *
 * \param[in, out] mc          pointer to the model context structure
 * \param[in]      abs_perm    property struct. for the absolute permeability
 * \param[in]      flag        optional metadata
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_init(cs_gwf_sspf_t       *mc,
                 cs_property_t       *abs_perm,
                 cs_flag_t            flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for saturated single-phase flows in a porous
 *        media. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in]      flag         optional settings for the module
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_init_setup(cs_flag_t           flag,
                       cs_gwf_sspf_t      *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of single-phase flows in a saturated
 *        porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_finalize_setup(const cs_cdo_connect_t        *connect,
                           const cs_cdo_quantities_t     *cdoq,
                           cs_gwf_sspf_t                 *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of single-phase flows in a
 *        saturated porous media
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  type of operation(s) to do
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_update(const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq,
                   const cs_time_step_t         *ts,
                   cs_flag_t                     update_flag,
                   cs_flag_t                     option_flag,
                   cs_gwf_sspf_t                *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the steady-state of the groundwater flows module in case of
 *        single-phase flows in a saturated porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      flag       type of additional treatment(s) to do
 * \param[in, out] mc         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_compute_steady_state(const cs_mesh_t             *mesh,
                                 const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *cdoq,
                                 const cs_time_step_t        *time_step,
                                 cs_flag_t                    flag,
                                 cs_gwf_sspf_t               *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new hydraulic state for the groundwater flows module.
 *        Case of single-phase flows in a saturated porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      flag       type of additional treatment(s) to do
 * \param[in, out] mc         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_compute(const cs_mesh_t              *mesh,
                    const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    const cs_time_step_t         *time_step,
                    cs_flag_t                     flag,
                    cs_gwf_sspf_t                *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of single phase flows in a saturated porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_extra_op(const cs_cdo_connect_t         *connect,
                     const cs_cdo_quantities_t      *cdoq,
                     cs_flag_t                       post_flag,
                     cs_gwf_sspf_t                  *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of single-phase flows in a saturated porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] abs_perm     property for the absolute permeability
 * \param[in] mc           pointer to the model context structure
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_extra_post(int                        mesh_id,
                       cs_lnum_t                  n_cells,
                       const cs_lnum_t            cell_ids[],
                       cs_flag_t                  post_flag,
                       const cs_property_t       *abs_perm,
                       const cs_gwf_sspf_t       *mc,
                       const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_SSPF_H__ */
