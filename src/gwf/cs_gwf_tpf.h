#ifndef __CS_GWF_TPF_H__
#define __CS_GWF_TPF_H__

/*============================================================================
 * Main functions dedicated to the modelling of two-phase flows in a porous
 * media. This media is always considered as unsaturated. Two sub-models are
 * considered: miscible (MTPF) or immiscible (ITPF)
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

#include "cs_gwf_hydraulic_model.h"

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
 * \brief Allocate and initialize the model context structure for two-phase
 *        flows in a porous media
 *
 * \param[in] model       type of physical modelling
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tpf_t *
cs_gwf_tpf_create(cs_gwf_model_type_t      model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the modelling of two-phase
 *        flows in a porous media
 *
 * \param[in, out] p_tpf   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_free(cs_gwf_tpf_t    **p_tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the relaxation property by value and set this value.
 *
 * \param[in] tpf  pointer to the model context structure
 * \param[in] val  reference value used to set the relaxation property
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_define_relax_pty_by_value(cs_gwf_tpf_t *tpf,
                                     double        val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the relaxation property by value and set this value.
 *
 * \param[in] tpf          pointer to the model context structure
 * \param[in] func         function pointer to a time function
 * \param[in] func_context context related to this function
 *
 * \return a pointer to the created definition (\ref cs_xdef_t structure)
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_gwf_tpf_define_relax_pty_by_time_func(cs_gwf_tpf_t   *tpf,
                                         cs_time_func_t *func,
                                         void           *func_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of two-phase flows.
 *        Common to the different sub-models relying on two-phase flows.
 *
 * \param[in] tpf      pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_log_setup(cs_gwf_tpf_t          *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a two-phase flows model in porous media
 *
 * \param[in, out] tpf         pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init(cs_gwf_tpf_t            *tpf,
                cs_property_type_t       perm_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for two-phase flows in porous media. At this
 *        stage, all soils have been defined and equation parameters are set.
 *        Case of a miscible or immiscible model.
 *
 * \param[in]      post_flag   optional postprocessing request(s)
 * \param[in, out] tpf         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_setup(cs_flag_t         post_flag,
                      cs_gwf_tpf_t     *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      flag      optional settings for the module
 * \param[in, out] tpf       pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_finalize_setup(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq,
                          cs_flag_t                    flag,
                          cs_gwf_tpf_t                *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_values(const cs_cdo_connect_t        *connect,
                       const cs_cdo_quantities_t     *cdoq,
                       cs_gwf_tpf_t                  *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new state for the groundwater flows module.
 *        Case of two-phase flows in porous media.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] tpf          pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_compute(const cs_mesh_t               *mesh,
                   const cs_cdo_connect_t        *connect,
                   const cs_cdo_quantities_t     *cdoq,
                   const cs_time_step_t          *time_step,
                   cs_flag_t                      option_flag,
                   cs_gwf_tpf_t                  *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Operate a "current to previous" step on fields or arrays which have
 *        at least a storage of the previous step (time t^n when computing
 *        t^{n+1})
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in, out] tpf      pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_current_to_previous(const cs_cdo_connect_t *connect,
                               cs_gwf_tpf_t           *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of a two-phase flow model in
 *        porous media (miscible or immiscible). To operate a "current to
 *        previous" step, one has to call the dedicated function \ref
 *        cs_gwf_tpf_current_to_previous()
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval    time at which properties are evaluated if needed
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] tpf          pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_update(const cs_mesh_t           *mesh,
                  const cs_cdo_connect_t    *connect,
                  const cs_cdo_quantities_t *cdoq,
                  double                     time_eval,
                  cs_flag_t                  option_flag,
                  cs_gwf_tpf_t              *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of miscible or immiscible two-phase flows in porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts         pointer to a cs_time_step_t struct.
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] tpf        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_op(const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    const cs_time_step_t         *ts,
                    cs_flag_t                     post_flag,
                    cs_gwf_tpf_t                 *tpf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of saturated two-phase flows (tpf) in porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] abs_perm     property for the absolute permeability
 * \param[in] tpf          pointer to the model context structure
 * \param[in] connect      pointer to additional connectivities for CDO
 * \param[in] cdoq         pointer to additional mesh quantities for CDO
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_post(int                         mesh_id,
                      cs_lnum_t                   n_cells,
                      const cs_lnum_t             cell_ids[],
                      cs_flag_t                   post_flag,
                      const cs_property_t        *abs_perm,
                      const cs_gwf_tpf_t         *tpf,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *cdoq,
                      const cs_time_step_t       *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TPF_H__ */
