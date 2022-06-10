#ifndef __CS_GWF_H__
#define __CS_GWF_H__

/*============================================================================
 * Set of main functions to handle the groundwater flow module with CDO
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

#include "cs_base.h"
#include "cs_equation.h"
#include "cs_gwf_param.h"
#include "cs_gwf_priv.h"
#include "cs_gwf_soil.h"
#include "cs_gwf_tracer.h"

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
 * \brief  Check if the groundwater flow module has been activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]   model           type of physical modelling
 * \param[in]   option_flag     optional flag to specify this module
 * \param[in]   post_flag       optional automatic postprocessing
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_activate(cs_gwf_model_type_t      model,
                cs_flag_t                option_flag,
                cs_flag_t                post_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to groundwater flows
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the main structure which manages a two-phase flow model
 *
 * \return a pointer to the structure cs_gwf_two_phase_t
 */
/*----------------------------------------------------------------------------*/

cs_gwf_two_phase_t *
cs_gwf_get_two_phase_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the numerical options related to the two phase flow models
 *
 * \param[in] use_coupled_solver         true/false
 * \param[in] use_incremental_solver     true/false
 * \param[in] use_properties_on_submesh  true/false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_two_phase_numerical_options(bool    use_coupled_solver,
                                       bool    use_incremental_solver,
                                       bool    use_properties_on_submesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters defining the two-phase flow model.
 *         Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] l_diffusivity_h  diffusivity of the main gas component in the
 *                             liquid phase
 * \param[in] w_molar_mass     molar mass of the main liquid component
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 * \param[in] henry_constant   constant in the Henry law
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_miscible_two_phase_model(cs_real_t       l_mass_density,
                                    cs_real_t       l_viscosity,
                                    cs_real_t       g_viscosity,
                                    cs_real_t       l_diffusivity_h,
                                    cs_real_t       w_molar_mass,
                                    cs_real_t       h_molar_mass,
                                    cs_real_t       ref_temperature,
                                    cs_real_t       henry_constant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters defining the immiscible two-phase flow model.
 *         Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_immiscible_two_phase_model(cs_real_t       l_mass_density,
                                      cs_real_t       l_viscosity,
                                      cs_real_t       g_viscosity,
                                      cs_real_t       h_molar_mass,
                                      cs_real_t       ref_temperature);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the flag dedicated to the post-processing of the GWF module
 *
 * \param[in]  post_flag             flag to set
 * \param[in]  reset                 reset post flag before
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_post_options(cs_flag_t       post_flag,
                        bool            reset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the advection field related to the Darcy flux in the liquid
 *         phase
 *
 * \return a pointer to a cs_adv_field_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_gwf_get_adv_field(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a new cs_gwf_soil_t structure. An initialization by
 *         default of all members is performed.
 *         Case of a soil with an isotropic absolute permeability
 *
 * \param[in]  z_name      name of the volume zone corresponding to the soil
 * \param[in]  density     value of the bulk mass density
 * \param[in]  k_abs       absolute (or intrisic) permeability (scalar-valued)
 * \param[in]  porosity    value of the porosity (saturated moisture content)
 * \param[in]  model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_iso_soil(const char                *z_name,
                    double                     density,
                    double                     k_abs,
                    double                     porosity,
                    cs_gwf_soil_model_t        model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a new cs_gwf_soil_t structure. An initialization by
 *         default of all members is performed.
 *
 * \param[in]  z_name      name of the volume zone corresponding to the soil
 * \param[in]  density     value of the bulk mass density
 * \param[in]  k_abs       absolute (or intrisic) permeability (tensor-valued)
 * \param[in]  porosity    value of the porosity (saturated moisture content)
 * \param[in]  model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_aniso_soil(const char                *z_name,
                      double                     density,
                      double                     k_abs[3][3],
                      double                     porosity,
                      cs_gwf_soil_model_t        model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module

 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation. Tracer is advected thanks to the darcian velocity
 *         and diffusion/reaction parameters result from a physical modelling.
 *         Terms solved in this equation are activated according to predefined
 *         settings. The advection field corresponds to that of the liquid
 *         phase.
 *
 * \param[in]  tr_model   physical modelling to consider (0 = default settings)
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_tracer(cs_gwf_tracer_model_t     tr_model,
                  const char               *eq_name,
                  const char               *var_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation.  Tracer is advected thanks to the darcian
 *         velocity and diffusion/reaction parameters result from a physical
 *         modelling. Terms are activated according to predefined settings.
 *         Modelling of the tracer parameters are left to the user
 *
 * \param[in] eq_name         name of the tracer equation
 * \param[in] var_name        name of the related variable
 * \param[in] init_setup      function pointer (predefined prototype)
 * \param[in] finalize_setup  function pointer (predefined prototype)
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_user_tracer(const char                       *eq_name,
                       const char                       *var_name,
                       cs_gwf_tracer_init_setup_t       *init_setup,
                       cs_gwf_tracer_finalize_setup_t   *finalize_setup);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context of the model after the activation of the
 *         module and a first settings of the model parameters (physical and
 *         numerical). At this stage, cs_user_parameters() has not been called
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_model_context(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the groundwater flow model and its related
 *         equations.
 *         At this stage, all soils have been defined and equation parameters
 *         are set (cs_user_parameters() has been called).
 *         Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last initialization step of the groundwater flow module. At this
 *        stage, the mesh quantities are defined.
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_finalize_setup(const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the GWF module (done after all the setup phase and after
 *         the initialization of all equations)
 *         One sets an initial value to all quantities related to this module.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_values(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the groundwater system (pressure head, head in law, moisture
 *         content, darcian velocity, soil capacity or permeability if needed)
 *
 * \param[in]  mesh         pointer to a cs_mesh_t structure
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts           pointer to a cs_time_step_t structure
 * \param[in]  update_flag  metadata associated to the status of the update
 *                          step to perform
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update(const cs_mesh_t             *mesh,
              const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              const cs_time_step_t        *ts,
              cs_flag_t                    update_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute_steady_state(const cs_mesh_t              *mesh,
                            const cs_time_step_t         *time_step,
                            const cs_cdo_connect_t       *connect,
                            const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t              *mesh,
               const cs_time_step_t         *time_step,
               const cs_cdo_connect_t       *connect,
               const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_op(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of saturated single-phase flows (sspf) in porous media.
 *         Prototype of this function is given since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post_sspf(void                   *input,
                       int                     mesh_id,
                       int                     cat_id,
                       int                     ent_flag[5],
                       cs_lnum_t               n_cells,
                       cs_lnum_t               n_i_faces,
                       cs_lnum_t               n_b_faces,
                       const cs_lnum_t         cell_ids[],
                       const cs_lnum_t         i_face_ids[],
                       const cs_lnum_t         b_face_ids[],
                       const cs_time_step_t   *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of unsaturated single-phase flows (uspf) in porous media.
 *         Prototype of this function is given since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post_uspf(void                   *input,
                       int                     mesh_id,
                       int                     cat_id,
                       int                     ent_flag[5],
                       cs_lnum_t               n_cells,
                       cs_lnum_t               n_i_faces,
                       cs_lnum_t               n_b_faces,
                       const cs_lnum_t         cell_ids[],
                       const cs_lnum_t         i_face_ids[],
                       const cs_lnum_t         b_face_ids[],
                       const cs_time_step_t   *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of miscible two-phase flows (mtpf) in porous media.
 *         Prototype of this function is given since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post_mtpf(void                      *input,
                       int                        mesh_id,
                       int                        cat_id,
                       int                        ent_flag[5],
                       cs_lnum_t                  n_cells,
                       cs_lnum_t                  n_i_faces,
                       cs_lnum_t                  n_b_faces,
                       const cs_lnum_t            cell_ids[],
                       const cs_lnum_t            i_face_ids[],
                       const cs_lnum_t            b_face_ids[],
                       const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_H__ */
