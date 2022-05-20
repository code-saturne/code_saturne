#ifndef __CS_GWF_PARAM_H__
#define __CS_GWF_PARAM_H__

/*============================================================================
 * Types related to the groundwater flow module
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @name Flags specifying the kind of post-processing to perform in
 *       the groundwater flow module
 * @{
 *
 * \def CS_GWF_POST_SOIL_CAPACITY
 * \brief Activate the post-processing of the soil capacity (property in front
 *        of the unsteady term in Richards equation)
 *
 * \def CS_GWF_POST_LIQUID_SATURATION
 * \brief Activate the post-processing of the liquid saturation (also nammed
 *        "moisture content" in case of single phase flow)
 *
 * \def CS_GWF_POST_PERMEABILITY
 * \brief Activate the post-processing of the permeability field
 *
 * \def CS_GWF_POST_DARCY_FLUX_BALANCE
 * \brief Compute the overall balance at the different boundaries of
 *        the Darcy flux
 *
 * \def CS_GWF_POST_DARCY_FLUX_DIVERGENCE
 * \brief Compute in each control volume (vertices or cells w.r.t the space
 *        scheme) the divergence of the Darcy flux
 *
 * \def CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY
 * \brief Define a field at boundary faces for the Darcy flux and activate the
 *        post-processing
 *
 * \def CS_GWF_POST_GAS_MASS_DENSITY

 * \brief Compute the mass density of the gas component for a miscible or
 *        immiscible two-phase flow model. One recalls that one assumes that
 *        there is no water in the gas phase and that the mass density is a
 *        function of the gas pressure through a perfect gas law.
 */

#define CS_GWF_POST_SOIL_CAPACITY              (1 << 0)
#define CS_GWF_POST_LIQUID_SATURATION          (1 << 1)
#define CS_GWF_POST_PERMEABILITY               (1 << 2)
#define CS_GWF_POST_DARCY_FLUX_BALANCE         (1 << 3)
#define CS_GWF_POST_DARCY_FLUX_DIVERGENCE      (1 << 4)
#define CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY     (1 << 5)
#define CS_GWF_POST_GAS_MASS_DENSITY           (1 << 6)

/*!
 * @}
 */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * @name Model definition for the groundwater flow module
 * @{
 *
 * \enum cs_gwf_model_type_t
 * \brief Type of system of equation(s) to consider for the physical modelling
 */

typedef enum {

  /*!
   * \brief Single phase (liquid phase) modelling in a porous media.
   *
   * All soils are assumed to be saturated. This yields several simplifications
   * in the Richards equation governing the water conservation. The Richards
   * equation is steady. The saturation is constant and there is no relative
   * permeability.
   */

  CS_GWF_MODEL_SATURATED_SINGLE_PHASE,

  /*!
   * \brief Single phase (liquid phase) modelling in a porous media.
   *
   * Some soils are not saturated and are described by a more complex model
   * such as the Van Genuchten-Mualen model. Simplifications made in the case
   * of \ref CS_GWF_MODEL_SATURATED_SINGLE_PHASE do not hold anymore. Richards
   * equation is unsteady and there may be a non-linearity to handle according
   * to the type of soil model. Soil properties such as permeability, soil
   * capacity and liquid saturation (also called moisture content) are neither
   * uniform nor steady.
   */

  CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE,

  /*!
   * \brief Miscible two phase flow modelling (gaseous and liquid phases) in
   *        porous media.
   *
   * A Richards-like equation is considered in each phase to take into account
   * the mass conservation of water and one other component. The component can
   * be disolved in the liquid phase. No water vapour is taken into
   * account. Please refer to \ref cs_gwf_two_phase_t for more details.
   */

  CS_GWF_MODEL_MISCIBLE_TWO_PHASE,

  /*!
   * \brief Immiscible two phase flow modelling (gaseous and liquid phases) in
   *        porous media.
   *
   * A Richards-like equation is considered in each phase to take into account
   * the mass conservation of water in the liquid phase and the conservation of
   * the other component in the gaseous phase. The model context is shared with
   * the miscible two-phase flow model. Please refer to \ref cs_gwf_two_phase_t
   * for more details.
   */

  CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE,

  CS_GWF_N_MODEL_TYPES     /*!< Number of predefined models (not a model) */

} cs_gwf_model_type_t;


/*!
 * @}
 * \enum cs_gwf_model_bit_t
 * \brief Additional modelling options either from the physical viewpoint or the
 *        numerical viewpoint
 */

typedef enum {

  /*!
   * @name Main physical modelling tags
   * @{
   *
   * \var CS_GWF_GRAVITATION
   * \brief Gravitation effects are taken into account in the Richards equation
   */

  CS_GWF_GRAVITATION                     = 1<< 0,  /* =   1 */

  /*!
   * @}
   * @name Main numerical flags
   * @{
   *
   * \var CS_GWF_FORCE_RICHARDS_ITERATIONS
   * \brief Even if the Richards equation is steady-state, this equation is
   *        solved at each iteration.
   */

  CS_GWF_FORCE_RICHARDS_ITERATIONS       = 1<< 6,  /* =   64 */

  /*!
   * \var CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE
   * \brief Compute the mean-value of the hydraulic head field and subtract
   *        this mean-value to get a field with zero mean-value. It's important
   *        to set this flag if no boundary condition is given.
   */

  CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE = 1<< 7,  /* =  128 */

  /*!
   * \var CS_GWF_ENFORCE_DIVERGENCE_FREE
   * \brief Activate a treatment to enforce a Darcy flux to be divergence-free
   */

  CS_GWF_ENFORCE_DIVERGENCE_FREE         = 1<< 8   /* =  256 */

} cs_gwf_model_bit_t;

/*!
 * @}
 * @name Soil modelling
 * @{
 *
 * \enum cs_gwf_soil_model_t
 * \brief Predefined hydraulic model of soils used in the groundwater flow
 *        module
 *
 * \var CS_GWF_SOIL_GENUCHTEN
 * Van Genuchten-Mualem laws defining the evolution of the effective liquid
 * saturne (also called dimensionless moisture content) and the relative
 * permeability
 *
 * The (effective) liquid saturation (also called moisture content) follows the
 * identity:
 * S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *         = (1 + |alpha . h|^n)^(-m)
 *
 * The isotropic relative permeability is defined as:
 * k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 * where m = 1 -  1/n
 *
 * \var CS_GWF_SOIL_SATURATED
 * Hydraulic model of soild where the soil is considered as saturated. In this
 * model, there no evolution taken into account. The liquid saturation and the
 * permeability are considered as constant.
 *
 * \var CS_GWF_SOIL_USER
 * User-defined model of soil
 */

typedef enum {

  CS_GWF_SOIL_GENUCHTEN,
  CS_GWF_SOIL_SATURATED,
  CS_GWF_SOIL_USER,

  CS_GWF_SOIL_N_HYDRAULIC_MODELS

} cs_gwf_soil_model_t;

/*!
 * @}
 * @name Tracer modelling
 * @{
 */

typedef cs_flag_t  cs_gwf_tracer_model_t;

/*!
 * \enum cs_gwf_tracer_model_bit_t
 * \brief Flags specifying the general behavior of a tracer associated to
 *        the groundwater flow module
 *
 * Elemental modelling choice either from the physical viewpoint or the
 * numerical viewpoint for the transport of a tracer
 *
 */

typedef enum {

  /*!
   * \brief User-defined tracer.
   *
   * All terms can be modified with user functions
   */

  CS_GWF_TRACER_USER                        = 1<< 0, /* =    1 */

  /* Physical phenomena to consider */
  /* ------------------------------ */

  /*!
   * \brief EK model with 3 parameters
   *
   * Add the sorption phenomena to the default tracer equation. Case of the EK
   * model with 3 parameters. Sorption is assumed to be infinite
   */

  CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS    = 1<< 1, /* =    2 */

  /*!
   * \brief EK model with 5 parameters.
   *
   * Add the sorption phenomena to the default tracer equation in the case of
   * the EK model with 5 parameters. Sorption is assumed to be finite.  An
   * additional equation related to the concentration of sorpted tracer in the
   * second kind of sites.
   */

  CS_GWF_TRACER_SORPTION_EK_5_PARAMETERS    = 1<< 2, /* =    4 */

  /*!
   * \brief Add the precipitation phenomena to the default tracer equation
   */

  CS_GWF_TRACER_PRECIPITATION               = 1<< 4, /* =    16 */

} cs_gwf_tracer_model_bit_t;

/*!
 * @}
 */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PARAM_H__ */
