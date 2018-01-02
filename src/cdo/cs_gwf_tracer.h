#ifndef __CS_GWF_TRACER_H__
#define __CS_GWF_TRACER_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function pointer prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the phisical properties related to a
 *         tracer modelling
 *
 * \param[in, out] input        pointer to a structure cast on-the-fly
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_update_t) (void                        *input,
                          const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to free the input of a tracer model
 *
 * \param[in, out] input     pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_free_input_t) (void      *input);

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Type of predefined modelling for the groundwater flows */
typedef enum {

  CS_GWF_TRACER_STANDARD, /* Default behavior of a tracer in groundwater flow
                            module */
  CS_GWF_TRACER_USER,     /* User-defined behavior */
  CS_GWF_N_TRACER_MODELS

} cs_gwf_tracer_model_t;

/* Set of parameters related to a tracer equation attached to a standard
   modelling */

typedef struct {

  /* Parameters (array of size: n_soils) */
  double    *rho_kd;        // Bulk density times the distribution coefficient
  double    *alpha_l;       // Longitudinal dispersivity
  double    *alpha_t;       // Transversal dispersivity
  double    *wmd;           // Water molecular diffusivity
  double    *reaction_rate; /* First order decay coefficient (related to the
                               reaction term) */

  /* Variables used for the update of physical properties */
  cs_field_t   *darcy_velocity_field;
  cs_field_t   *moisture_content;

} cs_gwf_std_tracer_input_t;

/* Set of parameters describing a tracer */
typedef struct {

  int                          id;       /* tracer id */
  cs_equation_t               *eq;       /* related equation */

  /* Physical modelling adopted for this tracer */
  cs_gwf_tracer_model_t        model;

  cs_field_t                  *diffusivity; /* NULL if no diffusion term is
                                               build in the tracer equation */
  int                          reaction_id; /* id related to the reaction
                                               term in the tracer equation */

  /* Pointer to an input structure according to the model */
  void                        *input;

  /* Pointers to functions */
  cs_gwf_tracer_update_t      *update_properties;
  cs_gwf_tracer_free_input_t  *free_input;

} cs_gwf_tracer_t;

/*============================================================================
 * Public function pointer prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set the parameters related to a tracer equation
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_setup_t) (const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the terms to build in the algebraic
 *         system for a tracer equation according to the settings
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_add_terms_t) (cs_gwf_tracer_t             *tracer);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new cs_gwf_tracer_t structure and initialize its members by
 *         default.
 *         Add a new equation related to the groundwater flow module.
 *         This equation is a specific transport equation.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in]   tracer_id   id number of the soil
 * \param[in]   eqname      name of the tracer equation
 * \param[in]   varname     name of the related variable
 * \param[in]   model       model related to this tracer
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_init(int                      tracer_id,
                   const char              *eq_name,
                   const char              *var_name,
                   cs_adv_field_t          *adv_field,
                   cs_gwf_tracer_model_t    model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_gwf_tracer_t structure
 *
 * \param[in, out]  tracer   pointer to a cs_gwf_tracer_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_free(cs_gwf_tracer_t     *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a tracer for a specified soil when the tracer is attached to
 *         the default model
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or NULL if all
 *                                 soils are selected)
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_standard_tracer(cs_gwf_tracer_t   *tracer,
                           const char        *soil_name,
                           double             wmd,
                           double             alpha_l,
                           double             alpha_t,
                           double             distrib_coef,
                           double             reaction_rate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add terms to the algebraic system related to a tracer equation
 *         according to the settings.
 *         Case of the standar tracer modelling
 *         Rely on the generic functinon: cs_gwf_tracer_add_terms_t
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_standard_add_terms(cs_gwf_tracer_t     *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_standard_setup(const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TRACER_H__ */
