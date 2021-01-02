#ifndef __CS_GWF_TRACER_H__
#define __CS_GWF_TRACER_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Typedef definition
 *============================================================================*/

typedef cs_flag_t  cs_gwf_tracer_model_t;
typedef struct _gwf_tracer_t  cs_gwf_tracer_t;

/*============================================================================
 * Public function pointer prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the phisical properties related to a
 *         tracer modelling
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_update_t) (cs_gwf_tracer_t             *tracer,
                          cs_real_t                    t_eval,
                          const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to free the input of a tracer model
 *
 * \param[in, out] tracer     pointer to a structure cs_gwf_tracer_t
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_free_input_t) (cs_gwf_tracer_t      *tracer);

/*============================================================================
 * Structure definitions
 *============================================================================*/


/*!
 * @name Flags specifying the general behavior of a tracer associated to
 *       the groundwater flow module
 * @{
 *
 * \enum cs_gwf_tracer_model_bit_t
 * \brief elemental modelling choice either from the physical viewpoint or the
 * numerical viewpoint for the transport of a tracer
 *
 * \def CS_GWF_TRACER_USER
 * \brief user-defined tracer. All terms can be modified with user functions
 *
 * \def CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS
 * \brief Add the sorption phenomena to the default tracer equation. Case of
 *        the EK model with 3 parameters. Sorption is assumed to be infinite
 *
 * \def CS_GWF_TRACER_SORPTION_EK_5_PARAMETERS
 * \brief Add the sorption phenomena to the default tracer equation. Case of the
 *        EK model with 5 parameters. Sorption is assumed to be finite.  An
 *        additional equation related to the concentration of sorpted tracer in
 *        the second kind of sites.
 *
 * \def CS_GWF_TRACER_PRECIPITATION
 * \brief Add the precipitation phenomena to the default tracer equation
 */

/* Type of predefined modelling for the groundwater flows */
typedef enum {

     CS_GWF_TRACER_USER                        = 1<< 0, /* =    1 */

     /* Physical phenomena to consider */
     /* ------------------------------ */

     CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS    = 1<< 1, /* =    2 */
     CS_GWF_TRACER_SORPTION_EK_5_PARAMETERS    = 1<< 2, /* =    4 */


     CS_GWF_TRACER_PRECIPITATION               = 1<< 4, /* =    16 */

} cs_gwf_tracer_model_bit_t;

/*! @} */

/* Set of parameters related to a tracer equation attached to a standard
   modelling */

typedef struct {

  /* Parameters to determine the behavior in each soil
   * (array of size: n_soils)
   */

  /* Common settings shared by all physical modelling */
  /* ------------------------------------------------ */

  double    *rho_bulk;      /* bulk density (kg.m^-3) */
  double    *kd0;           /* reference value of the distribution coefficient
                               (m^".kg^-1) */
  double    *rho_kd;        /* Derived quantity: rho_bulk*kd0 */

  double    *alpha_l;       /* Longitudinal dispersivity */
  double    *alpha_t;       /* Transversal dispersivity */

  double    *wmd;           /* Water molecular diffusivity (m^2.s^-1) */

  double    *reaction_rate; /* First order decay coefficient (related to the
                               reaction term) */

  /* Precipitation members (set to NULL if not used) */
  /* ----------------------------------------------- */

  double       *conc_w_star;    /* maximal value of the concentraction of tracer
                                   in the liquid phase. Exceeded quantities are
                                   stored in the solid phase (->
                                   conc_precip). These values corresponds to the
                                   user settings */

  cs_real_t    *conc_satura;    /* array storing the value of the saturated
                                   concentration in the liquid phase at vertices
                                   (only used in case of CDOVB schemes or CDOVCB
                                   schemes */
  cs_real_t    *conc_precip;    /* array storing the concentration in the
                                   precipitation (solid) storage. The size may
                                   vary w.r.t. to the discrtization scheme.
                                */

  cs_field_t   *precip_field;

  /* Sorption members (set to NULL if not used) */
  /* ------------------------------------------ */

  double       *k0_plus;        /* kinetic coefficient towards site 2 locations
                                   (m^3.kg^-1.s^-1) */
  double       *k0_minus;       /* kinetic coefficient from site 2 locations
                                   (s^-1) */

  cs_real_t    *conc_site2;     /* array allocated to n_cells */

  /* Variables used for the update of physical properties (shared pointers) */

  cs_field_t   *darcy_velocity_field;
  cs_field_t   *moisture_content;    /* also called the saturation, denoted by
                                        \theta (-) */

} cs_gwf_tracer_input_t;

/* Set of parameters describing a tracer structure */
/* ----------------------------------------------- */

struct _gwf_tracer_t{

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
  cs_gwf_tracer_update_t      *update_diff_tensor;
  cs_gwf_tracer_update_t      *update_precipitation;
  cs_gwf_tracer_free_input_t  *free_input;

};

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
 * \param[in]   eq_name     name of the tracer equation
 * \param[in]   var_name    name of the related variable
 * \param[in]   adv_field   pointer to a cs_adv_field_t structure
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
 * \brief  For a specified soil set the main parameters corresponding to a
 *         default modelling of a tracer transport
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
cs_gwf_set_main_tracer_param(cs_gwf_tracer_t   *tracer,
                             const char        *soil_name,
                             double             wmd,
                             double             alpha_l,
                             double             alpha_t,
                             double             distrib_coef,
                             double             reaction_rate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a specified soil set the parameters corresponding to a
 *         precipitation modelling of a tracer transport
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or NULL if all
 *                                 soils are selected)
 * \param[in]      conc_w_star     value of the saturated concentration in the
 *                                 liquid phase
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_precip_tracer_param(cs_gwf_tracer_t   *tracer,
                               const char        *soil_name,
                               double             conc_w_star);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add terms to the algebraic system related to a tracer equation
 *         according to the settings.
 *         Case of the standard tracer modelling
 *         Rely on the generic function: cs_gwf_tracer_add_terms_t
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_add_terms(cs_gwf_tracer_t     *tracer);

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
cs_gwf_tracer_setup(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant,
                    cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display the main features related to a tracer
 *
 * \param[in]  tracer   pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_log_setup(const cs_gwf_tracer_t     *tracer);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TRACER_H__ */
