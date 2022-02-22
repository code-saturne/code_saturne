#ifndef __CS_GWF_TRACER_H__
#define __CS_GWF_TRACER_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
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
#include "cs_base.h"
#include "cs_equation.h"
#include "cs_gwf_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Typedef definition
 *============================================================================*/

typedef struct _gwf_tracer_t  cs_gwf_tracer_t;

/*============================================================================
 * Public function pointer prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function to update the first setup stage (the one done before
 *        building mesh and its related quantities) for a tracer equation
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_init_setup_t) (cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function to finalize the setup of parameters related to a
 *        tracer equation. At this stage, mesh and its related quantities have
 *        been built.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      adv           pointer to an advection field structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_finalize_setup_t) (const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  const cs_adv_field_t        *adv,
                                  cs_gwf_tracer_t             *tracer);

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
(cs_gwf_tracer_free_context_t) (cs_gwf_tracer_t      *tracer);

/*============================================================================
 * Structure definitions
 *============================================================================*/

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

  const cs_field_t  *darcy_velocity_field;

} cs_gwf_tracer_default_context_t;


/* Set of parameters describing a tracer structure */
/* ----------------------------------------------- */

struct _gwf_tracer_t{

  cs_gwf_model_type_t          hydraulic_model;

  /* Physical modelling adopted for this tracer */

  cs_gwf_tracer_model_t        model;

  cs_field_t                  *diffusivity; /* NULL if no diffusion term is
                                               build in the tracer equation */
  int                          reaction_id; /* id related to the reaction
                                               term in the tracer equation */

  /*! \var eq
   *  \brief pointer to the related equation structure
   */

  cs_equation_t                  *equation;

  /* Pointer to a context structure according to the model */

  void                           *context;

  /* Pointers to functions */

  /*!
   * \var update_diff_tensor
   *      Function used to update the diffusion tensor (dispersion + diffusion)
   *
   * \var update_precipitation
   *      Function used to update the quantities related to the precipitation
   * model
   *
   * \var finalize_setup
   *       This is a function pointer to finalize the setup of a tracer
   * equation. There is a function pointer by default but this can be
   * overloaded by a user-defined function in the case of a user-defined
   * tracer.
   *
   * \var init_setup
   *      This is a function pointer to initialize the setup (adding terms in
   * an equation). At this stage, the mesh has not been loaded.  There is a
   * function pointer by default but this can be overloaded by a user-defined
   * function in the case of a user-defined tracer.
   *
   * \var free_context
   *      Function to free quantities or structure associated to the context
   * structure of a tracer.
   */

  cs_gwf_tracer_update_t           *update_diff_tensor;
  cs_gwf_tracer_update_t           *update_precipitation;
  cs_gwf_tracer_finalize_setup_t   *finalize_setup;
  cs_gwf_tracer_init_setup_t       *init_setup;
  cs_gwf_tracer_free_context_t     *free_context;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the cs_gwf_tracer_t structure associated to
 *         the name given as parameter
 *
 * \param[in]  eq_name    name of the tracer equation
 *
 * \return the pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_by_name(const char   *eq_name);

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
 * \param[in]   tr_model        model related to this tracer
 * \param[in]   gwf_model       main model for the GWF module
 * \param[in]   eq_name         name of the tracer equation
 * \param[in]   var_name        name of the related variable
 * \param[in]   adv_field       pointer to a cs_adv_field_t structure
 * \param[in]   init_setup      function pointer (predefined prototype)
 * \param[in]   finalize_setup  function pointer (predefined prototype)
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_add(cs_gwf_tracer_model_t            tr_model,
                  cs_gwf_model_type_t              gwf_model,
                  const char                      *eq_name,
                  const char                      *var_name,
                  cs_adv_field_t                  *adv_field,
                  cs_gwf_tracer_init_setup_t      *init_setup,
                  cs_gwf_tracer_finalize_setup_t  *finalize_setup);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all tracers
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the max. value of the theta parameter associated to a time
 *        scheme. Loop on all tracer equations.
 *
 * \return the computed value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tracer_get_time_theta_max(void);

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
cs_gwf_tracer_set_main_param(cs_gwf_tracer_t   *tracer,
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
cs_gwf_tracer_set_precip_param(cs_gwf_tracer_t   *tracer,
                               const char        *soil_name,
                               double             conc_w_star);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup step for tracer equations. Soils and equation
 *        parameters are defined at this stage.
 *        Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the tracer setup
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_finalize_setup(const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the diffusion tensor related to each tracer equation
 *
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_update_diff_tensor(cs_real_t                    t_eval,
                                 const cs_mesh_t             *mesh,
                                 const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display the main features related to each tracer
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_log_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state for all tracer equations.
 *         Nothing is done if all equations are unsteady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_compute_steady_all(const cs_mesh_t              *mesh,
                                 const cs_time_step_t         *time_step,
                                 const cs_cdo_connect_t       *connect,
                                 const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for all tracer equations.
 *         Nothing is done if all equations are steady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_compute_all(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add terms to the algebraic system related to a tracer equation
 *         according to the settings.
 *         Case of the default tracer modelling
 *         Rely on the generic function: cs_gwf_tracer_add_terms_t
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_default_init_setup(cs_gwf_tracer_t     *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation case of
 *         a fully saturated flow model
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      adv           pointer to an advection field structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_sat_finalize_setup(const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *quant,
                                 const cs_adv_field_t        *adv,
                                 cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation in case of
 *         an unsaturated flow model
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      adv           pointer to an advection field structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_unsat_finalize_setup(const cs_cdo_connect_t      *connect,
                                   const cs_cdo_quantities_t   *quant,
                                   const cs_adv_field_t        *adv,
                                   cs_gwf_tracer_t             *tracer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a given set of cells of the field related
 *         to a tracer equation. This integral turns out to be exact for linear
 *         functions.
 *
 * \param[in]    connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]    cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]    tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]    z_name    name of the volumic zone where the integral is done
 *                         (if NULL or "" all cells are considered)
 *
 * \return the value of the integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tracer_integrate(const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_gwf_tracer_t      *tracer,
                        const char                 *z_name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TRACER_H__ */
