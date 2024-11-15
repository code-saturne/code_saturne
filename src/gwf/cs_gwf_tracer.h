#ifndef __CS_GWF_TRACER_H__
#define __CS_GWF_TRACER_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
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
#include "cs_base.h"
#include "cs_equation.h"
#include "cs_gwf_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Typedef definition
 *============================================================================*/

/*! \struct cs_gwf_tracer_t
 *
 *  \brief Set of parameters describing a tracer structure
 */

typedef struct _gwf_tracer_t  cs_gwf_tracer_t;

/* \struct cs_gwf_tracer_default_context_t
 *
 * \brief Set of parameters related to a tracer equation attached to a standard
 *        modelling
 */

typedef struct _gwf_tracer_default_context_t  cs_gwf_tracer_default_context_t;

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
 * \brief  Generic function to update the properties related to a tracer.
 *         This function depends on a numerical scheme and a physical model.
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    null or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_update_t) (cs_gwf_tracer_t             *tracer,
                          void                        *context,
                          const cs_time_step_t        *ts,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the quantity of tracer using an integral of the tracer
 *        concentration over a given set of cells (the cells of the volume
 *        zone). Several terms can be computed.
 *        The quantity is given in moles. A parallel operation (a sum
 *        reduction) is performed inside this function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq        equation related to a tracer
 * \param[in]      zone      pointer to a volume zone structure
 * \param[in]      context   pointer to a context structure for a tracer
 * \param[in, out] results   resulting array of values
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_tracer_integrate_t)(const cs_cdo_connect_t          *connect,
                            const cs_cdo_quantities_t       *cdoq,
                            const cs_equation_t             *eq,
                            const cs_zone_t                 *zone,
                            void                            *context,
                            double                           results[]);

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Structure to handle of radioactive decay chain  */

typedef struct {

  char                  *name;       /* Name of the decay chain */
  int                    id;         /* Position in the array of chain decays */

  int                    n_tracers;  /* Number of tracers in the decay chain */
  cs_gwf_tracer_unit_t   unit;       /* Type of unit measure:
                                      *  - mole (classical tracer)
                                      *  - Becquerel (activity)
                                      */

  cs_gwf_tracer_t      **tracers;   /* Array of tracers associated to the
                                       decay chain. tracers[0] is the parent
                                       without any ancestor. */

  cs_xdef_t            **st_defs;   /* List of source term definitions */

} cs_gwf_tracer_decay_chain_t;

/* Set of parameters related to a tracer equation attached to a standard
   modelling */

struct _gwf_tracer_default_context_t {

  /* Common settings shared by all physical modelling */
  /* ------------------------------------------------ */

  double     decay_coef;   /* First order decay coefficient (related to the
                              reaction term). This value is intrinsic to the
                              component (a radioactive element for instance)
                              and not to the soil. */

  /* The following parameters are defined for each soil (arrays of size equal
   * to n_soils) */

  double    *rho_bulk;      /* bulk density (kg.m^-3) */
  double    *kd0;           /* reference value of the distribution coefficient
                               (m^".kg^-1) */
  double    *rho_kd;        /* Derived quantity: rho_bulk*kd0 */

  double    *alpha_l;       /* Longitudinal dispersivity */
  double    *alpha_t;       /* Transversal dispersivity */

  double    *wmd;           /* Water molecular diffusivity (m^2.s^-1) */

  /* Precipitation members (set to null if not used) */
  /* ----------------------------------------------- */

  double       *conc_l_star;    /* maximal value of the concentration of
                                 * tracer in the liquid phase in mol/m^3. There
                                 * is one user-defined value for each soil. The
                                 * exceeded quantities are stored in the solid
                                 * phase (-> precip_mass). These values
                                 * corresponds to the user settings
                                 */

  cs_real_t    *precip_mass;    /* array storing the mass of precipitate
                                 * (solid) in the dedicated auxiliary
                                 * storage. The size of this array may vary
                                 * w.r.t. to the discretization scheme.
                                 */

  cs_field_t   *precip_field;   /* field structure storing the (interpolated)
                                 * values of the concentration of precipitate
                                 * in mol/kg in each cell.
                                 */

  /* Sorption members (set to null if not used) */
  /* ------------------------------------------ */

  double       *k0_plus;        /* kinetic coefficient towards site 2 locations
                                   (m^3.kg^-1.s^-1) */
  double       *k0_minus;       /* kinetic coefficient from site 2 locations
                                   (s^-1) */

  cs_real_t    *conc_site2;     /* array allocated to n_cells whatever the
                                   space discretization is */

  /* Variables used for the update of physical properties (shared pointers) */

  const cs_field_t    *darcy_velocity_field;

};

/* Set of parameters describing a tracer structure */
/* ----------------------------------------------- */

struct _gwf_tracer_t {

  /*!
   * @name Physical modelling information for a tracer
   * @{
   */

  /*! \var hydraulic_model
   *       Type of hydraulic model to consider. This is an information shared
   *       with the main gorundwater flow structure
   */

  cs_gwf_model_type_t     hydraulic_model;

  /*! \var model
   *       Type of tracer model to consider. 0 corresponds to the default
   *       behavior.
   */

  cs_gwf_tracer_model_t   model;

  /*! \var diffusivity
   *       Field related to the property associated to the diffusion term.
   *       null if no diffusion term is build in the tracer equation.
   */

  cs_field_t             *diffusivity;

  /*! \var reaction_id
   *       Since there can be several reaction terms associated to an
   *       equation. One stores the id related to the reaction term wich is
   *       automatically added when a radioactive tracer is considered.
   */

  int                     reaction_id;

  /*!
   * @}
   * @name Other members
   * @{
   */

  /*! \var equation
   *       Pointer to the related equation structure
   */

  cs_equation_t          *equation;

  /*! \var context
   *       Pointer to a context structure cast on-the-fly according to the
   *       model.
   */

  void                   *context;

  /*! \var chain_position_id
   *       Position of the current tracer inside a decay chain. The default
   *       value is -1 meaning that this is not a tracer associated to a decay
   *       chain.
   */

  int                     chain_position_id;

  /*! \var chain_id
   *       id in the array of decay chains. (-1 if not  associated to a decay
   *       chain; this is the default behavior)
   */

  int                     chain_id;

  /*!
   * @}
   * @name Function pointers associated to a tracer
   * @{
   */

  /*!
   * \var update_diff_pty
   *      Function used to update the diffusion property which is a tensor in
   *      the most generic case (dispersion + diffusion)
   *
   * \var update_precipitation
   *      Function used to update the quantities related to the precipitation
   *      model
   *
   * \var update_decay_chain_st
   *      Function used to update the source term induced by the parent in a
   *      decay chain
   *
   * \var integrate
   *      Function to compute the quantity of tracer inside a volume. The way
   *      to compute this quantity may be optimized according to the hydraulic
   *      model or the tracer modelling.
   *
   * \var init_setup
   *      This is a function pointer to initialize the setup (adding terms in
   *      an equation). At this stage, the mesh has not been loaded.  There is
   *      a function pointer by default but this can be overloaded by a
   *      user-defined function in the case of a user-defined tracer.
   *
   * \var finalize_setup
   *      This is a function pointer to finalize the setup of a tracer
   *      equation. There is a function pointer by default but this can be
   *      overloaded by a user-defined function in the case of a user-defined
   *      tracer.
   *
   * \var free_context
   *      Function to free quantities or structure associated to the context
   *      structure of a tracer.
   */

  cs_gwf_tracer_update_t           *update_diff_pty;
  cs_gwf_tracer_update_t           *update_precipitation;
  cs_gwf_tracer_update_t           *update_decay_chain_st;

  cs_gwf_tracer_integrate_t        *integrate;

  cs_gwf_tracer_init_setup_t       *init_setup;
  cs_gwf_tracer_finalize_setup_t   *finalize_setup;
  cs_gwf_tracer_free_context_t     *free_context;

  /*!
   * @}
   */
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
 * \brief Create a new cs_gwf_tracer_t structure and initialize its members.
 *        This creation of a new tracer is fully done in the case of a default
 *        tracer. Additional settings has to be done in the case of a
 *        user-defined tracer.
 *
 *        Add a new equation related to the groundwater flow module. This
 *        equation is a specific transport equation. The tracer is advected
 *        thanks to the darcian velocity (in the liquid phase) which is given
 *        by the resolution of the Richards equation. Diffusion and reaction
 *        coefficients result from a physical modelling.
 *
 * \param[in] tr_model        model related to this tracer
 * \param[in] gwf_model       main model for the GWF module
 * \param[in] eq_name         name of the tracer equation
 * \param[in] var_name        name of the related variable
 * \param[in] adv_field       pointer to a cs_adv_field_t structure
 * \param[in] lambda          value of the first order decay coefficient
 * \param[in] chain_position  -1 if not used or the id in the chain position
 * \param[in] chain_id        -1 or id of the associated decay chain
 * \param[in] init_setup      function pointer (predefined prototype)
 * \param[in] finalize_setup  function pointer (predefined prototype)
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
                  double                           lambda,
                  int                              chain_position,
                  int                              chain_id,
                  cs_gwf_tracer_init_setup_t      *init_setup,
                  cs_gwf_tracer_finalize_setup_t  *finalize_setup);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all tracer structures and all decay chains
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
 * \brief Set the main parameters corresponding to a default modelling of a
 *        tracer transport equation for a specified soil
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or null if all
 *                                 soils are selected)
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_set_soil_param(cs_gwf_tracer_t   *tracer,
                             const char        *soil_name,
                             double             wmd,
                             double             alpha_l,
                             double             alpha_t,
                             double             distrib_coef);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a specified soil set the parameters corresponding to a
 *         precipitation modelling of a tracer transport
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or null if all
 *                                 soils are selected)
 * \param[in]      conc_l_star     value of the saturated concentration in the
 *                                 liquid phase
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_set_precip_param(cs_gwf_tracer_t   *tracer,
                               const char        *soil_name,
                               double             conc_l_star);

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
 * \brief Update the diffusion property related to each tracer equation
 *        The update strategy depends on the soil/tracer features and also
 *        on the hydraulic model.
 *
 * \param[in] ts        pointer to a cs_time_step_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 * \param[in] quant     pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_update_diff_pty(const cs_time_step_t        *ts,
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
 * \brief Compute the integral of the tracer concentration over a given set of
 *        cells. This integral gives the number of moles of tracer inside the
 *        related volume. Moreover, it turns out to be exact for linear
 *        functions. A parallel operation (a sum reduction) is performed
 *
 * \param[in]    connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]    cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]    tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]    z_name    name of the volumic zone where the integral is done
 *                         (if null or "" all cells are considered)
 *
 * \return the value of the integral (number of moles in the zone)
 *         parallel synchronization is done
 */
/*----------------------------------------------------------------------------*/

double
cs_gwf_tracer_integrate(const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_gwf_tracer_t      *tracer,
                        const char                 *z_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the quantity of tracer using the integral of the tracer
 *        concentration over a given set of cells. Two terms are computed: one
 *        for the quantity of moles inside the liquid phase and another one for
 *        the quantity of tracer inside the precipitation state (in moles). A
 *        parallel operation (a sum reduction) is performed inside this
 *        function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]      z_name    name of the volume zone where the integral is
 *                           done (if null or "" all cells are considered)
 * \param[in, out] results   array of values. [0]= the quantity of moles
 *                           in the liquid phase, [1]= the quantity of
 *                           moles inside the precipitation state
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_integrate_by_terms(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *cdoq,
                                 const cs_gwf_tracer_t      *tracer,
                                 const char                 *z_name,
                                 double                      results[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a decay chain structure to manage several linked tracers
 *
 * \param[in] n_tracers    number of tracers equations
 * \param[in] chain_name   name of the decay chain
 * \param[in] unit         type of unit used in the tracer equations
 *
 * \return a pointer to the new cs_gwf_tracer_decay_chain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_create_decay_chain(int                      n_tracers,
                                 const char              *chain_name,
                                 cs_gwf_tracer_unit_t     unit);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the decay chain structure associated to the given id
 *        If not found, it returns the null pointer.
 *
 * \param[in] id   id of the decay chain to retrieve
 *
 * \return a pointer to a new cs_gwf_tracer_decay_chain_t structure or null
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_decay_chain_by_id(int        id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the decay chain structure associated to the name given as
 *        parameter. If not found, it returns the null pointer.
 *
 * \param[in] chain_name   name of the decay chain
 *
 * \return a pointer to a new cs_gwf_tracer_decay_chain_t structure or null
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_decay_chain_by_name(const char      *chain_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the tracer structure for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a null
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_gwf_tracer_t structure or null
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_decay_chain_get_tracer(cs_gwf_tracer_decay_chain_t  *tdc,
                                     int                           id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the equation structure for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a null
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_equation_t structure or null
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_tracer_decay_chain_get_equation(cs_gwf_tracer_decay_chain_t  *tdc,
                                       int                           id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the equation parameters for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a null
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_equation_param_t structure or null
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_gwf_tracer_decay_chain_get_equation_param(cs_gwf_tracer_decay_chain_t  *tdc,
                                             int                           id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TRACER_H__ */
