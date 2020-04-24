#ifndef __CS_SOLIDIFICATION_H__
#define __CS_SOLIDIFICATION_H__

/*============================================================================
 * Header to handle the solidification module with CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_equation.h"
#include "cs_navsto_param.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_solidification.h

  \brief Structure and routines handling the solidification module dedicated to
         the resolution of electro-magnetic equations

*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Structure and type definitions
 *============================================================================*/

typedef cs_flag_t  cs_solidification_model_t;

/*! \enum cs_navsto_param_model_bit_t
 *  \brief Bit values for physical modelling related to the Navier-Stokes system
 *  of equations
 *
 * \var CS_SOLIDIFICATION_MODEL_STOKES
 * Stokes equations (mass and momentum) with the classical choice of variables
 * i.e. velocity and pressure. Mass density is assumed to be constant.
 *
 * \var CS_SOLIDIFICATION_MODEL_NAVIER_STOKES
 * Navier-Stokes equations: mass and momentum with a constant mass density
 *
 * \var CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87
 * Modelling introduced in Voller and Prakash entitled:
 * "A fixed grid numerical modelling methodology for convection-diffusion mushy
 * region phase-change problems" Int. J. Heat Transfer, 30 (8), 1987.
 * No tracer. Only physical constants describing the solidification process are
 * used.
 *
 * \var CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 * The tracer is composed by an alloy with two chemical constituents
 */

typedef enum {

  /* Main modelling for the dynamic
     ------------------------------ */

  CS_SOLIDIFICATION_MODEL_STOKES                  = 1<<0, /* =   1 */
  CS_SOLIDIFICATION_MODEL_NAVIER_STOKES           = 1<<1, /* =   2 */

  /* Main modelling for the thermal system
     ------------------------------------- */

  CS_SOLIDIFICATION_MODEL_USE_TEMPERATURE         = 1<<2, /* =   4 */
  CS_SOLIDIFICATION_MODEL_USE_ENTHALPY            = 1<<3, /* =   8 */

  /* Solidification modelling
     ------------------------ */

  CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87       = 1<<4, /* =  16 */
  CS_SOLIDIFICATION_MODEL_BINARY_ALLOY            = 1<<5, /* =  32 */

} cs_solidification_model_bit_t;


/*!
 * @name Flags specifying automatic post-processing for the solidification
 *        module
 * @{
 *
 * \def CS_SOLIDIFICATION_POST_CBULK_ADIM
 *
 * \brief Compute and post-process (C_bulk - C_0)/C_0
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated
 * C_0 is the reference concentration
 *
 * \def CS_SOLIDIFICATION_POST_CLIQ
 * \brief Compute and post-process Cliq = C_l*g_l
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated.
 * g_l is the liquid fraction and C_l is the solute distribution (wt %)
 *
 * \def CS_SOLIDIFICATION_POST_CLIQ_ADIM
 * \brief Compute Cliq = C_l*g_l and post-process (Cliq - C_O)/C_0
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *       is activated.
 * g_l is the liquid fraction and C_l is the solute distribution (wt %)
 * C_0 is the reference concentration
 *
 * \def CS_SOLIDIFICATION_POST_CELL_STATE
 * \brief State related to each cell between (solid, mushy, liquid or eutectic)
 *
 * \def CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE
 * \brief Activate the (volumic) post-processing of the liquidus temperature
 *        in each cell
 *
 * \def CS_SOLIDIFICATION_POST_SEGREGATION_INDEX
 * \brief Activate the computation and output in the file solidification.dat
 *        for each time step of the segregation index defined by
 *        sqrt( 1/|Domaine| * \int_{Domain} ((C_bulk - C_0)/C_0)**2 )
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated
 *
 * \def CS_SOLIDIFICATION_POST_SOLID_FRACTION_PORTION
 * \brief Activate the computation and output in the file solidification.dat
 *        for each time step of the integral over the computational domain
 *        of the solid fraction divided by the volume of the domain.
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated
 *
 * \def CS_SOLIDIFICATION_ADVANCED_ANALYSIS
 * \brief Activate a set of post-processing (Advanced usage. Only for the
 * understanding of the solidification process)
*/

#define CS_SOLIDIFICATION_POST_CBULK_ADIM              (1 << 0) /* =   1 */
#define CS_SOLIDIFICATION_POST_CLIQ                    (1 << 1) /* =   2 */
#define CS_SOLIDIFICATION_POST_CLIQ_ADIM               (1 << 2) /* =   4 */
#define CS_SOLIDIFICATION_POST_CELL_STATE              (1 << 3) /* =   8 */
#define CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE    (1 << 4) /* =  16 */
#define CS_SOLIDIFICATION_POST_SEGREGATION_INDEX       (1 << 5) /* =  32 */
#define CS_SOLIDIFICATION_POST_SOLID_FRACTION_PORTION  (1 << 6) /* =  64 */
#define CS_SOLIDIFICATION_ADVANCED_ANALYSIS            (1 << 7) /* = 128 */

/*!
 * @name Flags specifying numerical options specific to the solidification
 *       module
 * @{
 *
 * \def CS_SOLIDIFICATION_SOLUTE_WITH_ADVECTIVE_SOURCE_TERM
 * \brief The solute equation related to the transport of the bulk concentration
 * is treated with a source term related to an explicit advection of the
 * quantity (C - Cl). The default behavior is to add a weighting coefficient
 * to the (implicit) advection term related to the liquid fraction
 *
 * \def CS_SOLIDIFICATION_UPDATE_GL_WITH_TAYLOR_EXPANSION
 * \brief The update of the liquid fraction using a Taylor expansion in time
 *   dgl/dt = dgl/dT*(dT/dt) + dgl/dC*(dC/dt)
 *
 * \def CS_SOLIDIFICATION_UPDATE_SOURCE_TERM_BY_STEP
 * \brief Update the source term related to the thermal equation considering
 * a path between the initial and final state. For each step, one considers to
 * add or not a source term.
 *
 * \def CS_SOLIDIFICATION_UPDATE_EUTECTIC_VOLLER
 * \brief Update the liquid fraction according to the model introduced in
 *  Voller & Prakash'89
*/

#define CS_SOLIDIFICATION_SOLUTE_WITH_ADVECTIVE_SOURCE_TERM  (1 << 0) /* =  1 */
#define CS_SOLIDIFICATION_UPDATE_GL_WITH_TAYLOR_EXPANSION    (1 << 1) /* =  2 */
#define CS_SOLIDIFICATION_UPDATE_SOURCE_TERM_BY_STEP         (1 << 2) /* =  4 */
#define CS_SOLIDIFICATION_UPDATE_EUTECTIC_VOLLER             (1 << 3) /* =  8 */

/*!
 * @}
 */

typedef struct _solidification_t  cs_solidification_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if solidification module is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_solidification_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the solidification module
 *
 * \param[in]  model            type of modelling
 * \param[in]  options          flag to handle optional parameters
 * \param[in]  post_flag        predefined post-processings
 * \param[in]  boundaries       pointer to the domain boundaries
 * \param[in]  algo_coupling    algorithm used for solving the NavSto system
 * \param[in]  ns_option        option flag for the Navier-Stokes system
 * \param[in]  ns_post_flag     predefined post-processings for Navier-Stokes
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_activate(cs_solidification_model_t      model,
                           cs_flag_t                      options,
                           cs_flag_t                      post_flag,
                           const cs_boundary_t           *boundaries,
                           cs_navsto_param_coupling_t     algo_coupling,
                           cs_flag_t                      ns_option,
                           cs_flag_t                      ns_post_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value of the epsilon parameter used in the forcing term
 *         of the momemtum equation
 *
 * \param[in]  forcing_eps    epsilon used in the penalization term to avoid a
 *                            division by zero
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_forcing_eps(cs_real_t    forcing_eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which described the Voller and
 *         Prakash modelling
 *
 * \param[in]  t_solidus      solidus temperature (in K)
 * \param[in]  t_liquidus     liquidus temperatur (in K)
 * \param[in]  latent_heat    latent heat
 * \param[in]  s_das          secondary dendrite space arms
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_voller_model(cs_real_t    t_solidus,
                                   cs_real_t    t_liquidus,
                                   cs_real_t    latent_heat,
                                   cs_real_t    s_das);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main numerical parameters which described a solidification
 *         process with a binary alloy (with component A and B)
 *
 * \param[in]  n_iter_max    max.number of iterations for the C/T equations
 * \param[in]  g_l_eps       tolerance requested between two iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_binary_alloy_param(int             n_iter_max,
                                         double          g_l_eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which described a solidification
 *         process with a binary alloy (with component A and B)
 *         Add a transport equation for the solute concentration to simulate
 *         the conv/diffusion of the alloy ratio between the two components of
 *         the alloy
 *
 * \param[in]  name          name of the binary alloy
 * \param[in]  varname       name of the unknown related to the tracer eq.
 * \param[in]  conc0         reference mixture concentration
 * \param[in]  beta          solutal dilatation coefficient
 * \param[in]  kp            value of the distribution coefficient
 * \param[in]  mliq          liquidus slope for the solute concentration
 * \param[in]  t_eutec       temperature at the eutectic point
 * \param[in]  t_melt        phase-change temperature for the pure material (A)
 * \param[in]  solute_diff   solutal diffusion coefficient in the liquid
 * \param[in]  latent_heat   latent heat
 * \param[in]  s_das         secondary dendrite arm spacing
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_binary_alloy_model(const char     *name,
                                         const char     *varname,
                                         cs_real_t       conc0,
                                         cs_real_t       beta,
                                         cs_real_t       kp,
                                         cs_real_t       mliq,
                                         cs_real_t       t_eutec,
                                         cs_real_t       t_melt,
                                         cs_real_t       solute_diff,
                                         cs_real_t       latent_heat,
                                         cs_real_t       s_das);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the solidification module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the Solidification module
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the solidification
 *         module
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_finalize_setup(const cs_cdo_connect_t       *connect,
                                 const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the solidification module in the log file dedicated to
 *         the setup
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context structure used to build the algebraic system
 *         This is done after the setup step.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_initialize(const cs_mesh_t              *mesh,
                             const cs_cdo_connect_t       *connect,
                             const cs_cdo_quantities_t    *quant,
                             const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the solidification module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_compute(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the solidification module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the solidification module.
 *         Prototype of this function is fixed since it is a function pointer
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
cs_solidification_extra_post(void                      *input,
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

#endif /* __CS_SOLIDIFICATION_H__ */
