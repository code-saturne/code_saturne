#ifndef __CS_CDO_TURBULENCE_H__
#define __CS_CDO_TURBULENCE_H__

/*============================================================================
 * Functions to handle the resolution of the turbulence modelling
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

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_property.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_NAVSTO_TOTAL_VISCOSITY    "total_viscosity"
#define CS_NAVSTO_LAM_VISCOSITY      "laminar_viscosity"
#define CS_NAVSTO_TURB_VISCOSITY     "turbulent_viscosity"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_turbulence_param_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         turbulence modelling. Several members are structures defined in
 *         cs_turbulence_model.h as a global variable. The purpose of this
 *         structure is to store all parameters in one place.
 */

typedef struct {

  /*! \var model
   * Main set of parameters to handle turbulence modelling. This
   * structure is shared with the legacy part.
   */

  cs_turb_model_t            *model;

  /*! \var rans_param
   * Main set of parameters to handle RANS modelling. This
   * structure is shared with the legacy part.
   * RANS means Reynolds Average Navier-Stokes
   */

  cs_turb_rans_model_t       *rans_param;

  /*! \var les_param
   * Main set of parameters to handle LES modelling. This
   * structure is shared with the legacy part.
   * LES means Large Eddy Simulation
   */

  cs_turb_les_model_t        *les_param;

  /*! \var reference_values
   *  Set of reference values associated to the turbulence modelling
   */

  cs_turb_ref_values_t       *reference_values;

} cs_turbulence_param_t;

/*! \struct cs_turbulence_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         turbulence modelling. Several members are structures defined in
 *         cs_turbulence_model.h
 */

typedef struct _cs_turbulence_t cs_turbulence_t;

/*============================================================================
 * Function pointer definition
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the context structure related to a given
 *         turbulence modelling
 *
 * \param[in]      tbm       pointer to a \ref cs_turb_model_t structure
 *
 * \return a pointer to a new allocated turbulence context structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_turb_init_context_t)(const cs_turb_model_t   *tbm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a given turbulence modelling
 *
 * \param[in, out]  tbc       pointer to a structure cast on-the-fly to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_turb_free_context_t)(void     *turb_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the turbulence
 *         model. This means that all related equations are built and then
 *         solved.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in, out] tbs        pointer to turbulence structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_turb_compute_t)(const cs_mesh_t            *mesh,
                    const cs_cdo_connect_t     *connect,
                    const cs_cdo_quantities_t  *quant,
                    const cs_time_step_t       *time_step,
                    cs_turbulence_t            *tbs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update properties, arrays related to the turbulent variables
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in]      tbs        pointer to a \ref cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_turb_update_t)(const cs_mesh_t                *mesh,
                   const cs_cdo_connect_t         *connect,
                   const cs_cdo_quantities_t      *quant,
                   const cs_time_step_t           *time_step,
                   const cs_turbulence_t          *tbs);

/*! \struct cs_turbulence_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         turbulence modelling. Several members are structures defined in
 *         cs_turbulence_model.h
 */

struct _cs_turbulence_t {

  /*!
   * @name Turbulence modelling
   * Set of parameters to handle turbulence modelling.
   * @{
   *
   * \var param
   * Main set of parameters to handle turbulence modelling. The members
   * of this structure are shared with the legacy part.
   */

  cs_turbulence_param_t     *param;

  /*! \var mom_eq
   * Pointer to the momentum equation structure.
   */
  cs_equation_t             *mom_eq;

  /*!
   * @}
   * @name Main properties related to the turbulence modelling
   */

  /*! \var rho
   *  mass density
   */

  cs_property_t              *rho;

  /*! \var mu_t
   *  total viscosity (dynamic turbulent + laminar)
   */

  cs_property_t              *mu_tot;

  /*! \var mu_l
   *  laminar viscosity
   */

  cs_property_t              *mu_l;

  /*! \var mu_t
   *  dynamic turbulent viscosity
   */

  cs_property_t              *mu_t;

  /*! \var mu_tot_array
   *  Array storing the value of the total viscosity in each cell
   */
  cs_real_t                  *mu_tot_array;

  /*!
   * @}
   * @name Main related fields
   * @{
   */

  /* \var mu_t_field
   * Field related to the turbulent viscosity
   */

  cs_field_t                 *mu_t_field;

  /*! \var rij
   *  Reynolds stress tensor
   */

  cs_field_t                 *rij;

  /*!
   * @}
   * @name Main structure and function pointers
   *
   * \brief The context structure is a structure cast on-the-fly and function
   * pointers are a limited set of functions associated to the main operations
   * for this structure: initialization, destruction, computation and update.
   *
   * @{
   */

  /*! \var context
   * Context related to the turbulence modelling. This structure is cast
   * on-the-fly according to the modelling choices
   */

  void                         *context;

  /*! \var init_context
   * Function pointer to initialize the context structure
   */

  cs_turb_init_context_t       *init_context;

  /*! \var free_context
   * Function pointer to de-allocate memory owned by the context structure
   */

  cs_turb_free_context_t       *free_context;

  /*! \var compute
   * Function pointer to compute all variables related to a turbulence model
   */

  cs_turb_compute_t            *compute;

  /*! \var compute_steady
   * Function pointer to compute in steady mode all variables related to a
   * turbulence model
   */

  cs_turb_compute_t            *compute_steady;

  /*! \var update
   * Function pointer to perform the update step (properties or arrays
   * associated to the variables of a turbulence model are updated)
   */

  cs_turb_update_t             *update;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure storing the set of parameters for the
 *         turbulence modelling
 *
 * \return a pointer to a new allocated cs_turbulence_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_turbulence_param_t *
cs_turbulence_param_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure managing the turbulence modelling
 *
 * \param[in] tbp       pointer to a \ref cs_turbulence_param_t structure
 *
 * \return a pointer to a new allocated cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

cs_turbulence_t *
cs_turbulence_create(cs_turbulence_param_t    *tbp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure managing the turbulence modelling
 *
 * \param[in, out]  p_turb_struct   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_free(cs_turbulence_t   **p_turb_struct);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the structure managing the turbulence modelling
 *
 * \param[in, out]  tbs     pointer to the structure to initialize
 * \param[in]       mom_eq  pointer to the structure mom_eq
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_setup(cs_turbulence_t     *tbs,
                         cs_equation_t       *mom_eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup of the turbulence modelling and especially the
 *         equations/properties and other related quantities
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] tbs        pointer to the turbulence main structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_finalize_setup(const cs_mesh_t            *mesh,
                             const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             const cs_time_step_t       *time_step,
                             cs_turbulence_t            *tbs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the values of quantities related to a turbulence model.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] tbs        pointer to the turbulence main structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_values(const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *time_step,
                          cs_turbulence_t             *tbs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the context structure related to the
 *         k-epsilon turbulence model
 *
 * \param[in]  tbm         structure which defines the turbulence model
 *
 * \return a pointer to a new allocated context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_turb_init_k_eps_context(const cs_turb_model_t      *tbm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to the k-epsilon turbulence model
 *
 * \param[in, out]  tbc   pointer to a structure context cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_turb_free_k_eps_context(void     *tbc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the turbulence
 *         model. This means that all related equations are built and then
 *         solved.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in, out] tbs        pointer to turbulence structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_k_eps(const cs_mesh_t            *mesh,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      const cs_time_step_t       *time_step,
                      cs_turbulence_t            *tbs);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update for the current time step the new state for the turbulence
 *         model. This is used to update the turbulent viscosity.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in]      tbs        pointer to a \ref cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_update_k_eps(const cs_mesh_t              *mesh,
                     const cs_cdo_connect_t       *connect,
                     const cs_cdo_quantities_t    *quant,
                     const cs_time_step_t         *time_step,
                     const cs_turbulence_t        *tbs);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_TURBULENCE_H__ */
