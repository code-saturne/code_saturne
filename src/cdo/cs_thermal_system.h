#ifndef __CS_THERMAL_SYSTEM_H__
#define __CS_THERMAL_SYSTEM_H__

/*============================================================================
 * Functions to handle cs_thermal_system_t structure
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

#include "cs_advection_field.h"
#include "cs_equation.h"
#include "cs_field.h"
#include "cs_property.h"
#include "cs_mesh.h"
#include "cs_source_term.h"
#include "cs_time_step.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Generic name given to fields and equations related to this module */

#define CS_THERMAL_EQNAME           "thermal_equation"
#define CS_THERMAL_CP_NAME          "thermal_capacity"
#define CS_THERMAL_LAMBDA_NAME      "thermal_conductivity"

/*!
 * @name Flags specifying automatic post-processing for the thermal module
 * @{
 *
 * \def CS_THERMAL_POST_ENTHALPY
 * \brief Create a field storing the enthalpy for post-processing and
 *        extra-operation
 *
 *
 */

#define CS_THERMAL_POST_ENTHALPY               (1 << 0) /* 1 */

/*!
 * @}
 */

/*=============================================================================
 * Structure and type definitions
 *============================================================================*/

typedef cs_flag_t  cs_thermal_model_type_t;

/*! \enum cs_thermal_model_type_bit_t
 *  \brief Bit values for physical modelling related to thermal system
 *
 * \def CS_THERMAL_MODEL_STEADY
 * \brief Disable the unsteady term of the thermal equation
 *
 * \def CS_THERMAL_MODEL_NAVSTO_ADVECTION
 * \brief Add an advection term arising from the velocity field solution
 *        of the Navier-Stokes equations
 *
 * \def CS_THERMAL_MODEL_USE_TEMPERATURE
 * \brief The thermal equation solved using the temperature variable
 *        (the default choice)
 *
 * \def CS_THERMAL_MODEL_USE_ENTHALPY
 * \brief The thermal equation solved using the enthalpy variable
 *        instead of the temperature (the default choice)
 *
 * \def CS_THERMAL_MODEL_USE_TOTAL_ENERGY
 * \brief The thermal equation solved using the total energy variable
 *        instead of the temperature (the default choice)
 *
 * \def CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY
 * \brief The thermal equation solved using an anisotropic conductivity
 *        instead of an isotropic one (the default choice)
 *
 */

typedef enum {

  CS_THERMAL_MODEL_STEADY                     = 1<<0,  /* =  1 */
  CS_THERMAL_MODEL_NAVSTO_ADVECTION           = 1<<1,  /* =  2 */

  /* Main variable to consider (by default the temperature in Kelvin)
     ---------------------------------------------------------------- */

  CS_THERMAL_MODEL_USE_TEMPERATURE            = 1<<2,  /* =  4 */
  CS_THERMAL_MODEL_USE_ENTHALPY               = 1<<3,  /* =  8 */
  CS_THERMAL_MODEL_USE_TOTAL_ENERGY           = 1<<4,  /* = 16 */

  /* Treatment of the diffusion term
     ------------------------------- */

  CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY   = 1<<5,  /* = 32 */

  /* Additional bit settings
     ----------------------- */

  CS_THERMAL_MODEL_IN_CELSIUS                 = 1<<7   /* = 128 */

} cs_thermal_model_type_bit_t;

/* Set of parameters related to the thermal module */

typedef struct {

  cs_flag_t   model;                  /* Choice of modelling */
  cs_flag_t   numeric;                /* General numerical options */
  cs_flag_t   post;                   /* Post-processing options */

  /* Equation associated to this module */
  /* ---------------------------------- */

  cs_equation_t  *thermal_eq;

  /* Properties associated to this module */
  /* ------------------------------------ */

  cs_property_t  *unsteady_property; /* Property of the unsteady term: rho.cp */
  cs_property_t  *lambda;            /* Thermal conductivity */
  cs_property_t  *cp;                /* Heat capacity */
  cs_property_t  *rho;               /* Mass density (may be shared) */
  cs_property_t  *kappa;             /* lambda/cp may be NULL*/

  /* value of lambda/cp in each cell (allocated if the related property kappa
   * is not NULL and defined by array) */
  cs_real_t      *kappa_array;

  /* Fields associated to this module */
  /* -------------------------------- */

  cs_field_t     *temperature;
  cs_field_t     *enthalpy;
  cs_field_t     *total_energy;

  /* Additional members */
  /* ------------------ */

  cs_real_t       ref_temperature;

  /* N.B.: Other reference values for properties are stored within each
   * property structure */

} cs_thermal_system_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the value of the reference temperature associated to a
 *        thermal system.
 *
 * \return the value of the reference temperature
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_system_get_reference_temperature(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the value of the reference temperature associated to the
 *        thermal system.
 *
 * \param[in]  ref     value of the reference temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_set_reference_temperature(cs_real_t    ref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the model flag related to a thermal system
 *
 * \return a flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_thermal_system_get_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Does the thermal system rely on the advection field associated to
 *        the Navier-Stokes equations ?
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_thermal_system_needs_navsto(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the thermal system has been activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_thermal_system_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the thermal system
 *
 * \param[in] model     model flag related to the thermal system
 * \param[in] numeric   (optional) numerical flag settings
 * \param[in] post      (optional) post-processing flag settings
 *
 * \return a pointer to a new allocated cs_thermal_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_thermal_system_t *
cs_thermal_system_activate(cs_flag_t         model,
                           cs_flag_t         numeric,
                           cs_flag_t         post);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the thermal system
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main equation related to the thermal system
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_thermal_system_get_equation(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the current temperature at face values
 *
 * \return the pointer to the array of face values.
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_thermal_system_get_face_temperature(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the thermal system
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the thermal system
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_finalize_setup(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_time_step_t       *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve a steady-state thermal system
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute_steady_state(const cs_mesh_t              *mesh,
                                       const cs_cdo_connect_t       *connect,
                                       const cs_cdo_quantities_t    *quant,
                                       const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the thermal system
 *
 * \param[in]  cur2prev   true="current to previous" operation is performed
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute(bool                          cur2prev,
                          const cs_mesh_t              *mesh,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set an initial values for all quantities related to this module
 *         This is done after the setup step.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_init_values(const cs_mesh_t             *mesh,
                              const cs_cdo_connect_t      *connect,
                              const cs_cdo_quantities_t   *quant,
                              const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the quantities related to the thermal module
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_update(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts,
                         bool                         cur2prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the thermal system
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the thermal system.
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_thermal_system_t structure)
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
cs_thermal_system_extra_post(void                      *input,
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
/*!
 * \brief  Summary of the main options related to cs_thermal_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_THERMAL_SYSTEM_H__ */
