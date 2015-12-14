#ifndef __CS_GROUNDWATER_H__
#define __CS_GROUNDWATER_H__

/*============================================================================
 * Compute the wall distance using the CDO framework
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Tag to build a flag dedicated to the groundwater module */
#define CS_GROUNDWATER_POST_MOISTURE  (1 <<  0) //    1: post the moisture content

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of predefined modelling for the groundwater flows */
typedef enum {

  CS_GROUNDWATER_MODEL_SATURATED, /* media is satured */
  CS_GROUNDWATER_MODEL_GENUCHTEN, /* Van Genuchten-Mualem laws for dimensionless
                                     moisture content and hydraulic conductivity
                                  */
  CS_GROUNDWATER_MODEL_TRACY,     /* Tracy model for unsaturated soils */
  CS_GROUNDWATER_MODEL_USER,      /* User-defined model */
  CS_GROUNDWATER_N_MODELS

} cs_groundwater_model_t;

typedef struct _groundwater_t  cs_groundwater_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \return a pointer to a new allocated cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_groundwater_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to groundwater flows
 *
 * \param[in, out]  gw     pointer to a cs_groundwater_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_groundwater_finalize(cs_groundwater_t   *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       keyname   name of key related to the member of adv to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_param(cs_groundwater_t    *gw,
                         const char          *keyname,
                         const char          *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_groundwater_t structure
 *
 * \param[in]  gw     pointer to a cs_groundwater_t struct. to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_summary(const cs_groundwater_t   *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]      connect          pointer to a cs_cdo_connect_t structure
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in]      model            keyword related to the model used
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_groundwater_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_init(const cs_cdo_connect_t  *connect,
                    int                      richards_eq_id,
                    const char              *model,
                    cs_property_t           *permeability,
                    cs_property_t           *soil_capacity,
                    cs_adv_field_t          *adv_field,
                    cs_groundwater_t        *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw               pointer to a cs_groundwater_t structure
 * \param[in]      tracer_eq_id     id related to the tracer equation
 * \param[in]      eqname           name of the equation
 * \param[in]      varname          name of the related variable
 * \param[in]      diff_property    pointer to a cs_property_t struct.
 * \param[in]      dispersivity     dispersivity for each axis (x, y, z]
 * \param[in]      bulk_density     value of the bulk density
 * \param[in]      distrib_coef     value of the distribution coefficient
 * \param[in]      reaction_rate    value of the first order rate of reaction
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_add_tracer(cs_groundwater_t   *gv,
                          int                 tracer_eq_id,
                          const char         *eq_name,
                          const char         *var_name,
                          cs_property_t      *diff_property,
                          cs_real_3_t         dispersivity,
                          double              bulk_density,
                          double              distrib_coef,
                          double              reaction_rate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the module dedicated to groundwater flows
 *
 * \param[in, out] equations    pointer to the array of cs_equation_t struct.
 * \param[in, out] gw           pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_automatic_settings(cs_equation_t      **equations,
                                  cs_groundwater_t    *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] eqs        array of pointers to cs_equation_t structures
 * \param[in, out] gw         pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_compute(const cs_mesh_t              *mesh,
                       const cs_time_step_t         *time_step,
                       double                        dt_cur,
                       const cs_cdo_connect_t       *connect,
                       const cs_cdo_quantities_t    *cdoq,
                       cs_equation_t                *eqs[],
                       cs_groundwater_t             *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined postprocessing for the groundwater module
 *
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 * \param[in]  gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_post(const cs_time_step_t      *time_step,
                    const cs_groundwater_t    *gw);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GROUNDWATER_H__ */
