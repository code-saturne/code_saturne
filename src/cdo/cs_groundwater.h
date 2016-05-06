#ifndef __CS_GROUNDWATER_H__
#define __CS_GROUNDWATER_H__

/*============================================================================
 * Compute the wall distance using the CDO framework
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/* List of available keys for setting the groundwater flow module */
typedef enum {

  CS_GWKEY_GRAVITATION,       // Take into acount gravitation effect
  CS_GWKEY_OUTPUT_MOISTURE,   // Activate post-processing for the moisture
  CS_GWKEY_N_KEYS

} cs_groundwater_key_t;

typedef enum {

  CS_SOILKEY_SAT_MOISTURE,  // Set the saturated moisture content
  CS_SOILKEY_RES_MOISTURE,  // Set the residual moisture content

  /* Keys specific to the Tracy model */
  CS_SOILKEY_TRACY_SAT_H,   // Head related to the saturated moisture content
  CS_SOILKEY_TRACY_RES_H,   // Head related to the residual moisture content

  CS_SOILKEY_N_KEYS

} cs_groundwater_soilkey_t;


/* Type of predefined modelling for the groundwater flows */
typedef enum {

  CS_GROUNDWATER_MODEL_COMPOSITE, /* Mixed of predefined groundwater model */
  CS_GROUNDWATER_MODEL_GENUCHTEN, /* Van Genuchten-Mualem laws for dimensionless
                                     moisture content and hydraulic conductivity
                                  */
  CS_GROUNDWATER_MODEL_SATURATED, /* media is satured */
  CS_GROUNDWATER_MODEL_TRACY,     /* Tracy model for unsaturated soils */
  CS_GROUNDWATER_MODEL_USER,      /* User-defined model */
  CS_GROUNDWATER_N_MODELS

} cs_groundwater_model_t;

/* Parameters defining the van Genuchten-Mualen law */
typedef struct {

  double  n;          // 1.25 < n < 6
  double  m;          // m = 1 - 1/n
  double  scale;      // scale parameter [m^-1]
  double  tortuosity; // tortuosity param. for saturated hydraulic conductivity

} cs_gw_genuchten_t;

typedef struct {

  double   h_r;
  double   h_s;

} cs_gw_tracy_t;

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
 * \brief  Get the number of requested soils
 *
 * \param[in]  gw        pointer to a cs_groundwater_t structure
 *
 * \return the number of requested soils
 */
/*----------------------------------------------------------------------------*/

int
cs_groundwater_get_n_soils(const cs_groundwater_t    *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       key       key related to the member of gw to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_param(cs_groundwater_t      *gw,
                         cs_groundwater_key_t   key,
                         const char            *keyval);

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
 * \param[in]      n_soils          number of soils to consider
 * \param[in]      n_tracers        number of tracers to consider
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_groundwater_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_initialize(const cs_cdo_connect_t  *connect,
                          int                      richards_eq_id,
                          int                      n_soils,
                          int                      n_tracer_eqs,
                          cs_property_t           *permeability,
                          cs_property_t           *soil_capacity,
                          cs_adv_field_t          *adv_field,
                          cs_groundwater_t        *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new type of soil to consider in the groundwater module
 *
 * \param[in, out] gw         pointer to a cs_groundwater_t structure
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      model_kw   keyword related to the model used
 * \param[in]      ks         value(s) of the saturated permeability
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_add_soil_by_value(cs_groundwater_t   *gw,
                                 const char         *ml_name,
                                 const char         *model_kw,
                                 const char         *pty_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       ml_name   name of the mesh location associated to this soil
 * \param[in]       key       key related to a member of the soil to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_soil_param(cs_groundwater_t          *gw,
                              const char                *ml_name,
                              cs_groundwater_soilkey_t   key,
                              const char                *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion and reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_groundwater_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      eqname          name of the equation
 * \param[in]      varname         name of the related variable
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_add_tracer(cs_groundwater_t    *gw,
                          int                  tracer_eq_id,
                          const char          *eqname,
                          const char          *varname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_groundwater_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      ml_name         name of the related mesh location
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_tracer_param(cs_groundwater_t    *gw,
                                int                  tracer_eq_id,
                                const char          *ml_name,
                                double               wmd,
                                double               alpha_l,
                                double               alpha_t,
                                double               bulk_density,
                                double               distrib_coef,
                                double               reaction_rate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the Richards equation
 *
 * \param[in, out] gw        pointer to a cs_groundwater_t structure
 * \param[in, out] richards  pointer to the related cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_richards_setup(cs_groundwater_t    *gw,
                              cs_equation_t       *richards);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a reaction term for a given tracer
 *
 * \param[in] gw         pointer to a cs_groundwater_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_groundwater_tracer_needs_reaction(const cs_groundwater_t    *gw,
                                     int                        eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a diffusion term for a given tracer
 *
 * \param[in] gw         pointer to a cs_groundwater_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_groundwater_tracer_needs_diffusion(const cs_groundwater_t    *gw,
                                      int                        eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for a tracer equation
 *
 * \param[in]      tracer_eq_id  id of the equation related to this tracer
 * \param[in, out] eq            pointer to the related cs_equation_t structure
 * \param[in, out] gw            pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_tracer_setup(int                  tracer_eq_id,
                            cs_equation_t       *eq,
                            cs_groundwater_t    *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      do_logcvg  output information on convergence or not
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
                       bool                          do_logcvg,
                       cs_equation_t                *eqs[],
                       cs_groundwater_t             *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_groundwater_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_list    list of cells (1 to n)
 * \param[in]      i_face_list  list of interior faces (1 to n)
 * \param[in]      b_face_list  list of boundary faces (1 to n)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_extra_post(void                      *input,
                          int                        mesh_id,
                          int                        cat_id,
                          int                        ent_flag[5],
                          cs_lnum_t                  n_cells,
                          cs_lnum_t                  n_i_faces,
                          cs_lnum_t                  n_b_faces,
                          const cs_lnum_t            cell_list[],
                          const cs_lnum_t            i_face_list[],
                          const cs_lnum_t            b_face_list[],
                          const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GROUNDWATER_H__ */
