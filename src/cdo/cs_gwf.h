#ifndef __CS_GWF_H__
#define __CS_GWF_H__

/*============================================================================
 * Set of main functions to handle the groundwater flow module with CDO
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

typedef enum {

  CS_SOILKEY_SAT_MOISTURE,  // Set the saturated moisture content
  CS_SOILKEY_RES_MOISTURE,  // Set the residual moisture content
  CS_SOILKEY_BULK_DENSITY,  // Set the bulk density of a soil

  /* Keys specific to the Tracy model */
  CS_SOILKEY_TRACY_SAT_H,   // Head related to the saturated moisture content
  CS_SOILKEY_TRACY_RES_H,   // Head related to the residual moisture content

  CS_SOILKEY_N_KEYS

} cs_gwf_soilkey_t;


/* Type of predefined modelling for the groundwater flows */
typedef enum {

  CS_GWF_HYDRAULIC_COMPOSITE, /* Mixed of predefined groundwater model */
  CS_GWF_HYDRAULIC_GENUCHTEN, /* Van Genuchten-Mualem laws for dimensionless
                                 moisture content and hydraulic conductivity */
  CS_GWF_HYDRAULIC_SATURATED, /* media is satured */
  CS_GWF_HYDRAULIC_TRACY,     /* Tracy model for unsaturated soils */
  CS_GWF_HYDRAULIC_USER,      /* User-defined model */
  CS_GWF_N_HYDRAULIC_MODELS

} cs_gwf_hydraulic_model_t;

/* Parameters defining the van Genuchten-Mualen law */
typedef struct {

  double  n;          // 1.25 < n < 6
  double  m;          // m = 1 - 1/n
  double  scale;      // scale parameter [m^-1]
  double  tortuosity; // tortuosity param. for saturated hydraulic conductivity

} cs_gwf_genuchten_t;

typedef struct {

  double   h_r;
  double   h_s;

} cs_gwf_tracy_t;

typedef struct _gwf_t  cs_gwf_t;

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

cs_gwf_t *
cs_gwf_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to groundwater flows
 *
 * \param[in, out]  gw     pointer to a cs_gwf_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_finalize(cs_gwf_t   *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of requested soils
 *
 * \param[in]  gw        pointer to a cs_gwf_t structure
 *
 * \return the number of requested soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(const cs_gwf_t    *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the gravity and set the gravitaty vector
 *
 * \param[in, out]  gw        pointer to a cs_gwf_t structure
 * \param[in]       gvec      values of the gravity vector
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_gravity_vector(cs_gwf_t              *gw,
                          const cs_real_3_t      gvec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Advanced setting: indicate where the darcian flux is stored
 *         cs_cdo_primal_cell is the default setting
 *         cs_cdo_dual_face_byc is a valid choice for vertex-based schemes
 *
 * \param[in, out]  gw              pointer to a cs_gwf_t structure
 * \param[in]       location_flag   where the flux is defined
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_darcian_flux_location(cs_gwf_t      *gw,
                                 cs_flag_t      location_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_gwf_t structure
 *
 * \param[in]  gw     pointer to a cs_gwf_t struct. to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_summary(const cs_gwf_t   *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in]      n_soils          number of soils to consider
 * \param[in]      n_tracers        number of tracers to consider
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_gwf_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_initialize(int                      richards_eq_id,
                  int                      n_soils,
                  int                      n_tracer_eqs,
                  cs_property_t           *permeability,
                  cs_property_t           *soil_capacity,
                  cs_adv_field_t          *adv_field,
                  cs_gwf_t                *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an isotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      k_s        value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_iso_soil_by_value(cs_gwf_t                   *gw,
                             cs_gwf_hydraulic_model_t    model,
                             const char                 *ml_name,
                             double                      k_s,
                             double                      theta_s,
                             double                      rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an orthotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      ks         value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_ortho_soil_by_value(cs_gwf_t                  *gw,
                               cs_gwf_hydraulic_model_t   model,
                               const char                *ml_name,
                               cs_real_t                 *ks,
                               double                     theta_s,
                               double                     rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an orthotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      k_s        value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_aniso_soil_by_value(cs_gwf_t                  *gw,
                               cs_gwf_hydraulic_model_t   model,
                               const char                *ml_name,
                               cs_real_t                 *ks,
                               double                     theta_s,
                               double                     rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_gwf_t structure
 *
 * \param[in, out]  gw        pointer to a cs_gwf_t structure
 * \param[in]       ml_name   name of the mesh location associated to this soil
 * \param[in]       key       key related to a member of the soil to set
 * \param[in]       val       value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_soil_param(cs_gwf_t           *gw,
                      const char         *ml_name,
                      cs_gwf_soilkey_t    key,
                      const cs_real_t     val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion and reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_gwf_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      eqname          name of the equation
 * \param[in]      varname         name of the related variable
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_add_tracer(cs_gwf_t       *gw,
                  int             tracer_eq_id,
                  const char     *eqname,
                  const char     *varname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_gwf_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      ml_name         name of the related mesh location
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_tracer_param(cs_gwf_t      *gw,
                        int            tracer_eq_id,
                        const char    *ml_name,
                        double         wmd,
                        double         alpha_l,
                        double         alpha_t,
                        double         distrib_coef,
                        double         reaction_rate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the Richards equation
 *
 * \param[in, out] gw        pointer to a cs_gwf_t structure
 * \param[in, out] richards  pointer to the related cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_richards_setup(cs_gwf_t        *gw,
                      cs_equation_t   *richards);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a reaction term for a given tracer
 *
 * \param[in] gw         pointer to a cs_gwf_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_tracer_needs_reaction(const cs_gwf_t    *gw,
                             int                eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a diffusion term for a given tracer
 *
 * \param[in] gw         pointer to a cs_gwf_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_tracer_needs_diffusion(const cs_gwf_t    *gw,
                              int                eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for a tracer equation
 *
 * \param[in]      tracer_eq_id  id of the equation related to this tracer
 * \param[in, out] eq            pointer to the related cs_equation_t structure
 * \param[in, out] gw            pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_setup(int               tracer_eq_id,
                    cs_equation_t    *eq,
                    cs_gwf_t         *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last initialization step of the groundwater flow module
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      n_equations  number of equations in the list
 * \param[in, out] equations    pointer to a list of cs_equation_t structures
 * \param[in, out] gw           pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_final_initialization(const cs_cdo_connect_t    *connect,
                            int                        n_equations,
                            cs_equation_t            **equations,
                            cs_gwf_t                  *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] eqs        array of pointers to cs_equation_t structures
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t              *mesh,
               const cs_time_step_t         *time_step,
               double                        dt_cur,
               const cs_cdo_connect_t       *connect,
               const cs_cdo_quantities_t    *cdoq,
               cs_equation_t                *eqs[],
               cs_gwf_t                     *gw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (cs_post_time_mesh_dep_output_t)
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
 * \param[in]      cell_list    list of cells (1 to n)
 * \param[in]      i_face_list  list of interior faces (1 to n)
 * \param[in]      b_face_list  list of boundary faces (1 to n)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post(void                      *input,
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

#endif /* __CS_GWF_H__ */
