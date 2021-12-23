#ifndef __CS_GWF_SOIL_H__
#define __CS_GWF_SOIL_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_mesh.h"
#include "cs_property.h"
#include "cs_volume_zone.h"

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
 * \brief  Generic function to update the physical properties related to a
 *         hydraulic model. The soil context depends on the type of soil.
 *
 * \param[in]      t_eval         time at which one performs the evaluation
 * \param[in]      mesh           pointer to a cs_mesh_t structure
 * \param[in]      connect        pointer to a cs_cdo_connect_t structure
 * \param[in]      quant          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone           pointer to a cs_zone_t
 * \param[in, out] soil_context   pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_update_t) (const cs_real_t              t_eval,
                        const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_zone_t             *zone,
                        void                        *soil_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set free the soil context of a soil structure
 *
 * \param[in, out] context      double pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_free_context_t)(void         **p_context);

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_gwf_soil_hydraulic_model_t
 *
 * \brief Predefined hydraulic model of soils used in the groundwater flow
 *        module
 *
 * \var CS_GWF_SOIL_GENUCHTEN
 * Van Genuchten-Mualem laws defining the evolution of the effective liquid
 * saturne (also called dimensionless moisture content) and the relative
 * permeability
 *
 * The (effective) liquid saturation (also called moisture content) follows the
 * identity:
 * S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *         = (1 + |alpha . h|^n)^(-m)
 *
 * The isotropic relative permeability is defined as:
 * k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 * where m = 1 -  1/n
 *
 * \var CS_GWF_SOIL_SATURATED
 * Hydraulic model of soild where the soil is considered as saturated. In this
 * model, there no evolution taken into account. The liquid saturation and the
 * permeability are considered as constant.
 *
 * \var CS_GWF_SOIL_USER
 * User-defined model of soil
 */

typedef enum {

  CS_GWF_SOIL_GENUCHTEN,
  CS_GWF_SOIL_SATURATED,
  CS_GWF_SOIL_USER,

  CS_GWF_SOIL_N_HYDRAULIC_MODELS

} cs_gwf_soil_hydraulic_model_t;

/*! \struct cs_gwf_soil_context_genuchten_t
 *
 * \brief Structure to handle the Van Genuchten-Mualen model of soil
 *
 *        See \ref CS_GWF_SOIL_GENUCHTEN). This structure stores the parameters
 *        defining the evolution laws for the liquid saturation and the
 *        relative permeability.
 */

typedef struct {

  /*!
   *@name Model parameters
   *@{
   *
   * \var saturated_moisture
   * Value of the liquid saturation (or moisture content)
   *
   * \var residual_moisture
   *      Also called residual liquid saturation
   *
   * \var saturated_permeability
   *      Permeability value (isotropic or anisotropic) when the soil is saturated
   *
   * \var n
   * Shape parameter. Should be 1.25 < n < 6
   *
   * \var m
   * Value depending of that of n
   * m = 1 - 1/n
   *
   * \var scale
   * Value of a scaling parameter in [m^-1]
   *
   * \var tortuosity
   * tortuosity param. for saturated hydraulic conductivity
   *
   *@}
   */

  double             residual_moisture;
  double             saturated_moisture;
  double             saturated_permeability[3][3];

  double             n;
  double             m;
  double             scale;
  double             tortuosity;

  /*!
   *@name  Array values
   *{
   */

  cs_real_t         *permeability_values;
  cs_real_t         *head_values;
  cs_real_t         *moisture_values;
  cs_real_t         *capacity_values;

  /*!
   *@}
   */

} cs_gwf_soil_context_genuchten_t;


/*! \struct cs_gwf_soil_context_saturated_t
 *
 * \brief Parameters defining a saturated soil in a given zone
 *
 * \var saturated_moisture
 * Value of the liquid saturation (or moisture content)
 *
 * \var saturated_permeability
 * Value of the permeability in the soil
 */

typedef struct {

  double       saturated_moisture;
  double       saturated_permeability[3][3];

} cs_gwf_soil_context_saturated_t;


/*! \struct cs_gwf_soil_t
 *
 * \brief Main structure to handle a soil in the groundawater flow module.
 *        Store a set of parameters and pointers describing a soil
 *
 */

typedef struct {

  /*!
 * \var id
 * id associated to a soil. Position in the array of soils.
 *
 * \var zone_id
 * id related to a volumic cs_zone_t structure (based on cells)
 *
 * \var bulk_density
 * Value of the mass density of the soil
 *
 * \var saturated_moisture
 * Value of the liquid saturation (or moisture content)
 *
 * \var model
 * Type of model describing the hydraulic behaviour of a soil (cf. \ref
 * \cs_gwf_soil_hydraulic_model_t for more details)
 *
 * \var context
 * Pointer to an context structure which depends on the model. This pointer has
 * to be casted on-the-fly.
 *
 * \var update_properties
 * Pointer to a function which manages the update of the properties describing
 * the porous media or used in associated equations (diffusion terms for
 * instance). These functions depend on the model of soil and the type of model
 * used in the groundwater flow module. May be set to NULL if there is no need
 * to update soil properties.
 *
 * \var free_context
 * Pointer to a function which free the context structure if needed. May be set
 * to NULL if there is nothing to free inside the context structure.
 */

  int                             id;
  int                             zone_id;

  double                          bulk_density;
  double                          saturated_moisture;

  cs_gwf_soil_hydraulic_model_t   model;

  void                           *context;

  /* Pointers to functions */

  cs_gwf_soil_update_t           *update_properties;
  cs_gwf_soil_free_context_t     *free_context;

} cs_gwf_soil_t;

/*============================================================================
 * User-defined function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of allocated soils
 *
 * \return the number of allocated soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its id
 *
 * \param[in]  id      id to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_id(int   id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its name
 *
 * \param[in]  name      name to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_name(const char    *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the saturated moisture for the given soil id
 *
 * \param[in]  soil_id     id of the requested soil
 *
 * \return the value of the saturated moisture
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_saturated_moisture(int   soil_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if all soils have been set as CS_GWF_SOIL_SATURATED
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_soil_all_saturated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check that at least one soil has been defined and the model of soil
 *         exists.
 *         Raise an error if a problem is encoutered.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_check(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new cs_gwf_soil_t structure and add it to the array of
 *         soils. An initialization by default of all members is performed.
 *
 * \param[in] zone                pointer to a volume zone structure
 * \param[in] model               type of modelling for the hydraulic behavior
 * \param[in] perm_type           type of permeability (iso/anisotropic)
 * \param[in] saturated_moisture  moisture content
 * \param[in] bulk_density        value of the mass density
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_create(const cs_zone_t                 *zone,
                   cs_gwf_soil_hydraulic_model_t    model,
                   cs_property_type_t               perm_type,
                   double                           saturated_moisture,
                   double                           bulk_density);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build an array storing the associated soil for each cell
 *
 * \param[in] n_cells      number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_build_cell2soil(cs_lnum_t    n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the associated soil for each cell
 */
/*----------------------------------------------------------------------------*/

const short int *
cs_gwf_get_cell2soil(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a saturated hydraulic model and attached to
 *         an isotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the saturated permeability
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_iso_saturated(cs_gwf_soil_t              *soil,
                              double                      k_s);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a saturated hydraulic model and attached to
 *         an anisotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the anisotropic saturated permeability
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_aniso_saturated(cs_gwf_soil_t              *soil,
                                double                      k_s[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a Van Genuchten-Mualen hydraulic model and
 *         attached to an isotropic saturated permeability
 *
 *         The (effective) liquid saturation (also called moisture content)
 *         follows the identity
 *         S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *                 = (1 + |alpha . h|^n)^(-m)
 *
 *         The isotropic relative permeability is defined as:
 *         k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 *         where m = 1 -  1/n
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the isotropic saturated permeability
 * \param[in]      theta_r    residual moisture
 * \param[in]      alpha      scale parameter (in m^-1)
 * \param[in]      n          shape parameter
 * \param[in]      L          turtuosity parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_iso_genuchten(cs_gwf_soil_t              *soil,
                              double                      k_s,
                              double                      theta_r,
                              double                      alpha,
                              double                      n,
                              double                      L);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a Van Genuchten-Mualen hydraulic model and
 *         attached to an anisotropic saturated permeability.
 *
 *         The (effective) liquid saturation (also called moisture content)
 *         follows the identity
 *         S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *                 = (1 + |alpha . h|^n)^(-m)
 *
 *         The isotropic relative permeability is defined as:
 *         k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 *         where m = 1 -  1/n
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the isotropic saturated permeability
 * \param[in]      theta_r    residual moisture/liquid saturation
 * \param[in]      alpha      scale parameter (in m^-1)
 * \param[in]      n          shape parameter
 * \param[in]      L          turtuosity parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_aniso_genuchten(cs_gwf_soil_t              *soil,
                                double                      k_s[3][3],
                                double                      theta_r,
                                double                      alpha,
                                double                      n,
                                double                      L);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a user-defined hydraulic model
 *
 * \param[in, out] soil               pointer to a cs_gwf_soil_t structure
 * \param[in]      context            pointer to a structure cast on-the-fly
 * \param[in]      update_func        function pointer to update propoerties
 * \param[in]      free_context_func  function pointer to free the context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_user(cs_gwf_soil_t                *soil,
                     void                         *context,
                     cs_gwf_soil_update_t         *update_func,
                     cs_gwf_soil_free_context_t   *free_context_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the properties of the groundwater flow module in the case where
 *         all soils are considered as saturated.
 *
 * \param[in, out]  permeability      pointer to a cs_property_t structure
 * \param[in, out]  moisture_content  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_saturated_set_properties(cs_property_t         *permeability,
                                     cs_property_t         *moisture_content);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the different arrays used in soil context for a GWF model set
 *         to unsaturated single-phase flows in a porous media.
 *
 * \param[in]  head              pointer to the current head values in cells
 * \param[in]  permeability      pointer to the current permeability values
 * \param[in]  moisture_content  pointer to the current moisture content values
 * \param[in]  capacity          pointer to the current soil capacity values
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_uspf_set_arrays(cs_real_t        head[],
                            cs_real_t        permeability[],
                            cs_real_t        moisture_content[],
                            cs_real_t        capacity[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the different arrays used in soil context for a GWF model set
 *         to miscible two-phase flows in a porous media.
 *
 * \param[in]  capillarity_p  current values of the capillarity pressure
 * \param[in]  l_saturation   current values of the liquid saturation
 * \param[in]  l_kr           current values of the relative liquid permeability
 * \param[in]  g_kr           current values of the relative gas permeability
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_mtpf_set_arrays(cs_real_t      capillarity_p[],
                            cs_real_t      l_saturation[],
                            cs_real_t      l_kr[],
                            cs_real_t      g_kr[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the soil properties
 *
 * \param[in]  time_eval         time at which one evaluates properties
 * \param[in]  mesh              pointer to the mesh structure
 * \param[in]  connect           pointer to the cdo connectivity
 * \param[in]  quant             pointer to the cdo quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_update(cs_real_t                     time_eval,
                   const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_SOIL_H__ */
