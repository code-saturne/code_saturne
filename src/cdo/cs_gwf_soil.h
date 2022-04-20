#ifndef __CS_GWF_SOIL_H__
#define __CS_GWF_SOIL_H__

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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_gwf_param.h"
#include "cs_gwf_priv.h"
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

typedef struct _gwf_soil_t cs_gwf_soil_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the properties related to a hydraulic
 *         model given the soil model. The soil parameters depend on the type
 *         of soil model.
 *
 * \param[in]      t_eval         time at which one performs the evaluation
 * \param[in]      mesh           pointer to a cs_mesh_t structure
 * \param[in]      connect        pointer to a cs_cdo_connect_t structure
 * \param[in]      quant          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone           pointer to a cs_zone_t
 * \param[in, out] soil           pointer to the soil structure to update
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_update_t) (const cs_real_t              t_eval,
                        const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_zone_t             *zone,
                        cs_gwf_soil_t               *soil);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set free the parameter structure associated to
 *         a soil
 *
 * \param[in, out] p_param     double pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_free_param_t)(void         **p_param);

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_gwf_soil_param_genuchten_t
 *
 * \brief Structure to handle the Van Genuchten-Mualen model of soil
 *
 *        See \ref CS_GWF_SOIL_GENUCHTEN. This structure stores the parameters
 *        defining the evolution laws for the liquid saturation and the
 *        relative permeability.
 */

typedef struct {

  /*!
   * \var residual_moisture
   *      Also called residual liquid saturation
   *
   * \var n
   *      Shape parameter. Should be 1.25 < n < 6
   *
   * \var m
   *      Value depending of that of n
   *      m = 1 - 1/n
   *
   * \var scale
   *      Value of a scaling parameter in [m^-1]
   *
   * \var tortuosity
   *      Tortuosity parameter for saturated hydraulic conductivity
   */

  double             residual_moisture;

  double             n;
  double             m;
  double             scale;
  double             tortuosity;

} cs_gwf_soil_param_genuchten_t;


/*! \struct _gwf_soil_t
 *
 * \brief Main structure to handle a soil in the groundwater flow module.
 *
 *        Store a set of parameters and pointers describing a soil and its
 *        related hydraulic model (shared with the main structure \ref
 *        cs_gwf_t)
 */

struct _gwf_soil_t {

  /*!
   * @name Metadata
   * @{

   * \var id
   *      id associated to a soil. Position in the array of soils.
   *
   * \var zone_id
   *      id related to a volumic cs_zone_t structure (based on cells)
   *
   * \var hydraulic_model
   *      Type of model use in the groundwater flow module to describe the
   *      hydraulic (see \ref cs_gwf_model_type_t for more details)
   *
   * \var hydraulic_context
   * Structure cast on-the-fly. This structure contains parameters, arrays,
   * properties and fields describing the hydraulic state. It depends on the
   * type of hydraulic model which is considered.
   *
   * @}
   * @name Soil features (whatever is the soil model)
   * @{
   *
   * \var bulk_density
   * Value of the mass density of the soil
   *
   * \var porosity
   *      Max. portion of volume in a soil where the liquid (or the gas) can
   *      be.  In a single-phase saturated model this corresponds to the
   *      saturated moisture or the max. liquid saturation.
   *
   * \var abs_permeability_dim
   *      =1 if isotropic or =3 if orthotropic or =9 if anisotropic
   *
   * \var abs_permeability
   *      Value of the intrisic permeability in the soil when all the porous
   *      media is filled with water (= saturated permeability)
   *
   * @}
   * @name Soil modelling
   * @{
   *
   * \var model
   *      Type of model describing the hydraulic behaviour of a soil (cf. \ref
   *      \cs_gwf_soil_model_t for more details)
   *
   * \var model_param
   *      Pointer to a structure cast on-the-fly (it depends on the type of
   *      soil model). This structure contains the set of parameters describing
   *      a soil. Can be set to NULL in the case of a saturated single-phase
   *      model.
   *
   * @}
   * @name Function pointers
   * @{
   *
   * \var update_properties
   *      Pointer to a function which manages the update of the properties
   *      describing the porous media or used in associated equations
   *      (diffusion terms for instance). These functions depend on the model
   *      of soil and the type of model used in the groundwater flow
   *      module. May be set to NULL if there is no need to update soil
   *      properties.
   *
   * \var free_model_param
   *      Pointer to a function which free the param structure if needed. May
   *      be set to NULL if there is nothing to free inside the param
   *      structure.
   *
   * @}
   */

  int                             id;
  int                             zone_id;

  cs_gwf_model_type_t             hydraulic_model;
  void                           *hydraulic_context;


  double                          bulk_density;
  double                          porosity;
  int                             abs_permeability_dim;
  double                          abs_permeability[3][3];

  cs_gwf_soil_model_t             model;
  void                           *model_param;

  /* Pointers to functions */

  cs_gwf_soil_update_t           *update_properties;
  cs_gwf_soil_free_param_t       *free_model_param;

};

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
 * \brief  Retrieve the max dim (aniso=9; iso=1) for the absolute permeability
 *         associated to each soil
 *
 * \return the associated max. dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_soil_get_permeability_max_dim(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if all soils have been set as CS_GWF_SOIL_SATURATED
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_soil_all_are_saturated(void);

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
 * \param[in] hydraulic_model     main hydraulic model for the module
 * \param[in] model               type of model for the soil behavior
 * \param[in] perm_type           type of permeability (iso/anisotropic)
 * \param[in] k_abs               absolute (intrisic) permeability
 * \param[in] porosity            porosity or max. moisture content
 * \param[in] bulk_density        value of the mass density
 * \param[in] hydraulic_context   pointer to the context structure
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_create(const cs_zone_t                 *zone,
                   cs_gwf_model_type_t              hydraulic_model,
                   cs_gwf_soil_model_t              model,
                   cs_property_type_t               perm_type,
                   double                           k_abs[3][3],
                   double                           porosity,
                   double                           bulk_density,
                   void                            *hydraulic_context);

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
 * \brief  Set a soil defined by a Van Genuchten-Mualen model
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
 * \param[in]      theta_r    residual moisture
 * \param[in]      alpha      scale parameter (in m^-1)
 * \param[in]      n          shape parameter
 * \param[in]      L          turtuosity parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_genuchten_param(cs_gwf_soil_t              *soil,
                                double                      theta_r,
                                double                      alpha,
                                double                      n,
                                double                      L);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a user-defined hydraulic model
 *
 * \param[in, out] soil              pointer to a cs_gwf_soil_t structure
 * \param[in]      context           pointer to a structure cast on-the-fly
 * \param[in]      update_func       function pointer to update propoerties
 * \param[in]      free_param_func   function pointer to free the context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_user(cs_gwf_soil_t               *soil,
                     void                        *context,
                     cs_gwf_soil_update_t        *update_func,
                     cs_gwf_soil_free_param_t    *free_param_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the definition of the soil porosity and absolute porosity (which
 *         are properties always defined). This relies on the definition of
 *         each soil.
 *
 * \param[in, out]  abs_permeability    pointer to a cs_property_t structure
 * \param[in, out]  soil_porosity       pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_shared_properties(cs_property_t      *abs_permeability,
                                  cs_property_t      *soil_porosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the definition of the soil porosity and absolute porosity (which
 *         are properties always defined). This relies on the definition of
 *         each soil.
 *
 * \param[in, out]  moisture_content  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_saturated_set_property(cs_property_t   *moisture_content);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the definition of some property(ies) in specific situations for
 *         the two-phase flow models
 *         This relies on the definition of each soil.
 *
 * \param[in, out]  mc  pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_tpf_set_property(cs_gwf_two_phase_t     *mc);

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
/*!
 * \brief  Update arrays associated to the definition of terms involved in the
 *         miscible two-phase flow model.
 *         Case of an isotropic absolute permeability.
 *
 * \param[in, out] mc            pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_iso_update_mtpf_terms(cs_gwf_two_phase_t    *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update arrays associated to the definition of terms involved in the
 *         immiscible two-phase flow model.
 *         Case of an isotropic absolute permeability.
 *
 * \param[in, out] mc            pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_iso_update_itpf_terms(cs_gwf_two_phase_t     *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update arrays associated to the definition of terms involved in the
 *         immiscible two-phase flow model.
 *         Case of an isotropic absolute permeability and an incremental solve
 *
 * \param[in]      ts          pointer to a cs_time_step_t structure
 * \param[in, out] mc          pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_iso_update_itpf_terms_incr(const cs_time_step_t    *ts,
                                       cs_gwf_two_phase_t      *mc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update arrays associated to the definition of terms involved in the
 *         immiscible two-phase flow model.
 *         Case of an isotropic absolute permeability with an incremental solve
 *         and a liquid saturation defined on a submesh.
 *
 * \param[in]      ts      pointer to a cs_time_step_t structure
 * \param[in]      connect pointer to a cs_cdo_connect_t structure
 * \param[in, out] mc      pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_iso_update_itpf_terms_incr_submesh(const cs_time_step_t    *ts,
                                               const cs_cdo_connect_t  *connect,
                                               cs_gwf_two_phase_t      *mc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_SOIL_H__ */
