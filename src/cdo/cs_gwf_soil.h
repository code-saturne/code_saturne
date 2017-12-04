#ifndef __CS_GWF_SOIL_H__
#define __CS_GWF_SOIL_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_time_step.h"
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
 *         hydraulic model. At least, moisture content and permeability are
 *         updated. If the simulation is time-depedent, then the soil capacity
 *         can be updated also.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      head_in_law  array of values for head used in law
 * \param[in]      zone         pointer to a cs_volume_zone_t
 * \param[in, out] input        pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_update_t) (const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *ts,
                        const cs_real_t             *head_values,
                        const cs_volume_zone_t      *zone,
                        void                        *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set free the input of a soil structure
 *
 * \param[in, out] input      pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_finalize_t)(void         *input);


/*============================================================================
 * Type definitions
 *============================================================================*/

/* Predefined hydraulic model of soil
 *=================================== */

/* Type of predefined modelling for the groundwater flows */
typedef enum {

  CS_GWF_SOIL_GENUCHTEN, /* Van Genuchten-Mualem laws for dimensionless
                            moisture content and hydraulic conductivity */
  CS_GWF_SOIL_SATURATED, /* media is satured */
  CS_GWF_SOIL_USER,      /* User-defined model */
  CS_GWF_SOIL_N_HYDRAULIC_MODELS

} cs_gwf_soil_hydraulic_model_t;

/* Structures used to handle soils with a Van Genuchten-Mualen modelling */
/* --------------------------------------------------------------------- */

/* Parameters defining a van Genuchten-Mualen hydraulic model */
typedef struct {

  double         bulk_density;
  double         residual_moisture;
  double         saturated_moisture;
  cs_real_33_t   saturated_permeability;

  /* Advanced parameters */
  double  n;          // 1.25 < n < 6
  double  m;          // m = 1 - 1/n
  double  scale;      // scale parameter [m^-1]
  double  tortuosity; // tortuosity param. for saturated hydraulic conductivity

} cs_gwf_soil_genuchten_param_t;

/* Input structure used to update the physical properties */
typedef struct {

  cs_real_t         *permeability_values;
  cs_real_t         *head_values;
  cs_real_t         *moisture_values;
  cs_real_t         *capacity_values;

} cs_gwf_genuchten_update_input_t;

/* Structures used to handle soils with a saturated modelling */
/* ---------------------------------------------------------- */

/* Parameters defining a saturated hydraulic model */
typedef struct {

  double         bulk_density;
  double         saturated_moisture;
  cs_real_33_t   saturated_permeability;

} cs_gwf_soil_saturated_param_t;

typedef struct {

  cs_real_t         *permeability_values;
  cs_real_t         *moisture_values;

} cs_gwf_saturated_update_input_t;

/* Set of parameters describing a soil */
typedef struct {

  int    id;       /* soil id */
  int    zone_id;  /* id related to a cs_volume_zone_t structure (based on
                      cells) */

  /* Physical modelling adopted for this soil */
  cs_gwf_soil_hydraulic_model_t   model;

  /* Pointer to an input structure according to the hydraulic model */
  void                           *input;

  /* Pointers to functions */
  cs_gwf_soil_update_t           *update_properties;
  cs_gwf_soil_finalize_t         *free_input;

} cs_gwf_soil_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a new cs_gwf_soil_t structure. A first initialization
 *         of all members by default is performed.
 *
 * \param[in]   z_name      name of the volume zone corresponding to the soil
 * \param[in]   model       type of modeling for the hydraulic behavior
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_add(const char                      *z_name,
                cs_gwf_soil_hydraulic_model_t    model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_free_all(void);

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
 * \brief  Retrieve the bulk density associated to the given soil structure
 *
 * \param[in] soil       pointer to a cs_gwf_soil_t structure
 *
 * \return
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_bulk_density(const cs_gwf_soil_t  *soil);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a saturated hydraulic model and attached to
 *         an isotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_iso_saturated_soil(cs_gwf_soil_t              *soil,
                              double                      k_s,
                              double                      theta_s,
                              double                      rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a saturated hydraulic model and attached to
 *         an isotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the anisotropic saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_aniso_saturated_soil(cs_gwf_soil_t              *soil,
                                double                     *k_s,
                                double                      theta_s,
                                double                      rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a Van Genuchten hydraulic model and attached
 *         to an anisotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the isotropic saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      theta_r    residual moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_iso_genuchten_soil(cs_gwf_soil_t              *soil,
                              double                      k_s,
                              double                      theta_s,
                              double                      theta_r,
                              double                      rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a Van Genuchten hydraulic model and attached
 *         to an anisotropic permeability
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      k_s        value of the anisotropic saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      theta_r    residual moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_aniso_genuchten_soil(cs_gwf_soil_t              *soil,
                                double                     *k_s,
                                double                      theta_s,
                                double                      theta_r,
                                double                      rho);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a user-defined hydraulic model and attached to
 *         an anisotropic permeability
 *
 * \param[in, out] soil           pointer to a cs_gwf_soil_t structure
 * \param[in]      input          pointer to a structure cast on-the-fly
 * \param[in]      update_func    pointer to the function used for updating
 * \param[in]      free_func      pointer to the function used for finalizing
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_user_soil(cs_gwf_soil_t              *soil,
                     void                       *input,
                     cs_gwf_soil_update_t       *update_func,
                     cs_gwf_soil_finalize_t     *free_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the properties of the groundwater flow module all soils are
 *         considered as saturated.
 *
 * \param[in, out]  permeability      pointer to a cs_property_t structure
 * \param[in, out]  moisture content  pointer to a cs_property_t structure
 * \param[in, out]  moisture field    pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_all_saturated(cs_property_t         *permeability,
                              cs_property_t         *moisture_content,
                              cs_field_t            *moisture_field);

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
 * \brief  Set the properties of the groundwater flow module thanks to
 *         cs_field_t structure. The consequence is that the related
 *         cs_property_t structure relies on only one definition (i.e. for the
 *         whole mesh). Fields are updated by using the update function pointer
 *         associated to each soil.
 *
 * \param[in, out]  permeability      pointer to a cs_property_t structure
 * \param[in]       permea_field      pointer to a cs_field_t structure
 * \param[in, out]  moisture content  pointer to a cs_property_t structure
 * \param[in]       moisture field    pointer to a cs_field_t structure
 * \param[in, out]  soil_capacity     pointer to a cs_property_t structure
 * \param[in]       capacity_field    pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_by_field(cs_property_t     *permeability,
                         cs_field_t        *permea_field,
                         cs_property_t     *moisture_content,
                         cs_field_t        *moisture_field,
                         cs_property_t     *soil_capacity,
                         cs_field_t        *capacity_field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_SOIL_H__ */
