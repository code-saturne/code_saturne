/*============================================================================
 * Main functions dedicated to soil management in groundwater flows
 *============================================================================*/

/* VERS */

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_field.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_reco.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_soil.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_soil.c

  \brief Main functions dedicated to soil management in groundwater flows when
         using CDO schemes

*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_SOIL_DBG 0

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_soil[] =
  " Stop execution. The structure related to a soil is empty.\n"
  " Please check your settings.\n";

static int  _n_soils = 0;
static cs_gwf_soil_t  **_soils = NULL;

/* The following array enables to get the soil id related to each cell.
   The array size is equal to n_cells */
static short int *_cell2soil_ids = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function that set free the given saturated soil
 *
 * \param[in, out] input      pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_free_saturated_soil(void       *input)
{
  cs_gwf_soil_saturated_param_t  *soil_param =
    (cs_gwf_soil_saturated_param_t *)input;

  BFT_FREE(soil_param);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function that compute the new values of the properties related to
 *         a soil with a saturated.
 *         Case of an isotropic permeability and an unsteady Richards eq.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      head_values  array of values for head used in law
 * \param[in]      zone         pointer to a cs_zone_t
 * \param[in, out] input        pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_saturated_iso_soil(const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts,
                           const cs_real_t             *head_values,
                           const cs_zone_t             *zone,
                           void                        *input)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);
  CS_UNUSED(head_values);

  /* Retrieve field values associated to properties to update */
  cs_field_t  *f = NULL;
  f = cs_field_by_name("permeability");
  cs_real_t  *permeability_values = f->val;
  f = cs_field_by_name("moisture_content");
  cs_real_t  *moisture_values = f->val;

  const cs_gwf_soil_saturated_param_t  *law =
    (cs_gwf_soil_saturated_param_t *)input;

  const  double  iso_satval = law->saturated_permeability[0][0];

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN) default(none)     \
  shared(zone, law, permeability_values, moisture_values)                  \
  firstprivate(iso_satval)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    /* Set the permeability value to the saturated values */
    permeability_values[c_id] = iso_satval;

    /* Set the moisture content (Se = 1 in this case)*/
    moisture_values[c_id] = law->saturated_moisture;

  } /* Loop on selected cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function that compute the new values of the properties related to
 *         a soil with a saturated.
 *         Case of an anisotropic permeability and an unsteady Richards eq.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      head_values  array of values for head used in law
 * \param[in]      zone         pointer to a cs_zone_t
 * \param[in, out] input        pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_saturated_aniso_soil(const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_time_step_t        *ts,
                             const cs_real_t             *head_values,
                             const cs_zone_t             *zone,
                             void                        *input)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);
  CS_UNUSED(head_values);

  /* Retrieve field values associated to properties to update */
  cs_field_t  *f = NULL;
  f = cs_field_by_name("permeability");
  cs_real_t  *permeability_values = f->val;
  f = cs_field_by_name("moisture_content");
  cs_real_t  *moisture_values = f->val;

  const cs_gwf_soil_saturated_param_t  *law =
    (cs_gwf_soil_saturated_param_t *)input;

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN) default(none)     \
  shared(zone, law, permeability_values, moisture_values)
  for (cs_lnum_t id = 0; id < zone->n_elts; id++) {

    const cs_lnum_t  c_id = zone->elt_ids[id];

    /* Set the permeability value to the saturated values */
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        permeability_values[9*c_id+3*i+j] = law->saturated_permeability[i][j];

    /* Set the moisture content (Se = 1 in this case)*/
    moisture_values[c_id] = law->saturated_moisture;

  } /* Loop on selected cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function that set free the given Van Genuchten soil
 *
 * \param[in, out] input      pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_free_genuchten_soil(void          *input)
{
  cs_gwf_soil_genuchten_param_t  *soil_param =
    (cs_gwf_soil_genuchten_param_t *)input;

  BFT_FREE(soil_param);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function that compute the new values of the properties related to
 *         a soil with a Van Genuchten-Mualen.
 *         Case of an isotropic permeability and an unsteady Richards eq.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      head_values  array of values for head used in law
 * \param[in]      zone         pointer to a cs_zone_t
 * \param[in, out] input        pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_genuchten_iso_soil(const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts,
                           const cs_real_t             *head_values,
                           const cs_zone_t             *zone,
                           void                        *input)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);

  /* Sanity checks */
  assert(head_values != NULL);

  /* Retrieve field values associated to properties to update */
  cs_field_t  *f = NULL;
  f = cs_field_by_name("permeability");
  cs_real_t  *permeability_values = f->val;
  f = cs_field_by_name("moisture_content");
  cs_real_t  *moisture_values = f->val;
  f = cs_field_by_name_try("soil_capacity");
  cs_real_t  *capacity_values = NULL;
  if (f != NULL)
    capacity_values = f->val;

  const cs_gwf_soil_genuchten_param_t  *law =
    (cs_gwf_soil_genuchten_param_t *)input;

  /* Up to now, only isotropic values are considered */
  const  double  iso_satval = law->saturated_permeability[0][0];
  const  double  delta_moisture =
    law->saturated_moisture - law->residual_moisture;

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN) default(none)   \
  shared(head_values, zone, law, permeability_values, moisture_values,   \
         capacity_values)                                                \
  firstprivate(iso_satval, delta_moisture)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head_values[c_id];

    if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */

      const double  coef = pow(fabs(law->scale * h), law->n);
      const double  se = pow(1 + coef, -law->m);
      const double  se_pow_overm = pow(se, 1/law->m);
      const double  coef_base = 1 - pow(1 - se_pow_overm, law->m);

      /* Set the permeability value */
      permeability_values[c_id] = iso_satval* pow(se, law->tortuosity) *
        coef_base*coef_base;

      /* Set the moisture content */
      moisture_values[c_id] = se*delta_moisture + law->residual_moisture;

      /* Set the soil capacity */
      if (capacity_values != NULL) {
        const double  ccoef = -law->n * law->m * delta_moisture;
        const double  se_m1 = se/(1. + coef);

        capacity_values[c_id] = ccoef * coef/h * se_m1;
      }

    }
    else {

      /* Set the permeability value to the saturated values */
      permeability_values[c_id] = iso_satval;

      /* Set the moisture content (Se = 1 in this case)*/
      moisture_values[c_id] = delta_moisture + law->residual_moisture;

      /* Set the soil capacity */
      if (capacity_values != NULL)
        capacity_values[c_id] = 0.;

    }

  } /* Loop on selected cells */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
                cs_gwf_soil_hydraulic_model_t    model)
{
  cs_property_t  *permeability = cs_property_by_name("permeability");
  cs_gwf_soil_t  *soil = NULL;

  BFT_MALLOC(soil, 1, cs_gwf_soil_t);

  int  soil_id = _n_soils;

  soil->id = soil_id;

  /* Attached a volume zone to the current soil */
  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);
  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Zone %s related to the same soil is not defined.\n"
              " Stop adding a new soil.", z_name);

  soil->zone_id = zone->id;
  soil->model = model;
  soil->input = NULL;

  /* Set function pointers */
  switch (model) {

  case CS_GWF_SOIL_SATURATED:

    switch (permeability->type) {
    case CS_PROPERTY_ISO:
      soil->update_properties = _update_saturated_iso_soil;
      break;
    case CS_PROPERTY_ANISO:
      soil->update_properties = _update_saturated_aniso_soil;
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid type of property for the permeability.\n"
                " Please check your settings.");
    }
    soil->free_input = _free_saturated_soil;
    break;

  case CS_GWF_SOIL_GENUCHTEN:
    switch (permeability->type) {
    case CS_PROPERTY_ISO:
      soil->update_properties = _update_genuchten_iso_soil;
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid type of property for the permeability.\n"
                " Please check your settings.");
    }
    soil->free_input = _free_genuchten_soil;
    break;

  case CS_GWF_SOIL_USER:
    /* cs_user_cdo_gwf_init_soil(z_name); */
    break;

  default:
    break; /* Nothing to do */

  } /* Switch on soil modeling */

  /* Store the new soils in the soil array */
  _n_soils++;
  BFT_REALLOC(_soils, _n_soils, cs_gwf_soil_t *);
  _soils[soil_id] = soil;

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_free_all(void)
{
  if (_n_soils < 1)
    return;
  assert(_soils != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->free_input != NULL)
      soil->free_input(soil->input);

    BFT_FREE(soil);

  } /* Loop on soils */

  BFT_FREE(_soils);
  BFT_FREE(_cell2soil_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of allocated soils
 *
 * \return the number of allocated soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(void)
{
  return _n_soils;
}

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
cs_gwf_soil_by_id(int   id)
{
  if (id > -1 && id < _n_soils)
    return _soils[id];
  else
    return NULL;
}

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
cs_gwf_soil_by_name(const char    *name)
{
  if (name == NULL)
    return NULL;

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *s = _soils[i];
    const cs_zone_t  *zone = cs_volume_zone_by_id(s->zone_id);

    if (strcmp(zone->name, name) == 0)
      return s;
  }

  /* Not found among the list */
  return NULL;
}

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
cs_gwf_soil_get_bulk_density(const cs_gwf_soil_t  *soil)
{
  assert(soil != NULL);

  cs_real_t  bulk_density = 1.0;

  switch (soil->model) {

  case CS_GWF_SOIL_GENUCHTEN:
    {
      const cs_gwf_soil_genuchten_param_t *param =
        (const cs_gwf_soil_genuchten_param_t *)soil->input;
      bulk_density = param->bulk_density;
    }
    break;

  case CS_GWF_SOIL_SATURATED:
    {
      const cs_gwf_soil_saturated_param_t *param =
        (const cs_gwf_soil_saturated_param_t *)soil->input;
      bulk_density = param->bulk_density;
    }
    break;

  case CS_GWF_SOIL_USER:
    cs_user_gwf_get_soil_density(soil, &bulk_density);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid model of soil.");

  }

  return bulk_density;
}

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
                              double                      rho)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_SATURATED)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not saturated\n", __func__);

  cs_gwf_soil_saturated_param_t  *soil_param = NULL;

  BFT_MALLOC(soil_param, 1, cs_gwf_soil_saturated_param_t);

  soil_param->bulk_density = rho;
  soil_param->saturated_moisture = theta_s;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      soil_param->saturated_permeability[i][j] = 0;
  for (int i = 0; i < 3; i++)
    soil_param->saturated_permeability[i][i] = k_s;

  soil->input = soil_param;
}

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
                                double                      rho)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_SATURATED)
    bft_error(__FILE__, __LINE__, 0,
              " %s : soil model is not saturated\n", __func__);

  cs_gwf_soil_saturated_param_t  *soil_param = NULL;

  BFT_MALLOC(soil_param, 1, cs_gwf_soil_saturated_param_t);

  soil_param->bulk_density = rho;
  soil_param->saturated_moisture = theta_s;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      soil_param->saturated_permeability[i][j] =  k_s[3*i+j];

  soil->input = soil_param;
}

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
                              double                      rho)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_GENUCHTEN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not Van Genuchten\n", __func__);

  cs_gwf_soil_genuchten_param_t  *soil_param = NULL;

  BFT_MALLOC(soil_param, 1, cs_gwf_soil_genuchten_param_t);

  soil_param->bulk_density = rho;
  soil_param->saturated_moisture = theta_s;
  soil_param->residual_moisture = theta_r;

  /* Additional advanced settings */
  soil_param->n = 1.56;
  soil_param->m = 1 - 1/soil_param->n;
  soil_param->scale = 0.036;
  soil_param->tortuosity = 0.5;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      soil_param->saturated_permeability[i][j] = 0;
  for (int i = 0; i < 3; i++)
    soil_param->saturated_permeability[i][i] = k_s;

  soil->input = soil_param;
}

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
                                double                      rho)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_GENUCHTEN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not Van Genuchten\n", __func__);

  cs_gwf_soil_genuchten_param_t  *soil_param = NULL;

  BFT_MALLOC(soil_param, 1, cs_gwf_soil_genuchten_param_t);

  soil_param->bulk_density = rho;
  soil_param->saturated_moisture = theta_s;
  soil_param->residual_moisture = theta_r;

  /* Additional advanced settings */
  soil_param->n = 1.56;
  soil_param->m = 1 - 1.56;
  soil_param->scale = 0.036;
  soil_param->tortuosity = 0.5;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      soil_param->saturated_permeability[i][j] = k_s[3*i+j];

  soil->input = soil_param;
}

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
                     cs_gwf_soil_finalize_t     *free_func)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_USER)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not user-defined.\n", __func__);

  /* Set function pointers */
  soil->input = input;
  soil->update_properties = update_func;
  soil->free_input = free_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the properties of the groundwater flow module all soils are
 *         considered as saturated.
 *
 * \param[in, out]  permeability      pointer to a cs_property_t structure
 * \param[in, out]  moisture_content  pointer to a cs_property_t structure
 * \param[in, out]  moisture_field    pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_all_saturated(cs_property_t         *permeability,
                              cs_property_t         *moisture_content,
                              cs_field_t            *moisture_field)
{
  CS_UNUSED(moisture_field);
  assert(permeability != NULL && moisture_content != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->model != CS_GWF_SOIL_SATURATED)
      bft_error(__FILE__, __LINE__, 0,
                " Invalid way of setting soil parameter.\n"
                " All soils are not considered as saturated.");

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_gwf_soil_saturated_param_t  *param =
      (cs_gwf_soil_saturated_param_t  *)soil->input;

    /* Set the permeability */
    switch (permeability->type) {
    case CS_PROPERTY_ISO:
      cs_property_def_iso_by_value(permeability,
                                   z->name,
                                   param->saturated_permeability[0][0]);
      break;

    case CS_PROPERTY_ORTHO:
      {
        cs_real_3_t  val = {param->saturated_permeability[0][0],
                            param->saturated_permeability[1][1],
                            param->saturated_permeability[2][2]};

        cs_property_def_ortho_by_value(permeability,
                                       z->name,
                                       val);
      }
      break;

    case CS_PROPERTY_ANISO:
      cs_property_def_aniso_by_value(permeability,
                                     z->name,
                      (double (*)[3])param->saturated_permeability);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of property.\n", __func__);
    }

    /* Set the moisture content */
    cs_property_def_iso_by_value(moisture_content,
                                 z->name,
                                 param->saturated_moisture);

  } /* Loop on soils */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build an array storing the associated soil for each cell
 *
 * \param[in] n_cells      number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_build_cell2soil(cs_lnum_t    n_cells)
{
  BFT_MALLOC(_cell2soil_ids, n_cells, short int);

  if (_n_soils == 1) {

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t j = 0; j < n_cells; j++)
      _cell2soil_ids[j] = 0;

  }
  else {

    assert(_n_soils > 1);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t j = 0; j < n_cells; j++)
      _cell2soil_ids[j] = -1; /* unset by default */

    for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

      const cs_gwf_soil_t  *soil = _soils[soil_id];
      const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      assert(z != NULL);

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < z->n_elts; j++)
        _cell2soil_ids[z->elt_ids[j]] = soil_id;

    } /* Loop on soils */

    /* Chcek if every cells is associated to a soil */
    for (cs_lnum_t j = 0; j < n_cells; j++)
      if (_cell2soil_ids[j] == -1)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: At least cell%d has no related soil.\n",
                  __func__, j);

  } /* n_soils > 1 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the associated soil for each cell
 */
/*----------------------------------------------------------------------------*/

const short int *
cs_gwf_get_cell2soil(void)
{
  return _cell2soil_ids;
}

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
 * \param[in, out]  moisture_content  pointer to a cs_property_t structure
 * \param[in]       moisture_field    pointer to a cs_field_t structure
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
                         cs_field_t        *capacity_field)
{
  /* All soils are considered all at once when property is defined by field */

  cs_property_def_by_field(permeability, permea_field);
  cs_property_def_by_field(moisture_content, moisture_field);

  if (soil_capacity != NULL)
    cs_property_def_by_field(soil_capacity, capacity_field);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void)
{
  const char  *meta = "  <GWF/Hydraulic Model>";
  cs_log_printf(CS_LOG_SETUP, "  <GWF/Soils>  n_soils %d", _n_soils);

  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_log_printf(CS_LOG_SETUP, "\n  <GWF/Soil %d> %s\n", soil->id, z->name);
    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
      {
        const cs_gwf_soil_genuchten_param_t  *si =
          (cs_gwf_soil_genuchten_param_t *)soil->input;

        cs_log_printf(CS_LOG_SETUP, "%s VanGenuchten-Mualen\n", meta);
        cs_log_printf(CS_LOG_SETUP, "    <Soil parameters>");
        cs_log_printf(CS_LOG_SETUP,
                      " residual_moisture %5.3e", si->residual_moisture);
        cs_log_printf(CS_LOG_SETUP,
                      " saturated_moisture %5.3e\n", si->saturated_moisture);
        cs_log_printf(CS_LOG_SETUP, "    <Soil parameters> n= %f, scale= %f,"
                      "tortuosity= %f\n", si->n, si->scale, si->tortuosity);
        cs_log_printf(CS_LOG_SETUP, "    <Soil saturated permeability>");
        cs_log_printf(CS_LOG_SETUP,
                      " [%-4.2e %4.2e %4.2e; %-4.2e %4.2e %4.2e;"
                      " %-4.2e %4.2e %4.2e]",
                      si->saturated_permeability[0][0],
                      si->saturated_permeability[0][1],
                      si->saturated_permeability[0][2],
                      si->saturated_permeability[1][0],
                      si->saturated_permeability[1][1],
                      si->saturated_permeability[1][2],
                      si->saturated_permeability[2][0],
                      si->saturated_permeability[2][1],
                      si->saturated_permeability[2][2]);
      }
      break;

    case CS_GWF_SOIL_SATURATED:
      {
        const cs_gwf_soil_saturated_param_t  *si =
          (cs_gwf_soil_saturated_param_t *)soil->input;

        cs_log_printf(CS_LOG_SETUP, "%s saturated\n", meta);
        cs_log_printf(CS_LOG_SETUP, "    <Soil parameters>");
        cs_log_printf(CS_LOG_SETUP,
                      " saturated_moisture %5.3e\n", si->saturated_moisture);
        cs_log_printf(CS_LOG_SETUP, "    <Soil saturated permeability>");
        cs_log_printf(CS_LOG_SETUP,
                      " [%-4.2e %4.2e %4.2e; %-4.2e %4.2e %4.2e;"
                      " %-4.2e %4.2e %4.2e]",
                      si->saturated_permeability[0][0],
                      si->saturated_permeability[0][1],
                      si->saturated_permeability[0][2],
                      si->saturated_permeability[1][0],
                      si->saturated_permeability[1][1],
                      si->saturated_permeability[1][2],
                      si->saturated_permeability[2][0],
                      si->saturated_permeability[2][1],
                      si->saturated_permeability[2][2]);
      }
      break;

    case CS_GWF_SOIL_USER:
      cs_log_printf(CS_LOG_SETUP, "%s user-defined\n", meta);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid model for groundwater module.\n"
                " Please check your settings.");

    } /* Switch model */

  } /* Loop on soils */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
