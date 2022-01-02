/*============================================================================
 * Main functions dedicated to soil management in groundwater flows
 *============================================================================*/

/* VERS */

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
#include <bft_printf.h>

#include "cs_field.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
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
 * \brief  Function that compute the new values of the properties related to
 *         a soil with a Van Genuchten-Mualen.
 *         Case of an isotropic permeability and an unsteady Richards eq.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      head_values   array of values for head used in law
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_soil_genuchten_iso(const cs_real_t              t_eval,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_zone_t             *zone,
                           void                        *soil_context)
{
  CS_UNUSED(t_eval);
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_gwf_soil_context_genuchten_t  *sc = soil_context;

  /* Only isotropic values are considered in this case */

  const double  iso_satval = sc->saturated_permeability[0][0];
  const double  delta_moisture = sc->saturated_moisture - sc->residual_moisture;
  const cs_real_t  *head = sc->head_values;

  /* Retrieve field values associated to properties to update */

  cs_real_t  *permeability = sc->permeability_values;
  cs_real_t  *moisture = sc->moisture_values;
  cs_real_t  *capacity = sc->capacity_values;

  assert(capacity != NULL && permeability != NULL && moisture != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)                \
  shared(head, zone, sc, permeability, moisture, capacity)              \
  firstprivate(iso_satval, delta_moisture)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head[c_id];

    if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */

      const double  coef = pow(fabs(sc->scale * h), sc->n);
      const double  se = pow(1 + coef, -sc->m);
      const double  se_pow_overm = pow(se, 1/sc->m);
      const double  coef_base = 1 - pow(1 - se_pow_overm, sc->m);

      /* Set the permeability value */

      permeability[c_id] =
        iso_satval* pow(se, sc->tortuosity) * coef_base*coef_base;

      /* Set the moisture content */

      moisture[c_id] = se*delta_moisture + sc->residual_moisture;

      /* Set the soil capacity */

      const double  ccoef = -sc->n * sc->m * delta_moisture;
      const double  se_m1 = se/(1. + coef);

      capacity[c_id] = ccoef * coef/h * se_m1;

    }
    else {

      /* Set the permeability value to the saturated values */

      permeability[c_id] = iso_satval;

      /* Set the moisture content (Se = 1 in this case)*/

      moisture[c_id] = delta_moisture + sc->residual_moisture;

      /* Set the soil capacity */

      capacity[c_id] = 0.;

    }

  } /* Loop on selected cells */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
 * \brief  Get the saturated moisture for the given soil id
 *
 * \param[in]  soil_id     id of the requested soil
 *
 * \return the value of the saturated moisture
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_saturated_moisture(int   soil_id)
{
  cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

  if (soil == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty soil.\n", __func__);

  return soil->saturated_moisture;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if all soils have been set as CS_GWF_SOIL_SATURATED
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_soil_all_saturated(void)
{
  for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = _soils[soil_id];
    if (soil->model != CS_GWF_SOIL_SATURATED)
      return false;

  }

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check that at least one soil has been defined and the model of soil
 *         exists.
 *         Raise an error if a problem is encoutered.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_check(void)
{
  if (_n_soils < 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is activated but no soil is defined.",
              __func__);
  if (_soils == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The soil structure is not allocated whereas %d soils"
              " have been added.\n", __func__, _n_soils);

  for (int i = 0; i < _n_soils; i++) {

    if (_soils[i]->model == CS_GWF_SOIL_N_HYDRAULIC_MODELS) {
      const cs_zone_t  *z = cs_volume_zone_by_id(_soils[i]->zone_id);
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid model of soil attached to zone %s\n",
                __func__, z->name);
    }

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new cs_gwf_soil_t structure and add it to the array of
 *         soils. An initialization by default of all members is performed.
 *
 * \param[in]   zone          pointer to a volume zone structure
 * \param[in]   model         type of modelling for the hydraulic behavior
 * \param[in]   perm_type     type of permeability (iso/anisotropic)
 * \param[in]   sat_moisture  value of the saturated moisture content
 * \param[in]   bulk_density  value of the mass density
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_create(const cs_zone_t                 *zone,
                   cs_gwf_soil_hydraulic_model_t    model,
                   cs_property_type_t               perm_type,
                   double                           sat_moisture,
                   double                           bulk_density)
{
  cs_gwf_soil_t  *soil = NULL;

  BFT_MALLOC(soil, 1, cs_gwf_soil_t);

  soil->id = _n_soils;
  soil->bulk_density = bulk_density;
  soil->saturated_moisture = sat_moisture;
  soil->model = model;

  /* Attached a volume zone to the current soil */

  assert(zone != NULL);
  soil->zone_id = zone->id;
  soil->context = NULL;

  soil->update_properties = NULL;
  soil->free_context = NULL;

  switch (model) {

  case CS_GWF_SOIL_SATURATED:
    {
      cs_gwf_soil_context_saturated_t  *sc = NULL;

      BFT_MALLOC(sc, 1, cs_gwf_soil_context_saturated_t);

      /* Default initialization */

      sc->saturated_moisture = sat_moisture;

      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          sc->saturated_permeability[ki][kj] = 0.0;

      sc->saturated_permeability[0][0] = 1.0;
      sc->saturated_permeability[1][1] = 1.0;
      sc->saturated_permeability[2][2] = 1.0;

      soil->context = sc;
    }
    break;

  case CS_GWF_SOIL_GENUCHTEN:
    {
      cs_gwf_soil_context_genuchten_t  *sc = NULL;

      BFT_MALLOC(sc, 1, cs_gwf_soil_context_genuchten_t);

      sc->residual_moisture = 0.;
      sc->saturated_moisture = sat_moisture;

      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          sc->saturated_permeability[ki][kj] = 0.0;

      sc->saturated_permeability[0][0] = 1.0;
      sc->saturated_permeability[1][1] = 1.0;
      sc->saturated_permeability[2][2] = 1.0;

      sc->n = 1.25;
      sc->m = 1 - 1./sc->n;
      sc->scale = 1.;
      sc->tortuosity = 1.;

      /* Pointer to property values will be set after */

      sc->permeability_values = NULL;
      sc->head_values = NULL;
      sc->moisture_values = NULL;
      sc->capacity_values = NULL;

      soil->context = sc;

      if (perm_type & CS_PROPERTY_ISO)
        soil->update_properties = _update_soil_genuchten_iso;
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of property for the permeability.\n"
                  " Please check your settings.", __func__);
    }
    break;

  case CS_GWF_SOIL_USER:
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of soil model\n", __func__);
    break; /* Nothing to do */

  } /* Switch on soil modeling */

  /* Store the new soils in the soil array */

  _n_soils++;
  BFT_REALLOC(_soils, _n_soils, cs_gwf_soil_t *);
  _soils[soil->id] = soil;

  return soil;
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

  if (_n_soils == 1)
    memset(_cell2soil_ids, 0, sizeof(short int)*n_cells);

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
                  " %s: At least cell %ld has no related soil.\n",
                  __func__, (long)j);

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

    if (soil->free_context != NULL)
      soil->free_context(&(soil->context));

    if (soil->context != NULL) {

      switch (soil->model) {

      case CS_GWF_SOIL_SATURATED:
        {
          cs_gwf_soil_context_saturated_t  *sc = soil->context;

          BFT_FREE(sc);
          sc = NULL;
        }
        break;

      case CS_GWF_SOIL_GENUCHTEN:
        {
          cs_gwf_soil_context_genuchten_t  *sc = soil->context;

          BFT_FREE(sc);
          sc = NULL;
        }
        break;

      default:
        cs_base_warn(__FILE__, __LINE__);
        bft_printf("%s: The context structure of a soil may not be freed.\n",
                   __func__);
        break;

      } /* Switch on predefined soil context */

    }
    BFT_FREE(soil);

  } /* Loop on soils */

  BFT_FREE(_soils);
  BFT_FREE(_cell2soil_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP, "  * GWF | Number of soils: %d\n", _n_soils);

  char  meta[64];
  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_log_printf(CS_LOG_SETUP, "\n        Soil.%d | Zone: %s\n",
                  soil->id, z->name);
    cs_log_printf(CS_LOG_SETUP, "\n        Soil.%d | Bulk.density: %6.3e\n",
                  soil->id, soil->bulk_density);

    sprintf(meta, "        Soil.%d |", soil->id);

    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
      {
        const cs_gwf_soil_context_genuchten_t  *sc = soil->context;

        cs_log_printf(CS_LOG_SETUP, "%s Model: VanGenuchten-Mualen\n", meta);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", meta);
        cs_log_printf(CS_LOG_SETUP,
                      " residual_moisture %5.3e", sc->residual_moisture);
        cs_log_printf(CS_LOG_SETUP,
                      " saturated_moisture %5.3e\n", sc->saturated_moisture);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", meta);
        cs_log_printf(CS_LOG_SETUP, " n= %f, scale= %f, tortuosity= %f\n",
                      sc->n, sc->scale, sc->tortuosity);
        cs_log_printf(CS_LOG_SETUP, "%s Saturated permeability\n", meta);
        cs_log_printf(CS_LOG_SETUP, "%s [%-4.2e %4.2e %4.2e;\n", meta,
                      sc->saturated_permeability[0][0],
                      sc->saturated_permeability[0][1],
                      sc->saturated_permeability[0][2]);
        cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e;\n", meta,
                      sc->saturated_permeability[1][0],
                      sc->saturated_permeability[1][1],
                      sc->saturated_permeability[1][2]);
        cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e]\n", meta,
                      sc->saturated_permeability[2][0],
                      sc->saturated_permeability[2][1],
                      sc->saturated_permeability[2][2]);
      }
      break;

    case CS_GWF_SOIL_SATURATED:
      {
        const cs_gwf_soil_context_saturated_t  *sc = soil->context;

        cs_log_printf(CS_LOG_SETUP, "%s Model: Saturated\n", meta);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters", meta);
        cs_log_printf(CS_LOG_SETUP,
                      " saturated_moisture %5.3e\n", sc->saturated_moisture);
        cs_log_printf(CS_LOG_SETUP, "%s Saturated permeability\n", meta);
        cs_log_printf(CS_LOG_SETUP, "%s [%-4.2e %4.2e %4.2e;\n", meta,
                      sc->saturated_permeability[0][0],
                      sc->saturated_permeability[0][1],
                      sc->saturated_permeability[0][2]);
        cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e;\n", meta,
                      sc->saturated_permeability[1][0],
                      sc->saturated_permeability[1][1],
                      sc->saturated_permeability[1][2]);
        cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e]\n", meta,
                      sc->saturated_permeability[2][0],
                      sc->saturated_permeability[2][1],
                      sc->saturated_permeability[2][2]);
      }
      break;

    case CS_GWF_SOIL_USER:
      cs_log_printf(CS_LOG_SETUP, "%s Model: User-defined\n", meta);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid model for groundwater module.\n"
                " Please check your settings.");

    } /* Switch model */

  } /* Loop on soils */

  cs_log_printf(CS_LOG_SETUP, "\n");
}

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
                              double                      k_s)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_context_saturated_t  *sc = soil->context;

  if (soil->model != CS_GWF_SOIL_SATURATED)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not saturated\n", __func__);
  if (sc == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);

  /* Default initialization is the identity matrix */

  for (int i = 0; i < 3; i++)
    sc->saturated_permeability[i][i] = k_s;
}

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
                                double                      k_s[3][3])
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_context_saturated_t  *sc = soil->context;

  if (soil->model != CS_GWF_SOIL_SATURATED)
    bft_error(__FILE__, __LINE__, 0,
              "%s : soil model is not saturated\n", __func__);
  if (sc == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      sc->saturated_permeability[i][j] =  k_s[i][j];
}

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
                              double                      L)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_context_genuchten_t  *sc = soil->context;

  if (soil->model != CS_GWF_SOIL_GENUCHTEN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not Van Genuchten\n", __func__);
  if (sc == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);
  if (n <= FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n = %6.4e (the shape parameter).\n"
              "This value should be > 0.\n", __func__, n);

  sc->residual_moisture = theta_r;

  /* Default initialization is the identity matrix */

  for (int i = 0; i < 3; i++)
    sc->saturated_permeability[i][i] = k_s;

  /* Additional advanced settings */

  sc->n = n;
  sc->m = 1 - 1/sc->n;
  sc->scale = alpha;
  sc->tortuosity = L;
}

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
                                double                      L)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_context_genuchten_t  *sc = soil->context;

  if (soil->model != CS_GWF_SOIL_GENUCHTEN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not Van Genuchten\n", __func__);
  if (sc == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);
  if (n <= FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n = %6.4e (the shape parameter).\n"
              "This value should be > 0.\n", __func__, n);

  sc->residual_moisture = theta_r;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      sc->saturated_permeability[i][j] = k_s[i][j];

  /* Additional advanced settings */

  sc->n = n;
  sc->m = 1 - 1/sc->n;
  sc->scale = alpha;
  sc->tortuosity = L;
}

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
                     cs_gwf_soil_free_context_t   *free_context_func)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_USER)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not user-defined.\n", __func__);

  /* Set function pointers */
  soil->context = context;
  soil->update_properties = update_func;
  soil->free_context = free_context_func;
}

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
                                     cs_property_t         *moisture_content)
{
  assert(permeability != NULL && moisture_content != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->model != CS_GWF_SOIL_SATURATED)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid way of setting soil parameter.\n"
                " All soils are not considered as saturated.", __func__);

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_gwf_soil_context_saturated_t  *sc = soil->context;

    /* Set the permeability */

    if (permeability->type & CS_PROPERTY_ISO)
      cs_property_def_iso_by_value(permeability,
                                   z->name,
                                   sc->saturated_permeability[0][0]);

    else if (permeability->type & CS_PROPERTY_ORTHO) {

      cs_real_3_t  val = {sc->saturated_permeability[0][0],
                          sc->saturated_permeability[1][1],
                          sc->saturated_permeability[2][2]};

      cs_property_def_ortho_by_value(permeability, z->name, val);

    }
    else if (permeability->type & CS_PROPERTY_ANISO) {

      cs_property_def_aniso_by_value(permeability,
                                     z->name,
                      (double (*)[3])sc->saturated_permeability);

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of property.\n", __func__);

    /* Set the moisture content */

    assert(fabs(sc->saturated_moisture - soil->saturated_moisture) < FLT_MIN);

    cs_property_def_iso_by_value(moisture_content,
                                 z->name,
                                 sc->saturated_moisture);

  } /* Loop on soils */

}

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
                            cs_real_t        capacity[])
{
  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
      {
        cs_gwf_soil_context_genuchten_t  *sc = soil->context;

        sc->permeability_values = permeability;
        sc->head_values = head;
        sc->moisture_values = moisture_content;
        sc->capacity_values = capacity;
      }
      break;

    default:
      break; /* Do nothing. For user-defined soils, one has to do similar
                things in cs_user_gwf.c */

    }

  } /* Loop on soils */
}

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
                   const cs_cdo_quantities_t    *quant)
{
  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];
    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);

    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
    case CS_GWF_SOIL_USER:
      soil->update_properties(time_eval,
                              mesh, connect, quant,
                              zone,
                              soil->context);
      break;

    default:
      break; /* Do nothing (for instance in the case of a saturated soil which
                is constant (steady and uniform) */

    }

  } /* Loop on soils */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
