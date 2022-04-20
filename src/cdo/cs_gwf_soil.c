/*============================================================================
 * Main functions dedicated to soil management in groundwater flows
 *============================================================================*/

/* VERS */

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
#include "cs_gwf_priv.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_constants.h"
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

#define CS_GWF_SOIL_DBG  0

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
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_soil_genuchten_iso(const cs_real_t              t_eval,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_zone_t             *zone,
                           cs_gwf_soil_t               *soil)
{
  CS_UNUSED(t_eval);
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  if (soil == NULL)
    return;

  assert(soil->hydraulic_model ==  CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE);

  /* Retrieve the soil parameters */

  cs_gwf_soil_param_genuchten_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_unsaturated_single_phase_t  *hc = soil->hydraulic_context;

  /* Only isotropic values are considered in this case */

  const double  iso_satval = soil->abs_permeability[0][0];
  const double  delta_m = soil->porosity - sp->residual_moisture;
  const cs_real_t  *head = hc->head_in_law;

  /* Retrieve field values associated to properties to update */

  cs_real_t  *permeability = hc->permeability_field->val;
  cs_real_t  *moisture = hc->moisture_field->val;
  cs_real_t  *capacity = hc->capacity_field->val;

  assert(capacity != NULL && permeability != NULL && moisture != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)                \
  shared(head, zone, sp, permeability, moisture, capacity)              \
  firstprivate(iso_satval, delta_m)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head[c_id];

    if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */

      const double  coef = pow(fabs(sp->scale * h), sp->n);
      const double  se = pow(1 + coef, -sp->m);
      const double  se_pow_overm = pow(se, 1/sp->m);
      const double  coef_base = 1 - pow(1 - se_pow_overm, sp->m);

      /* Set the permeability value : abs_perm * rel_perm */

      permeability[c_id] =
        iso_satval * pow(se, sp->tortuosity) * coef_base*coef_base;

      /* Set the moisture content (or liquid saturation) */

      moisture[c_id] = se*delta_m + sp->residual_moisture;

      /* Set the soil capacity = \frac{\partial S_l}{partial h} */

      const double  ccoef = -sp->n * sp->m * delta_m;
      const double  se_m1 = se/(1. + coef);

      capacity[c_id] = ccoef * coef/h * se_m1;

    }
    else {

      /* Set the permeability value to the saturated values */

      permeability[c_id] = iso_satval;

      /* Set the moisture content (Sle = 1 in this case)*/

      moisture[c_id] = delta_m + sp->residual_moisture;

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

  return soil->porosity;  /* = saturated moisture or max. liquid satiration */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the max dim (aniso=9; iso=1) for the absolute permeability
 *         associated to each soil
 *
 * \return the associated max. dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_soil_get_permeability_max_dim(void)
{
  int dim = 0;

  if (_n_soils < 1)
    return dim;

  if (_soils == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The soil structure is not allocated whereas %d soils"
              " have been added.\n", __func__, _n_soils);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    dim = CS_MAX(dim, soil->abs_permeability_dim);

  }

  return dim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if all soils have been set as CS_GWF_SOIL_SATURATED
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_soil_all_are_saturated(void)
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
                   void                            *hydraulic_context)
{
  cs_gwf_soil_t  *soil = NULL;

  BFT_MALLOC(soil, 1, cs_gwf_soil_t);

  soil->id = _n_soils;

  /* Attached a volume zone to the current soil */

  assert(zone != NULL);
  soil->zone_id = zone->id;

  /* Members related to the hydraulic model */

  soil->hydraulic_model = hydraulic_model;
  soil->hydraulic_context = hydraulic_context;

  /* Members related to the soil parameters/model */

  soil->model = model;
  soil->model_param = NULL;

  soil->bulk_density = bulk_density;
  soil->porosity = porosity;

  for (int ki = 0; ki < 3; ki++)
    for (int kj = 0; kj < 3; kj++)
      soil->abs_permeability[ki][kj] = k_abs[ki][kj];

  switch (perm_type) {

  case CS_PROPERTY_ISO:
    soil->abs_permeability_dim = 1;
    break;

  case CS_PROPERTY_ANISO:
    soil->abs_permeability_dim = 9;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of absolute permeability.\n", __func__);

  } /* Switch on the type of absolute permeability */

  /* Initialize function pointers */

  soil->update_properties = NULL;
  soil->free_model_param = NULL;

  /* Initialization which are specific to a soil model */

  switch (model) {

  case CS_GWF_SOIL_SATURATED:
    if (hydraulic_model != CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid type of soil with the general hydraulic model.\n"
                " In a saturated single-phase model, all soils have to be"
                " of type CS_GWF_SOIL_SATURATED.\n", __func__);
    break;

  case CS_GWF_SOIL_GENUCHTEN:
    {
      cs_gwf_soil_param_genuchten_t  *sp = NULL;

      BFT_MALLOC(sp, 1, cs_gwf_soil_param_genuchten_t);

      sp->residual_moisture = 0.;

      sp->n = 1.25;
      sp->m = 1 - 1./sp->n;
      sp->scale = 1.;
      sp->tortuosity = 1.;

      soil->model_param = sp;

      if (perm_type & CS_PROPERTY_ISO)
        if (hydraulic_model == CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE)
          soil->update_properties = _update_soil_genuchten_iso;
        else
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid type of hydraulic model.\n"
                    " Please check your settings.", __func__);
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of property for the permeability.\n"
                  " Please check your settings.", __func__);
    }
    break;

  case CS_GWF_SOIL_USER:
    break; /* All has to be done by the user */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of soil model\n", __func__);
    break; /* Nothing to do */

  } /* Switch on the soil model */

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

    if (soil->free_model_param != NULL)
      soil->free_model_param(&(soil->model_param));

    if (soil->model_param != NULL) {

      switch (soil->model) {

      case CS_GWF_SOIL_GENUCHTEN:
        {
          cs_gwf_soil_param_genuchten_t  *sp = soil->model_param;

          BFT_FREE(sp);
          sp = NULL;
        }
        break;

      default:
        cs_base_warn(__FILE__, __LINE__);
        bft_printf("%s: The context structure of a soil may not be freed.\n",
                   __func__);
        break;

      } /* Switch on predefined soil context */

    }

    /* The hydraulic context is shared and thus is freed during the free of the
       cs_gwf_t structure */

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

  char  id[64];
  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    sprintf(id, "        Soil.%d |", soil->id);

    cs_log_printf(CS_LOG_SETUP, "\n%s Zone: %s\n", id, z->name);
    cs_log_printf(CS_LOG_SETUP, "%s Bulk.density: %.1e\n",
                  id, soil->bulk_density);
    cs_log_printf(CS_LOG_SETUP, "%s Max.Porosity: %.3e (=saturated_moisture)\n",
                  id, soil->porosity);
    cs_log_printf(CS_LOG_SETUP, "%s Absolute permeability\n", id);
    cs_log_printf(CS_LOG_SETUP, "%s [%-4.2e %4.2e %4.2e;\n", id,
                  soil->abs_permeability[0][0],
                  soil->abs_permeability[0][1],
                  soil->abs_permeability[0][2]);
    cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e;\n", id,
                  soil->abs_permeability[1][0],
                  soil->abs_permeability[1][1],
                  soil->abs_permeability[1][2]);
    cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e]\n", id,
                  soil->abs_permeability[2][0],
                  soil->abs_permeability[2][1],
                  soil->abs_permeability[2][2]);

    /* Display the model parameters */

    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
      {
        const cs_gwf_soil_param_genuchten_t  *sp = soil->model_param;

        cs_log_printf(CS_LOG_SETUP, "%s Model: *VanGenuchten-Mualen*\n", id);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP,
                      " residual_moisture %5.3e\n", sp->residual_moisture);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP, " n= %f, scale= %f, tortuosity= %f\n",
                      sp->n, sp->scale, sp->tortuosity);
      }
      break;

    case CS_GWF_SOIL_SATURATED:
        cs_log_printf(CS_LOG_SETUP, "%s Model: *Saturated*\n", id);
      break;

    case CS_GWF_SOIL_USER:
      cs_log_printf(CS_LOG_SETUP, "%s Model: *User-defined*\n", id);
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
                                double                      L)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_param_genuchten_t  *sp = soil->model_param;

  if (soil->model != CS_GWF_SOIL_GENUCHTEN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not Van Genuchten\n", __func__);
  if (sp == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);
  if (n <= FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n = %6.4e (the shape parameter).\n"
              "This value should be > 0.\n", __func__, n);

  sp->residual_moisture = theta_r;

  /* Additional advanced settings */

  sp->n = n;
  sp->m = 1 - 1/sp->n;
  sp->scale = alpha;
  sp->tortuosity = L;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a soil defined by a user-defined model
 *
 * \param[in, out] soil              pointer to a cs_gwf_soil_t structure
 * \param[in]      param             pointer to a structure cast on-the-fly
 * \param[in]      update_func       function pointer to update propoerties
 * \param[in]      free_param_func   function pointer to free the param struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_user(cs_gwf_soil_t                *soil,
                     void                         *param,
                     cs_gwf_soil_update_t         *update_func,
                     cs_gwf_soil_free_param_t     *free_param_func)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_USER)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not user-defined.\n", __func__);

  /* Set pointers */

  soil->model_param = param;
  soil->update_properties = update_func;
  soil->free_model_param = free_param_func;
}

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
                                  cs_property_t      *soil_porosity)
{
  assert(abs_permeability != NULL && soil_porosity != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    /* Define the absolute permeability */

    if (abs_permeability->type & CS_PROPERTY_ISO) {

      assert(soil->abs_permeability_dim == 1);
      cs_property_def_iso_by_value(abs_permeability,
                                   z->name,
                                   soil->abs_permeability[0][0]);

    }
    else if (abs_permeability->type & CS_PROPERTY_ANISO) {

      cs_property_def_aniso_by_value(abs_permeability,
                                     z->name,
                                     soil->abs_permeability);

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of property.\n", __func__);

    /* Set the soil porosity */

    cs_property_def_iso_by_value(soil_porosity,
                                 z->name,
                                 soil->porosity);

  } /* Loop on soils */
}

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
cs_gwf_soil_saturated_set_property(cs_property_t   *moisture_content)
{
  assert(moisture_content != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->model != CS_GWF_SOIL_SATURATED)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid way of setting soil parameter.\n"
                " All soils are not considered as saturated.", __func__);

    /* Set the moisture content. In this case, one set the moisture content to
       the soil porosity since one considers that the soil is fully
       saturated */

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_property_def_iso_by_value(moisture_content,
                                 z->name,
                                 soil->porosity);

  } /* Loop on soils */
}

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
cs_gwf_soil_tpf_set_property(cs_gwf_two_phase_t     *mc)
{
  if (mc == NULL)
    return;
  assert(mc->time_hg_pty != NULL);

  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];
    assert(soil != NULL);

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    /* Set the value in each soil of the unsteady term for the Hydrogen
       conservation equation */

    cs_property_def_iso_by_value(mc->time_hg_pty,
                                 z->name,
                                 soil->porosity * mh_ov_rt);

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

    cs_gwf_soil_t  *soil = _soils[i];
    assert(soil != NULL);

    switch (soil->model) {

    case CS_GWF_SOIL_GENUCHTEN:
    case CS_GWF_SOIL_USER:
      {
        assert(soil->update_properties != NULL);

        const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);

        soil->update_properties(time_eval,
                                mesh, connect, quant,
                                zone,
                                soil);
      }
      break;

    default:
      break; /* Do nothing (for instance in the case of a saturated soil which
                is constant (steady and uniform) */

    } /* Switch on the soil model */

  } /* Loop on soils */

}

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
cs_gwf_soil_iso_update_mtpf_terms(cs_gwf_two_phase_t     *mc)
{
  if (mc == NULL)
    return;

  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const cs_real_t  *g_cell_pr = mc->g_cell_pressure;

  /* In the immiscible case, mc->l_diffusivity_h should be set to 0 */

  const double  h_diff_const = (mc->l_diffusivity_h > 0) ?
    hmh * mc->l_mass_density * mc->l_diffusivity_h / mc->w_molar_mass : 0.;

  const cs_real_t  *l_sat = mc->l_saturation->val;
  const cs_real_t  *l_cap = mc->l_capacity;

  for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = _soils[soil_id];
    assert(soil != NULL);
    assert(soil->hydraulic_model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != NULL);

    const double  w_time_coef = soil->porosity * mc->l_mass_density;
    const double  h_time_coefa = soil->porosity * mh_ov_rt;
    const double  h_time_coefb = soil->porosity*hmh - h_time_coefa;
    const double  wl_diff_coef = soil->abs_permeability[0][0]/mc->l_viscosity;
    const double  hg_diff_coef = soil->abs_permeability[0][0]/mc->g_viscosity;
    const double  h_diff_coef = soil->porosity * h_diff_const;

    /* Main loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

      const double  l_diff_coef = wl_diff_coef * mc->l_rel_permeability[c_id];

      /* Water conservation equation. Updates arrays linked to properties which
         define computed terms */

      mc->time_wg_array[c_id] = w_time_coef * l_cap[c_id];

      mc->time_wl_array[c_id] = -mc->time_wg_array[c_id];

      mc->diff_wl_array[c_id] = mc->l_mass_density * l_diff_coef;

      /* Hydrogen conservation equation. Updates arrays linked to properties
         which define computed terms */

      mc->time_hg_array[c_id] =
        h_time_coefa + h_time_coefb*(l_sat[c_id] + l_cap[c_id]*g_cell_pr[c_id]);

      mc->diff_hg_array[c_id] = mh_ov_rt * g_cell_pr[c_id] /* g_rho */
        * mc->g_rel_permeability[c_id] * hg_diff_coef;
      if (h_diff_coef > 0) /* If not = immiscible case */
        mc->diff_hg_array[c_id] += h_diff_coef * l_sat[c_id];

      mc->time_hl_array[c_id] = -h_time_coefb * g_cell_pr[c_id] * l_cap[c_id];

      mc->diff_hl_array[c_id] = hmh * l_diff_coef * g_cell_pr[c_id];

    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */
}

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
cs_gwf_soil_iso_update_itpf_terms(cs_gwf_two_phase_t     *mc)
{
  if (mc == NULL)
    return;

  /* One assumes that there is no diffusion of the gas component in the liquid
     phase --> This removes some terms w.r.t. the miscible model */

  assert(mc->henry_constant < 1e-15); /* Should be very low */

  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const cs_real_t  *g_cell_pr = mc->g_cell_pressure;

  /* In the immiscible case, mc->l_diffusivity_h should be set to 0 */

  const cs_real_t  *l_sat = mc->l_saturation->val;
  const cs_real_t  *l_cap = mc->l_capacity;
  const cs_real_t  *krl =  mc->l_rel_permeability;
  const cs_real_t  *krg =  mc->g_rel_permeability;

  for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = _soils[soil_id];
    assert(soil != NULL);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != NULL);

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  phi_rhol = phi * mc->l_mass_density;

    const double  l_diff_coef = k_abs/mc->l_viscosity;
    const double  wl_diff_coef = mc->l_mass_density * l_diff_coef;
    const double  hg_diff_coef = k_abs/mc->g_viscosity;

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
    cs_lnum_t c_min_id = 0, c_max_id = 0, c_mid_id = 0;

    if (zone->n_elts > 0) {

      c_min_id = zone->n_elts - 1, c_max_id = 0;
      c_mid_id = zone->elt_ids[zone->n_elts/2];

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];
        if (c_id < c_min_id) c_min_id = c_id;
        if (c_id > c_max_id) c_max_id = c_id;

      }

    } /* n_elts > 0 */
#endif

    /* Main loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

      const double  rhog = mh_ov_rt * g_cell_pr[c_id];
      const double  rho_hl = hmh *  g_cell_pr[c_id];
      const double  sl = l_sat[c_id];
      const double  sg = 1 - sl;
      const double  dsl_dpc = l_cap[c_id];

      /* Water conservation equation. Updates arrays linked to properties which
         define computed terms */

      mc->time_wg_array[c_id] = phi_rhol * dsl_dpc;
      mc->time_wl_array[c_id] = -mc->time_wg_array[c_id];

      mc->diff_wl_array[c_id] = wl_diff_coef * krl[c_id];

      /* Hydrogen conservation equation. Updates arrays linked to properties
         which define computed terms */

      mc->time_hg_array[c_id] = phi *
        ( mh_ov_rt*sg  + (rho_hl - rhog)*dsl_dpc + hmh*sl );
      mc->time_hl_array[c_id] = phi * (rhog - rho_hl) * dsl_dpc;

      mc->diff_hg_array[c_id] = rhog * hg_diff_coef * krg[c_id];
      mc->diff_hl_array[c_id] = rho_hl * l_diff_coef * krl[c_id];

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
      if (c_id == c_min_id || c_id == c_mid_id || c_id == c_max_id) {

        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d |wl block| time_pty % 6.4e, diff_pty % 6.4e\n",
                      c_id, mc->time_wl_array[c_id], mc->diff_wl_array[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         |wg block| time_pty % 6.4e\n",
                      mc->time_wg_array[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         |hl block| time_pty % 6.4e, diff_pty % 6.4e\n",
                      mc->time_hl_array[c_id], mc->diff_hl_array[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         |hg block| time_pty % 6.4e, diff_pty % 6.4e\n",
                      mc->time_hg_array[c_id], mc->diff_hg_array[c_id]);
#if CS_GWF_SOIL_DBG > 1
        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d >> rhog %6.4e, Pg_cell % 6.4e, Pc_cell %6.4e\n",
                      c_id, rhog, g_cell_pr[c_id],
                      mc->capillarity_cell_pressure[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         >> sl=%6.4e,"
                      " dsl_dpc=% 6.4e, krl=%6.4e, krg=%6.4e\n",
                      sl, dsl_dpc, krl[c_id], krg[c_id]);
#endif
      }
#endif
    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */
}

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
                                       cs_gwf_two_phase_t      *mc)
{
  if (mc == NULL)
    return;

  /* One assumes that there is no diffusion of the gas component in the liquid
     phase --> This removes some terms w.r.t. the miscible model */

  assert(mc->henry_constant < 1e-15); /* Should be very low */

  const cs_real_t  inv_dtcur = 1./ts->dt[0];
  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const cs_real_t  *g_cell_pr = mc->g_cell_pressure;

  /* In the immiscible case, mc->l_diffusivity_h is set to 0 */

  const cs_real_t  *l_sat = mc->l_saturation->val;
  const cs_real_t  *l_sat_pre = mc->l_saturation->val_pre;
  const cs_real_t  *krl =  mc->l_rel_permeability;
  const cs_real_t  *krg =  mc->g_rel_permeability;

  for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = _soils[soil_id];
    assert(soil != NULL);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != NULL);

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  phi_rhol = phi * mc->l_mass_density;

    const double  l_diff_coef = k_abs/mc->l_viscosity;
    const double  wl_diff_coef = mc->l_mass_density * l_diff_coef;
    const double  hg_diff_coef = k_abs/mc->g_viscosity;

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
    cs_lnum_t c_min_id = 0, c_max_id = 0, c_mid_id = 0;

    if (zone->n_elts > 0) {

      c_min_id = zone->n_elts - 1, c_max_id = 0;
      c_mid_id = zone->elt_ids[zone->n_elts/2];

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];
        if (c_id < c_min_id) c_min_id = c_id;
        if (c_id > c_max_id) c_max_id = c_id;

      }

    } /* n_elts > 0 */
#endif

    /* Main loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

      const double  rhog = mh_ov_rt * g_cell_pr[c_id];
      const double  rho_hl = hmh *  g_cell_pr[c_id];
      const double  sl = l_sat[c_id];         /* last (k,n+1) computed value */
      const double  sl_pre = l_sat_pre[c_id]; /* previous (n) computed value */
      const double  sg = 1 - sl;
      const double  dsl_dt = (sl - sl_pre) * inv_dtcur;

      /* Water conservation equation. Updates arrays linked to properties which
         define computed terms */

      mc->diff_wl_array[c_id] = wl_diff_coef * krl[c_id];
      mc->srct_wl_array[c_id] = -phi_rhol * dsl_dt;

      /* Hydrogen conservation equation. Updates arrays linked to properties
         which define computed terms */

      mc->time_hg_array[c_id] = phi * ( mh_ov_rt*sg  + hmh*sl );
      mc->diff_hg_array[c_id] = rhog * hg_diff_coef * krg[c_id];
      mc->reac_hg_array[c_id] = phi * (hmh - mh_ov_rt) * dsl_dt;
      mc->diff_hl_array[c_id] = rho_hl * l_diff_coef * krl[c_id];

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
      if (c_id == c_min_id || c_id == c_mid_id || c_id == c_max_id) {

        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d |wl block| diff_pty % 6.4e\n",
                      c_id, mc->diff_wl_array[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         |hg block| time_pty % 6.4e, diff_pty % 6.4e,"
                      " reac_pty % 6.4e\n",
                      mc->time_hg_array[c_id], mc->diff_hg_array[c_id],
                      mc->reac_hg_array[c_id]);
#if CS_GWF_SOIL_DBG > 1
        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d >> rhog %6.4e, Pg_cell % 6.4e, Pc_cell %6.4e\n",
                      c_id, rhog, g_cell_pr[c_id],
                      mc->capillarity_cell_pressure[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         >> sl=%6.4e, sl_pre=%6.4e, krl=%6.4e,"
                      " krg=%6.4e\n", sl, sl_pre, krl[c_id], krg[c_id]);
#endif
      }
#endif
    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */
}

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
                                               cs_gwf_two_phase_t      *mc)
{
  if (mc == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 1
  cs_iter_algo_t  *algo = mc->nl_algo;
  cs_real_t  *wl_time = NULL, *wl_src = NULL;
  cs_real_t  *hg_src = NULL, *hg_time = NULL;
  cs_real_t  *slc = NULL, *hg_reac = NULL;

  BFT_MALLOC(wl_time, connect->n_cells, cs_real_t);
  BFT_MALLOC(wl_src, connect->n_cells, cs_real_t);
  BFT_MALLOC(hg_src, connect->n_cells, cs_real_t);
  BFT_MALLOC(hg_time, connect->n_cells, cs_real_t);
  BFT_MALLOC(hg_reac, connect->n_cells, cs_real_t);
  BFT_MALLOC(slc, connect->n_cells, cs_real_t);
#endif

  /* One assumes that there is no diffusion of the gas component in the liquid
     phase --> This removes some terms w.r.t. the miscible model */

  assert(mc->henry_constant < 1e-15); /* Should be very low */

  const cs_real_t  inv_dtcur = 1./ts->dt[0];
  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const double  delta_h = hmh - mh_ov_rt;
  assert(delta_h < 0);

  const cs_real_t  *g_pr = mc->g_pressure->val;
  const cs_real_t  *g_pr_pre = mc->g_pressure->val_pre;
  const cs_real_t  *l_pr = mc->l_pressure->val;
  const cs_real_t  *l_pr_pre = mc->l_pressure->val_pre;
  const cs_real_t  *g_cell_pr = mc->g_cell_pressure;

  /* In the immiscible case, mc->l_diffusivity_h is set to 0 */

  const cs_real_t  *l_sat = mc->l_saturation_submesh;
  const cs_real_t  *l_sat_pre = mc->l_saturation_submesh_pre;
  const cs_real_t  *krl =  mc->l_rel_permeability;
  const cs_real_t  *krg =  mc->g_rel_permeability;

  for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = _soils[soil_id];
    assert(soil != NULL);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != NULL);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
    cs_lnum_t c_min_id = 0, c_max_id = 0, c_mid_id = 0;

    if (zone->n_elts > 0) {

      c_min_id = zone->n_elts - 1, c_max_id = 0;
      c_mid_id = zone->elt_ids[zone->n_elts/2];

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];
        if (c_id < c_min_id) c_min_id = c_id;
        if (c_id > c_max_id) c_max_id = c_id;

      }

    } /* n_elts > 0 */
#endif

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi_rhol = soil->porosity * mc->l_mass_density;
    const double  phi_ov_dt = soil->porosity * inv_dtcur;
    const double  phi_rhol_ov_dt = phi_rhol * inv_dtcur;
    const double  l_diff_coef = k_abs/mc->l_viscosity;
    const double  wl_diff_coef = mc->l_mass_density * l_diff_coef;
    const double  hg_diff_coef = k_abs/mc->g_viscosity;

    /* Main loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 1
      wl_time[c_id] = 0, hg_time[c_id] = 0, hg_reac[c_id] = 0;
      wl_src[c_id] = 0, hg_src[c_id] = 0., slc[c_id] = 0;
#endif

      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

        const cs_lnum_t  v = c2v->ids[j];
        const double  pg = g_pr[v];            /* last (k,n+1) value */
        const double  pg_pre = g_pr_pre[v];    /* last (n) value */
        const double  pl = l_pr[v];            /* last (k,n+1) value */
        const double  pl_pre = l_pr_pre[v];    /* last (n) value */
        const double  lc = mc->l_capacity[j];

        const double  sl = l_sat[j];           /* last (k,n+1) value */

        /* Water equation */

        if (lc > 0)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid dsl/dpc > 0\n", __func__);

        mc->time_wl_array[j] = -lc*phi_rhol;
        mc->srct_wl_array[j] = phi_rhol_ov_dt * lc * (pg_pre - pg);

        /* Hydrogen equation */

        mc->time_hg_array[j] = soil->porosity * (mh_ov_rt*(1-sl) + hmh*sl);
        mc->srct_hg_array[j] = 0;
        mc->reac_hg_array[j] = 0;

        const double  coef = phi_ov_dt * lc * delta_h;

        if (pg > 0)
          mc->time_hg_array[j] += lc * delta_h * pg * soil->porosity;

        else { /* pg < 0 */

          /* lc * delta_h > 0 ==> lc * delta_h * pg < 0 */

          if (pg > pg_pre)
            mc->reac_hg_array[j] += coef * (pg - pg_pre);

          else {  /* pg < pg_pre and pg < 0 */

            mc->reac_hg_array[j] += -coef * pg_pre;
            mc->srct_hg_array[j] += -coef * pg*pg;

          }

        } /* sign of pg */

        if (pl < pl_pre)
          mc->reac_hg_array[j] += coef * (pl_pre - pl);

        else { /* pl > pl_pre */

          if (pg > 0)
            mc->srct_hg_array[j] += coef * (pl - pl_pre) * pg;

          else {  /* pg < 0 */

            mc->reac_hg_array[j] += -coef * pl;
            mc->srct_hg_array[j] += -coef * pl_pre * pg;

          }

        }

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 1
        wl_time[c_id] += mc->time_wl_array[j];
        wl_src[c_id] += mc->srct_wl_array[j];
        hg_time[c_id] += mc->time_hg_array[j];
        hg_src[c_id] += mc->srct_hg_array[j];
        hg_reac[c_id] += mc->srct_hg_array[j];
        slc[c_id] += sl;
#endif
      } /* Loop on cell vertices */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 1
      const double  ratio = 1./(c2v->idx[c_id+1] - c2v->idx[c_id]);

      slc[c_id] *= ratio;
      hg_time[c_id] *= ratio;
      hg_reac[c_id] *= ratio;
      hg_src[c_id] *= ratio;
      wl_src[c_id] *= ratio;
      wl_time[c_id] *= ratio;
#endif

      /* Water conservation equation. Updates arrays linked to properties which
         define computed terms */

      mc->diff_wl_array[c_id] = wl_diff_coef * krl[c_id];

      /* Hydrogen conservation equation. Updates arrays linked to properties
         which define computed terms */

      const double  rhog = mh_ov_rt * g_cell_pr[c_id];
      const double  rho_hl = hmh * g_cell_pr[c_id];

      mc->diff_hg_array[c_id] = rhog * hg_diff_coef * krg[c_id];
      mc->diff_hl_array[c_id] = rho_hl * l_diff_coef * krl[c_id];

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 0
      if (c_id == c_min_id || c_id == c_mid_id || c_id == c_max_id) {

        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d |wl block| diff_pty %6.4e srct",
                      c_id, mc->diff_wl_array[c_id]);
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          cs_log_printf(CS_LOG_DEFAULT, " % 6.4e", mc->srct_wl_array[j]);
        cs_log_printf(CS_LOG_DEFAULT, "\n");

        cs_log_printf(CS_LOG_DEFAULT, "         |hg block|"
                      "diff_pty %6.4e (time_pty; reac; srct)",
                      mc->diff_hg_array[c_id]);
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          cs_log_printf(CS_LOG_DEFAULT, " (% 6.4e; % 6.4e; % 6.4e)",
                        mc->time_hg_array[j], mc->reac_hg_array[j],
                        mc->srct_hg_array[j]);
        cs_log_printf(CS_LOG_DEFAULT, "\n");
#if CS_GWF_SOIL_DBG > 1
        cs_log_printf(CS_LOG_DEFAULT,
                      "c_id%4d >> rhog %6.4e, Pg_cell % 6.4e, Pc_cell %6.4e\n",
                      c_id, rhog, g_cell_pr[c_id],
                      mc->capillarity_cell_pressure[c_id]);
        cs_log_printf(CS_LOG_DEFAULT,
                      "         >> krl=%6.4e, krg=%6.4e\n",
                      krl[c_id], krg[c_id]);
#endif
      }
#endif
    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SOIL_DBG > 1
  char  label[32];

  sprintf(label, "wl_time_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    wl_time, NULL, NULL,
                    ts);
  sprintf(label, "wl_src_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    wl_src, NULL, NULL,
                    ts);
  sprintf(label, "hg_src_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    hg_src, NULL, NULL,
                    ts);
  sprintf(label, "hg_time_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    hg_time, NULL, NULL,
                    ts);
  sprintf(label, "hg_reac_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    hg_reac, NULL, NULL,
                    ts);
  sprintf(label, "sl_iter%02d", algo != NULL ? algo->n_algo_iter : -1);
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_DEFAULT,
                    label,
                    1,
                    false,
                    false,
                    CS_POST_TYPE_cs_real_t,
                    slc, NULL, NULL,
                    ts);

  BFT_FREE(wl_time);
  BFT_FREE(wl_src);
  BFT_FREE(hg_src);
  BFT_FREE(hg_time);
  BFT_FREE(hg_reac);
  BFT_FREE(slc);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
