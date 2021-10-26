/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_gwf-example.c
 *
 * \brief Main user function for setting of a calculation with CDO for the
 *        groundwater flow module
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_tracy_struct] */

/* Parameters defining the saturated hydraulic Tracy model */

typedef struct {

  double    bulk_density;
  double    k_s;
  double    h_s;
  double    h_r;
  double    theta_s;
  double    theta_r;

} cs_tracy_param_t;
/*! [param_cdo_gwf_tracy_struct] */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the physical properties related to a
 *         hydraulic model. At least, moisture content and permeability are
 *         updated. If the simulation is time-depedent, then the soil capacity
 *         can be updated also.
 *
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      head_values  array of values for head used in law
 * \param[in]      zone         pointer to a cs_zone_t
 * \param[in, out] input        pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_set_user_update_soil] */
static void
_tracy_update(const cs_real_t              t_eval,
              const cs_mesh_t             *mesh,
              const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              const cs_real_t             *head_values,
              const cs_zone_t             *zone,
              void                        *input)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  /* Retrieve field values associated to properties to update */
  cs_field_t  *f = NULL;
  f = cs_field_by_name("permeability");
  cs_real_t  *permeability_values = f->val;
  f = cs_field_by_name("moisture_content");
  cs_real_t  *moisture_values = f->val;
  f = cs_field_by_name_try("soil_capacity");
  cs_real_t  *capacity_values = NULL;
  assert(f != NULL);
  capacity_values = f->val;

  const cs_tracy_param_t  *law = (cs_tracy_param_t *)input;
  const double  delta_moisture = law->theta_s - law->theta_r;

  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head_values[c_id];
    const cs_real_t  k_r = (h - law->h_r)/(law->h_s - law->h_r);
    const double  isoval = law->k_s * k_r;

    /* Set the permeability value to the saturated values */
    permeability_values[c_id] = isoval;

    /* Set the moisture content (Se = 1 in this case)*/
    moisture_values[c_id] = law->theta_r + k_r * delta_moisture;

    /* Set the capacity values */
    capacity_values[c_id] = delta_moisture /(law->h_s - law->h_r);

  } /* Loop on selected cells */

}
/*! [param_cdo_gwf_set_user_update_soil] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set free the input of a soil structure
 *
 * \param[in, out] input      pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_set_user_free_soil] */
static void
_tracy_free(void         *input)
{
  cs_tracy_param_t  *tp = (cs_tracy_param_t  *)input;

  BFT_FREE(tp);
}
/*! [param_cdo_gwf_set_user_free_soil] */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for each soil and tracer how is defined each term of the
 *         the tracer equation. Soils and tracer equations have to be added
 *         previously in \ref cs_user_model
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_gwf_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* =========
     Set soils
     ========= */

  /*! [param_cdo_gwf_set_user_soil] */
  {
    /* User-defined structure to keep what is needed to define the behavior of
       the soil */

    cs_tracy_param_t  *tp = NULL;

    BFT_MALLOC(tp, 1, cs_tracy_param_t);

    /* Set of parameters for this soil */

    tp->k_s = 1.15741e-4;
    tp->h_s = 0.;
    tp->h_r = -100;
    tp->theta_s = 0.45;
    tp->theta_r = 0.15;
    tp->bulk_density = 1.0;

    /* Retrieve the soil by its name */

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_name("cells");

    /* Set the structure and the function pointers */

    cs_gwf_set_user_soil(soil,          /* soil structure */
                         tp,            /* soil context structure */
                         _tracy_update, /* function to update the soil */
                         _tracy_free);  /* function to free the context */

    /* The function to free the context as well as the pointer to the soil
       context can be set to NULL if there is no need to consider it */
  }
  /*! [param_cdo_gwf_set_user_soil] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
