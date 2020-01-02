/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* Permeability in each subdomain */
static const double k1 = 1e5, k2 = 1;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the explicit definition of the problem for the Richards eq.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time      when ?
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      pt_ids    list of elements ids (to access coords and fill)
 * \param[in]      coords    where ?
 * \param[in]      compact   true:no indirection, false:indirection for filling
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval    result of the function
 */
/*----------------------------------------------------------------------------*/

static inline void
get_tracer_sol(cs_real_t          time,
               cs_lnum_t          n_points,
               const cs_lnum_t   *pt_ids,
               const cs_real_t   *xyz,
               bool               compact,
               void              *input,
               cs_real_t         *retval)
{
  CS_UNUSED(input);

  /* Physical parameters */
  const double  magnitude = 2*k1/(k1 + k2);
  const double  x_front = magnitude * time;

  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t  i = 0; i < n_points; i++) {

      const cs_lnum_t  id = pt_ids[i];
      const double  x = xyz[3*id];
      if (x <= x_front)
        retval[id] = 1;
      else
        retval[id] = 0;
    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t  i = 0; i < n_points; i++) {
      const double  x = xyz[3*pt_ids[i]];
      if (x <= x_front)
        retval[i] = 1;
      else
        retval[i] = 0;
    }

  }
  else {

    for (cs_lnum_t  i = 0; i < n_points; i++) {
      const double  x = xyz[3*i];
      if (x <= x_front)
        retval[i] = 1;
      else
        retval[i] = 0;
    }

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for each soil and tracer how is defined each term of the
 *         the tracer equation. Soils and tracer equations have to be added
 *         previously
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

  /*! [param_cdo_gwf_set_soil] */
  {
    cs_gwf_soil_t  *s1 = cs_gwf_soil_by_name("soil1");
    cs_gwf_set_iso_saturated_soil(s1,
                                  k1,     // saturated permeability
                                  1.0,    // saturated moisture
                                  1.0);   // bulk density (useless)

    cs_gwf_soil_t  *s2 = cs_gwf_soil_by_name("soil2");
    cs_gwf_set_iso_saturated_soil(s2,
                                  k2,     // saturated permeability
                                  1.0,    // saturated moisture
                                  1.0);   // bulk density (useless)
  }
  /*! [param_cdo_gwf_set_soil] */

  /* ===========
     Set tracers (soil have to be defined first)
     ===========

     Set parameters related to each tracer equation in each soil */

  /*! [param_cdo_gwf_set_tracer] */
  {
    cs_gwf_tracer_t *tr = cs_gwf_tracer_by_name("Tracer_01");
    cs_gwf_set_standard_tracer(tr,
                               NULL,        // zone name or NULL for all
                               0.,          // water molecular diffusivity
                               0., 0.,      // alpha (longi. and transvesal)
                               0.,          // distribution coef.
                               0.);         // 1st order decay coef.
  }
  /*! [param_cdo_gwf_set_tracer] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
