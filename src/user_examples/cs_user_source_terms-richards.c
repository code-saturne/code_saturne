/*============================================================================
 * User source terms example for the Richards equation.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms-richards.c
 *
 * \brief User source terms example for the Richards equation.
 *
 * See the reference \ref cs_user_source_terms.c for documentation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define source terms.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(st_imp);

  /*! [richards_leaching_init] */
  /* local variables */

  const cs_field_t  *f = cs_field_by_id(f_id);
  const cs_real_t  t = cs_glob_time_step->t_cur;

  const cs_real_t  *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  /* Map saturation and delay */

  char f_name[64];
  snprintf(f_name, 63, "%s_delay", f->name); f_name[63] = '\0';
  const cs_real_t  *delay = cs_field_by_name(f_name)->val;

  const cs_real_t  *saturation = cs_field_by_name("saturation")->val;
  /*! [richards_leaching_init] */

  /* logging */
  /*! [richards_leaching_log] */
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if (eqp->verbosity >= 1)
    bft_printf(" User source terms for variable %s\n",
               cs_field_get_label(f));
  /*! [richards_leaching_log] */

  /* Leaching from volume zone labeled as "LEACHING_ZONE" */

  /*! [richards_leaching] */
  cs_real_t lambda = 1.e-2;        /* first order decay coefficient */
  cs_real_t leaching_time = 1.e2;  /* leaching duration */

  const cs_zone_t *z = cs_volume_zone_by_name_try("LEACHING_ZONE");

  /* progressive leaching */

  if (z != NULL && t < leaching_time) {
    cs_real_t leaching_volume = z->f_measure;

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      cs_lnum_t c_id = z->elt_ids[i];

      st_exp[c_id] =   cell_f_vol[i] * exp(-lambda*t)
                     / (  leaching_time * leaching_volume
                        * saturation[i] * delay[i]);
    }
  }
  /*! [richards_leaching] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
