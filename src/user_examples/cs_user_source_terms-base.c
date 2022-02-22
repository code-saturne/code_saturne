/*============================================================================
 * Base examples for additional right-hand side source terms for
 * variable equations (momentum, scalars, turbulence...).
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
 * \file cs_user_source_terms-base.c
 *
 * \brief Base examples for additional right-hand side source terms for
 *   variable equations (momentum, scalars, turbulence...).
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
  /*! [st_meta] */
  /* field structure */
  const cs_field_t  *f = cs_field_by_id(f_id);

  /* mesh quantities */
  const cs_lnum_t  n_cells = domain->mesh->n_cells;
  const cs_real_t  *cell_f_vol = domain->mesh_quantities->cell_vol;
  /*! [st_meta] */

  /* Scalar variance indicator: if var_f_id > -1, the field
     is a variance of field var_f_id */

  /*! [field_is_variance] */
  int var_f_id = cs_field_get_key_int(f, cs_field_key_id("first_moment_id"));
  /*! [field_is_variance] */

  /* Density */

  /*! [density_2] */
  const cs_real_t  *cpro_rom = CS_F_(rho)->val;
  /*! [density_2] */

  /* Example of arbitrary source term for the scalar f, name "scalar_2"
   *
   *                       S = A * f + B
   *
   *        appearing in the equation under the form
   *
   *                  rho*df/dt = S (+ regular terms in the equation)
   *
   * In the following example:
   *   A = - rho / tauf
   *   B =   rho * prodf
   *
   * with
   *   tauf   = 10.0   [ s  ] (dissipation time for f)
   *   prodf  = 100.0  [ [f]/s ] (production of f by unit of time)
   *
   * which yields
   *   st_imp[i] = cell_f_vol[i]* A = - cell_f_vol[i]*rho/tauf
   *   st_exp[i] = cell_f_vol[i]* B =   cell_f_vol[i]*rho*prodf
   */

  /*! [src_term_applied] */
  if (strcmp(f->name, "scalar_2") == 0) {

    /* logging */

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

    if (eqp->verbosity >= 1)
      bft_printf(" User source terms for variable %s\n",
                 cs_field_get_label(f));

    /* apply source terms to all cells */

    const cs_real_t tauf  = 10.0;
    const cs_real_t prodf = 100.0;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      st_imp[i] = - cell_f_vol[i]*cpro_rom[i]/tauf;
      st_exp[i] =   cell_f_vol[i]*cpro_rom[i]*prodf;
    }

  }
  /*! [src_term_applied] */

  /* Example of arbitrary volumic heat term in the equation for enthalpy h.
   *
   * In the considered example, a uniform volumic source of heating is imposed
   * in the cells of a volume zone named "heater".
   *
   * The global heating power if Pwatt (in W); the total fluid volume of the
   * selected zone is available in the zone structure, and copied to
   * a local value voltf (in m^3).
   *
   * This yields
   *    st_imp[i] = 0
   *    st_exp[i] = cell_f_vol[i]*pwatt/voltf;
   */

  /* Warning :
   * It is assumed here that the thermal scalar is an enthalpy.
   * If the scalar is a temperature, PWatt does not need to be divided
   * by Cp because Cp is put outside the diffusion term and multiplied
   * in the temperature equation as follows:
   *
   * rho*Cp*cell_f_vol*dT/dt + .... =  cell_f_vol[i]* Pwatt/voltf
   */

  /*! [ex_3_apply] */
  if (f == CS_F_(h)) { /* enthalpy */

    /* logging */

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

    if (eqp->verbosity >= 1)
      bft_printf(" User source terms for variable %s\n",
                 cs_field_get_label(f));

    /* apply source terms in zone cells */

    const cs_zone_t *z = cs_volume_zone_by_name_try("heater");

    if (z != NULL) {

      cs_real_t pwatt = 100.0;
      cs_real_t voltf = z->f_measure;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        cs_lnum_t c_id = z->elt_ids[i];
        st_exp[c_id] = cell_f_vol[c_id]*pwatt/voltf;
      }

    }

  }
  /*! [ex_3_apply] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
