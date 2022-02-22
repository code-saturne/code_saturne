/*============================================================================
 * User functions for input of calculation parameters.
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
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

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

/*============================================================================
 * Local User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the values associated with a mass injection.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       NULL or indirection list of elements
 * \param[in]      xyz           where ?
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res           result of the function. Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [inlet_cal_analytic_func] */
static void
_define_injection(cs_real_t           time,
                  cs_lnum_t           n_elts,
                  const cs_lnum_t    *elt_ids,
                  const cs_real_t    *xyz,
                  bool                dense_output,
                  void               *input,
                  cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(xyz);

  cs_field_t *f = input; /* Current field for which this injection is called */

  const cs_zone_t *z = cs_volume_zone_by_name("mass_injection");

  /* Injection mass flow rate */
  cs_real_t xgamma = 30000;

  if (z->measure <= 0)
    bft_error(__FILE__, __LINE__, 0,
              "In function %s: volume of zone %d (%s) is %g.",
              __func__, z->id, z->name, z->measure);

  /* Assume output directly to main arrays, so not handle dense case */
  if (dense_output == false)
    bft_error(__FILE__, __LINE__, 0,
              "Function %s currently only handles dense output.",
              __func__);

  /* For pressure
     ------------ */

  if (f == CS_F_(p)) {

    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
    const cs_real_t *cell_f_vol = mq->cell_f_vol;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      res[i] = xgamma;
    }

    /* Logging: compute/check the mass flow rate in the volume */

    static bool _logged_mass_flux = false;
    if (_logged_mass_flux == false && cs_log_default_is_active()) {

      const cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      if (eqp->verbosity >= 1 || _logged_mass_flux == false) {
        cs_real_t flucel = 0.;
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t c_id = elt_ids[i];
          flucel += cell_f_vol[c_id] * res[i];
        }
        cs_parall_sum(1, CS_REAL_TYPE, &flucel);

        cs_log_printf(CS_LOG_DEFAULT,
                      ("      Mass rate generated in the domain: %14.5e\n"
                       "      ----------------------------------\n"),
                      flucel);

        _logged_mass_flux = true;
      }

    }

  }
}
/*! [inlet_cal_analytic_func] */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* Example 1: simulation of an inlet condition by mass source terms */

  /*! [inlet_cal] */
  cs_real_t wind = 0.1;
  cs_real_t dh = 0.5;

  const char z_name[] = "mass_injection";

  /* Pressure */

  double mass_in[1] = {30000};

  cs_equation_add_volume_mass_injection_by_value
    (cs_field_get_equation_param(CS_F_(p)), z_name, mass_in);

  /* Velocity */

  double vel_in[3] = {0, wind, 0};

  cs_equation_add_volume_mass_injection_by_value
    (cs_field_get_equation_param(CS_F_(vel)), z_name, vel_in);

  /* Turbulence values */

  cs_turbulence_inflow_volume_mass_injection_ke_hyd_diam
    (z_name,
     cs_math_sq(wind),
     dh,
     cs_glob_fluid_properties->ro0,
     cs_glob_fluid_properties->viscl0);

  /* Scalars */

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);
    int scalar_id = (f->type & CS_FIELD_VARIABLE) ?
      cs_field_get_key_int(f, cs_field_key_id("scalar_id")) : -1;

    if (scalar_id >= 0) {
      double val = 1;
      cs_equation_add_volume_mass_injection_by_value
        (cs_field_get_equation_param(f), z_name, &val);
    }

  }
  /*! [inlet_cal] */

  /* Variant where an analytic function is used for the pressure;
     the function used here also computes and logs the total mass rate */

  {
    /*! [inlet_cal_analytic] */
    cs_field_t *f = CS_F_(p);

    cs_equation_add_volume_mass_injection_by_analytic
      (cs_field_get_equation_param(f),
       z_name,
       _define_injection,    // associated function
       f);                   // input structure (use field to access info)
    /*! [inlet_cal_analytic] */
  }

  /* Example 2 : simulation of a suction (by a pump for instance) with a
   *             total rate of 80 000 kg/s.
   *             The suction rate is supposed to be uniformly distributed
   *             on all the cells in the "suction_pump" zone. */

  /*! [suction_pump] */
  double mass_out[1] = {-80000};

  cs_equation_add_volume_mass_injection_by_qov
    (cs_field_get_equation_param(CS_F_(p)), "suction_pump", mass_out);
  /*! [suction_pump] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
