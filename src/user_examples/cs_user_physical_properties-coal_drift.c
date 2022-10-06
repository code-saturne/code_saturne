/*============================================================================
 * User definition of physical properties.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties-coal_drift.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  /*![init]*/
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_lnum_t n_cells_ext = domain->mesh->n_cells_with_ghosts;

  /* Number of fields */
  const int n_fields = cs_field_n_fields();

  /* Key id for drift scalar */
  const int keydri = cs_field_key_id("drift_scalar_model");

  /* Key id of the coal scalar class */
  const int keyccl = cs_field_key_id("scalar_class");

  const int *iym1 = cs_glob_lagr_coal_comb->iym1;

  cs_real_t *cpro_ym1_3 = cs_field_by_id(iym1[2])->val;
  cs_real_t *cpro_ym1_5 = cs_field_by_id(iym1[4])->val;
  cs_real_t *cpro_ym1_7 = cs_field_by_id(iym1[6])->val;
  cs_real_t *cpro_ym1_8 = cs_field_by_id(iym1[7])->val;
  cs_real_t *cpro_ym1_9 = cs_field_by_id(iym1[8])->val;
  cs_real_t *cpro_ym1_11 = cs_field_by_id(iym1[10])->val;
  cs_real_t *cpro_ym1_12 = cs_field_by_id(iym1[11])->val;

  cs_real_t *visco;
  BFT_MALLOC(visco, n_cells_ext, cs_real_t);
  /*![init]*/

  /* The following examples should be adapted by the user
   * ---------------------------------------------------- */

  /* Example
   * ------- */

  /*![example_1]*/

  /* Temperature */
  const cs_real_t *cpro_temp = CS_F_(t)->val;

  /* Gas density */
  cs_real_t *cpro_rom1 = cs_field_by_name("rho_gas")->val;

  /* First initialization */
  if (cs_glob_time_step->nt_cur <= 1) {
    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
    const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

    const int n_class = cs_glob_combustion_model->coal.nclacp;
    const double *rho20  = cs_glob_combustion_model->coal.rho20;
    const double *diam20 = cs_glob_combustion_model->coal.diam20;

    cs_array_set_value_real(n_cells, 1, viscl0, visco);
    cs_array_set_value_real(n_cells, 1, ro0, cpro_rom1);

    for (int icla = 0; icla < n_class; icla++) {
      char rho_name[64], diam_name[64];

      snprintf(rho_name, 63, "rho_p_%02d", icla+1);
      snprintf(diam_name, 63, "diam_p_%02d", icla+1);

      cs_real_t *cpro_rom2 = cs_field_by_name(rho_name)->val;
      cs_real_t *cpro_diam2 = cs_field_by_name(diam_name)->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_rom2[c_id] = rho20[icla];
        cpro_diam2[c_id] = diam20[icla];
      }
    }
  }

  /*--- Gas viscosity function of temperature ---*/

  /*---1-O2 2-CO 3-H2 4-N2 5-SO2 6-NH3 7-CO2-----*/

  const cs_real_t aa1 = 4.0495e-6, bb1 = 6.2200e-8;
  const cs_real_t cc1 = -2.3032e-11, dd1 = 4.4077e-15;

  const cs_real_t aa2 = 9.9987e-6, bb2 = 5.1578e-8;
  const cs_real_t cc2 = -1.8383e-11, dd2 = 3.33307e-15;

  const cs_real_t aa3 = 2.8940e-6, bb3 = 2.22508e-8;
  const cs_real_t cc3 = -8.041e-12, dd3 = 1.4619e-15;

  const cs_real_t aa4 = 4.3093e-6, bb4 = 5.0516e-8;
  const cs_real_t cc4 = -1.7869e-11, dd4 = 3.2136e-15;

  const cs_real_t aa5 = -1.9889e-6, bb5 = 5.365e-8;
  const cs_real_t cc5 = -1.4286e-11, dd5 = 2.1639e-15;

  const cs_real_t aa6 = -1.293e-6, bb6 = 4.1194e-8;
  const cs_real_t  cc6 = -1.7720e-11, dd6 = 1.8699e-15;

  const cs_real_t aa7 = 4.4822e-7, bb7 = 5.4327e-8;
  const cs_real_t  cc7 = -1.7581e-11, dd7 = 2.9979e-15;

  /*--------------------------------------------------------------------
   *  law                    mu   = a + b T + c T**2 + d T**3
   *      so      cpro_viscl(iel) = a +b*xvart+c*xvart**2 + d*xvart**3
   *--------------------------------------------------------------------*/

  if (cs_glob_time_step->nt_cur > 1) {

     for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

       const cs_real_t xvart = cpro_temp[c_id];
       const cs_real_t xvart2 = cs_math_pow2(cpro_temp[c_id]);
       const cs_real_t xvart3 = cs_math_pow3(cpro_temp[c_id]);

       const cs_real_t visco_O2  = aa1 + xvart*bb1 + cc1*xvart2 + dd1*xvart3;
       const cs_real_t visco_CO  = aa2 + xvart*bb2 + cc2*xvart2 + dd2*xvart3;
       const cs_real_t visco_H2  = aa3 + xvart*bb3 + cc3*xvart2 + dd3*xvart3;
       const cs_real_t visco_N2  = aa4 + xvart*bb4 + cc4*xvart2 + dd4*xvart3;
       const cs_real_t visco_SO2 = aa5 + xvart*bb5 + cc5*xvart2 + dd5*xvart3;
       const cs_real_t visco_NH3 = aa6 + xvart*bb6 + cc6*xvart2 + dd6*xvart3;
       const cs_real_t visco_CO2 = aa7 + xvart*bb7 + cc7*xvart2 + dd7*xvart3;

       /* Viscosity of the mixing */
       visco[c_id] =   (  cpro_ym1_8[c_id] * visco_O2
                        + cpro_ym1_3[c_id] * visco_CO
                        + cpro_ym1_5[c_id] * visco_H2
                        + cpro_ym1_12[c_id]* visco_N2
                        + cpro_ym1_11[c_id]* visco_SO2
                        + cpro_ym1_7[c_id] * visco_NH3
                        + cpro_ym1_9[c_id] * visco_CO2)
                     / (  cpro_ym1_8[c_id] + cpro_ym1_3[c_id]
                        + cpro_ym1_5[c_id] + cpro_ym1_12[c_id]
                        + cpro_ym1_11[c_id] + cpro_ym1_7[c_id]
                        + cpro_ym1_9[c_id]);

     }

  }

  /* get x1 = 1 - sum cpro_x2 */
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;

  /* All gas scalars have the same drift as if1m(1)
   *-----------------------------------------------*/

  /* Loop on fields */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *fld =  cs_field_by_id(f_id);

    /* Index of the scalar class (<0 if the scalar belongs to the gas phase) */
    int icla =  cs_field_get_key_int(fld, keyccl);
    int iscdri = cs_field_get_key_int(fld, keydri);

    /*  We only handle here one scalar with a drift per gas class */
    if ((icla < 0) && (iscdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX)) {

      /* Name of the drift scalar */
      cs_real_t *cpro_taupg
        = cs_field_by_composite_name("drift_tau", fld->name)->val;

      /* Initialize to 0 */
      cs_array_set_value_real(n_cells, 1, 0., cpro_taupg);
    }
  }

  /* Loop over coal particle classes
   * We only handle here coal class with a drift
   * ------------------------------------------- */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *fld =  cs_field_by_id(f_id);

    /* Index of the scalar class (<0 if the scalar belongs to the gas phase) */
    int icla =  cs_field_get_key_int(fld, keyccl);
    int iscdri = cs_field_get_key_int(fld, keydri);

    /* We only handle here one scalar with a drift per particle class */
    if ((icla > -1) && (iscdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX)) {

      char rho_name[64], diam_name[64], x2_name[64];
      snprintf(rho_name, 63, "rho_p_%02d", icla+1);
      snprintf(diam_name, 63, "diam_p_%02d", icla+1);
      snprintf(x2_name, 63, "x_p_%02d", icla+1);

      cs_real_t *cpro_x2 = cs_field_by_name(x2_name)->val;
      cs_real_t *cpro_rom2 = cs_field_by_name(rho_name)->val;
      cs_real_t *cpro_diam2 = cs_field_by_name(diam_name)->val;

      /* Name of the drift scalar */
      cs_real_t *cpro_taup
        = cs_field_by_composite_name("drift_tau", fld->name)->val;

      /* Computation of the relaxation time of the particles
       * the drift is therefore v_g = tau_p * g
       * --------------------------------------------------- */

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        /* Simple model for Low Reynolds Numbers */
        cpro_taup[c_id] =    cpro_x1[c_id] * cpro_rom2[c_id]
                           * cs_math_pow2(cpro_diam2[c_id])
                           / (18.0*visco[c_id]);
      }

      /* Drift for the gas:
       * tau_pg = - Sum_i X2_i v_gi */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_taup[c_id] += -(cpro_taup[c_id]*cpro_x2[c_id]);
    }
  }

  /* Free memory */
  BFT_FREE(visco);

  /*![example_1]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
