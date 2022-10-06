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
 * \file cs_user_physical_properties-richards_unsat.c
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
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  const int kivisl = cs_field_key_id("diffusivity_id");

  /* Norm of gravity vector */
  const cs_real_t gravn = cs_math_3_norm(cs_glob_physical_constants->gravity);

  /* Hydraulic head (H=h+z) and darcian velocity */
  const cs_real_t *cvar_pr = CS_F_(p)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  /* Field for the capacity table (C = grad(theta)/grad(h)) */
  cs_real_t *capacity = cs_field_by_name("capacity")->val;

  /* Field for the saturation table (theta) */
  cs_real_t *saturation = cs_field_by_name("saturation")->val;

  /* Field for the soil density table (It is bulk density) */
  cs_real_t *soil_density = cs_field_by_name("soil_density")->val;

  /* Field for tensorial dipersion */
  cs_real_6_t *visten
    = (cs_real_6_t *)cs_field_by_name("anisotropic_turbulent_viscosity")->val;

  /* Field for diffusion (dipersion and molecular diffusion)
   *  for the transport part */
  cs_real_t *cpro_vscalt = NULL;
  const cs_field_t *fld = cs_field_by_name("scalar1");
  const int ifcvsl = cs_field_get_key_int(fld, kivisl);
  if (ifcvsl >= 0)
    cpro_vscalt = cs_field_by_id(ifcvsl)->val;

  /* Field for permeability for the flow part */
  cs_real_t *permeability = NULL;
  cs_real_6_t *tensor_permeability = NULL;
  cs_field_t *f_perm = cs_field_by_name("permeability");

  if (f_perm->dim == 1)
    permeability = f_perm->val;
  else if (f_perm->dim == 6)
    tensor_permeability = (cs_real_6_t *)f_perm->val;

  /* Example physical properties for single variably saturated soil
   * -------------------------------------------------------------- */

  /* Flow part
   * ---------*/

  /*![richards_set_genuch]*/
  /* Set intrinsic permeability (only depends on soil) */

  cs_real_t ki = 1;
  cs_real_t ki_xyz[3] = {1., 1., 1e-1};

  /* Set values of the Van Genuchten model parameters */
  const cs_real_t ks_param = 0.3;
  const cs_real_t thetar_param = 0.078;
  const cs_real_t thetas_param = 0.3;
  const cs_real_t n_param = 1.56;
  const cs_real_t m_param = 1-1/n_param; /* (Mualem condition) */
  const cs_real_t l_param = 0.5;
  const cs_real_t alpha_param = 0.036;
  /*![richards_set_genuch]*/

  /* Loop on all cells */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /*![richards_set_press]*/
    /* Switch from hydraulic head (H=h+z) to pressure head (h) */
    cs_real_t darcy_h = cvar_pr[c_id];
    if (gravn > cs_math_epzero)
      darcy_h = cvar_pr[c_id] - cell_cen[c_id][2];
     /*![richards_set_press]*/

    /*![richards_sat_part]*/
    /* Saturated part (h<=0) */
    if (darcy_h >= 0) {

      capacity[c_id] = 0.;
      saturation[c_id] = thetas_param;

      if (f_perm->dim == 1) {
        permeability[c_id] = ks_param*ki;
      }
      else if (f_perm->dim == 6) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          tensor_permeability[c_id][ii] = ks_param*ki_xyz[ii];
        for (cs_lnum_t ii = 3; ii < 6; ii++)
          tensor_permeability[c_id][ii] = 0.0;
      }
    }
    /*![richards_sat_part]*/

    /*![richards_unsat_part]*/
    /* Unsaturated part (h<0) */
    else {

      const cs_real_t tmp_1 = pow(cs_math_fabs(alpha_param*darcy_h), n_param);
      const cs_real_t tmp_2 = 1.0 / (1.0 + tmp_1);
      const cs_real_t se_param = pow(tmp_2, m_param);

      capacity[c_id] =   -m_param * n_param * tmp_1 * (thetas_param-thetar_param)
                       * se_param * tmp_2 / darcy_h;
      saturation[c_id] = thetar_param + se_param*(thetas_param-thetar_param);

      if (f_perm->dim == 1) {
        permeability[c_id] = ks_param*ki*pow(se_param, l_param)*
                             cs_math_pow2(1.0-pow(1.0-tmp_2, m_param));
      }
      else if (f_perm->dim == 6) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          tensor_permeability[c_id][ii]
            = ks_param * ki_xyz[ii] * pow(se_param, l_param)
                       * cs_math_pow2(1.0 - pow(1.0-tmp_2, m_param));
        for (cs_lnum_t ii = 3; ii < 6; ii++)
          tensor_permeability[c_id][ii] = 0.0;
      }
    }
    /*![richards_unsat_part]*/
  } /* End of loop on cells */

  /* Transport part for one solute with anisotropic dispersion and sorption
   * ====================================================================== */

  /*![richards_unsat_trpt_init]*/

  /* Set values of the longitudinal and transversal dirpersivity */
  const cs_real_t darcy_anisotropic_dispersion_l = 2.0;
  const cs_real_t darcy_anisotropic_dispersion_t = 1.e-1;

  /* Set value of the molecular diffusion */
  const cs_real_t molecular_diffusion = 1.e-3;
  /*![richards_unsat_trpt_init]*/

  /*![richards_unsat_mol_diff]*/
  if (cpro_vscalt != NULL) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_vscalt[c_id] = saturation[c_id]*molecular_diffusion;
  }
  /*[richards_unsat_mol_diff]*/

  /*![richards_unsat_aniso_disp]*/

  /* Computation of the isotropic dispersivity */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Computation of the norm of the velocity */
    const cs_real_t velocity_norm = cs_math_3_norm(vel[c_id]);
    const cs_real_t velocity_norm_pe = velocity_norm + 1.e-15;

    /* Tensorial dispersion is stored in visten*/
    const cs_real_t tmp_lt
      = darcy_anisotropic_dispersion_l-darcy_anisotropic_dispersion_t;

    visten[c_id][0] =   darcy_anisotropic_dispersion_t*velocity_norm
                      + tmp_lt*cs_math_pow2(vel[c_id][0])/(velocity_norm_pe);
    visten[c_id][1] =   darcy_anisotropic_dispersion_t*velocity_norm
                      + tmp_lt*cs_math_pow2(vel[c_id][1])/(velocity_norm_pe);
    visten[c_id][2] =   darcy_anisotropic_dispersion_t*velocity_norm
                      + tmp_lt*cs_math_pow2(vel[c_id][2])/(velocity_norm_pe);

    visten[c_id][3] = tmp_lt*vel[c_id][1]*vel[c_id][0]/(velocity_norm_pe);
    visten[c_id][4] = tmp_lt*vel[c_id][1]*vel[c_id][2]/(velocity_norm_pe);
    visten[c_id][5] = tmp_lt*vel[c_id][2]*vel[c_id][0]/(velocity_norm_pe);

  }
  /*![richards_unsat_aniso_disp]*/

  /*![richards_unsat_soilwater_partition]*/

  /* Set soil density (bulk density!) for delay computation
     (delay = 1 + soil_density * K_d / saturation)*/

  cs_array_set_value_real(n_cells, 1, 1.5, soil_density);

  /* Get soil-water partition structure */
  cs_gwf_soilwater_partition_t sorption_scal;
  int key_part = cs_field_key_id("gwf_soilwater_partition");
  cs_field_get_key_struct(fld, key_part, &sorption_scal);

  /* Field for kd */
  cs_real_t *cpro_kd = cs_field_by_id(sorption_scal.ikd)->val;

  /* Field for EK model parameters (kplus and kminus) */
  cs_real_t *cpro_kplus = cs_field_by_id(sorption_scal.ikp)->val;
  cs_real_t *cpro_kminus = cs_field_by_id(sorption_scal.ikm)->val;

  /* Set sorption parameters */
  cs_array_set_value_real(n_cells, 1, 5., cpro_kd);

  /* if EK model is chosen, set specific parameters */
  cs_array_set_value_real(n_cells, 1, 1e-3, cpro_kplus);
  cs_array_set_value_real(n_cells, 1, 1e-4, cpro_kminus);

  /* Field for cpro_mxsol index (if precipitation option is activated) */
  cs_real_t *cpro_mxsol = cs_field_by_id(sorption_scal.imxsol)->val;
  cs_array_set_value_real(n_cells, 1, 10., cpro_mxsol);

  /*![richards_unsat_soilwater_partition]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
