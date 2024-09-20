/*============================================================================
 * Base thermal model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_array.h"
#include "cs_dispatch.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_location.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"


#include "cs_xdef.h"
#include "cs_cf_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_thermal_model.c
        base thermal model data.

  \struct cs_thermal_model_t

  \brief Thermal model descriptor.

  Members of this thermal model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_thermal_model_t::thermal_variable
        Thermal variable solved for this physical model.
           - 0: no thermal model
           - 1: temperature
           - 2: enthalpy
           - 3: total energy (only for compressible module)
           - 4: internal energy
  \var  cs_thermal_model_t::itherm
        \deprecated alias/old name for thermal_variable

  \var  cs_thermal_model_t::temperature_scale
        Temperature scale
        The specification of the temperature scale in a consistent
        manner with the values used (initial and boundary conditions)
        is especially important in case of radiation modelling.
        - 0: none
        - 1: Kelvin
        - 2: Celsius
  \var  cs_thermal_model_t::itpscl
        \deprecated alias/old name for temperature_scale
  \var  cs_thermal_model_t::has_kinetic_st
        Take kinetic source term in energy equation into account
        (see Amino, Flageul, Carissimo, Tiselj, Benhamadouche, Ferrand 2022)
        - 0 (default)
        - 1
  \var  cs_thermal_model_t::cflt
        Take kinetic source term in energy equation into account
        (see Amino, Flageul, Carissimo, Tiselj, Benhamadouche, Ferrand 2022)
        - false (default)
        - true
  \var  cs_thermal_model_t::cflp
        Take kinetic source term in energy equation into account
        (see Amino, Flageul, Carissimo, Tiselj, Benhamadouche, Ferrand 2022)
        - false (default)
        - true
  \var  cs_thermal_model_t::has_pdivu
        Add to the right hand side the term equal to -p div(u)
  \var  cs_thermal_model_t::has_dissipation
        Add to the right hand side the thermal dissipation term
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main thermal model structure and associated pointer */

static cs_thermal_model_t  _thermal_model = {
  {.thermal_variable = CS_THERMAL_MODEL_INIT},
  {.temperature_scale = CS_TEMPERATURE_SCALE_KELVIN},
  .has_kinetic_st = 0,
  .cflt = false,
  .cflp = false,
  .has_pdivu = 0,
  .has_dissipation = 0
};

const cs_thermal_model_t  *cs_glob_thermal_model = &_thermal_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_thermal_model_get_pointers(int     **itherm,
                                int     **itpscl);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to members of the global thermal model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   itherm  --> pointer to cs_glob_thermal_model->thermal_variable
 *   itpscl  --> pointer to cs_glob_thermal_model->temperature_scale
 *----------------------------------------------------------------------------*/

void
cs_f_thermal_model_get_pointers(int  **itherm,
                                int  **itpscl)
{
  *itherm = (int *) &(_thermal_model.thermal_variable);
  *itpscl = (int *) &(_thermal_model.temperature_scale);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Return thermal field (temperature, enthalpy, total energy according to
 *        thermal model).
 *
 * \return   pointer to thermal field
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_thermal_model_field(void)
{
  cs_field_t *th_f;
  switch (_thermal_model.thermal_variable) {
  case CS_THERMAL_MODEL_TEMPERATURE:
    th_f = CS_F_(t);
    break;
  case CS_THERMAL_MODEL_ENTHALPY:
    th_f = CS_F_(h);
    break;
  case CS_THERMAL_MODEL_TOTAL_ENERGY:
    th_f = CS_F_(e_tot);
    break;
  default:
    th_f = nullptr;
  }

  return th_f;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cs_glob_thermal_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_thermal_model_t *
cs_get_glob_thermal_model(void)
{
  return &_thermal_model;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Print the thermal model structure to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_thermal_model_log_setup(void)
{
  int thermal_variable = cs_glob_thermal_model->thermal_variable;
  int temperature_scale = cs_glob_thermal_model->temperature_scale;

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Thermal model options\n"
                 "---------------------\n\n"
                 "  Continuous phase:\n\n"));

  const char *thermal_variable_value_str[]
    = {N_("no thermal model"),
       N_("temperature)"),
       N_("enthalpy"),
       N_("total energy")};

  const char *temperature_scale_value_str[]
    = {N_("none"),
       N_("temperature in Kelvin"),
       N_("temperature in Celsius")};

  cs_log_printf(CS_LOG_SETUP,
                ("    Thermal model\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    thermal_variable:   %d (%s)\n"),
                thermal_variable, _(thermal_variable_value_str[thermal_variable]));

  cs_log_printf(CS_LOG_SETUP,
                ("    Temperature scale\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    temperature_scale:  %d (%s)\n"),
                temperature_scale, _(temperature_scale_value_str[temperature_scale]));

  cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    Thermal variable solved: %s (field id %d)\n"),
       tf->name, tf->id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize thermal variables if needed
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_init(void)
{
  cs_field_t *f_cv = cs_field_by_name_try("isochoric_heat_capacity");
  if (f_cv != nullptr)
    cs_thermal_model_cv(f_cv->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of the square of sound velocity multiplied
 *        by gamma.
 *
 * \param[in]      cp      array of isobaric specific heat values for dry air
 * \param[in]      temp    array of temperature values
 * \param[in]      pres    array of pressure values
 * \param[in,out]  fracv   array of volume fraction values
 * \param[in,out]  fracm   array of mass fraction values
 * \param[in,out]  frace   array of energy fraction values
 * \param[out]     dc2     array of the values of the square of sound velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_c_square(const cs_real_t  cp[],
                          const cs_real_t  temp[],
                          const cs_real_t  pres[],
                          const cs_real_t  fracv[],
                          const cs_real_t  fracm[],
                          const cs_real_t  frace[],
                          cs_real_t        dc2[])
{
  CS_NO_WARN_IF_UNUSED(cp);
  CS_NO_WARN_IF_UNUSED(fracm);

  /*  Local variables */

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int ieos = cs_glob_cf_model->ieos;
  const cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  const cs_fluid_properties_t *phys_prop = cs_glob_fluid_properties;

  cs_dispatch_context ctx;

  /* no specific eos : the pressure equation is a Poisson equation */
  cs_field_t *fhyd = cs_field_by_name_try("H2");

  /* Ideal gas */
  if (ieos == CS_EOS_GAS_MIX && fhyd != nullptr) {
    /* WIP : only available in this function for hydrogen and air */
    cs_real_t rh = 4157.; /* R/MH2 */
    cs_real_t *yhyd = cs_field_by_name("H2")->val;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      dc2[c_id] = 1. / (temp[c_id]*((1. - yhyd[c_id])*rair + yhyd[c_id]*rh));
    });
  }
  else if (ieos == CS_EOS_IDEAL_GAS) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      dc2[c_id] = 1. / (rair * temp[c_id]);
    });
  }

  /* Ideal gas mixture (only water accounted for). TODO : other gases */
  else if (ieos == CS_EOS_MOIST_AIR) {
    /* B, C are the Antoine's law constants */
    const cs_real_t B = 17.438;
    const cs_real_t C = 239.78;
    const cs_real_t cvv = phys_prop->cpv0 - phys_prop->r_v_cnst;

    const cs_real_t p0 = phys_prop->p0;
    const cs_real_t rvsra = phys_prop->rvsra;
    const cs_real_t cva = phys_prop->cp0 - phys_prop->r_pg_cnst;
    const cs_real_t cvl = phys_prop->cvl;
    const cs_real_t l00 = phys_prop->l00;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      if (fracv[c_id] < frace[c_id]) {

        cs_real_t prest = pres[c_id] + p0;

        // Saturation pressure
        cs_real_t tc = temp[c_id] - cs_physical_constants_celsius_to_kelvin;
        cs_real_t ps = cs_air_pwv_sat(tc);

        // partial rho / partial p
        cs_real_t drhodp = -prest / (rair * cs_math_pow2(temp[c_id])
            * (1. - frace[c_id] + fracv[c_id] * rvsra))
          + (1. /ps)*(prest)*B*C
          /(rair*temp[c_id]*cs_math_pow2(prest*(1. /ps)
                - (1. - 1. /rvsra))
              *cs_math_pow2(1. - frace[c_id] + fracv[c_id]*rvsra)
              *cs_math_pow2(C + temp[c_id]
                - cs_physical_constants_celsius_to_kelvin));

        // partial e / partial p
        cs_real_t dedp = -(1. /ps)*(1. /rvsra)*(l00+temp[c_id]*(cvv - cvl))
          /cs_math_pow2(prest*(1. /ps) - (1. - 1. /rvsra));

        // partial rho / partial T
        cs_real_t drhodt = rair*((1. - frace[c_id] + fracv[c_id] *rvsra)
            + temp[c_id] *B *C *prest *(1./ps)
              /cs_math_pow2(prest*(1. /ps) - (1. - 1. /rvsra))
            *cs_math_pow2(C + temp[c_id]
              - cs_physical_constants_celsius_to_kelvin));

        // partial e/ partial T
        cs_real_t dedt = cs_thermal_model_demdt(prest,
                                                temp[c_id],
                                                frace[c_id],
                                                rvsra, cva, cvv, cvl, l00);

        // compute drho/dp (at constant internal energy)
        dc2[c_id] = -drhodt * dedp /dedt + drhodp;
      }
      else {
        dc2[c_id] = 1. / (rair * temp[c_id] * (1. - frace[c_id] + fracv[c_id] * rvsra));
      }
    });
  }
  else {
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., dc2);
  }

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the
 *        temperature at constant pressure.
 *
 * \param[in]  pres  array of pressure values
 * \param[in]  temp  array of temperature values (in Kelvin)
 * \param[in]  yw    array of the total water mass fraction
 * \param[in]  rvsra  ratio gas constant h2o / dry air
 * \param[in]  cva    difference between heat capacity of the dry air
 *                    and r_pg_const
 * \param[in]  cvv    difference beteween heat capacity of the water
 *                    in its gaseous phase and r_v_cnst
 * \param[in]  cpl    heat capacity of the water in its liquid phase
 * \param[in]  l00    water latent heat
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE cs_real_t
cs_thermal_model_demdt(const cs_real_t  pres,
                       const cs_real_t  temp,
                       const cs_real_t  yw,
                       const cs_real_t  rvsra,
                       const cs_real_t  cva,
                       const cs_real_t  cvv,
                       const cs_real_t  cpl,
                       const cs_real_t  l00)
{
  /*  Local variables */
  /* sat = A + B*t / (C + t) */
  cs_real_t sat = 6.4147
                + 17.438 * (temp - cs_physical_constants_celsius_to_kelvin)
                / (239.78 + temp - cs_physical_constants_celsius_to_kelvin);

  cs_real_t f = l00 - cpl*temp;

  cs_real_t d = cva*(1.0 - yw) + cpl * yw;
  //cs_real_t dem1 = exp(-sat);
  //0.622 * (cvv - cpl)/ (pres *  exp(-sat) - 0.378);
  cs_real_t demdt = d + (1. /rvsra) * (cvv - cpl)
    / (pres *  exp(-sat) - (1. - 1. /rvsra))
    + (1. /rvsra) * 17.438 * 239.78 * pres * (f + cvv * temp)
    * exp(-sat)
    / (cs_math_pow2(239.78 + temp - cs_physical_constants_celsius_to_kelvin)
        * cs_math_pow2(pres * exp(-sat) - (1. - 1. /rvsra)));

  return demdt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the derivative of the internal energy related to the
 *        temperature at constant internal energy.
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     temp    array of temperature values (in Kelvin)
 * \param[in]     yw      array of the total water mass fraction
 * \param[in]     rvsra   ratio gas constant h2o / dry air
 * \param[in]     cva     difference between heat capacity of the dry air
 *                        and r_pg_const
 * \param[in]     cvv     difference beteween heat capacity of the water
 *                        in its gaseous phase and r_v_cnst
 * \param[in]     cpl     heat capacity of the water in its liquid phase
 * \param[in]     l00     water latent heat
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE cs_real_t
cs_thermal_model_demdt_ecsnt(const cs_real_t  pres,
                             const cs_real_t  temp,
                             const cs_real_t  yw,
                             const cs_real_t  rvsra,
                             const cs_real_t  cva,
                             const cs_real_t  cvv,
                             const cs_real_t  cpl,
                             const cs_real_t  l00)
{
  // dedt at constant pressure:

  cs_real_t dedt = cs_thermal_model_demdt(pres, temp, yw,
                                          rvsra, cva, cvv, cpl, l00);

  cs_real_t sat =   6.4147 + 17.438 * (temp
                  - cs_physical_constants_celsius_to_kelvin)
              / (239.78 + temp - cs_physical_constants_celsius_to_kelvin);
  cs_real_t F = l00 - cpl *temp;

  cs_real_t D = cva *(1.0 - yw) + cpl * yw;
  // cs_real_t dem1 = exp(-sat);
  // 0.622 * (cvv - cpl)/ (pres *  exp(-sat) - 0.378);
  // dedp at constant temperature
  cs_real_t dedp =   -exp(-sat) *(1. /rvsra) *(l00 + temp  *(cvv - cpl))
                   / cs_math_pow2(pres *exp(-sat) - (1. - 1. /rvsra));
  cs_real_t dpdt = dedt /dedp;
  cs_real_t demdt = D + (1. /rvsra) * (cvv - cpl)
    / (pres *  exp(-sat) - (1. - 1. /rvsra))
    - (1. / rvsra) * (F + cvv * temp) * exp(-sat)
    / cs_math_pow2(pres * exp(-sat) - (1. - 1. /rvsra)) *
    (dpdt - 17.438 * 239.78 * exp(-sat)
     /cs_math_pow2(239.78 + temp - cs_physical_constants_celsius_to_kelvin));

  return demdt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief First pass to compute the contribution of the kinetic energy based
 *        source term from the prediction step
 *
 * \param[in]       i_massflux    inner mass flux used in the momentum equation
 * \param[in]       b_massflux    boundary mass flux used in the momentum equation
 * \param[in]       vela      velocity at previous time step
 * \param[in]       vel       velocity at iteration k
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_kinetic_st_prepare(const cs_real_t  i_massflux[],
                                    const cs_real_t  b_massflux[],
                                    const cs_real_t  vela[][3],
                                    const cs_real_t  vel[][3])
{
  if (cs_glob_thermal_model->has_kinetic_st != 1)
    return;

  const cs_real_t *restrict dt = CS_F_(dt)->val;
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_t *restrict cell_f_vol = mq->cell_f_vol;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;

  /* The source term is composed of two parts
   * - one can be fully computed in the prediction step
   * - the second has to be corrected by the ratio of density after
   *   the correction step, so it is stored in val_pre */
  cs_real_t *sk = cs_field_by_name("kinetic_energy_thermal_st")->val;
  cs_real_t *sk_pre = cs_field_by_name("kinetic_energy_thermal_st")->val_pre;

  const cs_field_t *f_vel = CS_F_(vel);
  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);
  const cs_real_t thetv = eqp_u->theta;

  cs_real_t *croma = CS_F_(rho)->val_pre;
  cs_real_t *cromaa = CS_F_(rho)->vals[2];

  cs_mesh_adjacencies_update_cell_i_faces();

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;
  const short int *cell_i_faces_sgn = ma->cell_i_faces_sgn;

  /* Get useful arrays */
  cs_real_3_t *i_uk
    = (cs_real_3_t *)cs_field_by_name("inner_face_velocity")->val;
  cs_real_3_t *b_uk
    = (cs_real_3_t *)cs_field_by_name("boundary_face_velocity")->val;

  /* Unsteady part */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t rhoa_theta = thetv * croma[c_id] + (1. - thetv) * cromaa[c_id];
    cs_real_t cnorm_v = cs_math_3_square_norm(vela[c_id]);

    /* rho^(n-1/2,k-1) |u^(n-1/2)|^2/2 */
    sk[c_id] = 0.5 * cell_f_vol[c_id] * rhoa_theta
             * cs_math_3_square_norm(vela[c_id]) / dt[c_id];

    cs_real_t cnorm_dv = cs_math_3_square_distance(vel[c_id], vela[c_id]);
    /* rho^(n-1/2,k-1) (|du^(n-1/2)|^2 - |u^(n-1/2)|^2)/2 */
    sk_pre[c_id] = 0.5 * cell_f_vol[c_id] * rhoa_theta * (cnorm_dv - cnorm_v)
                 / dt[c_id];

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      const short int sign = cell_i_faces_sgn[cidx];

      cs_real_t norm_v = cs_math_3_square_norm(i_uk[f_id]);
      cs_real_t norm_dv = cs_math_3_square_distance(i_uk[f_id], vel[c_id]);

      sk[c_id] -= sign * 0.5 * i_massflux[f_id] * norm_v;
      sk_pre[c_id] -= sign * 0.5 * i_massflux[f_id] * (norm_dv - norm_v);
    }

    /* Loop on boundary faces */
    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t f_id = cell_b_faces[cidx];

      cs_real_t norm_v = cs_math_3_square_norm(b_uk[f_id]);
      cs_real_t norm_dv = cs_math_3_square_distance(b_uk[f_id], vel[c_id]);

      sk[c_id] -= 0.5 * b_massflux[f_id] * norm_v;
      sk_pre[c_id] -= 0.5 * b_massflux[f_id] * (norm_dv - norm_v);
    }

  });

  ctx.wait(); // needed for the following synchronization

  if (m->halo != nullptr) {
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, sk);
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, sk_pre);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the computation of the kinetic energy based source term
 *
 * \param[in]       cromk1    density values at time n+1/2,k-1
 * \param[in]       cromk     density values at time n+1/2,k
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_kinetic_st_finalize(const cs_real_t  cromk1[],
                                     const cs_real_t  cromk[])
{

  if (cs_glob_thermal_model->has_kinetic_st != 1)
    return;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  /* The source term is composed of two parts
   * - one can be fully computed in the prediction step
   * - the second has to be corrected by the ratio of density after
   *   the correction step, so it is stored in val_pre */
  cs_real_t *sk = cs_field_by_name("kinetic_energy_thermal_st")->val;
  cs_real_t *sk_pre = cs_field_by_name("kinetic_energy_thermal_st")->val_pre;

  /* finalize the computation */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    sk[c_id] += cromk1[c_id]/cromk[c_id] * sk_pre[c_id];
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the kinetic source term if needed
 *
 * \param[in, out]  smbrs  RHS of the thermal equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_add_kst(cs_real_t  smbrs[])
{
  if (cs_glob_thermal_model->has_kinetic_st == 1) {

    const cs_field_t *f_sk = cs_field_by_name_try("kinetic_energy_thermal_st");

    if (f_sk != nullptr) {
      cs_dispatch_context ctx;
      const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
      const cs_real_t *kst = f_sk->val;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        smbrs[c_id] += kst[c_id];
      });

      ctx.wait(); // needed for the next cs_solve_equation CPU function
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the pressure equation.
 *
 * \param[in]       croma     density values at the last time iteration
 * \param[in]       trav2     predicted velocity
 * \param[in]       cvara_pr  pressure values at the last time iteration
 * \param[in]       i_massflux    face mass fluxes
 * \param[in, out]  cflp      CFL condition related to the pressure equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cflp(const cs_real_t  croma[],
                      const cs_real_t  trav2[][3],
                      const cs_real_t  cvara_pr[],
                      const cs_real_t  i_massflux[],
                      cs_real_t        cflp[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *restrict dt = CS_F_(dt)->val;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *alphafij = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;

  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_field_t *f_vel = CS_F_(vel);

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);
  const cs_real_3_t *i_face_u_normal = (const cs_real_3_t *)fvq->i_face_u_normal;
  const cs_real_3_t *b_face_u_normal = (const cs_real_3_t *)fvq->b_face_u_normal;

  cs_real_t theta = eqp_u->theta;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  if (eqp_u->blencv > 0 && eqp_u->ischcv == 1) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      cs_real_t surf = i_face_surf[f_id];
      cs_real_t dist = i_dist[f_id];

      cs_real_t i_dot = cs_math_3_dot_product(trav2[c_id0], i_face_u_normal[f_id]);
      cs_real_t j_dot = cs_math_3_dot_product(trav2[c_id1], i_face_u_normal[f_id]);

      cs_real_t dp = (cvara_pr[c_id0] - cvara_pr[c_id1]);

      cs_real_t fluxi = surf * dt[c_id0] / (croma[c_id0] * cell_f_vol[c_id0])
        * (  (alphafij[f_id] * i_dot + (1. - alphafij[f_id]) * j_dot)
           + (1. - theta) * dt[c_id0] * dp / dist);

      cs_real_t fluxj = - surf * dt[c_id1] / (croma[c_id1] * cell_f_vol[c_id1])
        * (  (alphafij[f_id] * i_dot + (1 - alphafij[f_id]) * j_dot)
           + (1. - theta) * dt[c_id1] * (-dp) / dist);

      if (c_id0 < n_cells)
        cs_dispatch_sum(&cflp[c_id0], fluxi, i_sum_type);
      if (c_id1 < n_cells)
        cs_dispatch_sum(&cflp[c_id1], fluxj, i_sum_type);

    });
  }
  else if (eqp_u->blencv <= 0 && eqp_u->ischcv == 1) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      cs_real_t surf = i_face_surf[f_id];
      cs_real_t dist = i_dist[f_id];

      cs_real_t dp = (cvara_pr[c_id0] - cvara_pr[c_id1]);

      cs_real_t dot = (i_massflux[f_id] > 0) ?
                       cs_math_3_dot_product(trav2[c_id0], i_face_u_normal[f_id]):
                       cs_math_3_dot_product(trav2[c_id1], i_face_u_normal[f_id]);

      cs_real_t fluxi = surf * dt[c_id0] / (croma[c_id0] * cell_f_vol[c_id0])
        * (dot + (1. - theta) * dt[c_id0] * dp / dist);


      cs_real_t fluxj = - surf * dt[c_id1] / (croma[c_id1] * cell_f_vol[c_id1])
        * (dot + (1. - theta) * dt[c_id1] * (-dp) / dist);

      if (c_id0 < n_cells)
        cs_dispatch_sum(&cflp[c_id0], fluxi, i_sum_type);
      if (c_id1 < n_cells)
        cs_dispatch_sum(&cflp[c_id1], fluxj, i_sum_type);

    });
  }

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t surf = b_face_surf[f_id];

      cs_real_t fluxi = surf * dt[c_id] / (croma[c_id] * cell_f_vol[c_id])
        * cs_math_3_dot_product(trav2[c_id], b_face_u_normal[f_id]);

      cs_dispatch_sum(&cflp[c_id], fluxi, b_sum_type);

  });

  ctx.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the isochoric heat capacity
 *
 * \param[in]     xcvv      isobaric heat capacity
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cv(cs_real_t  *xcvv)
{
  /* Get global data */
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;

  cs_dispatch_context ctx;

  if (cs_glob_cf_model->ieos == CS_EOS_MOIST_AIR) {
    /* get useful arrays and constants */
    cs_real_t *yw = cs_field_by_name("yw")->val;
    cs_real_t *yv = cs_field_by_name("yv")->val;
    cs_real_t cva = phys_pro->cp0 - phys_pro->r_pg_cnst;
    cs_real_t cvv = phys_pro->cpv0 - phys_pro->r_v_cnst;
    cs_real_t cvl = phys_pro->cvl;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      xcvv[c_id] =   cva * (1. - yw[c_id]) + cvv * yv[c_id]
                   + (yw[c_id] - yv[c_id]) * cvl;
    });
  }
  else if (cs_glob_cf_model->ieos == CS_EOS_IDEAL_GAS) {
    if (phys_pro->icp > 0) {
      cs_real_t *cp = CS_F_(cp)->val;
      cs_real_t r_pg_cnst = phys_pro->r_pg_cnst;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        xcvv[c_id] = cp[c_id] - r_pg_cnst;
      });
    }
    else {
      cs_real_t cp0 = phys_pro->cp0;
      cs_real_t r_pg_cnst = phys_pro->r_pg_cnst;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        xcvv[c_id] = cp0 - r_pg_cnst;
      });
    }
  }
  else { /* quid when ieos = CS_EOS_MOIST_AIR */
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 1.0, xcvv);
  }

  ctx.wait(); // needed for CPU cs_solve_equation.c
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and add the dissipation term of the thermal equation to
 *        its right hand side.
 *
 * \param[in]      vistot  array for the total viscosity
 * \param[in]      gradv   tensor for the velocity gradient
 * \param[in,out]  smbrs   array of equation right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_dissipation(const cs_real_t  vistot[],
                             const cs_real_t  gradv[][3][3],
                             cs_real_t        smbrs[])
{
  if (cs_glob_thermal_model->has_dissipation != 0) {

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
    const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;

    const cs_lnum_t n_cells = m->n_cells;

    cs_dispatch_context ctx;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrs[c_id]
        +=   2. * cell_f_vol[c_id] * vistot[c_id]
           * (  cs_math_pow2(gradv[c_id][0][0])
              + cs_math_pow2(gradv[c_id][1][1])
              + cs_math_pow2(gradv[c_id][2][2]))
           + 0.5 * (  cs_math_pow2(gradv[c_id][1][0] + gradv[c_id][0][1])
                    + cs_math_pow2(gradv[c_id][2][0] + gradv[c_id][0][2])
                    + cs_math_pow2(gradv[c_id][2][1] + gradv[c_id][1][2]))
           - 1./3. * cs_math_pow2(  gradv[c_id][0][0]
                                  + gradv[c_id][1][1]
                                  + gradv[c_id][2][2]);
    });

    ctx.wait(); // needed for CPU cs_solve_equation
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the Newton method to compute the temperature from the
 *        internal energy
 *
 * \param[in]       method    method used to compute the temperature
 * \param[in]       pk1       pressure values at the last inner iteration
 * \param[in]       th_scal   internal energy values
 * \param[in]       cvar_pr   pressure values
 * \param[in]       cvara_pr  pressure values at the last time iteration
 * \param[in]       yw        total water mass fraction
 * \param[in, out]  yv        vapor of water mass fraction
 * \param[in, out]  temp      temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_newton_t(int               method,
                          const cs_real_t  *pk1,
                          const cs_real_t   th_scal[],
                          const cs_real_t   cvar_pr[],
                          const cs_real_t   cvara_pr[],
                          const cs_real_t   yw[],
                          cs_real_t         yv[],
                          cs_real_t         temp[])
{
  /* Get global data */
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_field_t *f_vel = CS_F_(vel);

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;
  const cs_physical_constants_t *pc = cs_glob_physical_constants;

  cs_dispatch_context ctx;

  /* Newton method error threshold */
  cs_real_t epsy = 1.e-7;

  /* Thermodynamic constants */
  const cs_real_t cva = phys_pro->cp0 - phys_pro->r_pg_cnst;
  const cs_real_t cvv = phys_pro->cpv0 - phys_pro->r_v_cnst;
  const cs_real_t cvl = phys_pro->cvl;
  const cs_real_t l00 = phys_pro->l00;
  const cs_real_t p0 = phys_pro->p0;
  const cs_real_t ro0 = phys_pro->ro0;
  const cs_real_t rvsra = phys_pro->rvsra;

  /* Two methods can be used to correct yv
   * The first, yv_cor = 1, performs the newton method
   * using both temperature and pressure at n+1,k.
   * The second, yv_cor = 2, performs a increment on
   * the vapour mass fraction following dp */

  const cs_real_t *xyzp0 = phys_pro->xyzp0;
  const cs_real_t *gxyz = pc->gravity;
#if defined(HAVE_ACCEL)
  cs_real_t *_gxyz = nullptr, *_xyzp0 = nullptr;
  if (cs_get_device_id() > -1) {
    CS_MALLOC_HD(_gxyz, 3, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(_xyzp0, 3, cs_real_t, cs_alloc_mode);
    for (int i = 0; i < 3; i++) {
      _gxyz[i] = cs_glob_physical_constants->gravity[i];
      _xyzp0[i] = phys_pro->xyzp0[i];
    }

    cs_mem_advise_set_read_mostly(_gxyz);
    cs_mem_advise_set_read_mostly(_xyzp0);

    xyzp0 = _xyzp0;
    gxyz = _gxyz;
  }
#endif

  if (method == 1) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t xcvv =   cva * (1 - yw[c_id]) + cvv * yv[c_id]
                       + cvl * (yw[c_id] - yv[c_id]);

      temp[c_id] = th_scal[c_id] / xcvv - l00 * yv[c_id] / xcvv;

      // yvsat
      cs_real_t pres =   cvar_pr[c_id] + p0
                       + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                              cell_cen[c_id],
                                                              gxyz);

      cs_real_t tc = temp[c_id] - cs_physical_constants_celsius_to_kelvin;

      // After fortran migration to C, cs_air_yw_sat could be call
      cs_real_t xsat = cs_air_x_sat(tc, pres);
      cs_real_t ysat = xsat/(1. + xsat);
      //cs_real_t ysat = cs_air_yw_sat(tc, pres);

      cs_real_t yv_ = 0.;

      if (yv[c_id] < yw[c_id]) {
        yv_ = ysat;
        xcvv = cva * (1. - yw[c_id]) + cvv *yv_ + cvl *(yw[c_id] - yv_);
        // Estimate new temperature
        cs_real_t em_ = temp[c_id] * xcvv + l00 *yv_;
        cs_real_t errort = cs_math_fabs(th_scal[c_id] - em_);

        while (errort > epsy) {
          cs_real_t demdt = cs_thermal_model_demdt_ecsnt(pres,
                                                         temp[c_id],
                                                         yw[c_id],
                                                         rvsra,
                                                         cva,
                                                         cvv,
                                                         cvl,
                                                         l00);

          temp[c_id] = (th_scal[c_id] - em_) /demdt + temp[c_id];
          cs_real_t tc_ = temp[c_id] - cs_physical_constants_celsius_to_kelvin;
          // After fortran migration to C, cs_air_yw_sat could be call
          cs_real_t xv_ = cs_air_x_sat(tc_, pres);
          yv_ = xv_/(1. + xv_);
          //yv_ = cs_air_yw_sat(tc_, pres);
          xcvv = cva * (1 - yw[c_id]) + cvv *yv_ + cvl *(yw[c_id] - yv_);
          em_ = temp[c_id] *xcvv + l00 *yv_;
          errort = cs_math_fabs(th_scal[c_id] - em_);
        }

        if (yv_ > yw[c_id]) {
          yv[c_id] = yw[c_id];
          xcvv = cva * (1 - yw[c_id]) + cvv * yv[c_id];
          temp[c_id] = th_scal[c_id] /xcvv - yv[c_id] * l00 / xcvv;
        }
        else {
          yv[c_id] = yv_;
        }

      }
      else { //previous iteration not saturation
        // verifiy if saturation is reached:
        if (yw[c_id] > ysat) {
          yv_ = ysat;
          xcvv = cva * (1. - yw[c_id]) + cvv * yv_ + cvl * (yw[c_id] - yv_);
          // Estimate new temperature
          cs_real_t em_ = temp[c_id] * xcvv + l00 * yv_;
          cs_real_t errort = cs_math_fabs(th_scal[c_id] - em_);

          while (errort > epsy) {
            cs_real_t demdt = cs_thermal_model_demdt_ecsnt(pres,
                                                           temp[c_id],
                                                           yw[c_id],
                                                           rvsra,
                                                           cva,
                                                           cvv,
                                                           cvl,
                                                           l00);

            temp[c_id] = (th_scal[c_id] - em_) / demdt + temp[c_id];

            cs_real_t tc_ = temp[c_id] - cs_physical_constants_celsius_to_kelvin;
            // After fortran migration to C, cs_air_yw_sat could be call
            cs_real_t xv_ = cs_air_x_sat(tc_, pres);
            yv_ = xv_/(1. + xv_);
            //yv_ = cs_air_yw_sat(tc_, pres);
            xcvv = cva * (1. - yw[c_id]) + cvv *yv_ + cvl * (yw[c_id] - yv_);
            em_ = temp[c_id] *xcvv + l00 *yv_;
            errort = cs_math_fabs(th_scal[c_id] - em_);
          }

          if (yv_ > yw[c_id]) {
            yv[c_id] = yw[c_id];
            xcvv = cva * (1. - yw[c_id]) + cvv * yv[c_id];
            temp[c_id] = th_scal[c_id] / xcvv - yv[c_id] * l00 / xcvv;
          }
          else {
            yv[c_id] = yv_;
          }
        }
        else {
          yv_ = yw[c_id];
        }
      }

      if (yv_ > yw[c_id]) {
        yv[c_id] = yw[c_id];
        xcvv = cva * (1. - yw[c_id]) + cvv * yv[c_id];
        temp[c_id] = th_scal[c_id] /xcvv - yv[c_id] * l00 / xcvv;
      }
      else {
        yv[c_id] = yv_;
      }
      //tempk[c_id] = temp[c_id];

    }); /* End of loop on cells */

  }
  else {  /* if (method != 1) */

    cs_real_t _coef = (eqp_u->theta >= 1) ?  1. :  2.;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      if (yv[c_id] < yw[c_id]) {
        cs_real_t xcvv =   cva * (1 - yw[c_id])
                         + cvv * yv[c_id]
                         + cvl * (yw[c_id] - yv[c_id]);

        // Saturation pressure
        cs_real_t tc = temp[c_id] - cs_physical_constants_celsius_to_kelvin;
        cs_real_t ps = cs_air_pwv_sat(tc);

        // dyv/dp;
        cs_real_t pres =   cvar_pr[c_id] + phys_pro->p0
                         + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                cell_cen[c_id],
                                                                gxyz);

        cs_real_t dyvdp
          = - (1. /rvsra) * ps / cs_math_pow2(pres - (1. - 1./rvsra) * ps);

        // correction of yv
        cs_real_t drop = (  _coef * cvar_pr[c_id] - (_coef - 1.)
                          *  cvara_pr[c_id] - pk1[c_id]) *dyvdp;

        yv[c_id] = yv[c_id] + drop;
        temp[c_id] = th_scal[c_id] / xcvv - l00 * yv[c_id] /xcvv;

        if (yv[c_id] > yw[c_id]) {
          yv[c_id] = yw[c_id];
          xcvv = cva * (1. - yw[c_id]) + cvv * yv[c_id];
        }
        else {
          xcvv =   cva * (1. - yw[c_id])
                 + cvv * yv[c_id] + cvl * (yw[c_id] - yv[c_id]);
        }
        temp[c_id] = th_scal[c_id] / xcvv - yv[c_id] * l00 / xcvv;
        //tempk[c_id] = temp[c_id];
      }

    }); /* End of loop on cells */

  }

  ctx.wait(); // needed for the cs_solve_equation CPU function

#if defined(HAVE_ACCEL)
  CS_FREE_HD(_gxyz);
  CS_FREE_HD(_xyzp0);
#endif

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the term pdivu to the thermal equation rhs.
 *
 * \param[in, out]  smbrs     array of the right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_pdivu(cs_real_t  smbrs[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /*  Local variables */
  int thermal_variable = cs_glob_thermal_model->thermal_variable;
  int has_pdivu = cs_glob_thermal_model->has_pdivu;
  cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  cs_real_t p0 = cs_glob_fluid_properties->p0;

  if (!(has_pdivu))
    return;

  cs_real_t *pdivu;
  CS_MALLOC_HD(pdivu, n_cells_ext, cs_real_t, cs_alloc_mode);
  cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, 0., pdivu);

  cs_real_t *i_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;
  cs_real_t *b_massflux
    =  cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  /* Get pressure gradient, including the pressure increment gradient */

  const cs_real_3_t *gradp = nullptr, *gradphi = nullptr;

  cs_field_t *f_pg = cs_field_by_name_try("algo:pressure_gradient");
  if (f_pg != nullptr)
    gradp = (const cs_real_3_t *)f_pg->val;
  cs_field_t *f_pig = cs_field_by_name_try("algo:pressure_increment_gradient");
  if (f_pig != nullptr)
    gradphi = (const cs_real_3_t *)f_pig->val;

  const cs_real_t *cvara_pr = CS_F_(p)->val_pre;
  const cs_real_t *cvar_rho = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;
  const cs_real_t *temp = CS_F_(t)->val;
  const cs_real_t *tempa = CS_F_(t)->val_pre;

  const cs_field_t *f_vel = CS_F_(vel);
  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);
  const cs_real_t thetv = eqp_u->theta;
  //cs_real_t _coef = 1. + 2. * (1. - thetv);

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Case for humid air */

  /* A discuter */
  /* Devient obsolète? */

  /*if (cpro_yv != nullptr && cpro_yw != nullptr) {
    if (cpro_yva == nullptr) cpro_yva = cpro_yv;
    if (cpro_ywa == nullptr) cpro_ywa = cpro_yw;

    _pdivu_humid(thetv,
    temp, tempa,
    cvar_var, cvara_var,
    vel,
    xcvv,
    cpro_yw, cpro_ywa,
    cpro_yv, cpro_yva,
    gradp,
    gradphi,
    smbrs);
    return;
    }*/

  /* General case */
  int method = 2;

  if (   thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
      || thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {
    if (method == 1) { // deprecated
      // Interior faces contribution
      ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
        cs_lnum_t c_id0 = i_face_cells[f_id][0];
        cs_lnum_t c_id1 = i_face_cells[f_id][1];

        if (i_massflux[f_id] >= 0.) {
          //pres = (_coef *cvar_pr[ii] - _coef *cvar_pr[ii]) + p0;
          cs_real_t pres = p0;
          cs_real_t presa = cvara_pr[c_id0] + p0;

          cs_real_t flux = i_massflux[f_id] / cvar_rho[c_id0]
            * (thetv * pres + (1.-thetv) * presa);

          cs_dispatch_sum(&pdivu[c_id0],  flux, i_sum_type);
          cs_dispatch_sum(&pdivu[c_id1], -flux, i_sum_type);

          //pdivu[c_id0] += thetv *rair *temp[c_id0] *i_massflux[f_id]
          //           + (1-thetv) *rair *tempa[c_id0] *i_massflux[f_id];

          //pdivu[c_id1] -= thetv *rair *temp[c_id0] *i_massflux[f_id]
          //           + (1-thetv) *rair *tempa[c_id0] *i_massflux[f_id];

        }
        else {
          cs_real_t pres = p0;
          cs_real_t presa = cvara_pr[c_id1] + p0;

          cs_real_t fluxi = i_massflux[f_id] / cvar_rho[c_id1]
            * (thetv * pres + (1.-thetv) * presa);

          cs_dispatch_sum(&pdivu[c_id0],  fluxi, i_sum_type);
          cs_dispatch_sum(&pdivu[c_id1], -fluxi, i_sum_type);

          //pdivu[c_id0] += thetv *rair *temp[c_id1] *i_massflux[f_id]
          //           + (1-thetv) *rair *tempa[c_id1] *i_massflux[f_id];
          //pdivu[c_id1] -= thetv *rair *temp[c_id1] *i_massflux[f_id]
          //           + (1-thetv) *rair *tempa[c_id1] *i_massflux[f_id];

        }
      });

      /* Boundary faces contribution */
      ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
        cs_lnum_t c_id = b_face_cells[f_id];
        cs_real_t pres = p0;
        cs_real_t presa = cvara_pr[c_id] + p0;
        cs_real_t flux = b_massflux[f_id] / cvar_rho[c_id]
          * (thetv * pres + (1.-thetv) * presa);

        cs_dispatch_sum(&pdivu[c_id], flux, b_sum_type);

        //pdivu[c_id] += thetv *rair *temp[c_id] *i_massflux[f_id]
        //             + (1-thetv) *rair *tempa[c_id] *i_massflux[f_id];

      });
    } // method 2
    else {
      // Interior faces contribution
      ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
        cs_lnum_t c_id0 = i_face_cells[f_id][0];
        cs_lnum_t c_id1 = i_face_cells[f_id][1];

        if (i_massflux[f_id] >= 0.) {
          cs_real_t flux = rair * i_massflux[f_id]
            * (thetv * temp[c_id0] + (1. - thetv) * tempa[c_id0]);

          cs_dispatch_sum(&pdivu[c_id0],  flux, i_sum_type);
          cs_dispatch_sum(&pdivu[c_id1], -flux, i_sum_type);
        }
        else {
          cs_real_t flux = rair * i_massflux[f_id]
            * (thetv * temp[c_id1] + (1. - thetv) * tempa[c_id1]);

          cs_dispatch_sum(&pdivu[c_id0],  flux, i_sum_type);
          cs_dispatch_sum(&pdivu[c_id1], -flux, i_sum_type);
        }
      });

      /* Boundary faces contribution */
      ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
        cs_lnum_t c_id = b_face_cells[f_id];

        cs_real_t flux = rair * b_massflux[f_id]
          * (thetv * temp[c_id] + (1. - thetv) * tempa[c_id]);

        cs_dispatch_sum(&pdivu[c_id], flux, b_sum_type);
      });
    }

    ctx.wait(); // needed for the following synchro
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, pdivu);

    /* pdiv(u) = div(pu) - u.grad p */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      pdivu[c_id] -= cell_f_vol[c_id]
                       * (  vel[c_id][0] *(gradp[c_id][0] + gradphi[c_id][0])
                          + vel[c_id][1] *(gradp[c_id][1] + gradphi[c_id][1])
                          + vel[c_id][2] *(gradp[c_id][2] + gradphi[c_id][2]));
      smbrs[c_id] -= pdivu[c_id];
    });

    ctx.wait(); // needed for the cs_solve_equation CPU function

  } /* End test on thermal_variable */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the CFL number related to the thermal equation
 *
 * \param[in]     croma       array of density values at the last time iteration
 * \param[in]     tempk       array of the temperature
 * \param[in]     tempka      array of the temperature at the previous time step
 * \param[in]     xcvv        array of the isochoric heat capacity
 * \param[in]     vel         array of the velocity
 * \param[in]     i_massflux  array of the inner faces mass fluxes
 * \param[in]     b_massflux  array of the boundary faces mass fluxes
 * \param[in]     cflt        CFL condition related to thermal equation
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_model_cflt(const cs_real_t  croma[],
                      const cs_real_t  tempk[],
                      const cs_real_t  tempka[],
                      const cs_real_t  xcvv[],
                      const cs_real_t  vel[][3],
                      const cs_real_t  i_massflux[],
                      const cs_real_t  b_massflux[],
                      cs_real_t        cflt[])
{
  // TODO: make it compatible for others EOS
  /* Get global data */
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  cs_real_t *restrict dt = CS_F_(dt)->val;

  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;
  const cs_field_t *f_vel = CS_F_(vel);

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();

  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;
  const short int *cell_i_faces_sgn = ma->cell_i_faces_sgn;

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);
  const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;

  cs_real_t thetv = eqp_u->theta;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;

  cs_real_3_t *gradp
    = (cs_real_3_t *)cs_field_by_name("algo:pressure_gradient")-> val;
  cs_real_3_t *gradphi
    = (cs_real_3_t *)cs_field_by_name("algo:pressure_increment_gradient")-> val;
  cs_real_t gammagp = phys_pro->cp0/(phys_pro->cp0  - phys_pro->r_pg_cnst);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Initialization */
    cflt[c_id] =   dt[c_id] * (gammagp - 1.0)
                 * (  vel[c_id][0] *(gradp[c_id][0] + gradphi[c_id][0])
                    + vel[c_id][1] *(gradp[c_id][1] + gradphi[c_id][1])
                    + vel[c_id][2] *(gradp[c_id][2] + gradphi[c_id][2]))
                 / (croma[c_id] * tempka[c_id] * xcvv[c_id]);

    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    /* Loop on interior faces of cell c_id */
    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      const short int sign = cell_i_faces_sgn[cidx];

      if (i_massflux[f_id] > 0)
        cflt[c_id] += sign * (dt[c_id] / (croma[c_id] *cell_f_vol[c_id]))
          * (i_massflux[f_id]) * (thetv*(gammagp - 1.) * tempk[c_id]/tempka[c_id]
                              + (1. - thetv) *(2. - gammagp));
    }

    /* Loop on boundary faces of cell c_id */
    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t f_id = cell_b_faces[cidx];

      cflt[c_id] += (dt[c_id] / (croma[c_id] * cell_f_vol[c_id]))
        * b_massflux[f_id] * (  thetv*(gammagp - 1.) * tempk[c_id]/tempka[c_id]
                              + (1. - thetv) * (2. - gammagp));
    }

  });

  ctx.wait(); // needed for the following synchronization
  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cflt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
