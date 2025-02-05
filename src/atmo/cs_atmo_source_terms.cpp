/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"

#include "alge/cs_gradient.h"
#include "atmo/cs_atmo.h"
#include "atmo/cs_atmo_profile_std.h"
#include "atmo/cs_intprf.h"
#include "base/cs_array.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_measures_util.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_thermal_model.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "pprt/cs_physical_model.h"
#include "rayt/cs_rad_transfer.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_source_terms.h"

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atmo_source_terms.cpp compute source terms
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static local variables
 *============================================================================*/

static bool r3_is_defined = false;
static int treated_scalars = 0;

static cs_real_t qliqmax, r3max;

static cs_real_3_t *grad1 = nullptr, *grad2 = nullptr;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_atr1vf(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field with homogeneous Neumann BC's.
 *----------------------------------------------------------------------------*/

static void
_gradient_homogeneous_neumann_sca(cs_real_t         pvar[],
                                  cs_real_3_t       grad[])
{

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

    // options for gradient calculation
  const cs_field_t *f = cs_field_by_name("ym_water");
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  int imrgra = eqp->imrgra;
  int n_r_sweeps = eqp->nswrgr;
  int iwarnp = eqp->verbosity;
  int imligp = eqp->imligr;
  cs_real_t epsrgp = eqp->epsrgr;
  cs_real_t climgp = eqp->climgr;
  char var_name[32];

  strcpy(var_name, "Work array");
  var_name[31] = '\0';

  /* Choose gradient type */

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  /* Compute gradient */

  cs_gradient_scalar(var_name,
                     gradient_type,
                     halo_type,
                     1,             /* inc */
                     n_r_sweeps,
                     0,             /* iphydp */
                     1,             /* w_stride */
                     iwarnp,
                     (cs_gradient_limit_t)imligp,
                     epsrgp,
                     climgp,
                     nullptr,          /* f_ext */
                     nullptr,          /* bc_coeffs */
                     pvar,
                     nullptr,          /* c_weight */
                     nullptr,          /* coupling */
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the mean volumic radius of the droplets
 *
 * \param[in]     cpro_rho       density
 * \param[in]     cpro_liqwt     liquid water
 * \param[in]     cvar_ntdrp     number of droplets
 * \param[out]    r3             droplet equation radius
 */
/*----------------------------------------------------------------------------*/

static void
_compute_radius_volumic_droplets(const cs_lnum_t  n_cells,
                                 const cs_real_t  *cpro_rho,
                                 const cs_real_t  *cpro_liqwt,
                                 const cs_real_t  *cvar_ntdrp,
                                 cs_real_t        *r3)
{
  const cs_real_t rho_water = 1.0e3;
  const cs_real_t conversion = 1e6; // passing from 1/cm**3 to 1/m**3

  r3max = 0.0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t nc = cvar_ntdrp[c_id];
    const cs_real_t rho  = cpro_rho[c_id];
    const cs_real_t qliq = cpro_liqwt[c_id];
    if (qliq >= 1e-8) {
      nc = fmax(nc, 1.0);
      r3[c_id] = pow(0.75/cs_math_pi*(rho*qliq)/(rho_water*nc*conversion), (1.0/3));
    }
    else {
      r3[c_id] = 0.0;
    }
    r3max = fmax(r3[c_id], r3max);
  }

  cs_parall_max(1, CS_REAL_TYPE, &r3max);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute deposition velocity
 *
 * \param[in]       tempc         Temperature in Celsius
 * \param[in]       rom           Density
 * \param[in]       pres          Pressure
 * \param[in]       cfnns         Non neutral correction coefficient for scalars
 * \param[in]       ustar         Boundary ustar
 * \param[in]       rugt          Boundary thermal roughness
 * \param[in]       rcloudvolmoy  droplet equation radius
 * \param[in]       wg            sedimentation velocity
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_compute_deposition_vel(const cs_real_t tempc,
                        const cs_real_t rom,
                        const cs_real_t pres,
                        const cs_real_t cfnns,
                        const cs_real_t ustar,
                        const cs_real_t rugt,
                        const cs_real_t rcloudvolmoy,
                        const cs_real_t wg)
{
  // deposition velocity
  cs_real_t depo_vel = 0.0;

  // deposition is computed only for first level
  const cs_real_t ckarm  = 0.4;
  const cs_real_t eps0   = 3.0;
  const cs_real_t dzmin  = 4.0;
  const cs_real_t alpha  = 1.5;
  const cs_real_t gamma  = 0.56;
  const cs_real_t arecep = 0.01;
  const cs_real_t cbolz  = 1.38e-23;

  const cs_real_t dp   = dzmin/2.0;
  const cs_real_t temp = tempc + cs_physical_constants_celsius_to_kelvin;

  const cs_real_t muair = 1.83e-5*(416.16/(temp+120.0))
                        * (pow(temp/296.16, 1.50));

  const cs_real_t nuair = muair/rom;
  const cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  const cs_real_t lpm = (2.0*muair/pres)*(sqrt(0.125*cs_math_pi*rair*temp));
  const cs_real_t ccunning = 1.0
                           + (lpm/rcloudvolmoy)
                           * (1.257 + 0.4*exp(-1.1*rcloudvolmoy/lpm));

  const cs_real_t dbrow = cbolz*temp*ccunning/(6*cs_math_pi*muair*rcloudvolmoy);
  const cs_real_t cebro = pow(nuair, (-1.0)*gamma)/dbrow;

  const cs_real_t ather = ckarm/log((dp+rugt)/rugt);
  const cs_real_t raero = 1.0 / (ather * ustar * cfnns);

  const cs_real_t st    = wg*ustar/(9.81*arecep);
  const cs_real_t ceimp = pow(st/(st+alpha), 2.);
  const cs_real_t ceint = 2.0*(pow(rcloudvolmoy/arecep, 2.));

  // Fog or cloud droplet deposition
  if (ustar > 0.0) {
    const cs_real_t rsurf = 1.0/(eps0*ustar*(ceimp+ceint+cebro)*exp(-sqrt(st)));
    depo_vel = 1.0 / (raero + rsurf);
  }
  else
    depo_vel = 0.0;

  return depo_vel;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the gradient of the two following quantities
 *        1) grad1 = rho*qliq*V(r3)*exp(5*sc^2)
 *        2) grad2 = nc*V(r3)*exp(-sc^2)
 *
 * \param[in]    m               pointer to mesh structure
 * \param[in]    mq              pointer to mesh quantities structure
 * \param[in]    at_opt          pointer to atmo options structure
 * \param[in]    phys_pro        pointer to  fluid properties structure
 * \param[in]    r3              droplet equation radius
 * \param[in]    cpro_rho        density
 * \param[in]    cpro_met_p      meteo pressure
 * \param[in]    cpro_tempc      real temperature
 * \param[in]    cpro_liqwt      liquid water
 * \param[in]    cvar_ntdrp      number of droplets
 *
 */
/*----------------------------------------------------------------------------*/

static void
_compute_gradient(const cs_mesh_t                *m,
                  const cs_mesh_quantities_t     *mq,
                  const cs_atmo_option_t         *at_opt,
                  const cs_fluid_properties_t    *phys_pro,
                  const cs_real_t                *r3,
                  const cs_real_t                *cpro_rho,
                  const cs_real_t                *cpro_met_p,
                  const cs_real_t                *cpro_tempc,
                  const cs_real_t                *cpro_liqwt,
                  const cs_real_t                *cvar_ntdrp)
{

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  if (r3max < 1.0e-10) {
    cs_array_real_fill_zero(3*n_cells, (cs_real_t*)grad1);
    cs_array_real_fill_zero(3*n_cells, (cs_real_t*)grad2);
  }

  // compute sedimentation velocity
  cs_real_t *sed_vel = nullptr;
  BFT_MALLOC(sed_vel, n_cells, cs_real_t);

  //taup g, with taup = cuning * d^2 * rhop / (18 * mu) ...
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    sed_vel[c_id] = 1.19e8 * pow(r3[c_id], 2.0);

  // take into account deposition if enabled
  if (at_opt->deposition_model > 0) {

    cs_real_t *pres = nullptr;
    BFT_MALLOC(pres, n_cells, cs_real_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (at_opt->meteo_profile == 0) {
        cs_real_t dum = 0.0;
        cs_atmo_profile_std(0.0, phys_pro->p0, phys_pro->t0,
                            cell_cen[c_id][2], &pres[c_id], &dum, &dum);
      }
      else if (at_opt->meteo_profile == 1) {
        pres[c_id] = cs_intprf(at_opt->met_1d_nlevels_t,
                               at_opt->met_1d_ntimes,
                               at_opt->z_temp_met,
                               at_opt->time_met,
                               at_opt->hyd_p_met,
                               cell_cen[c_id][2],
                               cs_glob_time_step->t_cur);
      }
      else {
        pres[c_id] = cpro_met_p[c_id];
      }
    }

    cs_real_t *ustar  = cs_field_by_name("boundary_ustar")->val;
    cs_real_t *rugd   = cs_field_by_name("boundary_roughness")->val;
    cs_real_t *rugt   = cs_field_by_name("boundary_thermal_roughness")->val;
    cs_real_t *bcfnns = cs_field_by_name("non_neutral_scalar_correction")->val;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id ++) {
      if (   cs_glob_bc_type[face_id] != CS_ROUGHWALL
          && cs_glob_bc_type[face_id] != CS_SMOOTHWALL)
        continue;
      const cs_lnum_t c_id = b_face_cells[face_id];
      if (r3[c_id] <= 0.0 || rugd[face_id] <= 0.0)
        continue;
      const cs_real_t depo_vel = _compute_deposition_vel(cpro_tempc[c_id],
                                                         cpro_rho[c_id],
                                                         pres[c_id],
                                                         bcfnns[face_id],
                                                         ustar[face_id],
                                                         rugt[face_id],
                                                         r3[c_id],
                                                         sed_vel[c_id]);
      sed_vel[c_id] += depo_vel;
    }

    BFT_FREE(pres);
  } // end deposition_model > 0

  cs_real_t *local_field = nullptr;

  /* Computation of the gradient of rho*qliq*V(r3)*exp(5*sc^2)
   * it corresponds to div(qliq* exp() rho vel_d) */
  BFT_MALLOC(local_field, n_cells_ext, cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    local_field[c_id] = cpro_rho[c_id]                 // mass density of the air kg/m3
                      * cpro_liqwt[c_id]               // total liquid water content kg/kg
                      * sed_vel[c_id]                  // deposition velocity m/s
                      * exp(5*pow(at_opt->sigc, 2.0)); // coefficient coming from log-norm
                                                       // law of the droplet spectrum

  _gradient_homogeneous_neumann_sca(local_field, grad1);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    local_field[c_id] = cpro_rho[c_id]               // mass density of the air kg/m3
                      * cvar_ntdrp[c_id]             // total liquid water content kg/kg
                      * sed_vel[c_id]                // deposition velocity m/s
                      * exp(-pow(at_opt->sigc, 2.0)); // coefficient coming from log-norm
                                                     // law of the droplet spectrum
  BFT_FREE(sed_vel);

  _gradient_homogeneous_neumann_sca(local_field, grad2);

  BFT_FREE(local_field);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[])
{
  /* Microphysics options */
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  bool rain = at_opt->rain;
  bool autoconversion = at_opt->autoconversion;
  bool accretion = at_opt->accretion;
  bool rupture = at_opt->rupture;
  bool autocollection_cloud = at_opt->autocollection_cloud;
  bool autocollection_rain = at_opt->autocollection_rain;
  int cloud_type =at_opt->cloud_type;

  /* Mesh options */
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  /* Physical properties */
  cs_real_t *rho_m = (cs_real_t *)CS_F_(rho)->val; /* Mixture density */
  cs_real_t *rho_h = cs_field_by_name_try("rho_humid_air")->val;
  cs_real_t rho_l = 1.e3;
  cs_real_t *yc = cs_field_by_name_try("liquid_water")->val;
  cs_field_t *f_ym_w = CS_F_(ym_w);     /* Water mass fraction*/
  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val;     /* Water mass fraction*/
  cs_real_t *yr = cs_field_by_name_try("ym_l_r")->val;
  cs_real_t *nc = cs_field_by_name_try("number_of_droplets")->val;
  cs_real_t *nr = cs_field_by_name_try("n_r")->val;

  /* Field ids */
  int nc_id = cs_field_by_name_try("number_of_droplets")->id;
  int yr_id = cs_field_by_name_try("ym_l_r")->id;
  int nr_id = cs_field_by_name_try("n_r")->id;

  /* Add source terms for the rain option */
  if (rain == true) {

    /* Initialize properties based on cloud type */
    cs_real_t rs;
    cs_real_t sigma;
    if (cloud_type == 0){
      rs = 9.59e-6;
      sigma = 0.15;
    }
    if (cloud_type == 1){
      rs = 7.47e-6;
      sigma = 0.28;
    }

    /* Go through each cell */
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

      /* Initialize variables to store source terms for each fields */

      /* Compute r3 */
      cs_real_t r3;
      if (nc[cell_id] >= 1){
        r3 = pow((0.75*rho_h[cell_id]*yc[cell_id])
                          / (cs_math_pi*rho_l*nc[cell_id]), 1./3.);
      }
      else {
        r3 = 0.;
      }

      /* Autoconversion */
      if (autoconversion == true && r3 >= rs) {

        cs_real_t exp_sig2 = exp( 9. * sigma * sigma);

        //FIXME find the formula with the good dimensions...
        cs_real_t a = 0.00729*(cs_math_pow4(1.e5 * r3) * sqrt(exp_sig2-1.)-0.4)
          *( 1e6*r3*pow(exp_sig2 - 1.,1./6.)-7.5);

        //TODO check a >0
        if (a < 0.){
          a= 0.;
        }

        /* Total water mass fraction (except rain) */
        if (f_id == f_ym_w->id) {
          cs_real_t _exp = a * rho_m[cell_id] * cs_math_pow2(yc[cell_id])
            * cell_f_vol[cell_id];
          imp_st[cell_id] -= _exp / fmax(ym_w[cell_id], cs_math_epzero);
          exp_st[cell_id] -= _exp;
        }
        if (f_id == yr_id) {
          //TODO should be multiplied by ym_w^{n+1}/ym_w^{n}
          exp_st[cell_id] += a * rho_m[cell_id] * cs_math_pow2(yc[cell_id])
            * cell_f_vol[cell_id];
        }
        if (f_id == nc_id) {
          cs_real_t _exp =  a * rho_m[cell_id] * cs_math_pow2(yc[cell_id])
            * cell_f_vol[cell_id];
          imp_st[cell_id] -= _exp / nc[cell_id];
          exp_st[cell_id] -= _exp;
        }
        if (f_id == nr_id) {
          exp_st[cell_id] += 3. * a * cell_f_vol[cell_id]
            / (4. * cs_math_pi * cs_math_pow3(fmax(41e-6, r3)) * rho_l);
        }
      }

      /* Accretion */
      if (accretion == true) {
        cs_real_t tau = yr[cell_id] /fmax(yc[cell_id] + yr[cell_id], cs_math_epzero);
        cs_real_t phi = tau / (tau + 5e-4);
        cs_real_t kr = 5.78;

        if (f_id == f_ym_w->id) {
          cs_real_t _exp = kr * rho_m[cell_id] * cs_math_pow2(yc[cell_id])
            * phi * cell_f_vol[cell_id];
          imp_st[cell_id] -= _exp / fmax(ym_w[cell_id], cs_math_epzero);
          exp_st[cell_id] -= _exp;
        }
        if (f_id == yr_id) {
          exp_st[cell_id] += kr * rho_m[cell_id]*yc[cell_id]*yr[cell_id]*phi*cell_f_vol[cell_id];
        }
        if (f_id == nc_id) {
          exp_st[cell_id] += kr * rho_m[cell_id]*nc[cell_id]*yr[cell_id]*phi*cell_f_vol[cell_id];
        }
      }

      /* Self Collection Cloud */
      if (autocollection_cloud == true) {
        if (f_id == nc_id) {
          cs_real_t _exp = 9.44e9*cs_math_pow2(rho_m[cell_id]*yc[cell_id])
            * exp(9*sigma*sigma) * cell_f_vol[cell_id];
          imp_st[cell_id] -= _exp / fmax(nc[cell_id], 1.);
          exp_st[cell_id] -= _exp;
        }
      }

      /* Self Collection Rain */
      if (autocollection_rain == true) {
        if (f_id == nr_id) {
          cs_real_t kr = 5.78;
          cs_real_t e = 1.;
          if (rupture) {
            if (r3 >= 3.e-4 && r3 < 1.e-3)
              e = exp(-5.e-3 * (r3 - 3.e-4));
            else if (r3 >= 1.e-3)
              e = 0.;
          }
          cs_real_t _yr = fmax(yr[cell_id], 0.);
          imp_st[cell_id] -= e * kr * rho_m[cell_id] * _yr;
          exp_st[cell_id] -= e * kr * rho_m[cell_id] * _yr * nr[cell_id];
        }
      }

      /* Evaporation TODO */
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Additional right-hand side source terms for scalar equations
 *        taking into account dry and humid atmospheric variables.
 *        If 1D atmospheric radiative is used additional source terms for
 *        the thermal scalar equation to take into account the radiative forcing.
 *
 * \param[in]     f_id          field id
 * \param[in,out] st_exp        Explicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_scalar_source_term(int              f_id,
                           cs_real_t        st_exp[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  const cs_lnum_t n_cells =m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_vol = mq->cell_vol;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  cs_field_t *th_f = cs_thermal_model_field();
  const cs_field_t *fld = cs_field_by_id(f_id);

  const cs_real_t *cpro_met_p = nullptr;
  const cs_real_t *cpro_tempc = nullptr;
  const cs_real_t *cpro_rho   = CS_F_(rho)->val;

  const cs_real_t cp0 = phys_pro->cp0;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  if (at_opt->meteo_profile > 1)
    cpro_met_p = cs_field_by_name("meteo_pressure")->val;

  if (cs_field_by_name_try("real_temperature") != nullptr)
    cpro_tempc = cs_field_by_name_try("real_temperature")->val;

  /*-------------------------------------------------------------------
   * Taking into account radiative forcing for the 1D radiative module
   *    (if the 3D module is not activated)
   *-------------------------------------------------------------------*/

  if (   at_opt->radiative_model_1d > 0
      && cs_glob_rad_transfer_params->type == 0
      && cs_glob_physical_model_flag[CS_ATMOSPHERIC] > 0) {

    const cs_real_t *cvar_pottemp = th_f->val;

    // Source terms in the equation of the liquid potential temperature
    if (fld == th_f) {
      cs_real_t *ray3Di = nullptr, *ray3Dst = nullptr;
      BFT_MALLOC(ray3Di, n_cells, cs_real_t);
      BFT_MALLOC(ray3Dst, n_cells, cs_real_t);
      /* Call the 1D radiative model
       * Compute the divergence of the IR and solar radiative fluxes: */
      cs_f_atr1vf();

      /* Cressman interpolation of the 1D radiative fluxes on the 3D mesh:
       * Infra red */
      cs_measures_set_t *ms
        = cs_measures_set_by_id(at_opt->infrared_1D_profile);
      cs_cressman_interpol(ms,
                           ray3Di,
                           1);
      // Sun
      ms = cs_measures_set_by_id(at_opt->solar_1D_profile);
      cs_cressman_interpol(ms,
                           ray3Dst,
                           1);

      /* Store radiative fluxes for droplet
       *  nucleation model for humid atmosphere */
      if (cs_field_by_name_try("radiative_cooling") != nullptr) {
        cs_real_t *cpro_rad_cool = cs_field_by_name("radiative_cooling")->val;
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_rad_cool[c_id] = ray3Dst[c_id] - ray3Di[c_id];
      }

      // Explicit source term for the thermal scalar equation:
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t cp_rho = cp0*cell_vol[c_id]*cpro_rho[c_id];
       // Conversion Temperature -> Potential Temperature
        const cs_real_t pot_temp = cvar_pottemp[c_id] / (cpro_tempc[c_id] + tkelvi);
        st_exp[c_id] += cp_rho*(ray3Dst[c_id]-ray3Di[c_id])*pot_temp;
      }

      BFT_FREE(ray3Di);
      BFT_FREE(ray3Dst);
    }

  }

  /* Take into source terms for thetal,
     qw and nc due to sedimentation of drops */

  /* for humid atmo. physics only */
  if (   at_opt->sedimentation_model == 1
      && cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    cs_real_t *r3 = cs_field_by_name("droplet_eq_radius")->val;
    const cs_real_t *cpro_liqwt = cs_field_by_name("liquid_water")->val;

    // Test minimum liquid water to carry out drop sedimentation
    qliqmax = 0.0;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      qliqmax = fmax(cpro_liqwt[c_id], qliqmax);
    cs_parall_max(1, CS_REAL_TYPE, &qliqmax);

    if (!r3_is_defined) {
      const cs_real_t *cvar_ntdrp = cs_field_by_name("number_of_droplets")->val;

      // Compute the mean value: (<r^3>)**1/3
      _compute_radius_volumic_droplets(n_cells,
                                       cpro_rho,
                                       cpro_liqwt,
                                       cvar_ntdrp,
                                       r3);
      r3_is_defined = true;

      BFT_MALLOC(grad1, n_cells_ext, cs_real_3_t);
      BFT_MALLOC(grad2, n_cells_ext, cs_real_3_t);

      if (qliqmax > 1e-8)
        _compute_gradient(m,
                          mq,
                          at_opt,
                          phys_pro,
                          r3,
                          cpro_rho,
                          cpro_met_p,
                          cpro_tempc,
                          cpro_liqwt,
                          cvar_ntdrp);
    } // r3_not_defined

    // For thermal field
    if (th_f == fld) {
      if (qliqmax > 1e-8) {
        const cs_real_t clatev = phys_pro->clatev;
        const cs_real_t rair   = phys_pro->r_pg_cnst;
        const cs_real_t pref   = cs_glob_atmo_constants->ps;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t dum = 0.0, pp = 0.0;
          if (at_opt->meteo_profile == 0)
            cs_atmo_profile_std(0.0, phys_pro->p0, phys_pro->t0,
                                cell_cen[c_id][2], &pp, &dum, &dum);
          else if (at_opt->meteo_profile == 1)
            pp = cs_intprf(at_opt->met_1d_nlevels_t,
                           at_opt->met_1d_ntimes,
                           at_opt->z_temp_met,
                           at_opt->time_met,
                           at_opt->hyd_p_met,
                           cell_cen[c_id][2],
                           cs_glob_time_step->t_cur);
          else
            pp = cpro_met_p[c_id];

          st_exp[c_id] -= clatev*pow(pref/pp, (rair/phys_pro->cp0))
                          *cell_vol[c_id]*grad1[c_id][2];
        }
      }
      treated_scalars += 1;
    }

    // For ym water
    else if (fld == cs_field_by_name("ym_water")) {
      if (qliqmax > 1e-8) {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          st_exp[c_id] -= cell_vol[c_id]*grad1[c_id][2];
      }
      treated_scalars += 1;
    }

    // For Number of droplets
    else if (fld == cs_field_by_name("number_of_droplets")) {
      if (qliqmax > 1e-8) {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          st_exp[c_id] += cell_vol[c_id]*grad2[c_id][2];
      }
      treated_scalars += 1;
    }

    if (treated_scalars%3 == 0) {
      r3_is_defined = false;
      BFT_FREE(grad1);
      BFT_FREE(grad2);
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional right-hand side source terms
 *         for momentum equation in case of free inlet
 *
 * \param[in,out] exp_st        Explicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_source_term_for_inlet(cs_real_3_t        exp_st[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_vol = mq->cell_vol;
  const cs_real_3_t *cell_cen = mq->cell_cen;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *cpro_rho = CS_F_(rho)->val;
  const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

  cs_real_t *tot_vol = nullptr;
  cs_real_3_t *mom_a = nullptr;
  cs_real_3_t *mom_met_a = nullptr;

  cs_real_3_t *cpro_momst
    = (cs_real_3_t *)cs_field_by_name("momentum_source_terms")->val;

  // Bulk momentum
  int n_level = 1;
  // Variable in z
  if (at_opt->open_bcs_treatment != 1)
    n_level = fmax(at_opt->met_1d_nlevels_d, 1);

  CS_MALLOC(tot_vol, n_level, cs_real_t);
  CS_MALLOC(mom_a, n_level, cs_real_3_t);
  CS_MALLOC(mom_met_a, n_level, cs_real_3_t);

  cs_real_3_t *mom = at_opt->mom_cs;
  cs_real_3_t *mom_met = at_opt->mom_met;

  const cs_real_t uref = cs_glob_turb_ref_values->uref;

  // Save previous values
  cs_array_copy<cs_real_t>(3*n_level,
                           (const cs_real_t *)(mom),
                           (cs_real_t *)(mom_a));

  cs_array_copy<cs_real_t>(3*n_level,
                           (const cs_real_t *)(mom_met),
                           (cs_real_t *)(mom_met_a));

  for (int l_id =0; l_id < n_level; l_id++) {
    for (int ii = 0; ii < 3; ii++) {
      mom[l_id][ii] = 0.0;
      mom_met[l_id][ii] = 0.0;
    }
    tot_vol[l_id] = 0.0;
  }

  /* Computation of the target momentum bulk
     using the interpolated mean velocity field */

  cs_real_3_t *cpro_vel_target = nullptr;
  if (cs_field_by_name_try("meteo_velocity") != nullptr)
    cpro_vel_target
      = (cs_real_3_t *)cs_field_by_name_try("meteo_velocity")->val;
  else
    CS_MALLOC(cpro_vel_target, n_cells_ext, cs_real_3_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    int l_id = 0;
    cs_real_t xuent = 0.0, xvent = 0.0;
    const cs_real_t zent = cell_cen[c_id][2];

    // Get level id
    if (at_opt->met_1d_nlevels_d > 0) {
      cs_real_t dist_min = fabs(zent - at_opt->z_dyn_met[0]);
      for (int id = 1; id < n_level; id++) {
        const cs_real_t dist_ent = fabs(zent -at_opt->z_dyn_met[id]);
        if (dist_ent >= dist_min)
          continue;
        l_id = id;
        dist_min = dist_ent;
      }
    }

    if (at_opt->theo_interp == 1 || at_opt->meteo_profile > 1) {
      xuent = cpro_vel_target[c_id][0];
      xvent = cpro_vel_target[c_id][1];
    }
    else {
      xuent = cs_intprf(at_opt->met_1d_nlevels_d,
                        at_opt->met_1d_ntimes,
                        at_opt->z_dyn_met,
                        at_opt->time_met,
                        at_opt->u_met,
                        zent,
                        cs_glob_time_step->t_cur);

      xvent = cs_intprf(at_opt->met_1d_nlevels_d,
                        at_opt->met_1d_ntimes,
                        at_opt->z_dyn_met,
                        at_opt->time_met,
                        at_opt->v_met,
                        zent,
                        cs_glob_time_step->t_cur);
    }

    tot_vol[l_id] += cell_vol[c_id];
    mom_met[l_id][0] += cpro_rho[c_id] * cell_vol[c_id] * xuent;
    mom_met[l_id][1] += cpro_rho[c_id] * cell_vol[c_id] * xvent;

  }

  if (cs_field_by_name_try("meteo_velocity") == nullptr)
    CS_FREE(cpro_vel_target);

  cs_parall_sum(n_level, CS_REAL_TYPE, tot_vol);
  cs_parall_sum(3*n_level, CS_REAL_TYPE, (cs_real_t *)(mom_met));

  for (int l_id = 0; l_id < n_level; l_id++) {
    if (tot_vol[l_id] <= 0.0)
      continue;
    for (int ii = 0; ii < 3; ii++)
      mom_met[l_id][ii] = mom_met[l_id][ii]/tot_vol[l_id];
  }

  /* Computation of the momentum using the computed velocity
     ------------------------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    int l_id = 0;
    const cs_real_t zent = cell_cen[c_id][2];

    // Get level id
    if (at_opt->met_1d_nlevels_d > 0) {
      cs_real_t dist_min = fabs(zent - at_opt->z_dyn_met[0]);
      for (int id = 1; id < n_level; id++) {
        const cs_real_t dist_ent = fabs(zent - at_opt->z_dyn_met[id]);
        if (dist_ent >= dist_min)
          continue;
        l_id = id;
        dist_min = dist_ent;
      }
    }

    for (int ii = 0; ii < 3; ii++)
      mom[l_id][ii] += cpro_rho[c_id] * cell_vol[c_id] *cvar_vel[c_id][ii];
  }

  cs_parall_sum(3*n_level, CS_REAL_TYPE, (cs_real_t *)mom);

  for (int l_id = 0; l_id < n_level; l_id++) {
    if (tot_vol[l_id] <= 0.0)
      continue;
    for (int ii = 0; ii < 3; ii++)
      mom[l_id][ii] = mom[l_id][ii]/tot_vol[l_id];
  }

  CS_FREE(tot_vol);

  /* Computation of the momentum source term
     --------------------------------------- */

  // First pass, reset previous values
  if (   cs_glob_time_step->nt_cur < 2
      || cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1) {
    cs_array_copy<cs_real_t>(3*n_level,
                             (const cs_real_t *)mom,
                             (cs_real_t *)mom_a);

    cs_array_copy<cs_real_t>(3*n_level,
                             (const cs_real_t *)(mom_met),
                             (cs_real_t *)(mom_met_a));

    cs_array_real_fill_zero(n_level, at_opt->dpdt_met);
  }

  // Delta of pressure integrated over a time step for each level
  cs_real_t *dpdtx = nullptr;
  cs_real_t *dpdty = nullptr;
  CS_MALLOC(dpdtx, n_level, cs_real_t);
  CS_MALLOC(dpdty, n_level, cs_real_t);

  for (int l_id = 0; l_id < n_level; l_id++) {

    // Momentum of CS and of the target
    const cs_real_t mom_norm = cs_math_3_norm(mom[l_id]);
    const cs_real_t mom_norm_a = cs_math_3_norm(mom_a[l_id]);

    const cs_real_t mom_met_norm = cs_math_3_norm(mom_met[l_id]);
    const cs_real_t mom_met_norm_a = cs_math_3_norm(mom_met_a[l_id]);

    at_opt->dpdt_met[l_id] += 0.5*(2.0*(mom_norm - mom_met_norm)
                                   - (mom_norm_a - mom_met_norm_a));

    // target meteo directions (current and previous)
    cs_real_3_t dir_met = {0.0, 0.0, 0.0};
    cs_real_3_t dir_met_a = {0.0, 0.0, 0.0};

    if (mom_met_norm > cs_math_epzero*uref*phys_pro->ro0)
      for (int ii = 0; ii < 3; ii++)
        dir_met[ii] = mom_met[l_id][ii] / mom_met_norm;
    if (mom_met_norm_a > cs_math_epzero*uref*phys_pro->ro0)
      for (int ii = 0; ii < 3; ii++)
        dir_met_a[ii] = mom_met_a[l_id][ii] / mom_met_norm_a;

    /* Delta of pressure in the target direction
     * Steady state DP and Rotating term due to transient meteo */
    dpdtx[l_id] = at_opt->dpdt_met[l_id] * dir_met[0]
                - mom_met_norm*(dir_met[0] - dir_met_a[0]);    //FIXME use directly umet?
    dpdty[l_id] = at_opt->dpdt_met[l_id] * dir_met[1]
                - mom_met_norm*(dir_met[1] - dir_met_a[1]);

  }

  CS_FREE(mom_a);
  CS_FREE(mom_met_a);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t dpdtx_ent = dpdtx[0];
    cs_real_t dpdty_ent = dpdty[0];

    if (n_level > 1) {
      cs_intprz(at_opt->met_1d_nlevels_d,
                at_opt->z_dyn_met,
                dpdtx,
                cell_cen[c_id][2],
                nullptr,
                &dpdtx_ent);

      cs_intprz(at_opt->met_1d_nlevels_d,
                at_opt->z_dyn_met,
                dpdty,
                cell_cen[c_id][2],
                nullptr,
                &dpdty_ent);
    }

    cpro_momst[c_id][0] = - dpdtx_ent / dt[c_id];
    cpro_momst[c_id][1] = - dpdty_ent / dt[c_id];
    cpro_momst[c_id][2] = 0.0; //FIXME not ok

    for (int ii = 0; ii < 3; ii++)
      exp_st[c_id][ii] += cpro_momst[c_id][ii] * cell_vol[c_id];

  }

  CS_FREE(dpdtx);
  CS_FREE(dpdty);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS

