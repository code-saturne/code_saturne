/*============================================================================
 * Atmospheric radiative fluxes for 1D scheme
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo.h"
#include "atmo/cs_atmo_chemistry.h"
#include "atmo/cs_air_props.h"
#include "atmo/cs_atmo_profile_std.h"
#include "atmo/cs_intprf.h"
#include "base/cs_array.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_field.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_measures_util.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "base/cs_physical_properties.h"
#include "base/cs_prototypes.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "pprt/cs_physical_model.h"
#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_1d_rad.h"

static void
_dump_scalar_to_file(const char *filename, double heuray, cs_real_t val)
{
  FILE *f = fopen(filename, "a");
  if (f != nullptr) {
    fprintf(f, "%14.7e %14.7e\n", heuray, val);
    fclose(f);
  }
}

static void
_dump_array_to_file(const char *filename,
                    double      heuray,
                    const       cs_real_t *arr,
                    int         size)
{
  FILE *f = fopen(filename, "a");
  if (f != nullptr) {
    fprintf(f, "%14.7e", heuray);
    for (int i = 0; i < size; i++) {
      fprintf(f, " %14.7e", arr[i]);
    }
    fprintf(f, "\n");
    fclose(f);
  }
}

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* local atmo 1d radiative structure */
static cs_atmo_1d_rad_t _atmo_1d_rad = {
  .radiative_model_1d = 0,
  .nvert = 1,
  .nlevels = 20,
  .nlevels_max = 0,
  .frequency = 1,
  .has_aerosol = 1,
  .verbosity = 0,
  .tausup = -1.0,
  .aod_o3_tot = 0.2,
  .aod_h2o_tot = 0.1,
  .aod_ir = 0.1,
  .conco2 = 3.5e-2*44.0/29.0,
  .gaero_o3 = 0.66,
  .gaero_h2o = 0.64,
  .piaero_o3 = 0.84,
  .piaero_h2o = 0.84,
  .black_carbon_frac = 0.0,
  .zaero = 6000.0,
  .xy = nullptr,
  .z = nullptr,
  .acinfe = nullptr,
  .dacinfe = nullptr,
  .aco2 = nullptr,
  .aco2s = nullptr,
  .daco2 = nullptr,
  .daco2s = nullptr,
  .acsup = nullptr,
  .acsups = nullptr,
  .dacsup = nullptr,
  .dacsups = nullptr,
  .tauzq = nullptr,
  .tauz = nullptr,
  .zq = nullptr,
  .ir_div = nullptr,
  .sol_div = nullptr,
  .iru = nullptr,
  .ird = nullptr,
  .solu = nullptr,
  .sold = nullptr,
  .qw = nullptr,
  .ql = nullptr,
  .qv = nullptr,
  .nc = nullptr,
  .fn = nullptr,
  .aerosols = nullptr,
  .albedo0 = nullptr,
  .emissi0 = nullptr,
  .temp0   = nullptr,
  .theta0  = nullptr,
  .qw0     = nullptr,
  .p0      = nullptr,
  .rho0    = nullptr
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_1d_rad_t *cs_glob_atmo_1d_rad = &_atmo_1d_rad;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cloud fraction estimate (we take the maximum)
 *----------------------------------------------------------------------------*/

static inline void
_cf_estimate(cs_real_t        &fn,
             cs_real_t        &taul,
             const cs_real_t  eps,
             const cs_real_t  val_1,
             const cs_real_t  val_2,
             const cs_real_t  fnerir,
             const cs_real_t  qqqqql)
{
  fn = std::max(fn, fnerir);
  if (fn < eps) {
    fn = val_1;
    taul = val_2;
  }
  else
    taul = exp(-qqqqql / fn);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute reflection and transmission
 *
 * \param[in]   pioc        Albedo of simple diffusion for cloud (water)
 * \param[in]   piaero      Albedo of simple diffusion for aerosol
 * \param[in]   gasym       Asymmetry factor for clouds
 * \param[in]   gaero       Asymmetry factor for aerosols
 * \param[in]   tauc        Optical depth for clouds
 * \param[in]   taua        Optical depth for aersols
 * \param[out]  ref         Reflection
 * \param[out]  tra         Transmission
 * \param[in]   epsc        clipping threshold
 * \param[in]   dqqv        Optical depth for Water vapor
 * \param[in]   mui         mui (Gauss quadrature)
 * \param[in]   muzero_cor  mu0 corrected
 */
/*----------------------------------------------------------------------------*/

static void
_compute_reflection_transmission(const cs_real_t pioc,
                                 const cs_real_t piaero,
                                 const cs_real_t gasym,
                                 const cs_real_t gaero,
                                 const cs_real_t tauc,
                                 const cs_real_t taua,
                                 cs_real_t &ref,
                                 cs_real_t &tra,
                                 const cs_real_t epsc,
                                 const cs_real_t dqqv,
                                 const cs_real_t mui,
                                 const cs_real_t muzero_cor)
{
  cs_real_t gas, fas, kt, gama1, gama2, tln;
  cs_real_t drt, extlnp, extlnm;
  cs_real_t pic, tau;
  cs_real_t gama3, gama4, a1, a2, ktmu, expmuo, exnmuo, rnum, tnum;

  /* Optical depth computation
     -------------------------- */
  tau = tauc + taua + dqqv;
  if (tau < epsc) {
    ref = 0.0;
    tra = 1.0;
    return;
  }

  /* Global and Direct solution
     -------------------------- */
  const bool glob = true; /* Solve global and direct */
  if (glob) {

    /* Pure diffusion atmosphere (pioc = 1) */
    if (pioc >= (1.0 - epsc)) { // TODO check .and. (taua .le. epsc))
      gama1 = (1. - gasym) / (2.0 * mui);
      ref = gama1 * tau / (1.0 + gama1 * tau);
      tra = 1.0 / (1.0 + gama1 * tau);
      return;
    }

    pic = (pioc * tauc + piaero * taua) / tau;

    /* Pure absorbing atmosphere (pioc = 0) */
    if (pic < epsc) {
      gama1 = 1.0 / mui;
      ref = 0.0;
      tra = exp(-gama1 * tau);
      return;
    }

    gas = (pioc * tauc * gasym + piaero * taua * gaero) / (pic * tau);
    /* Only forward diffusion */
    if (gas >= (1.0 - epsc)) {
      gama1 = (1.0 - pic) / mui;
      ref = 0.0;
      tra = exp(-gama1 * tau);
      return;
    }

    /* Joseph correction */
    fas = gas * gas;
    tau = (1.0 - pic * fas) * tau;
    pic = pic * (1. - fas) / (1.0 - pic * fas);
    gas = (gas - fas) / (1.0 - fas);
    gama1 = (1.0 - pic * (1.0 + gas) / 2.0) / mui;
    gama2 = pic * (1.0 - gas) / (2.0 * mui);

    kt = sqrt(gama1 * gama1 - gama2 * gama2);
    tln = kt * tau;
    extlnp = exp(tln);
    extlnm = exp(-tln);

    drt = (kt + gama1) * extlnp + (kt - gama1) * extlnm;

    ref = gama2 * (extlnp - extlnm) / drt;
    tra = 2.0 * kt / drt;
    return;
  }

  /* Diffuse and Direct solution */
  if (pioc >= (1.0 - epsc)) {
    gama1 = 0.5 * (1.0 - gasym) / mui;
    gama3 = 0.5 * (1.0 - gasym);
    exnmuo = exp(-tau / muzero_cor);

    ref = (gama1 * tau + (gama3 - gama1 * muzero_cor) * (1.0 - exnmuo))
          / (1.0 + gama1 * tau);
    tra = 1.0 - ref;
    return;
  }

  pic = (pioc * tauc + piaero * taua) / tau;

  if (pic < epsc) {
    ref = 0.0;
    tra = exp(-tau / muzero_cor);
    return;
  }

  gas = (pioc * tauc * gasym + piaero * taua * gaero) / (pic * tau);

  /* Only forward diffusion */
  if (gas >= (1.0 - epsc)) {
    gama1 = (1.0 - pic) / mui;
    gama2 = 0.0;
    gama3 = 0.5 * (1.0 - gas);
    gama4 = 0.5 * (1.0 + gas);
    expmuo = exp(+tau / muzero_cor * (1.0 - gama1 * muzero_cor));
    exnmuo = exp(-tau / muzero_cor * (1.0 + gama1 * muzero_cor));

    ref = gama3 * pic * (1.0 - exnmuo) / (1.0 + gama1 * muzero_cor);
    tra = exp(-tau / muzero_cor)
          * (1.0 - pic * gama4 * expmuo / (1.0 - gama1 * muzero_cor));
    return;
 }

  /* Joseph correction (Diffuse + Direct) */
  fas = gas * gas;
  tau = (1.0 - pic * fas) * tau;
  pic = pic * (1.0 - fas) / (1.0 - pic * fas);
  gas = (gas - fas) / (1.0 - fas);

  gama1 = (1.0 - pic * (1.0 + gas) / 2.0) / mui;
  gama2 = pic * (1.0 - gas) / (2.0 * mui);
  gama3 = 0.5 * (1.0 - 3.0 * gas * muzero_cor * mui);
  gama4 = 0.5 * (1.0 + 3.0 * gas * muzero_cor * mui);

  a1 = gama1 * gama4 + gama2 * gama3;
  a2 = gama1 * gama3 + gama2 * gama4;

  kt = sqrt(gama1 * gama1 - gama2 * gama2);
  ktmu = kt * muzero_cor;
  tln = kt * tau;
  extlnp = exp(tln);
  extlnm = exp(-tln);
  expmuo = exp(tau / muzero_cor);
  exnmuo = exp(-tau / muzero_cor);

  drt = (1.0 - ktmu * ktmu) * ((kt + gama1) * extlnp + (kt - gama1) * extlnm);

  rnum = (1.0 - ktmu) * (a2 + kt * gama3) * extlnp
       - (1.0 + ktmu) * (a2 - kt * gama3) * extlnm
       - 2.0 * kt * (gama3 - a2 * muzero_cor) * exnmuo;

  ref = pic * rnum / drt;

  tnum = (1.0 + ktmu) * (a1 + kt * gama4) * extlnp
       - (1.0 - ktmu) * (a1 - kt * gama4) * extlnm
       - 2.0 * kt * (gama4 + a1 * muzero_cor) * expmuo;

  tra = exnmuo * (1.0 - pic * tnum / drt);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Computes ozone amount above a given altitude
 *
 * \param[in] zh altitude
 */
 /*----------------------------------------------------------------------------*/

static inline cs_real_t
ozone_amount(const cs_real_t zh)
{
  constexpr cs_real_t a = 0.4;
  constexpr cs_real_t b = 20000.0;
  constexpr cs_real_t c = 5000.0;

  // precompute constant
  const cs_real_t k = 1.0 + exp(-b / c);

  return a * k / (1.0 + exp((zh - b) / c));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Computes derivative of ozone amount (ozone gradient) w.r.t altitude
 *
 * \param[in] zh altitude
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
ozone_gradient(const cs_real_t zh)
{
  constexpr cs_real_t a = 0.4;
  constexpr cs_real_t b = 20000.0;
  constexpr cs_real_t c = 5000.0;

  const cs_real_t k = 1.0 + exp(-b / c);
  const cs_real_t exp_term = exp((zh - b) / c);
  const cs_real_t denom = (1.0 + exp_term) * (1.0 + exp_term);

  return a / c * k * exp_term / denom;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute absorption of solar radiation by water vapor
 *
 * \param[in]  y   optical depth (or specific humidity)
 *
 * \return absorption coefficient
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
ray_sve(const cs_real_t y)
{
  const cs_real_t denom = pow(1.0 + 14.15 * y, 0.635) + 0.5925 * y;
  return 0.29 * y / denom;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Absorption derivative-function of the solar radiation by water vapor
 *
 * \param[in]  y   optical depth
 * \param[in]  dy  todo ?
 *
 * \return derivative of absorption coefficient
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
ray_sve_derivative(const cs_real_t y,
                   const cs_real_t dy)
{
  constexpr cs_real_t a = 0.29;
  constexpr cs_real_t b = 14.15;
  constexpr cs_real_t c = 0.635;
  constexpr cs_real_t d = 0.5925;

  const cs_real_t denom = pow(1.0 + b * y, c) + d * y;
  const cs_real_t num   = b * c * dy * pow(1.0 + b * y, c - 1.0) + d * dy;

  return a * dy / denom - a * y * num / (denom * denom);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Absorption function of the solar radiation by ozone
 *
 * \param[in]  x   optical depth for ozone
 *
 * \return absorption coefficient
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
ray_ozone_absorption(const cs_real_t x)
{
  constexpr cs_real_t a = 0.02118;
  constexpr cs_real_t b = 0.042;
  constexpr cs_real_t c = 0.000323;
  constexpr cs_real_t d = 1.082;
  constexpr cs_real_t e = 138.6;
  constexpr cs_real_t f = 0.805;
  constexpr cs_real_t g = 0.0658;
  constexpr cs_real_t h = 103.6;

  const cs_real_t ao3vis = a * x / (1.0 + (b + c * x) * x);

  const cs_real_t ao3uv  = d * x / pow(1.0 + e * x, f)
                         + g * x / (1.0 + pow(h * x, 3.0));

  return ao3vis + ao3uv;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Derivative of the absorption function of the solar radiation by ozone
 *
 * \param[in]  x   optical depth for ozone
 * \param[in]  dx  small increment (derivative multiplier)
 *
 * \return derivative of absorption coefficient
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
d_ray_ozone_absorption(const cs_real_t x,
                       const cs_real_t dx)
{
  constexpr cs_real_t a = 0.02118;
  constexpr cs_real_t b = 0.042;
  constexpr cs_real_t c = 0.000323;
  constexpr cs_real_t d = 1.082;
  constexpr cs_real_t e = 138.6;
  constexpr cs_real_t f = 0.805;
  constexpr cs_real_t g = 0.0658;
  constexpr cs_real_t h = 103.6;

  const cs_real_t den1 = 1.0 + b * x + c * x * x;
  const cs_real_t num1 = a * (1.0 - c * x * x);
  const cs_real_t term1 = num1 / (den1 * den1);

  const cs_real_t hx3 = h * x;
  const cs_real_t hx3_cub = hx3 * hx3 * hx3;
  const cs_real_t den2 = 1.0 + hx3_cub;
  const cs_real_t num2 = g * (1.0 - 2.0 * hx3_cub);
  const cs_real_t term2 = num2 / (den2 * den2);

  const cs_real_t ex = 1.0 + e * x;
  const cs_real_t exf = pow(ex, f);
  const cs_real_t num3 = d * (exf - x * f * e * pow(ex, f - 1.0));
  const cs_real_t den3 = pow(ex, 2.0 * f);
  const cs_real_t term3 = num3 / den3;
  return dx * (term1 + term2 + term3);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Store the incident radiation of the 1D model
 *
 * \param[in]  idx          spectral index
 * \param[in]  dim          dimension
 * \param[in]  nprofz       total number of measure points
 * \param[in]  bc_type      boundaries types
 * \param[in]  n_b_faces    number of boundary faces
 * \param[in]  muzero       mu0
 * \param[in]  muzero_cor   mu0 corrected
 * \param[in]  profz        z coordinates of measure points
 * \param[in]  profv        measured values
 * \param[in]  b_face_cog   center of gravity of border faces
 * \param[out] bpro_rad_inc spectral radiative incident flux
 *
 */
/*----------------------------------------------------------------------------*/

static void
_interpolate_boundary_rad_incident_flux(const int          idx,
                                        const int          dim,
                                        const int          nprofz,
                                        const int          bc_type[],
                                        const cs_lnum_t    n_b_faces,
                                        const cs_real_t    muzero,
                                        const cs_real_t    muzero_cor,
                                        const cs_real_t    profz[],
                                        const cs_real_t    profv[],
                                        const cs_real_3_t  *b_face_cog,
                                        cs_real_t          *bpro_rad_inc)
{
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    const int id = idx + face_id * dim;
    bpro_rad_inc[id] = 0.0;
    if (   (muzero > cs_math_epzero)
           && (   bc_type[face_id] != CS_SMOOTHWALL
               && bc_type[face_id] != CS_ROUGHWALL)) {
      cs_real_t var;
      cs_intprz(nprofz,
                profz,
                profv,
                b_face_cog[face_id][2],
                nullptr,
                &var);
      bpro_rad_inc[id] = var / muzero_cor;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Store the absorption coefficient of the 1D model
 *
 * \param[in]  idx          spectral index
 * \param[in]  dim          dimension
 * \param[in]  nprofz       total number of measure points
 * \param[in]  n_cells      number of cells
 * \param[in]  profz        z coordinates of measure points
 * \param[in]  profv        measured values
 * \param[in]  cell_cen     cell center coordinates
 * \param[out] cpro_coef    absorption coefficient
 *
 */
/*----------------------------------------------------------------------------*/

static void
_interpolate_coeff(const int          idx,
                   const int          dim,
                   const int          nprofz,
                   const cs_lnum_t    n_cells,
                   const cs_real_t    profz[],
                   const cs_real_t    profv[],
                   const cs_real_3_t  *cell_cen,
                   cs_real_t          *cpro_coef)
{
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t var;
    cs_intprz(nprofz,
              profz,
              profv,
              cell_cen[c_id][2],
              nullptr,
              &var);
    cpro_coef[idx + c_id * dim] = var;
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  1D Radiative scheme - IR H20 and dimere absorption
 *
 * \param[out]      tauv      transmission function for water vapor and dimers
 * \param[out]      dtauv     d(tauv)/dz
 * \param[in]       qqqq      optical depth for water vapor (z,z')
 * \param[in]       xqx       effective concentration absorption by water vapor
 * \param[in]       qqqqc     idem qqqq for dimers
 * \param[in]       xqc       idem xqx for dimers
 * \param[in]       ro        air density
 *
 */
/*----------------------------------------------------------------------------*/

static void
_compute_ir_h2o_dimer_transmission(cs_real_t        &tauv,
                                   cs_real_t        &dtauv,
                                   const cs_real_t  qqqq,
                                   const cs_real_t  xqx,
                                   const cs_real_t  qqqqc,
                                   const cs_real_t  xqc,
                                   const cs_real_t  ro)
{

  static constexpr cs_real_t A[5]
    = {7.76192e-7, 1.33836e-3, 1.66649e-1, 2.17686, 2.69020};
  static constexpr cs_real_t B[5]
    = {7.79097e-7, 1.36832e-3, 1.79601e-1, 2.70573, 5.15119};
  static constexpr cs_real_t AS[3] = {1.5075e-2, -3.6185e-2, 1.9245e-2};
  static constexpr cs_real_t BS[3] = {1.5075e-2,  1.9547e-1, 7.5271e-1};

  const cs_real_t u   = qqqq  / 10.0;
  const cs_real_t us  = qqqqc / 10.0;
  const cs_real_t xu  = xqx   / 10.0;
  const cs_real_t xus = xqc   / 10.0;

  cs_real_t za, dza;

  // Emissivity of water vapor
  if (u >= 0.01) {
    za  = 0.24 * log10(u + 0.01) + 0.622;
    dza = 0.24 / (u + 0.01) / log(10.0);
  } else {
    const cs_real_t u_shift = u + 3.59e-5;
    za  = 0.846 * pow(u_shift, 0.243) - 6.9e-2;
    dza = 0.846 * 0.243 * pow(u_shift, 0.243 - 1.0);
  }

  // Emissivity of dimers

  const cs_real_t n1 = A[0] + u * (A[1] + u * (A[2] + u * (A[3] + u * A[4])));
  const cs_real_t dn1
    = A[1] + u * (2.0 * A[2] + u * (3.0 * A[3] + u * 4.0 * A[4]));

  const cs_real_t d1
    = B[0] + u * (B[1] + u * (B[2] + u * (B[3] + u * (B[4] + u))));
  const cs_real_t dd1
    = B[1] + u * (2.0 * B[2] + u * (3.0 * B[3] + u * (4.0 * B[4] + u * 5.0)));

  const cs_real_t t1  = n1 / d1;
  const cs_real_t dt1 = dn1 / d1 - n1 * dd1 / d1 / d1;

  cs_real_t t2  = 0.0;
  cs_real_t dt2 = 0.0;

  if (us <= 0.5) {
    const cs_real_t n2  = AS[0] + us * (AS[1] + us * AS[2]);
    const cs_real_t dn2 = AS[1] + us * 2.0 * AS[2];
    const cs_real_t d2  = BS[0] + us * (BS[1] + us * (BS[2] + us));
    const cs_real_t dd2 = BS[1] + us * (2.0 * BS[2] + 3.0 * us);

    t2  = n2 / d2;
    dt2 = dn2 / d2 - n2 * dd2 / (d2 * d2);
  }

  // -------- Total transmission --------
  const cs_real_t abs0  = za + 0.4614 * t1 * (1.0 - t2);
  const cs_real_t dabs0
    = (dza * xu + 0.4614 * (dt1 * xu * (1.0 - t2) - t1 * dt2 * xus)) * ro;

  dtauv = -dabs0;
  tauv  = 1.0 - abs0;
}

/*----------------------------------------------------------------------------*/
/*
 * brief Compute carbon dioxide (CO2) and ozone (O3) absorption in the infrared
 *        (1D radiative scheme).
 *
 * This routine evaluates the total infrared absorption due to CO2 and O3
 * between two altitudes zz and zzp, including their differential contributions.
 *
 * param[in]  zz     Height above ground level (m)
 * param[in]  zzp    Intermediate altitude for ozone (m)
 * param[out] xa     Total IR absorption (CO2 + O3)
 * param[out] xda    Vertical derivative of absorption (CO2 + O3)
 * param[in]  u      Water vapor optical depth between (zz, zzp)
 * param[in]  q      Effective concentration for H2O absorption
 * param[in]  uco2   CO2 optical depth between (zz, zzp)
 * param[in]  qco2   Effective concentration for CO2 absorption
 * param[in]  ro     Air density (kg/m³)
 *
 */
/*----------------------------------------------------------------------------*/

static void
_compute_ir_co2_o3_absorption(const cs_real_t zz,
                              const cs_real_t zzp,
                              cs_real_t &xa,
                              cs_real_t &xda,
                              const cs_real_t u,
                              const cs_real_t q,
                              const cs_real_t uco2,
                              const cs_real_t qco2,
                              const cs_real_t ro)
{

  constexpr cs_real_t y[2]  = {0.0581, 0.0546};
  constexpr cs_real_t yo[2] = {0.0749, 0.0212};
  constexpr cs_real_t x[4]  = {1.33, -0.4572, 0.26, 0.286};
  constexpr cs_real_t xo[4] = {0.209, 7e-5, 0.436, -0.00321};
  constexpr cs_real_t xc[4] = {-0.00982, 0.0676, 0.421, 0.01022};

  cs_real_t dtauv, tauv;

  /* 1. H2O absorption into 15 µm band */
  if (u <= 20.0) {
    tauv  = x[0] + x[1] * pow(u + x[3], x[2]);
    dtauv = ro * q * x[1] * x[2] * pow(u + x[3], x[2] - 1.0);
  }
  else {
    dtauv = -0.2754 / log(10.0) * ro * q / u;
    tauv  = 0.33 - 0.2754 * (log10(u) - 1.3011);
  }

  /* 2. CO2 optical depth */
  cs_real_t ao, dao, duco2 = ro * qco2;
  if (uco2 <= 1.0) {
    ao  = xc[0] + xc[1] * pow(uco2 + xc[3], xc[2]);
    dao = duco2 * xc[1] * xc[2] * pow(uco2 + xc[3], xc[2] - 1.0);
  }
  else {
    ao  = y[0] + y[1] * log10(uco2);
    dao = y[1] / log(10.0) * duco2 / uco2;
  }

  /* 3. O3 optical depth */
  cs_real_t ao3, dao3;
  cs_real_t uo3  = abs(ozone_amount(zz) - ozone_amount(zzp));
  cs_real_t duo3 = -ozone_gradient(zz);
  if (uo3 <= 0.01) {
    ao3  = xo[0] * pow(uo3 + xo[1], xo[2]) + xo[3];
    dao3 = duo3 * xo[0] * xo[2] * pow(uo3 + xo[1], xo[2] - 1.0);
  }
  else {
    ao3  = yo[0] + yo[1] * log10(uo3);
    dao3 = yo[1] * duo3 / (log(10.0) * uo3);
  }

  /* 4. Total absorption */
  xa  = tauv * ao  + ao3;
  xda = tauv * dao + dtauv * ao + dao3;

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_initialize(void)
{
  int n_vert = 0;
  int n_level = 0;

  if (_atmo_1d_rad.radiative_model_1d != 0) {
    n_level = std::max(1, _atmo_1d_rad.nlevels_max);
    n_vert = std::max(1, _atmo_1d_rad.nvert);
  }
  int n_lp1 = n_level+1;

  if (_atmo_1d_rad.xy == nullptr)
    CS_MALLOC(_atmo_1d_rad.xy , 3*n_vert, cs_real_t);
  if (_atmo_1d_rad.z == nullptr)
    CS_MALLOC(_atmo_1d_rad.z , n_level, cs_real_t);
  if (_atmo_1d_rad.acinfe == nullptr)
    CS_MALLOC(_atmo_1d_rad.acinfe , n_level, cs_real_t);
  if (_atmo_1d_rad.dacinfe == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacinfe , n_level, cs_real_t);
  if (_atmo_1d_rad.aco2 == nullptr)
    CS_MALLOC(_atmo_1d_rad.aco2, n_level*n_level, cs_real_t);
  if (_atmo_1d_rad.aco2s == nullptr)
    CS_MALLOC(_atmo_1d_rad.aco2s, n_level*n_level, cs_real_t);
  if (_atmo_1d_rad.daco2 == nullptr)
    CS_MALLOC(_atmo_1d_rad.daco2, n_level*n_level, cs_real_t);
  if (_atmo_1d_rad.daco2s == nullptr)
    CS_MALLOC(_atmo_1d_rad.daco2s, n_level*n_level, cs_real_t);
  if (_atmo_1d_rad.acsup == nullptr)
    CS_MALLOC(_atmo_1d_rad.acsup, n_level, cs_real_t);
  if (_atmo_1d_rad.acsups == nullptr)
    CS_MALLOC(_atmo_1d_rad.acsups, n_level, cs_real_t);
  if (_atmo_1d_rad.dacsup == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacsup, n_level, cs_real_t);
  if (_atmo_1d_rad.dacsups == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacsups, n_level, cs_real_t);
  if (_atmo_1d_rad.tauzq == nullptr)
    CS_MALLOC(_atmo_1d_rad.tauzq, n_level+1, cs_real_t);
  if (_atmo_1d_rad.tauz == nullptr)
    CS_MALLOC(_atmo_1d_rad.tauz, n_level+1, cs_real_t);
  if (_atmo_1d_rad.zq == nullptr)
    CS_MALLOC(_atmo_1d_rad.zq, n_level+1, cs_real_t);
  if (_atmo_1d_rad.ir_div == nullptr)
    CS_MALLOC(_atmo_1d_rad.ir_div, n_level * n_vert, cs_real_t);
  if (_atmo_1d_rad.sol_div == nullptr)
    CS_MALLOC(_atmo_1d_rad.sol_div, n_level * n_vert, cs_real_t);
  if (_atmo_1d_rad.iru == nullptr)
    CS_MALLOC(_atmo_1d_rad.iru, n_lp1 * n_vert, cs_real_t);
  if (_atmo_1d_rad.ird == nullptr) {
    CS_MALLOC(_atmo_1d_rad.ird, n_lp1 * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_lp1 * n_vert, _atmo_1d_rad.ird);
  }
  if (_atmo_1d_rad.solu == nullptr)
    CS_MALLOC(_atmo_1d_rad.solu, n_lp1 * n_vert, cs_real_t);
  if (_atmo_1d_rad.sold == nullptr) {
    CS_MALLOC(_atmo_1d_rad.sold, n_lp1 * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_lp1 * n_vert, _atmo_1d_rad.sold);
  }
  if (_atmo_1d_rad.qw == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qw, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.qw);
  }
  if (_atmo_1d_rad.ql == nullptr) {
    CS_MALLOC(_atmo_1d_rad.ql, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.ql);
  }
  if (_atmo_1d_rad.qv == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qv, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.qv);
  }
  if (_atmo_1d_rad.nc == nullptr) {
    CS_MALLOC(_atmo_1d_rad.nc, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.nc);
  }
  if (_atmo_1d_rad.fn == nullptr) {
    CS_MALLOC(_atmo_1d_rad.fn, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.fn);
  }
  if (_atmo_1d_rad.aerosols == nullptr) {
    CS_MALLOC(_atmo_1d_rad.aerosols, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.aerosols);
  }
  if (_atmo_1d_rad.albedo0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.albedo0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.albedo0);
  }
  if (_atmo_1d_rad.emissi0== nullptr) {
    CS_MALLOC(_atmo_1d_rad.emissi0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.emissi0);
  }
  if (_atmo_1d_rad.temp0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.temp0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.temp0);
  }
  if (_atmo_1d_rad.theta0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.theta0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.theta0);
  }
  if (_atmo_1d_rad.qw0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qw0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.qw0);
  }
  if (_atmo_1d_rad.p0  == nullptr) {
    CS_MALLOC(_atmo_1d_rad.p0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.p0);
  }
  if (_atmo_1d_rad.rho0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.rho0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.rho0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute solar fluxes for both clear and cloudy atmosphere following
 * Lacis and Hansen (1974). The multiple diffusion is taken into account by an
 * addition method and overlapping between water vapor and liquid water with k
 * distribution method.
 * Some improvements from original version concerns:
 * - introduction of cloud fraction with hazardous recovering
 * - introduction of aerosol diffusion in the same way as for cloud droplets
 *   but with specific optical properties for aerosols.
 *
 * \param[in]   ivertc      index of vertical profile
 * \param[in]   k1          index for ground level
 * \param[in]   kmray       vertical levels number for radiation
 * \param[in]   heuray      Universal time (Hour)
 * \param[in]   imer1       sea index
 * \param[in]   albe        albedo
 * \param[in]   qqv         optical depth for water vapor (0,zqq)
 * \param[in]   qqvinf      idem qqv but for altitude above 11000m
 * \param[in]   zqq         vertical levels (interfaces)
 * \param[in]   zray        vertical levels (volumes)
 * \param[in]   qvray       specific humidity for water vapor
 * \param[in]   qlray       specific humidity for liquid water
 * \param[in]   fneray      cloud fraction
 * \param[in]   romray      air density
 * \param[in]   preray      pressure
 * \param[in]   temray      temperature
 * \param[out]  rayst       flux divergence of solar radiation
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_compute_solar(const int       ivertc,
                             const int       k1,
                             const int       kmray,
                             const cs_real_t heuray,
                             const int       imer1,
                             cs_real_t       *albe_p,
                             cs_real_t       qqv[],
                             const cs_real_t qqvinf,
                             const cs_real_t zqq[],
                             const cs_real_t zray[],
                             const cs_real_t qvray[],
                             cs_real_t       qlray[],
                             cs_real_t       fneray[],
                             const cs_real_t romray[],
                             const cs_real_t preray[],
                             const cs_real_t temray[],
                             cs_real_t       rayst[],
                             const cs_real_t ncray[])
{
  /* Minimal local class to access some 2D spans with shifted
     arguments, replacing ref(i, j) with ref(i+1, j).

     This is used to use a consistent syntax, with an array that
     would require the first indes to go to -1. */

  class mdspan_1_0_shifted_accessor {
  public:

    // Simple constructor
    CS_F_HOST_DEVICE
    mdspan_1_0_shifted_accessor(int stride, cs_real_t *data):
      _stride{stride},
      _data{data}
    {};

    // Overloaded () operator to access the (i+1,j)-th value tuple.
    CS_F_HOST_DEVICE
    inline cs_real_t &
    operator()
      (
       int  i,
       int  j
       ) const
    {
      return _data[(i+1) * _stride + j];
    }

  private:
    int         _stride;
    cs_real_t  *_data;
  };

  /* Regular variables */

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_3_t *b_face_cog = mq->b_face_cog;

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  const cs_atmo_1d_rad_t *at_1d_rad = cs_glob_atmo_1d_rad;

  const int kmx = at_1d_rad->nlevels_max;
  const int kmxp1 = kmx + 1;
  const int k0 = k1 - 1;
  const int nbmett = at_opt->met_1d_nlevels_t;

  cs_real_t albe = *albe_p;
  const cs_real_t cp0  = cs_glob_fluid_properties->cp0;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t cpvcpa = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;

  // data for pkn and kn distribution
  static constexpr cs_real_t kn[8] = {
    4.e-5, 0.002, 0.035, 0.377, 1.95, 9.40, 44.6, 190.0
  };

  static const cs_real_t pkn[8] = {
    0.647, 0.0698, 0.1443, 0.0584, 0.0335, 0.0225, 0.0158, 0.0087
  };

  // Data from Chuang 2002 for calculation of
  // cloud SSA taking into account black carbon
  static constexpr cs_real_t beta1[12] = {
    0.2382803, 0.2400113, 0.2471480, 0.2489583, 0.2542476, 0.2588392,
    0.2659081, 0.2700860, 0.2783093, 0.2814346, 0.2822860, 0.1797007
  };

  static constexpr cs_real_t beta2[12] = {
    0.2940957, 0.2936845, 0.2880274, 0.2871209, 0.2824498, 0.2775943,
    0.2698008, 0.2652960, 0.2564840, 0.2535739, 0.2487382, 0.1464709
  };

  static constexpr cs_real_t beta3[12] = {
    61.83657, 58.25082, 52.79042, 50.06907, 45.75322, 42.43440,
    37.03823, 32.32349, 25.99426, 20.05043, 12.76966, 3.843661
  };

  static constexpr cs_real_t beta4[12] = {
    574.2111, 565.0809, 519.0711, 494.0088, 448.3519, 409.9063,
    348.9051, 297.9909, 233.7397, 175.4385, 112.8208, 39.24047
  };

  static constexpr cs_real_t omega0[12] = {
    5.189239e-5, 2.261712e-5, 1.264190e-5, 9.446845e-6, 6.090293e-6,
    3.794524e-6, 1.735499e-6, 1.136807e-6, 2.261422e-6, 1.858815e-5,
    5.551822e-7, 2.325124e-5
  };

  static constexpr cs_real_t coeff_E_o3[12] = {
    0.0, 0.0, 0.0, 0.0,
    0.029193795335, 0.045219606416, 0.16880411522,
    0.186215099078, 0.5705673839,
    0.0, 0.0, 0.0
  };

  static constexpr cs_real_t coeff_E_h2o[12] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.27068650951091927, 0.6844296138881737, 0.044883876600907306
  };

  // Transmission function data for minor gases (Tmg)
  static constexpr cs_real_t amg[5] = {
    0.0721, 0.0062, 0.0326, 0.0192, 0.0003
  };

  static constexpr cs_real_t bmg[5] = {
    377.890, 243.670, 107.413, 166.095, 476.934
  };

  static constexpr cs_real_t cmg[5] = {
    0.5855, 0.4246, 0.5501, 0.4221, 0.4892
  };

  static constexpr cs_real_t dmg[5] = {
    3.1709, 1.7222, 0.9093, 0.7186, 0.1261
  };

  static constexpr cs_real_t umg[5] = {
    390.0, 0.075, 0.28, 1.6, 209500.0
  };

  // ref_(0,:) and refb_(0,!) : albedo
  cs_array_2d<cs_real_t> ref_(kmx+2, 8);
  cs_array_2d<cs_real_t> refb_(kmx+2, 8);

  mdspan_1_0_shifted_accessor ref(8, ref_.data());
  mdspan_1_0_shifted_accessor refb(8, refb_.data());

  cs_array_2d<cs_real_t> tau(kmx+1, 8);
  cs_array_2d<cs_real_t> pic(kmx+1, 8);
  cs_array_2d<cs_real_t> reft(kmx+1, 8);
  cs_array_2d<cs_real_t> upwf(kmx+1, 8);
  cs_array_2d<cs_real_t> fabso3c(kmx+1, 2);
  cs_array_2d<cs_real_t> tra(kmx+1, 8);
  cs_array_2d<cs_real_t> trad(kmx+1, 8);
  cs_array_2d<cs_real_t> trat(kmx+1, 8);
  cs_array_2d<cs_real_t> trard(kmx+1, 8);
  cs_array_2d<cs_real_t> dow(kmx+1, 8);
  cs_array_2d<cs_real_t> dowf(kmx+1, 8);
  cs_array_2d<cs_real_t> dowd(kmx+1, 8);
  cs_array_2d<cs_real_t> atln(kmx+1, 8);
  cs_array_2d<cs_real_t> absn(kmx+1, 8);
  cs_array_2d<cs_real_t> ufso3c(kmx+1, 2);

  cs_array<cs_real_t> w0_sir(kmx);
  cs_array<cs_real_t> w0_suv(kmx);
  cs_array<cs_real_t> g_apc_sir(kmx);
  cs_array<cs_real_t> g_apc_suv(kmx);
  cs_array<cs_real_t> ckup_suv_f(kmx);
  cs_array<cs_real_t> ckup_sir_f(kmx);
  cs_array<cs_real_t> ckdown_sir_r(kmx);
  cs_array<cs_real_t> ckdown_sir_f(kmx);
  cs_array<cs_real_t> ckdown_suv_r(kmx);
  cs_array<cs_real_t> ckdown_suv_f(kmx);

  cs_array<cs_real_t> dfs(kmx + 1);
  cs_array<cs_real_t> ufs(kmx + 1);
  cs_array<cs_real_t> drfs(kmx + 1);
  cs_array<cs_real_t> ddfs(kmx + 1);
  cs_array<cs_real_t> tauc(kmx + 1);
  cs_array<cs_real_t> dfso3(kmx + 1);
  cs_array<cs_real_t> ufso3(kmx + 1);
  cs_array<cs_real_t> fneba(kmx + 1);
  cs_array<cs_real_t> fabso3(kmx + 1);
  cs_array<cs_real_t> dffso3(kmx + 1);
  cs_array<cs_real_t> dfsh2o(kmx + 1);
  cs_array<cs_real_t> ufsh2o(kmx + 1);
  cs_array<cs_real_t> ddfso3(kmx + 1);
  cs_array<cs_real_t> fabsh2o(kmx + 1);
  cs_array<cs_real_t> dddfso3(kmx + 1);
  cs_array<cs_real_t> dffsh2o(kmx + 1);
  cs_array<cs_real_t> fnebmax(kmx + 1);
  cs_array<cs_real_t> ddfsh2o(kmx + 1);
  cs_array<cs_real_t> dddfsh2o(kmx + 1);

  cs_array<cs_real_t> gco3(kmx + 1);
  cs_array<cs_real_t> gch2o(kmx + 1);
  cs_array<cs_real_t> pic_o3(kmx + 1);
  cs_array<cs_real_t> tauao3(kmx + 1);
  cs_array<cs_real_t> tauah2o(kmx + 1);
  cs_array<cs_real_t> pic_h2o(kmx + 1);

  /* local initializations
     ===================== */

  const bool has_aerosol = (at_1d_rad->has_aerosol > 0);
  int ibase = 0;

  constexpr cs_real_t epsc = 1.e-8;
  constexpr cs_real_t z_ref = 0.647;

  for (int k = 0; k <= kmray; k++) {
    w0_sir[k] = 0.;
    w0_suv[k] = 0.;
    g_apc_sir[k] = 0.;
    g_apc_suv[k] = 0.;

    /* TODO: useful ? */
    if (qlray[k] < epsc) {
      qlray[k] = 0.;
      fneray[k] = 0.;
    }
  }

  for (int k = 0; k < kmxp1; k++) {
    dfs[k] = 0.;
    ddfs[k] = 0.;
    drfs[k] = 0.;
    tauc[k] = 0.;
    dfso3[k] = 0.;
    ufso3[k] = 0.;
    dfsh2o[k] = 0.;
    ufsh2o[k] = 0.;
    dffso3[k] = 0.;
    fabso3[k] = 0.;
    ddfso3[k] = 0.;
    dddfso3[k] = 0.;
    fnebmax[k] = 0.;
    fabsh2o[k] = 0.;
    dffsh2o[k] = 0.;
    ddfsh2o[k] = 0.;
    dddfsh2o[k] = 0.;
    if (has_aerosol)
      fneba[k] = 1.;
    for (int n = 0; n < 2; n++)
      fabso3c(k, n) = 0.;
    for (int n = 0; n < 8; n++) {
      dow(k, n)   = 0.;
      tau(k, n)   = 0.;
      pic(k, n)   = 0.;
      reft(k, n)  = 0.;
      atln(k, n)  = 0.;
      absn(k, n)  = 0.;
      trat(k, n)  = 1.;
      upwf(k, n)  = 0.;
      dowf(k, n)  = 0.;
      trad(k, n)  = 1.;
      dowd(k, n)  = 1.;
      tra(k, n)   = 1.;
      ref(k, n)   = 0.;
      refb(k, n)  = 0.;
    }
    gco3[k] = 0.;
    gch2o[k] = 0.;
    pic_o3[k] = 0.;
    tauao3[k] = 0.;
    tauah2o[k] = 0.;
    pic_h2o[k] = 0.;
  }

  for (int n = 0; n < 8; n++) {
    ref(k0, n)   = 0.;
    refb(k0, n)  = 0.;
  }

  cs_real_t refx = 0.;
  cs_real_t trax = 0.;
  cs_real_t refx0 = 0.;
  cs_real_t trax0 = 0.;

  /* Leighton 1980 - aerosol layer configuration */
  int iaero_top = 0;     /* id of the top of the aerosol layer */

  /* Initialization for Chuang (2002) black carbon in droplets */
  constexpr cs_real_t dm0 = 20.0;    /* droplet diameter (micrometers) */
  constexpr cs_real_t nu0 = 1.0e-8;  /* volume fraction of black carbon */

  /* Initialization of SSA-related data (Chuang 2002) */
  cs_real_t piocv[12] = {0.0}, copioc[12] = {0.0}, copioc20[12] = {0.0};

  /*  2 - calculation for muzero and solar constant fo
   * -------------------------------------------------
   *        muzero = cosin of zenithal angle
   *  (corrected to take curvature of the Earth)
   *        fo = solar constant in watt/m2
   *
   *
   *  careful : 0. h < heuray < 24. h
   *  --------- */

  cs_real_t omega, fo, za, muzero_cor;
  cs_atmo_compute_solar_angles(at_opt->latitude,
                               at_opt->longitude,
                               (cs_real_t)at_opt->squant,
                               heuray,
                               imer1,
                               &albe,
                               &za,
                               &muzero_cor,
                               &omega,
                               &fo);

  // if muzero is negative, it is night and solar radiation is not
  // computed
  cs_real_t qqvtot, deltaz;
  cs_real_t muzero = cos(za);

  if (muzero > cs_math_epzero) {
    /* Optical air mass
       Corrected value for very low angle
       cf. Kasten, F., Young, A.T., 1989.
       Revised optical air mass tables and approximation formula. */
    const cs_real_t m = 1. / muzero_cor;

    // as we suppressed m and mbar in dzx and dzy we can
    // eventually keep 1.9 for O3 to tested
    constexpr cs_real_t mbar = 1.9;
    // For H2O mbar is 5/3
    constexpr cs_real_t mbarh2o = 5.0 / 3.0;

    // for LH74 quadrature method in two-stream
    // sqrt is not available at compile time.
    constexpr cs_real_t mui = 0.57735026918962576451; // 1.0 / sqrt(3.0);

    //  3 -  albedos for O3 and Rayleigh diffusion
    // Note LH74 equation (18)
    constexpr cs_real_t rabar2 = 0.1440;
    const cs_real_t rabar = 0.2190 / (1.0 + 0.8160 * muzero);
    cs_real_t rrbar2s = 0.06850;
    const cs_real_t rrbar = 0.280/(1.0 + 6.430*muzero);
    // Note LH74: eq (15);
    cs_real_t rbar
      = rabar + (1.0 - rabar)*(1.0 - rabar2)*albe/(1.0 - rabar2*albe);

    //  4 - addition of one level for solar radiation
    qqvtot = qqvinf + qqv[kmray];
    qqv[kmray + 1] = qqvtot - qqvinf;

    // Transmission for minor gases
    cs_real_t tmg = 1.;
    for (int i = 0; i < 5; i++) {
      const cs_real_t num = amg[i] * m * umg[i];
      const cs_real_t den =   pow(1.0 + bmg[i] * m * umg[i], cmg[i])
                            + dmg[i] * m * umg[i];
     tmg = tmg * (1. - num / den);
    }

    //introduction of absorption by minor gases
    fo = fo * tmg;

    /* 5 - Solar radiation calculation for cloudy sky
       In order to take into account cloud fraction,
       multiple diffusion is achieved
       for both cloudy (index 1) and clear (index 2) sky

       5.1 cloud level determination (top for the top of the higher cloud,
       base for the bottom of the lower cloud) */

    int itop  = -1;
    for (int i = kmray; i >= k1; i--) {
      if (qlray[i] > epsc) {
        // FIXME: to be coherent with 3D version
        // For now we do not force fneray[i] = 1.0 (fractional cloud kept)
        if (itop == -1)
          itop = i + 1;
        ibase = i;
      }
    }

    // If no cloud found: assign to base level
    if (itop == -1) {
      itop  = k1;
      ibase = k1;
    }

    /* Calculation for optical parameters of clouds and aerosols
       (single scattering albedo, optical depth, radius, asymmetry factor) */

    fnebmax[kmray + 1] = 0.;
    cs_real_t tauctot = 0.;
    const cs_real_t pi = cs_math_pi;
    const cs_real_t sigc = at_opt->sigc;

    for (int i = kmray; i >= k1; i--) {
      cs_real_t req = 0.;
      deltaz = 0.;

      if (i >= ibase && i < itop) {
         // Liquid water density [g/m³]
        const cs_real_t wh2ol = 1.0e3 * romray[i] * qlray[i];
        // Mean droplet radius [µm]
        cs_real_t rm = 30.0 * wh2ol + 2.0;
        // Limit mean radius by level
        rm = (i <= nbmett) ? cs::min(10.0, rm) : cs::min(2.0, rm);

        // Effective radius
        if (ncray[i] > epsc && qlray[i] > epsc) {
          req = 1.0e3 * pow(  (3.0 * romray[i] * qlray[i])
                              / (4.0 * pi * ncray[i]), 1.0 / 3.0)
            * exp(sigc * sigc);
        }
        else {
          req = 1.5 * rm;
        }

        // Clipping climatological limits
        req = (req < 1.0) ? 1.0 : (req > 20.0 ? 20.0 : req);
        // Cloud optical thickness
        deltaz   = zqq[i + 1] - zqq[i];
        tauc[i]  = 1.5 * wh2ol * deltaz / req;
        tauctot += tauc[i];
      }
      else {
        req     = 0.;
        tauc[i] = 0.;
      }

      // Aerosol optical depth AOD
      fneba[i] = 0.;

      if (has_aerosol && zray[i] <= at_1d_rad->zaero) {
        iaero_top = std::max(i + 1, iaero_top);
        fneba[i]  = 1.;
        deltaz    = zqq[i+1] - zqq[i];

        // Distribution of AOD on the vertical homogeneous between 0 and zaero
        // Note, we used a formula based on concentration before v6.2.
        tauao3[i] = at_1d_rad->aod_o3_tot * deltaz / zqq[iaero_top];
        tauah2o[i]= at_1d_rad->aod_h2o_tot * deltaz / zqq[iaero_top];
      }

      // Calculation of SSA and Asymmetry factor for clouds using- Nielsen 2014
      // Note only the first two bands are taken, the third only is
      // approximately 0.
      // Useful for gas assymmetry, for albedo, either Chuang if BC, or Nielsen.
      const cs_real_t pioco3_1 = 1. - 33e-9 * req;  // FIXME MF: bad idea for precision
      const cs_real_t pioco3_2 = 1. - 1e-7  * req;
      const cs_real_t pioch2o_1 = 0.99999 - 149e-7 * req;
      const cs_real_t pioch2o_2 = 0.9985  - 92e-5  * req;

      fnebmax[i] = std::max(fnebmax[i + 1], fneray[i]);

      if (at_1d_rad->black_carbon_frac > epsc) {
        // SSA with black carbon fraction (Chuang 2002)
        const cs_real_t dm = req * 4.0 / 3.0;
        cs_real_t pioco3C = 0.;
        cs_real_t pioch2oC = 0.;

        // 4 to 12 because we there is no energy in the
        // 4 first spectral bands defined by Chuang 2002.
        for (int k = 4; k < 12; ++k) {
          const cs_real_t coef
            = (1.0 - exp(-beta3[k] * (at_1d_rad->black_carbon_frac - nu0)));
          const cs_real_t coef1
            = (1.0 - exp(-beta4[k] * (at_1d_rad->black_carbon_frac - nu0)));
          copioc20[k] = omega0[k] + beta1[k] * coef + beta2[k] * coef1;

          copioc[k] = (copioc20[k] * dm / dm0)
                    / (1.0 + 1.8 * copioc20[k] * (dm / dm0 - 1.0));

          piocv[k] = 1.0 - copioc[k];

          pioco3C  += coeff_E_o3[k]  * piocv[k];
          pioch2oC += coeff_E_h2o[k] * piocv[k];
        }

        pic_o3[i]  = pioco3C;
        pic_h2o[i] = pioch2oC;
      }

      else {
        // SSA and asymmetry factor (Nielsen 2014)
        const cs_real_t pioco3 = pioco3_1 * 0.24 + pioco3_2 * 0.76;
        const cs_real_t pioch2o = 0.60 * pioch2o_1 + 0.40 * pioch2o_2;
        pic_o3[i]  = pioco3;
        pic_h2o[i] = pioch2o;
      }

      // Gas asymmetry (Nielsen)
      const cs_real_t gasymo3 =
        (0.868 + 14e-5 * req - 61e-4 * exp(-0.25 * req)) * pioco3_1 * 0.24
        + (0.868 + 25e-5 * req - 63e-4 *exp(-0.25 * req)) * pioco3_2 * 0.76;

      const cs_real_t gasymh2o =
        (0.867 + 31e-5 * req - 78e-4 * exp(-0.195 * req)) * 0.60 * pioch2o_1
        + (0.864 + 54e-5 * req - 0.133 * exp(-0.194 * req)) * 0.40 * pioch2o_2;

      gco3[i]  = gasymo3;
      gch2o[i] = gasymh2o;
    }

    tauc[kmray + 1] = 0.;

    // 5.3 O3 absorption in presence of clouds

    // Calculation of the different albedos for O3 (LH 74)

    // Asymmetry factor and SSA for liquid water
    cs_real_t gasym = gco3[itop];
    cs_real_t pioc  = pic_o3[itop];

    // --- Calculation for cloudy layers ---
    _compute_reflection_transmission(pioc, 0., gasym, 0.,
                                     tauctot, 0.,
                                     refx, trax, epsc, 0.,
                                     mui, muzero_cor);
    cs_real_t rabarc = refx;
    cs_real_t tabarc = trax;

    // LH74 equation (15)
    cs_real_t rbarc = rabarc + tabarc * tabarc * albe / (1. - rabarc * albe);

    // --- Calculation for aerosol layers ---
    _compute_reflection_transmission(0., at_1d_rad->piaero_o3, 0.,
                                     at_1d_rad->gaero_o3, 0.,
                                     at_1d_rad->aod_o3_tot,
                                     refx, trax, epsc, 0.,
                                     mui, muzero_cor);

    cs_real_t rabara = refx;
    cs_real_t tabara = trax;

    // Effective reflectance for aerosol layer
    cs_real_t rbara = rabara + tabara * tabara * albe / (1. - rabara * albe);

    // in case there is an aerosol layer above the cloud layer

    if (has_aerosol && (iaero_top > itop)) {
      itop   = iaero_top;
      rbar   = rbara;
      rrbar2s = rabara;
    }

    if (at_1d_rad->verbosity > 0) {
      bft_printf("1D Radiative Model: top of cloud/aerosol layer "
                 "height = %f m (itop = %d)\n",
                 zqq[itop], itop);
    }

    // Calculation above the top of the cloud or aerosol layer

    //calculation have to start at the first level
    for (int l = itop; l <= kmray + 1; l++) {
      const cs_real_t zq    = zqq[l];
      cs_real_t zqp1;
      if (l == kmray + 1)
        zqp1 = zq;
      else
        zqp1 = zqq[l + 1];

      // Ozone amount traversed by the direct solar beam
      cs_real_t x   = m * ozone_amount(zq);
      cs_real_t xp1 = m * ozone_amount(zqp1);

      // Calculation of heat and radiation fluxes during cloudy sky
      cs_real_t zbas = zqq[itop];
      cs_real_t xstar =
        m * ozone_amount(zbas) + mbar * (ozone_amount(zbas) - ozone_amount(zq));
      cs_real_t xstarp1 =
        m * ozone_amount(zbas) + mbar * (ozone_amount(zbas) - ozone_amount(zqp1));

      // --- Cloudy sky absorption (index 0) ---
      const cs_real_t dud1 = 1. / (1. - rrbar2s * albe);
      fabso3c(l, 0) =
        muzero * fo * ( ( ray_ozone_absorption(x)
                        - ray_ozone_absorption(xp1)) * dud1
                       + rbarc * (ray_ozone_absorption(xstarp1)
                                  - ray_ozone_absorption(xstar)));

      // Direct downward radiation
      ddfso3[l] = muzero * fo * (z_ref - rrbar - ray_ozone_absorption(x));
      // Diffuse downward radiation (Rayleigh factor)
      dddfso3[l] = ddfso3[l] * (dud1 - 1.);
      // Global downward radiation
      dfso3[l] = ddfso3[l] * dud1;
      // Upward diffuse radiation under cloudy sky
      ufso3c(l, 0) =
        muzero * fo * (z_ref - rrbar - ray_ozone_absorption(xstar)) * rbarc;

      // Calculation of heat and radiation fluxes during  Clear sky
      zbas = zqq[k1];
      xstarp1 =   m * ozone_amount(zbas)
                + mbar * (ozone_amount(zbas) - ozone_amount(zqp1));

      xstar =   m * ozone_amount(zbas)
              + mbar * (ozone_amount(zbas) - ozone_amount(zq));

      fabso3c(l, 1) =
        muzero * fo * (  (ray_ozone_absorption(x) - ray_ozone_absorption(xp1))
                       * dud1
                       + albe * dud1 * (ray_ozone_absorption(xstarp1)
                                        - ray_ozone_absorption(xstar)));

      // Upward diffuse radiation for clear sky
      ufso3c(l, 1) =   muzero * fo
                             * (z_ref - rrbar - ray_ozone_absorption(xstar))
                             * albe * dud1;

      // sum depending on cloud fraction
      fabso3[l] =
        fnebmax[k1] * fabso3c(l, 0) + (1. - fnebmax[k1]) * fabso3c(l, 1);

      ufso3[l] =
        fnebmax[k1] * ufso3c(l, 0) + (1. - fnebmax[k1]) * ufso3c(l, 1);

      const cs_real_t dzx = ozone_gradient(zq);

      const cs_real_t denom1 = (z_ref - rrbar - ray_ozone_absorption(x));
      const cs_real_t denom2 = (z_ref - rrbar - ray_ozone_absorption(xstar));

      if (l < kmray + 1) {
        ckdown_suv_r[l] = d_ray_ozone_absorption(x, dzx) / denom1;
        ckdown_suv_f[l] = d_ray_ozone_absorption(x, dzx) / denom1;
        ckup_suv_f[l]   = d_ray_ozone_absorption(xstar, dzx) / denom2;
      }
    }

    // Calculation under the top of the cloud or the aerosol layer, the adding
    // Method with multiple diffusion is used

    // Top boundary conditions
    tra(kmray+1, 0) = 1.;

    // value from 44km to 11km
    trard(kmray+1, 0) = 1.;
    trad(kmray+1, 0) = 1.;
    ref(kmray+1, 0) = 0.;
    trat(kmray+1, 0) = 1.;
    reft(kmray+1, 0) = ref(kmray+1, 0);

    // Bottom boundary conditions
    ref(k0, 0) = albe;

    for (int l = k1; l <= kmray; l++) {
      tau(l, 0) = tauc(l) + tauao3(l);

      gasym = gco3[l];
      pioc  = pic_o3[l];

      //In the cloudy layers
      _compute_reflection_transmission(pioc, at_1d_rad->piaero_o3, gasym,
                                       at_1d_rad->gaero_o3,
                                       tauc[l], tauao3[l],
                                       refx, trax,
                                       epsc, 0., mui, muzero_cor);

      ref(l, 0) = fneray[l] * refx;
      tra(l, 0) = fneray[l] * trax;

      //In the aerosol layers
      _compute_reflection_transmission(0., at_1d_rad->piaero_o3, 0.,
                                       at_1d_rad->gaero_o3, 0., tauao3[l],
                                       refx, trax,
                                       epsc, 0., mui, muzero_cor);

      ref(l, 0) += (1. - fneray[l]) * refx;
      tra(l, 0) += (1. - fneray[l]) * trax;

      trard(l, 0) =   (fneray[l] * exp(-m * tauc[l]) + (1. - fneray[l]))
                    *  exp(-m * tauao3[l]);
    }

    // Downward addition of layers
    for (int l = kmray; l >= k1; l--) {
      cs_real_t drtt1 = 1. / (1. - reft(l+1, 0) * ref(l, 0));

      // Note:
      // R(top->l) = R(top->l+1)
      //           + T(top->l+1) R(l) T*(top->l+1) / (1 - R*(top->l+1)- R(l))
      // Equations 34 of LH74
      // Note R*(top->l) = R(top->l)
      //      T*(top->l+1) = T(top->l+1)

      // Eq. (33) LH74: reflection from top to level l
      reft(l, 0) =   reft(l+1, 0)
                   + trat(l+1, 0) * ref(l, 0) * trat(l+1, 0) * drtt1;

      // Eq. (34) LH74: transmission from top to level l
      trat(l, 0) = trat(l+1, 0) * tra(l, 0) * drtt1;

      // Transmission for direct radiation
      trad(l, 0) = trad(l+1, 0) * trard(l, 0);
    }

    // Upward addition (from surface to top)

    refb(k0, 0) = ref(k0, 0);

    for (int l = k1; l <= kmray+1; l++) {
      const cs_real_t dtrb1 = 1. -ref(l, 0) * refb(l-1, 0);
      refb(l, 0) = ref(l, 0) + tra(l, 0) * refb(l-1, 0) * tra(l, 0) / dtrb1;

      // Calculation of upward and downward fluxes and absorption
      const cs_real_t dud1 = 1. -reft(l, 0) * refb(l-1, 0);
      // Direct
      dowd(l, 0) = trad(l, 0);
      // Diffuse
      dowf(l, 0) = trat(l, 0) / dud1 - trad(l, 0);
      // Global == dowf(l,n) + dowd(l,n)
      dow(l, 0) = trat(l, 0) / dud1;
      upwf(l, 0) = refb(l-1, 0) * dow(l, 0);
      // Absorption from top to level l
      atln(l, 0) = 1. - reft(k1, 0) + upwf(l, 0) - dow(l, 0);
    }

    // If there is a cloud
    if (itop > k1) {
      for (int l = k1; l < itop; l++) {
        // addition of ozone absorption for heating in the
        // layers when adding method is used
        const cs_real_t zq   = zqq[l];
        const cs_real_t zqp1 = zqq[l + 1];
        const cs_real_t zbas = zqq[itop];
        // Note LH74 compute absn as atln[idx_ln] - atln[idx_lp1_n]
        // but can be simplified as:
        // upwf[idx_ln] - dow[idx_ln] - upwf[idx_lp1_n] + dow[idx_lp1_n]
        // absn[idx_ln] =   dow[idx_ln] * (refb(l-1,n) - 1.0)
        //                - dow[idx_lp1_n] * (refb[idx_ln] - 1.0)
        absn(l, 0) =   dow(l, 0) * (refb(l-1, 0) - 1.)
                     - dow(l+1, 0) * (refb(l, 0) - 1.);

        const cs_real_t x =    m * ozone_amount(zbas)
                            + mbar * (ozone_amount(zq)-ozone_amount(zbas));
        const cs_real_t xp1 =   m * ozone_amount(zbas)
                              + mbar * (ozone_amount(zqp1) - ozone_amount(zbas));
        const cs_real_t xstar = m * ozone_amount(zbas)
          + mbar * (ozone_amount(zqq[k1]) - ozone_amount(zbas))
          + mbar * (ozone_amount(zqq[k1]) - ozone_amount(zq));
        const cs_real_t xstarp1 = m * ozone_amount(zbas)
          + mbar * (ozone_amount(zqq[k1]) - ozone_amount(zbas))
          + mbar * (ozone_amount(zqq[k1]) - ozone_amount(zqp1));

        // taking into account ozone absorption in
        // the layers with clouds or aerosols
        const cs_real_t bas_ozone
          = (z_ref - ray_ozone_absorption(m*ozone_amount(zbas)));
        const cs_real_t x_ozone
          = ray_ozone_absorption(x) - ray_ozone_absorption(xp1);
        const cs_real_t ozone
          = (ray_ozone_absorption(xstarp1) - ray_ozone_absorption(xstar) );
        fabso3[l]
          = muzero * fo * (x_ozone + bas_ozone * absn(l, 0) + rbar * ozone);
        // fluxes calculation taking into account ozone absorption
        // Direct
        ddfso3[l]
          = muzero * fo * (z_ref - rrbar - ray_ozone_absorption(x)) * dowd(l, 0);

        // Diffuse:
        // two contributions: 1) transform direct into diffuse (dowf)
        //                    2) diffuse coming from the top
        dddfso3[l] =   muzero * fo * (z_ref - rrbar - ray_ozone_absorption(x))
                     * (dowf(l, 0) +   dow(l, 0) * albe * rrbar2s
                                     / (1. - rrbar2s * albe));

        // Global: dfso3 = ddfso3 + dddfso3  (we compute it via dow and factor)
        dfso3[l] = muzero * fo * (z_ref - rrbar - ray_ozone_absorption(x) )
                 * dow(l, 0) / (1. - rrbar2s * albe);

        // Upward (diffuse) radiation
        ufso3[l] = muzero * fo
                 * (z_ref - rrbar - ray_ozone_absorption(xstar) ) * upwf(l, 0);

        // gradient for absorption coefficient computation
        const cs_real_t dzx = ozone_gradient(zq);

        // absorption coefficient ckup and ckdown useful for 3D simulation
        ckdown_suv_r[l] = d_ray_ozone_absorption(x, dzx)
                        / (z_ref - rrbar - ray_ozone_absorption(x) );
        // optical depths must be changed to take into account transformation
        // direct->diffuse under cloud top
        ckdown_suv_f[l] = d_ray_ozone_absorption(x, dzx)
                        / (z_ref - rrbar - ray_ozone_absorption(x) );
        ckup_suv_f[l] = d_ray_ozone_absorption(xstar, dzx)
                      / (z_ref - rrbar - ray_ozone_absorption(xstar) );
      }

      // Calculation of upward flux above cloud or
      // aerosol layers taking into account
      // the upward flux transmitted by cloud or aerosol layers
      // if there is no cloud and no aerosol (itop == k1)
      // this term must NOT be added
      for (int l = itop; l <= kmray + 1; l++) {
        const cs_real_t zq   = zqq[l];
        const cs_real_t zbas = zqq[k1];
        const cs_real_t xstar = m * ozone_amount(zbas)
                              + mbar * (ozone_amount(zbas) - ozone_amount(zq));
        ufso3[l] = muzero * fo
                 * (z_ref - rrbar - ray_ozone_absorption(xstar) ) * upwf(itop, 0);
      }
    } // endif (itop > k1)

    // 6.4 Absorption by water vapor and liquid water (H20 band, SIR)

    // In that case we have to solve multiple diffusion. This is achieved
    // by means of the adding method following Lacis et Hansen, 1974
    // calculation of reflexivity and transmissivity for each vertical layer.
    for (int n = 0; n < 8; ++n) {
      // From 44 km -> 11 km
      cs_real_t dqqv = kn[n] * (qqvtot - qqv[kmray]) / 10.0;

      // Tod and ground Bcs
      tau(kmray+1, n) = dqqv;
      tra(kmray+1, n) = exp(-m * dqqv);     // == exp(-m * tau(kmray+1, n))

      // For direct radiation
      trard(kmray+1, n) = tra(kmray+1, n);  // exp(-m * tau(kmray+1, n));
      trad(kmray+1, n) = trard(kmray+1, n);

      ref(kmray+1, n) = 0.;
      ref(k0, n) = albe;

      trat(kmray+1, n) = tra(kmray+1, n);
      reft(kmray+1, n) = ref(kmray+1, n);

      // For direct radiation
      trad(kmray+1, n) = trard(kmray+1, n);

      for (int l = k1; l <= kmray; l++) {
        gasym = gch2o[l];
        pioc  = pic_h2o[l];

        // /10 for unit conversion
        dqqv = kn[n] * (qqv[l + 1] - qqv[l]) / 10.0;

        // In the cloud layers
        tau(l, n) = tauc[l] + dqqv + tauah2o[l];

        if (qlray[l] >= epsc) {
          _compute_reflection_transmission(pioc, at_1d_rad->piaero_h2o, gasym,
                                           at_1d_rad->gaero_h2o, tauc[l],
                                           tauah2o[l], refx, trax, epsc,
                                           dqqv, mui, muzero_cor);
          ref(l, n) =   fneray[l] * refx;
          tra(l, n) =   fneray[l] * trax
                      + (1. - fneray[l]) * exp(-mbarh2o * dqqv);
          if (has_aerosol) {
            _compute_reflection_transmission(0., at_1d_rad->piaero_h2o,
                                             0., at_1d_rad->gaero_h2o,
                                             0., tauah2o[l], refx0, trax0,
                                             epsc, dqqv, mui, muzero_cor);
            ref(l, n) = fneray[l] * refx + (1. - fneray[l]) * refx0;
            tra(l, n) = fneray[l] * trax + (1. - fneray[l]) * trax0;
          }
          // trard transmissivity for direct radiation
          trard(l, n) =   (fneray[l] * exp(-m * tauc[l]) + (1. - fneray[l]))
                        * exp(-m * (dqqv + tauah2o[l]));
        }

        else {
          // Clear sky layers
          ref(l, n) = 0.;
          tra(l, n) = exp(-mbarh2o * tau(l, n));
          trard(l, n) = exp(-m * (dqqv + tauah2o[l]));

          if (l >= itop) tra(l, n) = exp(-m * tau(l, n));

          if (has_aerosol) {
            _compute_reflection_transmission(0., at_1d_rad->piaero_h2o,
                                             0., at_1d_rad->gaero_h2o, 0.,
                                             tauah2o[l], refx, trax, epsc,
                                             dqqv, mui, muzero_cor);
            ref(l, n) = fneba[l] * refx;
            tra(l, n) = fneba[l] * trax
                        + (1. - fneba[l]) * exp(-mbarh2o * dqqv);
          }
        }
      }

      // Downward addition of layers
      for (int l = kmray; l >= k1; l--) {
        cs_real_t drtt1 = 1. / (1. - reft(l+1, n) * ref(l, n));
        // Note R(l->top) = R*(l->top)
        // R*(l->top) = R(l) + T*(l)*T(l) * R*(l+1>top)/(1 - R*(l+1->top).R(l))
        reft(l, n) =   reft(l+1, n)
                     + trat(l+1, n) * ref(l, n) * trat(l+1, n) * drtt1;
        // T(l->top) = T(l+1->top).T(l) /(1 - R*(l+1->top).R(l))
        // Note also it is equal to
        // T*(l->top) = T*(l+1->top)*T(l) /(1 - R*(l+1->top).R(l))
        trat(l, n) = trat(l+1, n) * tra(l, n) * drtt1;

        // trad for direct radiation
        trad(l, n) = trad(l+1, n) * trard(l, n);
      }

      // upward layer addition
      refb(k0, n) = ref(k0, n);

      // Note LH74 equation (33)
      // T*(l) = T(l)
      for (int l = k1; l <= kmray; l++) {
        cs_real_t dtrb1 = 1. / (1. - refb(l-1, n) * ref(l, n));
        refb(l, n) = ref(l, n) + tra(l, n) * refb(l-1, n) * tra(l, n) * dtrb1;
      }

      // Downward and upward fluxes
      for (int l = kmray + 1; l >= k1; l--) {
        // downward fluxes for direct radiation
        dowd(l, n) = trad(l, n);

        cs_real_t dud1 = 1. - reft(l, n) * refb(l-1, n);
        // Diffuse
        dowf(l,n) = trat(l, n) / dud1 - trad(l,n);

        // Global
        dow(l, n) = trat(l, n) / dud1;

        // Diffuse up = global up
        upwf(l, n) = refb(l-1, n) * dow(l, n);

        // calculation of absorption from 16km to level l
        // Note can be factorized
        atln(l, n) = pkn[n] * (1. - reft(k1, n) + upwf(l, n) - dow(l, n));
      }

      // Absorption per layer
      for (int l = kmray; l >= k1; l--) {
        // Note LH74 compute absn as:
        // absn(l,n) = atln(l,n) - atln(l+1,n)
        // but can be simplified as
        // upwf(l,n)-dow(l,n) - upwf(l+1,n) + dow(l+1,n)

        absn(l, n) = pkn[n] * (  dow(l, n)   * (refb(l-1, n) - 1.)
                               - dow(l+1, n) * (refb(l, n)   - 1.));
      }
    } // end loop 6.4

    //  summation over frequencies and estimation of absorption integrated
    //  on the whole spectru
    for (int l = kmray; l >= k1; l--) {
      fabsh2o[l] = 0.;
      for (int n = 0; n < 8; ++n) {
        fabsh2o[l] += absn(l, n);
      }

      fabsh2o[l] = fabsh2o[l] * fo * muzero;
    }

    // In the case with no clouds and aerosol in order to have exactly
    // the same expressions in 1D and 3D for the heating rate
    if (itop == k1) {
      for (int l = k1; l <= kmray; l++) {
        const cs_real_t y   = m * (qqvtot - qqv[l]);  // Note qqv(1) = 0
        const cs_real_t yp1 = m * (qqvtot - qqv[l+1]);
        const cs_real_t ystar = m * qqvtot + mbarh2o * qqv[l];
        const cs_real_t ystarp1= m * qqvtot + mbarh2o * qqv[l+1];

        fabsh2o[l] = muzero * fo * (  ray_sve(y) - ray_sve(yp1)
                                    + albe * (  ray_sve(ystarp1)
                                              - ray_sve(ystar)));

      }
    }

    // 5.5 heating in the layers
    for (int i = k1; i <= kmray; i++) {
      deltaz = zqq[i+1] - zqq[i];
      const cs_real_t cphum = cp0 * (1.0 + (cpvcpa - 1.) * qvray[i]);

      rayst[i] = (fabsh2o[i] + fabso3[i]) / deltaz / romray[i] / cphum;
    }

    // 5.6 calculation of solar fluxes
    // for global radiation, fd for direct radiation for the water vapor band
    // H20: SIR band

    for (int i = k1; i <= kmray + 1; i++) {
      dfsh2o[i]   = 0.;
      ufsh2o[i]   = 0.;
      ddfsh2o[i]  = 0.;
      dddfsh2o[i] = 0.;
      for (int n = 1; n < 8; ++n) {
        ufsh2o[i]  += pkn[n] * upwf(i, n);
        ddfsh2o[i] += pkn[n] * dowd(i, n);
        dddfsh2o[i] += pkn[n] * dowf(i, n);
      }
      ufsh2o[i]  *= fo * muzero;
      ddfsh2o[i] *= fo * muzero;
      dddfsh2o[i] *= fo * muzero;

      // Global
      dfsh2o[i] = ddfsh2o[i] + dddfsh2o[i];

      // The optical depth for gases have to be changed to take into account
      // the transformation of direct radiation in diffuse radiation under
      // the top of the cloud.
      cs_real_t ystar = m *   (qqvtot - qqv[itop])
                            + mbarh2o * (qqv[itop] + qqv[i]);

      cs_real_t y;
      if (i >= itop)
        y = m * (qqvtot - qqv[i]);
      else
        y = m * (qqvtot - qqv[itop]) + mbarh2o * (qqv[itop] - qqv[i]);

      // to test in the case with no clouds and aerosol in order to have exactly
      // the same expressions in 1D and 3D for the fluxes
      if (itop == k1) {
        dddfsh2o[i] = 0.;
        dfsh2o[i] = fo * muzero * (0.353 - ray_sve(y));
        ufsh2o[i]
          = fo * muzero * (0.353 - ray_sve(ystar)) * albe;
        ddfsh2o[i] = fo * muzero * (0.353 - ray_sve(y));
      }

      if (i < kmray + 1) {
        // (p_i/p_k) * (p_k/p0) = p_i/p0
        cs_real_t corp = preray[i] / 101315.0;
        cs_real_t dy =
          romray[i] * (qvray[i] * corp * sqrt(tkelvi / (temray[i] + tkelvi)));

        // Absorption coefficients
        ckdown_sir_r[i] = ray_sve_derivative(y, dy) / (0.353 - ray_sve(y));
        ckdown_sir_f[i] = ray_sve_derivative(y, dy) / (0.353 - ray_sve(y));
        ckup_sir_f[i] = ray_sve_derivative(ystar, dy) / (0.353 - ray_sve(ystar));
      }
    }

    // 6. Calculation of solar fluxes For the whole spectrum
    for (int k = k1; k <= kmray + 1; k++) {
      // Global (down and up) fluxes
      dfs[k] = dfsh2o[k] + dfso3[k];
      ufs[k] = ufsh2o[k] + ufso3[k];

      // direct radiation and diffuse mod (sum of vapor water band and O3 band)
      drfs[k] = ddfsh2o[k] + ddfso3[k];
      ddfs[k] = dddfsh2o[k] + dddfso3[k];

      _atmo_1d_rad.solu[k + ivertc * kmx] = ufs[k];
      _atmo_1d_rad.sold[k + ivertc * kmx] = dfs[k];
    }

    if (at_1d_rad->verbosity > 1) {
      bft_printf("\n--- 1D Radiative Model Solar Fluxes ---\n");
      bft_printf("fo = %e, mu0 = %e, albe = %e\n", fo, muzero, albe);
      bft_printf("itop = %d (z = %f m), k1 = %d (z = %f m)\n",
                 itop, zqq[itop], k1, zqq[k1]);
      bft_printf("dfsh2o[0] = %e, dfso3[0] = %e, dfs[0] = %e\n",
                 dfsh2o[k1], dfso3[k1], dfs[k1]);
      bft_printf("----------------------------------------\n\n");
    }

    // Note: Multiplication by transmission function for minor gases
    // Tmg is now taken into account by fo=fo*Tmg

    // Solar heating of the ground surface by downward global flux
    // see ground model
    // fos = dfs[k1] * (1. - albe);

    // Calculation of absorption coefficient ckup and ckdown useful
    // for 3D simulation
    // For water vapor without clouds and aerosols downward is only direct,
    // upward is only diffuse.

    // SIR band
    for (int k = k1; k <= kmray; k++) {
      cs_real_t tauapc = tauah2o[k] + tauc[k];
      deltaz = zqq[k + 1] - zqq[k];

      cs_real_t ckapcd = 0.;
      cs_real_t ckapcf = 0.;
      cs_real_t ck_aero_h2of = 0.;
      cs_real_t ck_aero_h2od = 0.;

      if (tauapc < epsc) {
        w0_sir[k] = 0.;
        g_apc_sir[k] = 0.;
      }
      else {
        cs_real_t picapc
          = (pic_h2o[k] * tauc[k] + at_1d_rad->piaero_h2o * tauah2o[k]) / tauapc;

        /* if we take into account asymmetry factor for forward diffuse radiation
         * Note apc means aerosols+clouds */
        cs_real_t gapc = (  pic_h2o[k] * tauc[k] * gch2o[k]
                          + at_1d_rad->piaero_h2o * tauah2o[k] * at_1d_rad->gaero_h2o)
                       / (tauapc * picapc);

        // Save values without Joseph correction for 3D
        g_apc_sir[k] = gapc;
        w0_sir[k] = picapc;

        // Absorption and forward diffusion
        ckapcf = tauapc / deltaz;

        // Direct (no Joseph correction)
        ckapcd = tauapc / (deltaz * muzero_cor);
        // Aerosol contributions
        ck_aero_h2of = tauah2o[k] / deltaz;
        ck_aero_h2od = tauah2o[k] / (deltaz * muzero_cor);
      }

      // Final coefficients
      ckup_sir_f[k] = ckup_sir_f[k]
                    + ck_aero_h2of * (1. - fneray[k])
                    + ckapcf * fneray[k];

      ckdown_sir_r[k] = ckdown_sir_r[k]
                      + ck_aero_h2od * (1. - fneray[k])
                      + ckapcd       * fneray[k];

      ckdown_sir_f[k] = ckdown_sir_f[k]
                      + ck_aero_h2of * (1. - fneray[k])
                      + ckapcf       * fneray[k];
    }

    // SUV band
    for (int k = k1; k <= kmray; k++) {
      cs_real_t tauapc = tauao3[k] + tauc[k];
      deltaz = zqq[k + 1] - zqq[k];

      cs_real_t ckapcd = 0.;
      cs_real_t ckapcf = 0.;
      cs_real_t ck_aero_o3f = 0.;
      cs_real_t ck_aero_o3d = 0.;

      if (tauapc < epsc) {
        w0_suv[k] = 0.;
        g_apc_suv[k] = 0.;
      }
      else {
        const cs_real_t picapc
          = (pic_o3[k] * tauc[k] + at_1d_rad->piaero_o3 * tauao3[k]) / tauapc;
        // if we take into account asymmetry factor for forward diffuse radiation
        const cs_real_t gapc = (  pic_o3[k] * tauc[k] * gco3[k]
                                +   at_1d_rad->piaero_o3 * tauao3[k]
                                  * at_1d_rad->gaero_o3)
                             / (tauapc * picapc);
        // Direct (no Joseph correction)
        ckapcd = tauapc / (deltaz * muzero_cor);
        w0_suv[k] = picapc;
        g_apc_suv[k] = gapc;
        // Absorption and forward diffusion
        ckapcf = tauapc / deltaz;
        ck_aero_o3f = tauao3[k] / deltaz;
        ck_aero_o3d = tauao3[k] / (deltaz * muzero_cor);
      }

      ckup_suv_f[k]
        = ckup_suv_f[k] + ck_aero_o3f * (1.0 - fneray[k]) + ckapcf * fneray[k];
      ckdown_suv_r[k]
        = ckdown_suv_r[k] + ck_aero_o3d * (1.0 - fneray[k]) + ckapcd * fneray[k];
      ckdown_suv_f[k]
        = ckdown_suv_f[k] + ck_aero_o3f * (1.0 - fneray[k]) + ckapcf * fneray[k];
    }

    /* In addition a source term has to be added in 3D for diffuse radiation.

       if mui=1/sqrt(3) quadrature method as LH74 the source term is added for
       both downward and upward radiation
       where  g1= 31/2(1-wo) , g3=(1-31/2µo)/2, g4=(1+31/2µo)/2 and t the total
       optical depth (gas + aerosol + cloud) t = tg + ta + tc.
       If mui=muzero_cor delta method the source term is added only in the downward
       diffuse radiation.
       In that condition the two-stream equations can be written:
       where  g1=(1-wo)/µo, g4=1 and t the total optical depth
       (gas + aerosol + cloud) t = tg + ta + tc. */

  } // end if muzero >= cs_math_epzero

  else {

    muzero = 0.;
    for (int k = k1; k <= kmray; k++) {
      rayst[k] = 0.;

      _atmo_1d_rad.solu[k + ivertc * kmx] = 0.;
      _atmo_1d_rad.sold[k + ivertc * kmx] = 0.;

      ckup_sir_f[k] = 0.;
      ckdown_sir_r[k] = 0.;
      ckdown_sir_f[k] = 0.;

      ckup_suv_f[k] = 0.;
      ckdown_suv_r[k] = 0.;
      ckdown_suv_f[k] = 0.;

      w0_sir[k] = 0.;
      w0_suv[k] = 0.;

      g_apc_sir[k] = 0.;
      g_apc_suv[k] = 0.;
    }

  }

  // Compute Boundary conditions for the 3D
  // (Director diFfuse) Solar radiance
  // at the top of the CFD domain and the absorption coefficients
  if (   cs_rad_time_is_active() == true
      && cs_field_by_name_try("spectral_rad_incident_flux") != nullptr) {

    const int *bc_type = cs_glob_bc_type;

    /* dimension of radiative fields */
    int idx = -1;
    const int dim = cs_glob_rad_transfer_params->nwsgg;
    cs_real_t *cpro_gapc = cs_field_by_name("asymmetry_factor")->val;
    cs_real_t *cpro_w0 = cs_field_by_name("simple_diffusion_albedo")->val;
    cs_real_t *cpro_ck_up = cs_field_by_name("rad_absorption_coeff_up")->val;
    cs_real_t *cpro_ck_down = cs_field_by_name("rad_absorption_coeff_down")->val;
    cs_real_t *bpro_rad_inc = cs_field_by_name("spectral_rad_incident_flux")->val;

    int rad_atmo_model = cs_glob_rad_transfer_params->atmo_model;
    // Direct Solar (denoted by _r) (for Solar IR band absorbed by H20)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR) {
      idx += 1;
      _interpolate_boundary_rad_incident_flux(idx, dim, kmx,
                                              bc_type, n_b_faces,
                                              muzero, muzero_cor,
                                              zqq, ddfsh2o,
                                              b_face_cog, bpro_rad_inc);
      _interpolate_coeff(idx, dim, kmx,
                         n_cells, zray, ckdown_sir_r,
                         cell_cen, cpro_ck_down);

      _interpolate_coeff(idx, dim, kmx,
                         n_cells, zray, w0_sir,
                         cell_cen, cpro_w0);

    }

    // Direct Solar (denoted by _r) (for visible UV (SUV) band absorbed by O3)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND) {
      idx += 1;
      _interpolate_boundary_rad_incident_flux(idx, dim, kmray + 1,
                                              bc_type, n_b_faces,
                                              muzero, muzero_cor,
                                              zqq, ddfso3,
                                              b_face_cog, bpro_rad_inc);

      _interpolate_coeff(idx, dim, kmx,
                         n_cells, zray, ckdown_suv_r,
                         cell_cen, cpro_ck_down);

      _interpolate_coeff(idx, dim, kmx,
                         n_cells, zray, w0_sir,
                         cell_cen, cpro_w0);
    }

    // Direct Solar (_r) if O3 band not activated: add to total solar band
    else if (rad_atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR) {
      _interpolate_boundary_rad_incident_flux(idx, dim, kmx,
                                              bc_type, n_b_faces,
                                              muzero, muzero_cor,
                                              zqq, ddfso3,
                                              b_face_cog, bpro_rad_inc);

      _interpolate_coeff(idx, dim, kmx,
                         n_cells, zray, ckdown_suv_r,
                         cell_cen, cpro_ck_down);
    }

    // Diffuse solar radiation incident up and down (SIR band)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR) {
      idx += 1;
      _interpolate_boundary_rad_incident_flux(idx, dim, kmray + 1,
                                              bc_type, n_b_faces,
                                              muzero, 1.0,
                                              zqq, dddfsh2o,
                                              b_face_cog, bpro_rad_inc);

      // Downward
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckdown_sir_f,
                         cell_cen, cpro_ck_down);

      // Upward
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckup_sir_f,
                         cell_cen, cpro_ck_up);

      // Simple diffusion albedo w0
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, w0_sir,
                         cell_cen, cpro_w0);

      // Asymmetry factor
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, g_apc_sir,
                         cell_cen, cpro_gapc);
    }

    // Diffuse solar radiation incident up and down (SUV - O3 band)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND) {
      idx += 1;
      _interpolate_boundary_rad_incident_flux(idx, dim, kmray + 1,
                                              bc_type, n_b_faces,
                                              muzero, 1.0,
                                              zqq, dddfso3,
                                              b_face_cog, bpro_rad_inc);

      // Downward
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckdown_suv_f,
                         cell_cen, cpro_ck_down);

      // Upward
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckup_suv_f,
                         cell_cen, cpro_ck_up);

      // Simple diffusion albedo w0
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, w0_suv,
                         cell_cen, cpro_w0);

      // Asymmetry factor
      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, g_apc_suv,
                         cell_cen, cpro_gapc);
    }
    else if (rad_atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR) {
      _interpolate_boundary_rad_incident_flux(idx, dim, kmray + 1,
                                              bc_type, n_b_faces,
                                              muzero, 1.0,
                                              zqq, dddfso3,
                                              b_face_cog, bpro_rad_inc);

      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckdown_suv_f,
                         cell_cen, cpro_ck_down);

      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, ckup_suv_f,
                         cell_cen, cpro_ck_up);

      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, w0_suv,
                         cell_cen, cpro_w0);

      _interpolate_coeff(idx, dim, kmray + 1,
                         n_cells, zray, g_apc_suv,
                         cell_cen, cpro_gapc);
  }

 } // end if

  // Write solar parameter dumps
  if (at_1d_rad->verbosity >= 2) {
    cs_array<cs_real_t> rayst_h2o(kmx);
    cs_array<cs_real_t> rayst_o3(kmx);
    for (int i = k1; i <= kmray; i++) {
      deltaz = zqq[i+1] - zqq[i];
      const cs_real_t cphum = cp0 * (1.0 + (cpvcpa - 1.) * qvray[i]);
      rayst_h2o[i] = fabsh2o[i] / deltaz / romray[i] / cphum;
      rayst_o3[i]  = fabso3[i]  / deltaz / romray[i] / cphum;
    }

    _dump_scalar_to_file("direct.txt", heuray, drfs[k1]);
    _dump_scalar_to_file("global.txt", heuray, dfs[k1]);
    _dump_scalar_to_file("direct_h2o.txt", heuray, ddfsh2o[k1]);
    _dump_scalar_to_file("global_h2o.txt", heuray, dfsh2o[k1]);
    _dump_scalar_to_file("direct_o3.txt", heuray, ddfso3[k1]);
    _dump_scalar_to_file("global_o3.txt", heuray, dfso3[k1]);

    _dump_array_to_file("heat_o3.txt",
                        heuray, rayst_o3.data(), kmx);
    _dump_array_to_file("heat_h2o.txt",
                        heuray, rayst_h2o.data(), kmx);
    _dump_array_to_file("sold_h2o.txt",
                        heuray, dfsh2o.data(), kmx + 1);
    _dump_array_to_file("sold_o3.txt",
                        heuray, dfso3.data(), kmx + 1);
    _dump_array_to_file("sold_o3_direct.txt",
                        heuray, ddfso3.data(), kmx + 1);
    _dump_array_to_file("sold_o3_diffuse.txt",
                        heuray, dddfso3.data(), kmx + 1);
    _dump_array_to_file("solu_h2o.txt",
                        heuray, ufsh2o.data(), kmx + 1);
    _dump_array_to_file("solu_o3.txt",
                        heuray, ufso3.data(), kmx + 1);
    _dump_array_to_file("qlray.txt",
                        heuray, qlray, kmx);
    _dump_array_to_file("sold_h2o_diffuse.txt",
                        heuray, dddfsh2o.data(), kmx + 1);
    _dump_array_to_file("sold_h2o_direct.txt",
                        heuray, ddfsh2o.data(), kmx + 1);
    _dump_array_to_file("fneray.txt",
                        heuray, fneray, kmx);
    _dump_array_to_file("fnebmax.txt",
                        heuray, fnebmax.data(), kmx + 1);
    _dump_array_to_file("taua_h2o.txt",
                        heuray, tauah2o.data(), kmx);
    _dump_array_to_file("taua_o3.txt",
                        heuray, tauao3.data(), kmx);
  }

 *albe_p = albe;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute infrared flux divergence profile and downward flux at
 *        ground level relying on a 1D radiative scheme.
 *
 * \param[in]  ivertc      index of vertical profile
 * \param[in]  k1          index for ground level
 * \param[in]  kmray       vertical levels number for radiation
 * \param[in]  emis        ground surface emissivity
 * \param[out] qqv         optical depth for water vapor (0,zqq)
 * \param[out] qqqv        idem qqv but for intermediates vertical levels (zray)
 * \param[in]  qqvinf      idem qqv but for altitude above 11000m
 * \param[in]  zqq         vertical levels (interfaces)
 * \param[in]  zray        vertical levels (volumes)
 * \param[in]  temray      temperature in Celsius
 * \param[in]  qvray       specific humidity for water vapor
 * \param[in]  qlray       specific humidity for liquid water
 * \param[in]  fnerir      cloud fraction
 * \param[in]  romray      air density
 * \param[in]  preray      pressure
 * \param[in]  aeroso      aerosol concentration in micro-g/m3
 * \param[in]  t_surf      surface temperature
 * \param[in]  p_surf      surface pressure
 * \param[out] rayi        IR flux divergence
 * \param[in]  ncray       Number of droplets interpolated on vertical grid
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_compute_infrared(const int        ivertc,
                                const int        k1,
                                const int        kmray,
                                const cs_real_t  emis,
                                cs_real_t        qqv[],
                                cs_real_t        qqqv[],
                                cs_real_t        *qqvinf_p,
                                cs_real_t        zqq[],
                                cs_real_t        zray[],
                                const cs_real_t  temray[],
                                const cs_real_t  qvray[],
                                const cs_real_t  qlray[],
                                cs_real_t        fnerir[],
                                const cs_real_t  romray[],
                                const cs_real_t  preray[],
                                const cs_real_t  aeroso[],
                                const cs_real_t  t_surf,
                                const cs_real_t  p_surf,
                                cs_real_t        rayi[],
                                const cs_real_t  ncray[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  const cs_atmo_1d_rad_t *at_1d_rad = cs_glob_atmo_1d_rad;

  const int kmx = at_1d_rad->nlevels_max;

  // downward IR flux at the ground
  cs_real_t foir = 0.;
  cs_real_t qqvinf = *qqvinf_p;

  const cs_real_t cp0 = cs_glob_fluid_properties->cp0;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t cpvcpa = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;

  cs_real_t *ufir = nullptr;
  cs_real_t *qc  = nullptr, *rov = nullptr;
  cs_real_t *roc = nullptr, *rol = nullptr;
  cs_real_t *qv0 = nullptr, *qqc = nullptr;
  cs_real_t *qql = nullptr, *tqq = nullptr;
  cs_real_t *dt4 = nullptr, *dfir = nullptr;
  cs_real_t *dz0 = nullptr, *kliq = nullptr;
  cs_real_t *qqqc = nullptr, *qqql = nullptr;
  cs_real_t *qco2 = nullptr, *qqco2 = nullptr;
  cs_real_t *ckup = nullptr, *ckdown = nullptr;
  cs_real_t *pspo = nullptr, *pspoqq = nullptr;
  cs_real_t *qqqco2 = nullptr, *roco2 = nullptr;

  // indexes for presence of clouds, aerosols
  bool inua = false, iaer = false;
  int iaero_top = -1;
  cs_real_t dz_aero = 0., xqqlinf = 0., taul = 0.;

  CS_MALLOC(qc,     kmx, cs_real_t);
  CS_MALLOC(rov,    kmx, cs_real_t);
  CS_MALLOC(roc,    kmx, cs_real_t);
  CS_MALLOC(rol,    kmx, cs_real_t);
  CS_MALLOC(qv0,    kmx, cs_real_t);
  CS_MALLOC(dt4,    kmx, cs_real_t);
  CS_MALLOC(dz0,    kmx, cs_real_t);
  CS_MALLOC(ckup,   kmx, cs_real_t);
  CS_MALLOC(qqqc,   kmx, cs_real_t);
  CS_MALLOC(qqql,   kmx, cs_real_t);
  CS_MALLOC(pspo,   kmx, cs_real_t);
  CS_MALLOC(qco2,   kmx, cs_real_t);
  CS_MALLOC(roco2,  kmx, cs_real_t);
  CS_MALLOC(ckdown, kmx, cs_real_t);
  CS_MALLOC(qqqco2, kmx, cs_real_t);

  CS_MALLOC(tqq,    kmx + 1, cs_real_t);
  CS_MALLOC(qqc,    kmx + 1, cs_real_t);
  CS_MALLOC(qql,    kmx + 1, cs_real_t);
  CS_MALLOC(ufir,   kmx + 1, cs_real_t);
  CS_MALLOC(dfir,   kmx + 1, cs_real_t);
  CS_MALLOC(kliq,   kmx + 1, cs_real_t);
  CS_MALLOC(qqco2,  kmx + 1, cs_real_t);
  CS_MALLOC(pspoqq, kmx + 1, cs_real_t);

  for (int k = 0; k < kmx; k++) {
    rov[k]  = 0.;
    roc[k]  = 0.;
    rol[k]  = 0.;
    roco2[k]  = 0.;
    qv0[k]  = 0.;
    qc[k]   = 0.;
    qco2[k] = 0.;
    qqv[k]  = 0.;
    qqc[k]  = 0.;
    qqco2[k]  = 0.;
    qql[k]  = 0.;
    qqqv[k] = 0.;
    qqqc[k] = 0.;
    qqqco2[k] = 0.;
    qqql[k] = 0.;
    pspo[k] = 0.;
    pspoqq[k] = 0.;
    dt4[k]  = 0.;
    tqq[k]  = 0.;
    rayi[k] = 0.;
    if (at_1d_rad->has_aerosol > 0 && aeroso[k] > 1.e-8) iaer = true;
  }

  pspoqq[kmx] = 0.;
  tqq[kmx] = 0.;
  qqv[kmx] = 0.;
  qqc[kmx] = 0.;
  qqco2[kmx]  = 0.;
  qql[kmx] = 0.;
  qqqv[kmx] = 0.;

  // Upper layers contribution (11000-44000) to optical depth
  qqvinf = 0.0050;
  constexpr cs_real_t qqcinf = 3.28e-7;
  const cs_real_t sig = cs_physical_constants_stephan;
  const cs_real_t cetyp = cs_glob_fluid_properties->rvsra / 1.1340;

  // Diffusion coefficient integrated over all directions
  constexpr cs_real_t beta = 1.660;

  // Calculation for the liquid water absorption coefficient
  // (same equivalent radius as solar radiation)
  for (int k = k1; k <= kmray; k ++) {
    // Liquid water content in g/m3 in the layers
    const cs_real_t wh2ol = 1.e3 * romray[k] * qlray[k];

    // Absorption coefficient depending on the equivalent radius
    // of the cloud droplets
    const cs_real_t rm = 30.0 * wh2ol + 2.0;

    // in case, microphysics data is available
    cs_real_t req = 1.5 * rm;
    if (ncray[k] > cs_math_epzero && qlray[k] > cs_math_epzero) {
      req = 1e6 * cbrt(    (3.0 * romray[k] * qlray[k])
                        / (4.0 * cs_math_pi * 1000.0 * ncray[k] * 1e6))
          * exp(cs_math_pow2(at_opt->sigc));
    }

    kliq[k] = beta * 0.75 * 1e3 / req;
    // Estimation of the inverse of the aerosol layer
    if (iaer == true && zray[k] <= at_1d_rad->zaero) {
      fnerir[k] = 1.0;
      iaero_top = std::max(iaero_top, k+1);
    }
  }

  if (iaero_top != -1)
    dz_aero = 1.0 / zqq[iaero_top];

  for (int k = k1; k < kmray; ++k) {
    dz0[k] = zray[k + 1] - zray[k];
  }
  dz0[kmray] = 0.0;

  // Warning zqq[kmx] is not 16 000 for the moment...
  cs_real_t zqq0 = 44000.0;
  cs_real_t tvsups, dtvsups;

  // 1. Computation to estimate absorption for water vapor and its dimer

  for (int k = k1; k <= kmray; ++k) {

    if (qlray[k] > 1e-8) inua = true;

    const cs_real_t corp = preray[k] / 101315.0;
    const cs_real_t cort = tkelvi / (temray[k] + tkelvi);
    pspo[k] = preray[k] / p_surf;

    qv0[k] = qvray[k] * corp * sqrt(cort);
    rov[k] = romray[k] * qv0[k];
    qco2[k] = at_1d_rad->conco2 * pow(corp, 0.75) * pow(cort, 0.325);
    roco2[k] = romray[k] * qco2[k];
    qc[k] =   qvray[k] * qvray[k] * corp
            * exp(1800.0 * (1.0 / (temray[k] + tkelvi) - 1.0 / 296.0))
            * cetyp;
    roc[k] = romray[k] * qc[k];
  }

  // Same for cloud layers (aerosol only for clear sky condition)
  if (inua) {
    for (int k = k1; k <= kmray; k++) {
      if (zray[k] <= at_1d_rad->zaero)
        rol[k] =   kliq[k] * romray[k] * qlray[k]
                 + beta * at_1d_rad->aod_ir * dz_aero;
      else
        rol[k] = kliq[k] * romray[k] * qlray[k];
    }
  }

  else if (iaer) {
    for (int k = k1; k <= kmray; k++)
      if (zray[k] <= at_1d_rad->zaero)
        rol[k] = beta * at_1d_rad->aod_ir * dz_aero;
  }

  // 2. optical depth calculation for water vapor and its dimer
  // qqq corresponds to standard levels, qq to intermediate levels

  // surface temperature (Kelvin)
  tqq[k1] = t_surf + tkelvi;

  //  Top interface temperature (isothermal above the domain)
  tqq[kmray + 1] = temray[kmray] + tkelvi;
  pspoqq[kmray + 1] = pspo[kmray];

  for (int k = k1 + 1; k <= kmray; ++k) {
    const cs_real_t alpha_k = (zqq[k] - zray[k-1]) / (zray[k] - zray[k-1]);
    tqq[k] = alpha_k * temray[k] + (1.0 - alpha_k) * temray[k-1] + tkelvi;
    pspoqq[k] = alpha_k * pspo[k] + (1.0 - alpha_k) * pspo[k - 1];
    // dT^4/dz × dz = 4 T³ ΔT
    dt4[k] = 4.0 * cs_math_pow3(tqq[k]) * (temray[k] - temray[k-1]);
  }

  // Boundary dual cell: transition from soil surface to center of layer 1
  dt4[k1] = 4.0 * cs_math_pow3(tqq[k1]) * (temray[k1] - t_surf);

  for (int k = k1; k <= kmray; ++k) {
    qqqv[k] = qqv[k] + (zray[k] - zqq[k]) * rov[k];
    qqqc[k] = qqc[k] + (zray[k] - zqq[k]) * roc[k];
    qqqco2[k] = qqco2[k] + (zray[k] - zqq[k]) * roco2[k];
    qqv[k+1] = qqqv[k] + (zqq[k+1] - zray[k]) * rov[k];
    qqc[k+1] = qqqc[k] + (zqq[k+1] - zray[k]) * roc[k];
    qqco2[k+1] = qqqco2[k] + (zqq[k+1] - zray[k]) * roco2[k];
  }

  pspoqq[k1] = pspo[k1];
  // TODO should be qqv(kmray+1)
  const cs_real_t xqqvinf   = qqqv[kmray] + qqvinf;
  const cs_real_t xqqcinf   = qqqc[kmray] + qqcinf;
  // Global integral from 0 to 44km
  const cs_real_t xqqco2inf = qqco2[kmray + 1];

  // 3. Optical depth calculation for liquid water in cloudy cases

  if (inua) {
    for (int k = k1; k <= kmray; k++) {
      qqql[k] = qql[k] + (zray[k] - zqq[k]) * rol[k];
      qql[k + 1] = qqql[k] + (zqq[k + 1] - zray[k]) * rol[k];
    }
    xqqlinf = qql[kmray + 1];
  }

  // 4. IR downward flux computation for the ground level

  cs_real_t fo = 0.0;
  const cs_real_t t41 = cs_math_pow4(tqq[k1]);
  const cs_real_t t4zt = cs_math_pow4(temray[kmray] + tkelvi);

  cs_real_t abcsups, dabcsups;
  cs_real_t abco2, dabco2, tauv, dtauv;

  // For cloudy sky (but also valid for clear sky)
  cs_real_t fn = fnerir[k1];
  for (int k = k1; k <= kmray; k++) {
    // cloud fraction estimation (we take the maximum)
    _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fnerir[k], qql[k]);

    _compute_ir_h2o_dimer_transmission
      (tauv, dtauv, qqv[k], qv0[k1], qqc[k], qc[k1], romray[k]);

    // Compute absorption from zbas to zqq
    _compute_ir_co2_o3_absorption(zqq[k1], zqq[k], abco2,
                                  dabco2, qqv[k], qv0[k1],
                                  qqco2[k], qco2[k1], romray[k]);

    fo -= (1.0 - (1.0 + fn * (taul - 1.0))
               * (tauv - abco2)) * dt4[k];
  }

  // Last level: 11km to 44km, assumed to be isothermal
  _compute_ir_h2o_dimer_transmission(tauv, dtauv, xqqvinf, qv0[k1], xqqcinf,
                                     qc[k1], romray[kmray]);

  _compute_ir_co2_o3_absorption(zqq[k1], zqq0, abco2, dabco2, xqqvinf, qv0[k1],
                                xqqco2inf, qco2[k1], romray[kmray]);

  fo += t4zt * (1.0 - (1.0 + fn * (taul - 1.0)) * (tauv - abco2));
  foir = fo;

  // IR flux calculation in the vertical layers
  cs_real_t qqqqv, qqqqc, qqqqco2, qqqql;
  for (int i = k1; i <= kmray; i++) {
    // cloud fraction estimation (we take the maximum)
    dfir[i] = 0.0;
    ufir[i] = 0.0;
    fn = fnerir[i];
    for (int k = i + 1; k <= kmray; k++) {
      qqqqv = qqv[k] - qqqv[i];
      qqqqc = qqc[k] - qqqc[i];
      qqqql = qql[k] - qqql[i];
      qqqqco2 = qqco2[k] - qqqco2[i];

      // Cloud fraction estimation (we take the maximum)
      _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fnerir[k], qqqql);

      // Downward fluxes
      _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[i], qqqqc,
                                         qc[i], romray[i]);

      // Compute from i to level k
      _compute_ir_co2_o3_absorption(zqq[i], zqq[k], abco2, dabco2, qqqqv,
                                    qv0[i], qqqqco2, qco2[i], romray[i]);

      dfir[i] -= (1.0 - (1.0 + fn * (taul - 1.0)) * (tauv - abco2)) * dt4[k];
    }

    qqqqv = xqqvinf - qqqv[i];
    qqqqc = xqqcinf - qqqc[i];
    qqqqco2 = xqqco2inf - qqqco2[i];

    _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[i],
                                       qqqqc, qc[i], romray[i]);
    _compute_ir_co2_o3_absorption(zqq[k1], zqq0, abco2, dabco2, qqqqv,
                                  qv0[i], qqqqco2, qco2[i], romray[i]);

    dfir[i] += t4zt * (1.0 - (1.0 + fn * (taul - 1.0)) * (tauv - abco2));

    // upward fluxes
    if (i > k1) {
      fn = fnerir[i];
      for (int k = k1; k <= i; k++) {
        // Cloud fraction estimation (we take the maximum)
        qqqqv = qqqv[i] - qqv[k];
        qqqqc = qqqc[i] - qqc[k];
        qqqql = qqql[i] - qql[k];
        qqqqco2 = qqqco2[i] - qqco2[k];
        _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fnerir[k], qqqql);

        _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[i],
                                           qqqqc, qc[i], romray[i]);
        _compute_ir_co2_o3_absorption(zqq[i], zqq[k], abco2, dabco2, qqqqv,
                                      qv0[i], qqqqco2, qco2[i], romray[i]);

        ufir[i] += (1.0 - (1.0 + fn * (taul - 1.0)) * (tauv - abco2)) * dt4[k];
      }
    }

    // Contribution of the upward flux reflected at the ground.
    // The modification proposed by Ponnulakshmi and al. (2009) for the upward
    // flux reflected by the ground when emissivity is different from 1
    // requires to compute the integral below (a3).
    cs_real_t a3 = 0.0;
    for (int k = k1; k <= kmray; k++) {
      qqqqv = qqqv[i] + qqv[k];
      qqqqc = qqqc[i] + qqc[k];
      qqqql = qqql[i] + qql[k];
      qqqqco2 = qqqco2[i] + qqco2[k];
      _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[i],
                                         qqqqc, qc[i], romray[i]);
      // Integration from i to k
      _compute_ir_co2_o3_absorption(zqq[i], zqq[k], abco2, dabco2, qqqqv,
                                    qv0[i], qqqqco2, qco2[i], romray[i]);

      fn = fnerir[k1];
      for (int ineb = k1; ineb <= kmray; ineb++)
        fn = std::max(fn, fnerir[ineb]);

      _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fn, qqqql);

      a3 -= (1.0 - (1.0 + fn * (taul - 1.0)) * (tauv - abco2)) * dt4[k];
    }

    qqqqv = qqqv[i] + xqqvinf;
    qqqqc = qqqc[i] + xqqcinf;
    qqqqco2 = qqqco2[i] + xqqco2inf;
    _compute_ir_h2o_dimer_transmission(tvsups, dtvsups, qqqqv, qv0[i],
                                       qqqqc, qc[i], romray[i]);
    _compute_ir_co2_o3_absorption(zqq[k1], zqq[i], abcsups, dabcsups,
                                  qqqqv, qv0[i], qqqqco2, qco2[i], romray[i]);

    const cs_real_t foirs = a3
      + (1.0- (1.0 + fn * (taul - 1.0)) * (tvsups - abcsups)) * t4zt;

    if (i > k1)
      ufir[i] += emis * t41 + (1.0 - emis) * foirs;
    else
      ufir[k1] = emis * t41 + (1.0 - emis) * foir;

  } //end IR

  // 5. Cooling in the vertical layers

  // For cloudy sky (also valid for clear sky)
  for (int k = k1; k <= kmray; k++) {

    const cs_real_t ctray
      = sig / romray[k] / (cp0 * (1.0 + (cpvcpa - 1.0) * qvray[k]));
    // a1: contribution from 0 to z
    cs_real_t a1 = 0.0;
    cs_real_t dul = rol[k];

    for (int kk = k1; kk <= k; kk++) {
      qqqqv = qqqv[k] - qqv[kk];
      qqqqc = qqqc[k] - qqc[kk];
      qqqql = qqql[k] - qql[kk];
      qqqqco2 = qqqco2[k] - qqco2[kk];

      // Integration from k to kk
      _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[k],
                                         qqqqc, qc[k], romray[k]);
      _compute_ir_co2_o3_absorption(zqq[k], zqq[kk], abco2, dabco2, qqqqv,
                                    qv0[k], qqqqco2, qco2[k], romray[k]);

      fn = fnerir[kk];
      for (int ineb = kk; ineb <= k; ineb++)
        fn = std::max(fn, fnerir[ineb]);
      _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fn, qqqql);

      a1 += dt4[kk] * (  (1.0 + fn * (taul - 1.0)) * (dabco2 - dtauv)
                       + taul * dul * (tauv - abco2));
    }

    // a2: contribution from z to ztop
    cs_real_t a2 = 0.0;
    for (int kk = k + 1; kk <= kmray; kk++) {
      qqqqv = qqv[kk] - qqqv[k];
      qqqqc = qqc[kk] - qqqc[k];
      qqqql = qql[kk] - qqql[k];
      qqqqco2 = qqco2[kk] -  qqqco2[k];
      _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[k],
                                         qqqqc, qc[k], romray[k]);
      _compute_ir_co2_o3_absorption(zqq[k], zqq[kk], abco2, dabco2, qqqqv,
                                    qv0[k], qqqqco2, qco2[k], romray[k]);
      fn = fnerir[kk];
      for (int ineb = k; ineb <= kk; ineb++)
        fn = std::max(fn, fnerir[ineb]);
      _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fn, qqqql);

      a2 += dt4[kk]*(  (1.0 + fn * (taul - 1.0)) * (dabco2-dtauv)
                       + taul * dul * (tauv - abco2));
    }

    // a3: contribution of the upward flux reflected at the ground
    // Modification of Ponnulakshmi and al. (2009)
    cs_real_t a3 = 0.0;
    for (int kk = k1; kk <= kmray; kk++) {
      qqqqv = qqqv[k] + qqv[kk];
      qqqqc = qqqc[k] + qqc[kk];
      qqqql = qqql[k] + qql[kk];
      qqqqco2 = qqqco2[k] + qqco2[kk];
      _compute_ir_h2o_dimer_transmission(tauv, dtauv, qqqqv, qv0[k],
                                         qqqqc, qc[k], romray[k]);
      _compute_ir_co2_o3_absorption(zqq[k], zqq[kk], abco2, dabco2, qqqqv,
                                    qv0[k], qqqqco2, qco2[k], romray[k]);

      fn = fnerir[k1];
      for (int ineb = k1; ineb <= kmray; ineb++)
        fn = std::max(fn, fnerir[ineb]);

      _cf_estimate(fn, taul, 1e-3, 0.0, 1.0, fn, qqqql);

      a3 += dt4[kk] * ( (1.0 + fn * (taul - 1.0)) * (dabco2 - dtauv)
          + taul * dul * (tauv - abco2));
    }

    // Contribution from z to infinity
    qqqqv = xqqvinf - qqqv[k];
    qqqqc = xqqcinf - qqqc[k];
    qqqql = xqqlinf - qqql[k];
    qqqqco2 = xqqco2inf - qqqco2[k];
    cs_real_t tvsup, dtvsup, abcsup, dabcsup;
    _compute_ir_h2o_dimer_transmission(tvsup, dtvsup, qqqqv, qv0[k],
                                       qqqqc, qc[k], romray[k]);
    _compute_ir_co2_o3_absorption(zqq[k], zqq0, abcsup, dabcsup, qqqqv,
                                  qv0[k], qqqqco2, qco2[k], romray[k]);

    cs_real_t tlsup;
    cs_real_t fns = fnerir[k];
    for (int ineb = k; ineb <= kmray; ineb++)
      fns = std::max(fns, fnerir[ineb]);

    _cf_estimate(fns, tlsup, 1e-3, 0.0, 1.0, fns, qqqql);

    // ground contribution transmitted by lower layers (0-z)
    qqqqv = qqqv[k];
    qqqqc = qqqc[k];
    qqqql = qqql[k];
    qqqqco2 = qqqco2[k];
    cs_real_t tvinfe, dtvinfe, abcinfe, dabcinfe;
    _compute_ir_h2o_dimer_transmission(tvinfe, dtvinfe, qqqqv, qv0[k],
                                       qqqqc, qc[k], romray[k]);

    _compute_ir_co2_o3_absorption(zqq[k1], zqq[k], abcinfe, dabcinfe, qqqqv,
                                  qv0[k], qqqqco2, qco2[k], romray[k]);
    cs_real_t tlinfe;
    cs_real_t fni = fnerir[k1];
    for (int ineb = k1; ineb <= k; ineb++)
      fni = std::max(fni, fnerir[ineb]);
    _cf_estimate(fni, tlinfe, 1e-3, 0.0, 0.0, fni, qqqql);

    // contribution of the upward flux reflected by the ground
    qqqqv = qqqv[k] + xqqvinf;
    qqqqc = qqqc[k] + xqqcinf;
    qqqql = qqql[k] + xqqlinf;
    qqqqco2 = qqqco2[k] + xqqco2inf;
    _compute_ir_h2o_dimer_transmission(tvsups, dtvsups, qqqqv, qv0[k],
                                       qqqqc, qc[k], romray[k]);

    _compute_ir_co2_o3_absorption(zqq[k1], zqq0, abcsups, dabcsups, qqqqv,
                                  qv0[k], qqqqco2, qco2[k], romray[k]);

    cs_real_t fnss = fnerir[k];
    cs_real_t tlsups;
    for (int ineb = k1; ineb <= kmray; ineb++)
      fnss = std::max(fnss, fnerir[ineb]);
    _cf_estimate(fnss, tlsups, 1e-3, 0.0, 0.0, fnss, qqqql);

    const cs_real_t term_1
      =   ((1.0 + fns*(tlsup - 1.0))
        * (dabcsup - dtvsup) + dul*tlsup*(tvsup - abcsup));
    const cs_real_t term_2
      =   ((1.0 + fnss*(tlsups - 1.0))
        * (dtvsups - dabcsups) - dul * tlsups * (tvsups - abcsups));

    /*! cooling rate in cloudy conditions
       Formula follows Ponnulakshmi and al. (2009)
       TODO should give back previous formula for dul=0 and tlsup=1 and tlsups=1
     */

    rayi[k] =  ctray * ( a1 - a2 + t4zt * term_1
                        -(1.0 - emis) * (a3 + t4zt * term_2));

    ckup[k] = (dabcinfe - dtvinfe) / (tvinfe - abcinfe);
    ckdown[k] = (dabcsup - dtvsup) / (tvsup - abcsup);
  }

  // Finalization: multiplication by sig
  for (int k = 0; k <= kmray; k++) {
    const int idx = k + ivertc * kmx;
    _atmo_1d_rad.iru[idx] = sig * ufir[k];
    _atmo_1d_rad.ird[idx] = sig * dfir[k];
  }

  // Compute Boundary conditions for the 3D
  // (Director diFfuse) Solar radiance
  // at the top of the CFD domain and the absorption coefficients

  if (   cs_rad_time_is_active() == true
      && cs_field_by_name_try("spectral_rad_incident_flux") != nullptr) {

    const int *bc_type = cs_glob_bc_type;
    /* dimension of radiative fields */
    int idx = -1;
    const int dim = cs_glob_rad_transfer_params->nwsgg;
    cs_real_t *cpro_ck_up = cs_field("rad_absorption_coeff_up")->val;
    cs_real_t *cpro_ck_down = cs_field("rad_absorption_coeff_down")->val;
    cs_real_t *bpro_rad_inc = cs_field("spectral_rad_incident_flux")->val;

    int rad_atmo_model = cs_glob_rad_transfer_params->atmo_model;
    // Direct solar radiation incident (for H2O band)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR)
      idx++;
    //  Direct solar radiation incident (for O3 band)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND)
      idx++;
    // Diffuse solar radiation incident
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR)
      idx++;
    // Diffuse solar radiation incident (for SUV O3 band)
    if (rad_atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND)
      idx++;

    // Infra Red radiation incident
    if (rad_atmo_model & CS_RAD_ATMO_3D_INFRARED) {
      idx++;
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        const int id = idx + face_id * dim;
        bpro_rad_inc[id] = 0.0;
        if (   bc_type[face_id] != CS_SMOOTHWALL
            && bc_type[face_id] != CS_ROUGHWALL) {
          cs_real_t var;
          cs_intprz(kmray + 1,
                    zqq,
                    dfir,
                    b_face_cog[face_id][2],
                    nullptr,
                    &var);
          bpro_rad_inc[id] = var;
        }
      }

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const int id = idx + c_id * dim;
        cs_real_t var;
        cs_intprz(kmx,
                  zray,
                  ckdown,
                  cell_cen[c_id][2],
                  nullptr,
                  &var);
        cpro_ck_down[id] = var;

        cs_intprz(kmx,
                  zray,
                  ckup,
                  cell_cen[c_id][2],
                  nullptr,
                  &var);
        cpro_ck_up[id] = var;
      }
    }
  }

  *qqvinf_p = qqvinf;

  CS_FREE(qc);
  CS_FREE(rov);
  CS_FREE(roc);
  CS_FREE(rol);
  CS_FREE(qv0);
  CS_FREE(qqc);
  CS_FREE(qql);
  CS_FREE(dt4);
  CS_FREE(tqq);
  CS_FREE(dz0);
  CS_FREE(ckup);
  CS_FREE(qqqc);
  CS_FREE(kliq);
  CS_FREE(qqql);
  CS_FREE(ufir);
  CS_FREE(dfir);
  CS_FREE(pspo);
  CS_FREE(qco2);
  CS_FREE(roco2);
  CS_FREE(qqco2);
  CS_FREE(ckdown);
  CS_FREE(qqqco2);
  CS_FREE(pspoqq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief free arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_finalize(void)
{
  CS_FREE(_atmo_1d_rad.xy);
  CS_FREE(_atmo_1d_rad.z);
  CS_FREE(_atmo_1d_rad.acinfe);
  CS_FREE(_atmo_1d_rad.dacinfe);
  CS_FREE(_atmo_1d_rad.aco2);
  CS_FREE(_atmo_1d_rad.aco2s);
  CS_FREE(_atmo_1d_rad.daco2);
  CS_FREE(_atmo_1d_rad.daco2s);
  CS_FREE(_atmo_1d_rad.acsup);
  CS_FREE(_atmo_1d_rad.acsups);
  CS_FREE(_atmo_1d_rad.dacsup);
  CS_FREE(_atmo_1d_rad.dacsups);
  CS_FREE(_atmo_1d_rad.tauzq);
  CS_FREE(_atmo_1d_rad.tauz);
  CS_FREE(_atmo_1d_rad.zq);
  CS_FREE(_atmo_1d_rad.ir_div);
  CS_FREE(_atmo_1d_rad.sol_div);
  CS_FREE(_atmo_1d_rad.iru);
  CS_FREE(_atmo_1d_rad.ird);
  CS_FREE(_atmo_1d_rad.solu);
  CS_FREE(_atmo_1d_rad.sold);
  CS_FREE(_atmo_1d_rad.albedo0);
  CS_FREE(_atmo_1d_rad.emissi0);
  CS_FREE(_atmo_1d_rad.temp0);
  CS_FREE(_atmo_1d_rad.theta0);
  CS_FREE(_atmo_1d_rad.qw0);
  CS_FREE(_atmo_1d_rad.p0);
  CS_FREE(_atmo_1d_rad.rho0);
  CS_FREE(_atmo_1d_rad.aerosols);
  CS_FREE(_atmo_1d_rad.fn);
  CS_FREE(_atmo_1d_rad.nc);
  CS_FREE(_atmo_1d_rad.qv);
  CS_FREE(_atmo_1d_rad.ql);
  CS_FREE(_atmo_1d_rad.qw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute radiative fluxes for the atmospheric model.
 *
 * Computes the source term for scalar equations from radiative forcing
 * (UV and IR radiative fluxes) with a 1D scheme.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_source_term(void)
{
  static bool initialized = false;

  const cs_time_step_t *ts = cs_glob_time_step;
  cs_atmo_1d_rad_t *atmo_1d_rad = cs_glob_atmo_1d_rad;

  /* Radiative fluxes are computed at the first call
     and then for a gvien frequency. */

  if ((ts->nt_cur % atmo_1d_rad->frequency) != 0 && initialized)
    return;

  /* 1. Initializations
   * ================== */

  const int kmx = atmo_1d_rad->nlevels_max;
  const int nlevels = atmo_1d_rad->nlevels;
  const int nvert = atmo_1d_rad->nvert;

  const cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  const int model_flag = cs_glob_physical_model_flag[CS_ATMOSPHERIC];
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_array<cs_real_t> temray(kmx);
  cs_array<cs_real_t> qvray(kmx);
  cs_array<cs_real_t> qlray(kmx);
  cs_array<cs_real_t> ncray(kmx);
  cs_array<cs_real_t> fneray(kmx);
  cs_array<cs_real_t> romray(kmx);
  cs_array<cs_real_t> preray(kmx);
  cs_array<cs_real_t> zproj(kmx);
  cs_array<cs_real_t> aeroso(kmx);
  cs_array_2d<cs_real_t> ttvert(nvert, kmx);
  cs_array_2d<cs_real_t> romvert(nvert, kmx);
  cs_real_t *tauz  = atmo_1d_rad->tauz;
  cs_real_t *tauzq  = atmo_1d_rad->tauzq;

  cs_span_2d<cs_real_t> qwvert(atmo_1d_rad->qw, nvert, kmx);
  cs_span_2d<cs_real_t> qlvert(atmo_1d_rad->ql, nvert, kmx);
  cs_span_2d<cs_real_t> qvvert(atmo_1d_rad->qv, nvert, kmx);
  cs_span_2d<cs_real_t> ncvert(atmo_1d_rad->nc, nvert, kmx);
  cs_span_2d<cs_real_t> fnvert(atmo_1d_rad->fn, nvert, kmx);
  cs_span_2d<cs_real_t> aevert(atmo_1d_rad->aerosols, nvert, kmx);
  cs_span_2d<cs_real_t> xyvert(atmo_1d_rad->xy, 3, nvert);
  cs_span_2d<cs_real_t> ir_div(atmo_1d_rad->ir_div, nvert, kmx);
  cs_span_2d<cs_real_t> sol_div(atmo_1d_rad->sol_div, nvert, kmx);

  cs_array_3d<cs_real_t> coords(nvert, kmx+1, 3);

  initialized = true;

  double heuray =   (double)at_opt->shour
                  + (double)at_opt->smin/60.
                  + at_opt->ssec/3600.;

  {
    cs_time_step_type_t idtvar = cs_glob_time_step_options->idtvar;
    if (idtvar == CS_TIME_STEP_CONSTANT || idtvar == CS_TIME_STEP_ADAPTIVE)
      heuray += ts->t_cur/3600.;
  }

  // Initialization:
  for (int k = 1; k < nlevels; k++) {
    preray[k] = 0.;
    temray[k] = 0.;
    qvray[k] = 0.;
    romray[k] = 0.;
    qlray [k] = 0.;
    ncray [k] = 0.;
    fneray[k] = 0.;
  }

  const cs_real_t *crom = cs_field("density")->val;
  const cs_real_t *cpro_tempc = cs_field("real_temperature")->val;

  const cs_real_t *cpro_pcliq = nullptr, *cvara_totwt = nullptr;
  const cs_real_t *cvara_ntdrp = nullptr, *nebdia = nullptr;

  if (model_flag == CS_ATMO_HUMID) {
    cpro_pcliq = cs_field("liquid_water")->val;
    cvara_totwt = cs_field("ym_water")->val_pre;
    cvara_ntdrp = cs_field("number_of_droplets")->val_pre;
    nebdia = cs_field("nebulosity_diag")->val;
  }

  /* 2.  Computing long-wave and short-wave radiative fluxes
   * ======================================================= */

  // Index of the bottom level
  const int k1 = 0;

  cs_span<cs_real_t> zq(atmo_1d_rad->zq, atmo_1d_rad->nlevels_max + 1);
  cs_span<cs_real_t> zray(atmo_1d_rad->z, atmo_1d_rad->nlevels_max);

  for (int ii = 0; ii < nvert; ii++) {  // (ixj) index

    cs_real_t xvert = xyvert(0, ii);
    cs_real_t yvert = xyvert(1, ii);

    // Addition of one level for solar radiation
    // TODO merge with 44km
    zq[kmx] = 16000.;

    // Coords are levels (faces in 3D) whereas zray is slice (cells in 3D)
    for (int k = 0; k < kmx+1; k++) {
      coords(ii, k, 0) = xvert;
      coords(ii, k, 1) = yvert;
      coords(ii, k, 2) = zq[k];
      if (k < kmx)
        zray[k] = 0.5 * (zq[k] + zq[k+1]);
    }
  }

  cs_interpol_grid_t *ig = cs_interpol_grid_by_id(at_opt->profiles_grid_id);

  if (ts->nt_cur == ts->nt_prev + 1) {
    cs_interpol_grid_init(ig, nvert*kmx, coords.data());  // face grid
  }

  /*
    Grid interpolation is refurbished to get a P0 interpolation
    and not a point-point interpolation.
    (i.e. we need to compute the mean of all code_saturne cells for a layer)
    A way could be to tag all nodes of the mesh is they belong to a ijk cell
    and then consider that a cell (or a face) belongs to a ijk if one node is
    in it (i.e. there is a non-zero intersection between the two).
    Then we can renormalize by the ratio "Volume(ijk)/SUM(cell_f_vol)".
  */

  // Interpolation from 3D to 1D in the fluid domain
  cs_interpol_field_on_grid(ig, cpro_tempc, ttvert.data());
  cs_interpol_field_on_grid(ig, crom, romvert.data());

  if (model_flag == CS_ATMO_HUMID) {
    // liquid content interpolation
    cs_interpol_field_on_grid(ig, cpro_pcliq, qlvert.data());
    // total water content interpolation
    cs_interpol_field_on_grid(ig, cvara_totwt, qwvert.data());

    // deduce vapor content interpolation
    for (int ii = 0; ii < nvert; ii++) {
      for (int k = 0; k < nlevels; k++)
        qvvert(ii, k) = qwvert(ii, k) - qlvert(ii, k);
    }

    cs_interpol_field_on_grid(ig, cvara_ntdrp, ncvert.data());
    cs_interpol_field_on_grid(ig, nebdia, fnvert.data());
  }

  // Interpolate soil (P0 interpolation), same for all verticals for the moment.

  cs_real_t *ground_albedo = atmo_1d_rad->albedo0;
  cs_real_t *ground_emissi = atmo_1d_rad->emissi0;
  cs_real_t *ground_temp = atmo_1d_rad->temp0;
  cs_real_t *ground_pressure = atmo_1d_rad->p0;
  cs_real_t *ground_totwat = atmo_1d_rad->qw0;
  cs_real_t *ground_density = atmo_1d_rad->rho0;

  if (at_opt->ground_model >= 1) {

    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

    int nbrsol;
    cs_lnum_t nfmodsol;
    const cs_lnum_t *elt_ids;

    const cs_zone_t *z = cs_boundary_zone_by_id(at_opt->ground_zone_id);
    cs_atmo_get_ground_zone(&nfmodsol, &nbrsol, &elt_ids);

    // Note: we use previous values of the ground to be coherent with
    // previous datasetting and what have seen the fluid.
    const cs_real_t *bvar_tempp = cs_field("ground_pot_temperature")->val_pre;
    const cs_real_t *bvar_total_water = cs_field("ground_total_water")->val_pre;

    const cs_real_t *bpro_albedo = cs_field("boundary_albedo")->val;
    const cs_real_t *bpro_emissi = cs_field("emissivity")->val;

    cs_real_t ground_mean_albedo = 0.;
    cs_real_t ground_mean_emissi = 0.;

    cs_real_t ground_mean_density = 0.;
    cs_real_t ground_mean_ttground  = 0.;
    cs_real_t ground_mean_totwat  = 0.;

    for (cs_lnum_t isol = 0; isol < nfmodsol; isol++) {
      cs_lnum_t face_id = elt_ids[isol];
      cs_real_t f_surf = b_face_surf[face_id];

      // Density: property of the boundary face:
      ground_mean_density += f_surf * crom[b_face_cells[face_id]];

      // Potential temperature for consistency with the code before.
      //TODO is it correct?
      ground_mean_ttground  += f_surf * (bvar_tempp[isol] - tkelvi);
      ground_mean_totwat  += f_surf * bvar_total_water[isol];

      ground_mean_albedo += f_surf * bpro_albedo[face_id];
      ground_mean_emissi += f_surf * bpro_emissi[face_id];
    }

    cs::parall::sum(ground_mean_density,
                    ground_mean_ttground,
                    ground_mean_totwat,
                    ground_mean_albedo,
                    ground_mean_emissi);

    ground_mean_density /= z->measure;
    ground_mean_ttground  /= z->measure;
    ground_mean_totwat  /= z->measure;
    ground_mean_albedo  /= z->measure;
    ground_mean_emissi  /= z->measure;

    // For now, all verticals have the same value.
    // TODO: automatic treatment for pressure ?

    for (int ii = 0; ii < nvert; ii++) {
      ground_albedo[ii] = ground_mean_albedo;
      ground_emissi[ii] = ground_mean_emissi;
      ground_temp[ii] = ground_mean_ttground;
      ground_totwat[ii] = ground_mean_totwat;
      ground_density[ii] = ground_mean_density;
    }

  } // ground model

  // Loop on the vertical array:

  const int distribution_model = at_opt->distribution_model;
  const int meteo_profile = at_opt->meteo_profile;
  const int met_1d_nlevels_d = at_opt->met_1d_nlevels_d;
  const int met_1d_ntimes = at_opt->met_1d_ntimes;
  const int nbmaxt = at_opt->met_1d_nlevels_max_t;
  const int nbmett = at_opt->met_1d_nlevels_t;

  const cs_real_t *z_temp_met = at_opt->z_temp_met;
  const cs_real_t *time_met = at_opt->time_met;
  const cs_real_t *hyd_p_met = at_opt->hyd_p_met;
  const cs_real_t *ttmet = at_opt->temp_met;
  const cs_real_t *qvmet = at_opt->qw_met;

  cs_real_t tausup = atmo_1d_rad->tausup;

  const cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  const cs_real_t rvsra = cs_glob_fluid_properties->rvsra;
  const cs_real_t abs_gz = cs::abs(cs_glob_physical_constants->gravity[2]);

  const cs_real_t p0 = cs_glob_fluid_properties->p0;
  const cs_real_t t0 = cs_glob_fluid_properties->t0;
  const double t_cur = ts->t_cur;

  for (int ii = 0; ii < nvert; ii++) {

    // FIXME the x, y position plays no role...
    // interpolation must be reviewed
    // cs_real_t xvert = xyvert(0, ii);
    // cs_real_t yvert = xyvert(1, ii);

    // Ground constants
    cs_real_t albedo = ground_albedo[ii];
    cs_real_t emis   = ground_emissi[ii];

    int imer1 = 0;

    // 2.1 Profiles used for the computation of the radiative fluxes
    // -------------------------------------------------------------

    // Loop over the variable defined until the top of the full domain.
    for (int k = 0; k < kmx; k++) {
      aeroso[k] = aevert(ii, k);
      if (model_flag == CS_ATMO_HUMID && distribution_model == 2) {
        qlray[k]  = qlvert(ii, k);
        ncray[k]  = ncvert(ii, k);
        fneray[k] = fnvert(ii, k);
      }
      else {
        // default values
        qlray[k] = 0.;
        ncray[k] = 0.;
        fneray[k] = 0.;
      }
    }

    // Ground variables
    // temray(0) = ground_temp[ii]
    // qvray(0)  = ground_totwat[ii]
    // romray(0) = ground_density[ii]
    // preray(0) = ground_pressure[ii]

    for (int k = 0; k < nlevels; k++) {

      temray[k] = ttvert(ii, k);
      qvray[k]  = qvvert(ii, k);
      romray[k] = romvert(ii, k);

      cs_real_t dummy;
      if (meteo_profile == 0)
        cs_atmo_profile_std(0., p0, t0, zray[k],
                            &preray[k], &dummy, &dummy);
      else if (meteo_profile == 1)
        preray[k] = cs_intprf(met_1d_nlevels_d, met_1d_ntimes,
                              z_temp_met, time_met, hyd_p_met, zray[k], t_cur);
      else {
        // TODO would be more coherent with an averaging of
        // "meteo_pressure" field.
        cs_atmo_profile_std(0., p0, t0, zray[k],
                            &preray[k], &dummy, &dummy);
      }

    }

    // Filling the additional levels
    // TODO do it before

    for (int k = nlevels; k < kmx; k++) {
      // Initialize with standard atmosphere above the domain.
      cs_atmo_profile_std(zray[nlevels-1],
                          preray[nlevels-1],
                          temray[nlevels-1]+tkelvi,
                          zray[k],
                          &preray[k],
                          &temray[k],
                          &romray[k]);
      temray[k] -= tkelvi;
      qvray[k] = 0.;
    }

    cs_user_atmo_1d_rad_prf(preray.data(),
                            temray.data(),
                            romray.data(),
                            qvray.data(),
                            qlray.data(),
                            ncray.data(),
                            aeroso.data());

    // Smoothing the temperature and humidity profile in the damping zone

    int ktamp = 0;
    if (meteo_profile == 1) {
      ktamp = cs::min(6, nlevels);
      for (int k = nlevels - ktamp; k < kmx; k++) {
        temray[k] = cs_intprf(nbmaxt, met_1d_ntimes, z_temp_met, time_met,
                              ttmet, zray[k], t_cur);
        qvray[k] = cs_intprf(nbmaxt, met_1d_ntimes, z_temp_met, time_met,
                             qvmet, zray[k], t_cur);
      }

      int icompt = 0;
      for (int k = nlevels-1; k > 0; k--) {
        icompt++;
        if (icompt <= 6) {
          cs_real_t zrac = 2.0 * (zray[k] - zray[nbmett-ktamp+2])
                              / (zray[nbmett-1] - zray[nbmett-ktamp-1]);
          cs_real_t w = (1. + tanh(zrac)) * 0.5;
          temray[k] = ttvert(ii, k) * (1. - w) + temray[k] * w;
          qvray[k]  = qvvert(ii, k) * (1. - w) + qvray[k] * w;
          qlray[k]  = qlvert(ii, k) * (1. - w) + qlray[k] * w;
        }
      }
    }

    // Clipping the humidity

    for (int k = 0; k < kmx; k++) {
      qvray[k] = cs::max(0., qvray[k]);
      qlray[k] = cs::max(0., qlray[k]);
    }

    // Computing pressure and density according to temperature
    // and qv profiles.

    for (int k = nlevels-ktamp; k < kmx; k++) {
      cs_real_t tmoy = 0.5 * (temray[k-1]+temray[k]) + tkelvi;
      cs_real_t rhum = rair*(1.+(rvsra-1.)*qvray[k]);
      cs_real_t rap = -abs_gz*(zray[k]-zray(k-1)) / (rhum*tmoy);
      preray[k] = preray[k-1] * exp(rap);
      romray[k] = preray[k] / ((temray[k] + tkelvi) * rhum);
    }

    // 2.2 Computing the radiative fluxes for the vertical
    // ---------------------------------------------------

    // Long-wave: InfraRed
    cs_atmo_1d_rad_compute_infrared(ii, k1, kmx-1, emis,
                                    tauzq, tauz, &tausup, zq,
                                    zray, temray, qvray, qlray,
                                    fneray, romray, preray, aeroso,
                                    ground_temp[ii], ground_pressure[ii],
                                    ir_div.sub_view(ii), ncray);

    // Short-wave: Sun
    cs_atmo_1d_rad_compute_solar(ii, k1, kmx-1, heuray,
                                 imer1, &albedo, tauzq, tausup, zq,
                                 zray, qvray, qlray,
                                 fneray, romray, preray, temray,
                                 sol_div.sub_view(ii), ncray);

    // FIXME: albedo may be modified by previous call, but the value
    // is not modified in the ground_albedo (e.g. albedo0) array.
    // If the value should only be modified locall, a "pass by value"
    // approach to the above function should be preffered
    // (this may be an artifact from the Fortran code).

  }

  atmo_1d_rad->tausup = tausup;

  int kmx_nvert = kmx*nvert;

  cs_array<int> cressm(kmx_nvert);
  cs_array<int> interp(kmx_nvert);
  cs_array_2d<cs_real_t> infrad(kmx_nvert, 3);

  cs_real_t h_cressman_radius = 1./8500;  // horizontal
  cs_real_t v_cressman_radius = 4./1.;    // vertical

  for (int ii = 0; ii < kmx_nvert; ii++) {
    cressm[ii] = 1;
    interp[ii] = 1;
    infrad(ii, 0) = h_cressman_radius;
    infrad(ii, 1) = h_cressman_radius;
    infrad(ii, 2) = v_cressman_radius;
  }

  // Map Infra Red 1D (rayi) on the structure idrayi
  cs_measures_set_map_values(cs_measures_set_by_id(at_opt->infrared_1D_profile),
                             kmx_nvert,
                             cressm.data(),
                             interp.data(),
                             coords.data(),
                             ir_div.data(),
                             infrad.data());

  // Map Sun 1D (rayst) on the structure idrayst
  cs_measures_set_map_values(cs_measures_set_by_id(at_opt->solar_1D_profile),
                             kmx_nvert,
                             cressm.data(),
                             interp.data(),
                             coords.data(),
                             sol_div.data(),
                             infrad.data());
}

/*----------------------------------------------------------------------------*/
