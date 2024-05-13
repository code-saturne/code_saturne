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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_atmo_aerosol_ssh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Beta function: beta(x,y) = gamma(x)*gamma(y)/gamma(x+y)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_beta(const cs_real_t x,
      const cs_real_t y)

{
  cs_real_t beta;

  beta = tgamma(x)*tgamma(y)/tgamma(x+y);

  return beta;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute hypergeometric function for |x| < 1 using a series
 *        (cf. for the definition of this function, see for example:
 *        http://mathworld.wolfram.com/hypergeometricfunction.html)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_hypser(cs_real_t  a,
        cs_real_t  b,
        cs_real_t  c,
        cs_real_t  x)
{
  cs_real_t hypser;

  const int maxiter = 10000;
  const cs_real_t error = 1.e-08;

  if (cs_math_fabs(x) >= 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s\n"
              "The x parameter should verify |x| < 1,  x = %12.5e",
              __func__, x);

  cs_real_t fac = 1;
  cs_real_t temp = fac;
  cs_real_t aa  = a;
  cs_real_t bb  = b;
  cs_real_t cc  = c;

  for (int nn = 1; nn < maxiter; nn++) {
    fac = ((aa*bb)/cc)*fac;
    fac = fac*x/nn;
    hypser = fac + temp;

    if (cs_math_fabs(hypser - temp) <= error)
      return hypser;

    aa += 1;
    bb += 1;
    cc += 1;
    temp = hypser;
  }

  return hypser;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Hypergeometric function
 *        (see http://mathworld.wolfram.com/hypergeometricfunction.html
 *        for definition)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_hypgeo(cs_real_t  a,
        cs_real_t  b,
        cs_real_t  c,
        cs_real_t  x)
{
  cs_real_t hypgeo;

  const cs_real_t pp = 0.1;

  /* Initialization */

  const cs_real_t gammaa   = tgamma(a);
  const cs_real_t gammab   = tgamma(b);
  const cs_real_t gammac   = tgamma(c);
  const cs_real_t gammabma = tgamma(b-a);
  const cs_real_t gammacma = tgamma(c-a);
  const cs_real_t gammaamb = tgamma(a-b);
  const cs_real_t gammacmb = tgamma(c-b);

  /* Compute hypergeometric function by convergent series for |x|<1 */

  if (x >= (pp - 1)) {
    hypgeo = _hypser(a, b, c, x);
                     }
  else if (x >= -1-pp) {
    cs_real_t y1 = _hypser(a, a+1.-c, a+1.-b, 1./x);
    cs_real_t y2 = _hypser(b, b+1.-c, b+1.-a, 1./x);
    hypgeo = (gammac*gammabma*y1*pow(-1./x, a))/(gammab*gammacma)
           + (gammac*gammaamb*y2*pow(-1./x, b))/(gammaa*gammacmb);
  }
  else {
    cs_real_t y1 = _hypser(a, a+1.-c, a+1.-b, 1./(-1.-pp));
    cs_real_t y2 = _hypser(b, b+1.-c, b+1.-a, 1./(-1.-pp));
    cs_real_t hyp1 = (gammac*gammabma*y1*pow(-1./(-1.-pp), a))/(gammab*gammacma)
                   + (gammac*gammaamb*y2*pow(-1./(-1.-pp), b))/(gammaa*gammacmb);
    cs_real_t hyp2 = _hypser(a, b, c, -1.+pp);
    hypgeo = hyp1 + (x - (-1.-pp))*(hyp2 - hyp1)/(2.*pp);
  }

  return hypgeo;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_initialize();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_finalize();

  BFT_FREE(cs_glob_atmo_chemistry->aero_conc_file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with aerosol concentrations
 *        and numbers from the external aerosol code.
 *
 * \param[out]  array  aerosol concentrations and numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_aero(cs_real_t  *array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_get_aero(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 *
 * \param[out]  array  gas concentrations
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t  *array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_get_gas(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_time_advance();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute aerosol cloud droplets nucleation when using the atmospheric
 * humid model using a microphysical model.
 *
 * It is taken into account as an additional step split from advection-diffusion
 * equation, hence the droplet number is first clipped if necessary.
 *
 * \param[out]  nc      droplet number (scalar) in 1/cm**3
 * \param[in]   rom     density of air in kg/m**3
 * \param[in]   qldia   mass fraction of liquid water in kg/kg
 * \param[in]   pphy    true pressure in pascals
 * \param[in]   refrad  radiative cooling
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_nuclea(cs_real_t         *nc,
                       const cs_real_t   *rom,
                       const cs_real_t   *qldia,
                       const cs_real_t   *pphy,
                       const cs_real_t   *refrad)
{
  /* Initialization
   * -------------- */

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;

  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;
  const cs_real_t *tempc = cs_field_by_name("real_temperature")->val;

  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t cp = phys_pro->cp0;
  cs_real_t rvsra = phys_pro->rvsra;
  cs_real_t clatev = phys_pro->clatev;
  cs_real_t rair = phys_pro->r_pg_cnst;

  const cs_real_t rhowater = 1000; // kg/m^3
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t nuc = 0.;
  cs_real_t fbeta = 0.;
  cs_real_t constc = 0.;
  cs_real_t constk = 0.;
  cs_real_t constmu = 0;
  cs_real_t constbeta = 0;
  cs_real_t rvap = rair*rvsra;
  const int modnuc = cs_glob_atmo_option->nucleation_model;

  if (modnuc == 1) {
    /* Constants for the model of Pruppacher and Klett (1997)
     * Case of Buffalo, New-York, Kocmond [1965] */
    constc = 3500.;
    constk = 0.9;

    /* fbeta : beta function (constk/2, 3/2) */
    fbeta = _beta(constk/2., 1.5);
  }
  else if (modnuc == 2) {
    /* Constants for model of Cohard and Pinty (1998)
     * (general case) */
    constc    = 3270.;
    constk    = 1.56;
    constmu   = 0.7;
    constbeta = 136.;

    /* Polluted Air Case. See Hudson and Li (1995)
     * constc    = 1865. ! 8700.d0 if Cohard and Pinty used on PARISFOG
     * constk    = 0.86
     * constmu   = 1.50
     * constbeta = 6.80 */
    fbeta = _beta(constk/2., 1.5);
  }

  /* Computation new Nc field (in cm^-3)
   * ----------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (qldia[c_id] <= cs_math_epzero) {

      nc[c_id] = 0.;

    }
    else { /* qldia[c_id] > cs_math_epzero */

      /* New droplets are created by nucleation if vel_z > 0 or refrad < 0,
       * then if number of new droplets is superior to nc,
       * the difference is added */

      if ((vel[c_id][2] > cs_math_epzero) || (refrad[c_id] < cs_math_epzero)) {

        cs_real_t tempk = tempc[c_id] + tkelvi;
        cs_real_t esat = cs_air_pwv_sat(tempc[c_id]);
        cs_real_t kka = ((5.69 + 0.017*tempc[c_id])/0.239)*1.e-3;
        cs_real_t ddv = 0.211*pow(tempk/tkelvi, 1.94)*(101325./pphy[c_id])*1.e-4;
        cs_real_t aa1 =   0.622*clatev*9.81/(rair*cp*cs_math_pow2(tempk))
                        - 9.81/(rair*tempk);
        cs_real_t aa2 =   rair*tempk/(0.622*esat)
                        + (0.622*cs_math_pow2(clatev))/(tempk*pphy[c_id]*cp);
        cs_real_t aa3 = 1.0/((1000.*rvap*tempk)/(esat*ddv)
                        + clatev*1000.*(clatev/(tempk*rvap) -1.0)/(kka*tempk));
        cs_real_t aa4 = -0.622*clatev/(rair*cs_math_pow2(tempk));

        if ((aa1*vel[c_id][2]+aa4*refrad[c_id]) > cs_math_epzero) {

          // 1.1  Model of Pruppacher and Klett (1997)
          if (modnuc == 1) {
            const cs_real_t cons2 = pow(constc, 2./(constk+2.));
            const cs_real_t cons3 = pow(0.01*(  aa1*vel[c_id][2]
                                              + aa4*refrad[c_id]), 3./2.);
            const cs_real_t cons4
              = (2.*cs_math_pi*rhowater*aa2*pow(aa3, 3./2.)*constk*fbeta);
            nuc = cons2*pow(cons3/cons4, constk/(constk + 2.));
          }

          // 1.2  Modele of Cohard and Pinty (1998)
          else if (modnuc == 2) {
            /* compute max supersaturation: sursat */
            cs_real_t yy = 1.0;
            cs_real_t sursat = 0;
            cs_real_t tmpsur = 0;
            for (int ii = 0; ii < 20; ii++) {
              tmpsur = sursat;
              yy = _hypgeo(constmu, constk/2, constk/2 + 3./2,
                        -constbeta*cs_math_pow2(sursat));
              const cs_real_t cons1
                = pow(0.01*(aa1*vel[c_id][2] + aa4*refrad[c_id]), 3./2);
              const cs_real_t cons2
                = (2*constk*constc*cs_math_pi*rhowater*aa2*fbeta
                   *pow(aa3, 3.0/2));

              sursat = pow((cons1/cons2)/yy, 1./(constk + 2.));
            }

            if (cs_math_fabs(tmpsur - sursat) > 1.e-2)
              bft_error(__FILE__, __LINE__, 0,
              " WARNING: Maximum saturation has not converged\n"
                "Residue = %12.5e", cs_math_fabs(tmpsur - sursat));

            nuc = constc*pow(sursat, constk)*_hypgeo(constmu,
                                                     constk/2,
                                                     constk/2 + 1.,
                                                     -constbeta*
                                                     cs_math_pow2(sursat));
            if (nuc < 0.)
              bft_error(__FILE__, __LINE__, 0,
              _(" ERROR: Cohard and Pindy model (1998).\n"
                " The nucleation is negative."));
          }

          /* 1.3  Model of Abdul-Razzak and al. (1998)
           * This model requires a fine knowledge of aerosols chemistry. */
          else if (modnuc == 3) {

            /* Warning: this set of constants fits for PARISFOG case.
             * They were determined using fine data on aerosols measured
             * during PARISFOG POI 13 */
            cs_real_t sigaero3 = 1.3;
            cs_real_t sigaero2 = 1.59;
            cs_real_t sigaero1 = 1.570;
            cs_real_t ntotal1 = 8700.;
            cs_real_t ntotal2 = 8300.;
            cs_real_t ntotal3 = 1000.;
            cs_real_t rayonm3 = 0.4e-6;
            cs_real_t rayonm2 = 0.055e-6;
            cs_real_t rayonm1 = 0.0165e-6;

            /* surface tension water-steam [N/m], tempk [K] */
            const cs_real_t tauvw = 0.0750;

            /* fractions in %
             * fisrt component: sulfate
             * second component: nitrate
             * third component: black carbon
             * fourth component: organic carbon */
            const cs_real_t fraero[4] = {0.25, 0.2, 0.164, 0.386};

            /* number of ions produced by dissociation of a salt molecule
               in water (if NaCl, nbion=2) */
            const cs_real_t nbion[4] = {3., 2., 0., 1.};

            /* osmotic coefficient */
            const cs_real_t coefosm[4] = {1., 1., 1., 1.};

            /* mass fraction of of soluble
               substance in aerosol mix: msolu/mmix */
            const cs_real_t frasolu[4] = {1., 1., 0., 0.1};

            /* molecular mass of water */
            const cs_real_t mmh2o = 18.e-3;

            /* molecular mass of aerosol components */
            const cs_real_t mmaero[4] = {132e-3, 80e-3, 250e-3, 250e-3};

            /* density of aerosol !FIXME to check */
            const cs_real_t rhoaero[4] = {1.77e3, 1.77e3, 1.77e3, 1.77e3};

            /* coefficient for Kelvin effect: coefa [/m] */
            const cs_real_t coefa = 2*tauvw/(rhowater*rvap*tempk);

            /* FIXME
             * cf Pruppacher and Klett 1997 (reprinted correction 2000) (6-28)
             * coefa = 3.3d-7 / tempk */

            /* coefficient for Raoult effect: coefb [-] */
            cs_real_t numcb = 0;
            cs_real_t dencb = 0;
            for (int ii = 0; ii < 4; ii++) {
              numcb += fraero[ii]*nbion[ii]*coefosm[ii]*frasolu[ii]/mmaero[ii];
              dencb += fraero[ii]/rhoaero[ii];
            }
            const cs_real_t coefb = mmh2o*numcb/(dencb*rhowater);

            /* supersaturation [-] */
            cs_real_t sursat1 = (2./sqrt(coefb))*(pow(coefa/3./rayonm1, 1.5));
            cs_real_t sursat2 = (2./sqrt(coefb))*(pow(coefa/3./rayonm2, 1.5));
            cs_real_t sursat3 = (2./sqrt(coefb))*(pow(coefa/3./rayonm3, 1.5));

            /* standard deviation function */
            cs_real_t foncf1 = 0.5*exp(2.5*cs_math_pow2(log(sigaero1)));
            cs_real_t foncf2 = 0.5*exp(2.5*cs_math_pow2(log(sigaero2)));
            cs_real_t foncf3 = 0.5*exp(2.5*cs_math_pow2(log(sigaero3)));
            cs_real_t foncg1 = 1.0 + 0.25*(log(sigaero1));
            cs_real_t foncg2 = 1.0 + 0.25*(log(sigaero2));
            cs_real_t foncg3 = 1.0 + 0.25*(log(sigaero3));

            /* coefficient corresponding to vertical velocity |aa1| */
            aa1 =   0.622*clatev*9.81/(rair*cp*cs_math_pow2(tempk))
                  - 9.81/(rair*tempk);

            /* coefficient corresponding to liquid water condensation |aa2| */
            aa2 =   rair*tempk/(0.622*esat)+(0.622*cs_math_pow2(clatev))
                  / (tempk*pphy[c_id]*cp);

            /* coefficient corresponding to the droplet growth |aa3| */
            ddv = 0.211*pow(tempk/tkelvi, 1.94)*(101325./pphy[c_id])*1.e-4;
            rvap = rvsra*rair;
            aa3 =   1/((1000*rvap*tempk)/(esat*ddv)
                  + clatev*1000*(clatev/(tempk*rvap)-1.0)/(kka*tempk));

            /* coefficient corresponding to infrared radiation |aa4| */
            aa4 = -0.622*clatev/(rair*cs_math_pow2(tempk));

            /* alphaV/G */
            cs_real_t coefavg = (aa1*vel[c_id][2] + aa4*refrad[c_id])/aa3;

            /* coefficient zeta */
            cs_real_t coefzeta = 2*coefa/3*sqrt(coefavg);

            /* coefficient eta - Ntot [m^(-3)] */
            const cs_real_t coefeta1
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal1*1.e6);
            const cs_real_t coefeta2
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal2*1.e6);
            const cs_real_t coefeta3
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal3*1.e6);

            /* max supersaturation (smax) */
            cs_real_t x_tp1 = foncf1*pow(coefzeta/coefeta1, 1.5);
            cs_real_t x_tp2 = pow((sursat1*sursat1)/(coefeta1+3*coefzeta), 0.75);
            cs_real_t smax1 = (x_tp1 + foncg1*x_tp2)/sursat1/sursat1;

            x_tp1 = foncf2*pow(coefzeta/coefeta2, 1.5);
            x_tp2 = pow((sursat2*sursat2)/(coefeta2+3*coefzeta), 0.75);
            cs_real_t smax2 = (x_tp1 + foncg2*x_tp2)/sursat2/sursat2;

            x_tp1 = foncf3*(pow(coefzeta/coefeta3, 1.5));
            x_tp2 = pow((sursat3*sursat3)/(coefeta3+3*coefzeta), 0.75);
            cs_real_t smax3 = (x_tp1 + foncg3*x_tp2)/sursat3/sursat3;

            cs_real_t smax = 1.0/sqrt(smax3 + smax2 + smax1);

            if ((smax > 1.0) || (sursat1 > 1.0) ||
                (sursat2 > 1.0) || (sursat3 > 1.0)) {
              bft_error(__FILE__, __LINE__, 0,
                        _("%s: Abdul-Razzak and al. model (1998).\n"
                          " Negative supersaturation."), __func__);
            }

            cs_real_t f3_ov_sqrt2 = 3.*sqrt(2.);
            cs_real_t ugauss1 =   (2.*log(sursat1/smax))
                                / (f3_ov_sqrt2*log(sigaero1));
            cs_real_t ugauss2 =   (2.*log(sursat2/smax))
                                / (f3_ov_sqrt2*log(sigaero2));
            cs_real_t ugauss3 =   (2.*log(sursat3/smax))
                                / (f3_ov_sqrt2*log(sigaero3));

            const cs_real_t nuc1 = 0.5 * ntotal1 * (1.0-erf(ugauss1));
            const cs_real_t nuc2 = 0.5 * ntotal2 * (1.0-erf(ugauss2));
            const cs_real_t nuc3 = 0.5 * ntotal3 * (1.0-erf(ugauss3));

            nuc = nuc1 + nuc2 + nuc3;

          } /* end nucleation model */

          /* 2. Compute difference */
          nuc = cs_math_fmax(nuc - nc[c_id], 0.0);
        }
        else {
          nuc = 0.0;
        }

        nc[c_id] += nuc;
      }

      /* 3. if qc > 0, w <= 0 and nc = 0,
       * we impose nc > 0 so that the mean volume radius = 10 microns */
      else if (nc[c_id] < cs_math_epzero) {

        nc[c_id] =   1.e-6*(3.*rom[c_id]*qldia[c_id])
                   / (4.*cs_math_pi*rhowater*cs_math_pow3(10e-6));

      } /* end Vel_z > cs_math_epzero*/

    } /* end qldia[c_id] > cs_math_epzero */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
