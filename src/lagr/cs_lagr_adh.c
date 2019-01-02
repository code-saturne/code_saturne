/*============================================================================
 * Methods for particle adhesion forces
 *============================================================================*/

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

/*============================================================================
 * Functions dealing with lagrangian adhesion forces
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_math.h"

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_physical_constants.h"
#include "cs_random.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_clogging.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_adh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro declarations
 *============================================================================*/

/*============================================================================
 * Local structure declarations
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Cut-off distance for adhesion forces (assumed to be the Born distance) */
static const cs_real_t  _d_cut_off = 1.65e-10;

/* Boltzmann constant */
static const double _k_boltz = 1.38e-23;

/* Free space permittivity */
static const cs_real_t _free_space_permit = 8.854e-12;

/* Faraday constant */
static const cs_real_t _faraday_cst = 9.648e4;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* VDW interaction between a sphere and a plate                               */
/*----------------------------------------------------------------------------*/

static void
_vdwsp(cs_real_t  cstham,
       cs_real_t  lambda_vdw,
       cs_real_t *distp,
       cs_real_t *rpart,
       cs_real_t *var)
{
  if (*distp < lambda_vdw / 2.0 / cs_math_pi)
    *var = - cstham * *rpart / (6.0 * *distp)
           * (1.0 / (1.0 + 14.0 * *distp / lambda_vdw + 5.0 * cs_math_pi / 4.9 * pow (*distp, 3.0)
                   / lambda_vdw / pow (*rpart, 2.0)));
  else
    *var = cstham * (2.45 / 60.0 / cs_math_pi * lambda_vdw
                     * (  (*distp - *rpart) / pow (*distp, 2)
                        - (*distp + 3.0 * *rpart) / pow ((*distp + 2.0 * *rpart), 2.0))
                     - 2.17 / 720.0 / pow (cs_math_pi, 2) * pow (lambda_vdw, 2.0)
                       * ((*distp - 2.0* *rpart) / pow (*distp, 3.0) - (*distp + 4.0 * *rpart) / pow ((*distp + 2.0 * *rpart), 3.0))
                     + 0.59 / 5040.0 / pow (cs_math_pi, 3.0) * pow (lambda_vdw, 3.0)
                       * ((*distp - 3.0 * *rpart) / pow (*distp, 4.0) - (*distp + 5.0 * *rpart) / pow ((*distp + 2.0 * *rpart), 4.0)));

}

/*----------------------------------------------------------------------------*/
/* VDW interaction between two spheres                                        */
/*----------------------------------------------------------------------------*/

static void
_vdwss(cs_real_t  cstham,
       cs_real_t  lambda_vdw,
       cs_real_t *distcc,
       cs_real_t *rpart1,
       cs_real_t *rpart2,
       cs_real_t *var)
{

  *var = - cstham * *rpart1 * *rpart2
        / (6.0 * (*distcc - *rpart1 - *rpart2) * (*rpart1 + *rpart2))
        * (  1.0
           -   5.32 * (*distcc - *rpart1 - *rpart2) / lambda_vdw
             * log (1.0 + lambda_vdw / (*distcc - *rpart1 - *rpart2) / 5.32));

}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

static void
_edlsp(cs_real_t *distp,
       cs_real_t *rpart,
       cs_real_t  phi_p,
       cs_real_t  phi_s,
       cs_real_t  tempf,
       cs_real_t *var)
{
  cs_lagr_physico_chemical_t *lag_pc = cs_glob_lagr_physico_chemical;
  cs_real_t charge  = 1.6e-19;

  cs_real_t ldebye = pow(   2000.0 * pow (_faraday_cst,2) * lag_pc->fion
                / (lag_pc->epseau * _free_space_permit * cs_physical_constants_r * tempf), -0.5);

  /* Reduced zeta potential    */
  cs_real_t lphi1  = lag_pc->valen * charge * phi_p / _k_boltz / tempf;
  cs_real_t lphi2  = lag_pc->valen * charge * phi_s / _k_boltz / tempf;

  /* xtended reduced zeta potential */
  /*  (following the work from Ohshima et al, 1982, JCIS, 90, 17-26)   */
  /* or the particle */

  cs_real_t tau    = *rpart / ldebye;
  lphi1  =   8.0 * tanh (lphi1 / 4.0)
           / (1.0 + sqrt (1.0 - (2.0 * tau + 1.0) / pow(tau + 1.0, 2.0) * pow(tanh(lphi1 / 4.0), 2.0)));

  /* For the plate   */

  lphi2  = 4.0 * tanh (lphi2 / 4.0);

  /* alculation for the EDL force   */

  cs_real_t alpha = sqrt ((*distp + *rpart) / *rpart) + sqrt (*rpart / (*distp + *rpart));

  cs_real_t omega1 = pow (lphi1, 2.0) + pow (lphi2, 2.0) + alpha * lphi1 * lphi2;

  cs_real_t omega2 = pow (lphi1, 2.0) + pow (lphi2, 2.0) - alpha * lphi1 * lphi2;

  cs_real_t gamma  = sqrt (*rpart / (*distp + *rpart)) * exp ( -1.0 / ldebye * *distp);

  *var   =  2.0 * cs_math_pi * lag_pc->epseau * _free_space_permit
          * pow(_k_boltz * tempf / charge, 2.0) * *rpart
          * (*distp + *rpart) / (*distp + 2.0 * *rpart)
          * (omega1 * log (1.0 + gamma) + omega2 * log (1.0 - gamma));

}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

static void
_edlss(cs_real_t *distcc,
       cs_real_t *rpart1,
       cs_real_t *rpart2,
       cs_real_t  phi1,
       cs_real_t  phi2,
       cs_real_t  tempf,
       cs_real_t *var)
{
  cs_lagr_physico_chemical_t *lag_pc = cs_glob_lagr_physico_chemical;
  cs_real_t charge = 1.6e-19;

  cs_real_t ldebye = pow(   2000.0 * pow (_faraday_cst,2) * lag_pc->fion
                / (lag_pc->epseau * _free_space_permit * cs_physical_constants_r * tempf), -0.5);

  /* Reduced zeta potential    */

  cs_real_t lphi1  = lag_pc->valen * charge * phi1 / _k_boltz / tempf;
  cs_real_t lphi2  = lag_pc->valen * charge * phi2 / _k_boltz / tempf;

  /* xtended reduced zeta potential */
  /*  (following the work from Ohshima et al, 1982, JCIS, 90, 17-26)   */
  /* For the first particle    */

  cs_real_t tau = *rpart1 / ldebye;
  lphi1 =   8.0 * tanh (lphi1 / 4.0)
           / (1.0 + sqrt (1.0 - (2.0 * tau + 1.0) / pow(tau + 1, 2.0) * pow(tanh(lphi1 / 4.0), 2.0)));

  /* For the second particle   */

  tau = *rpart2 / ldebye;
  lphi2 =  8.0 * tanh (lphi2 / 4.0)
         / (1.0 + sqrt (1.0 - (2.0 * tau + 1.0) / pow ((tau + 1.0), 2.0) * pow (tanh (lphi2 / 4.0), 2.0)));

  /* Calculation of the EDL force   */

  cs_real_t alpha  =  sqrt(*rpart2 * (*distcc - *rpart2) / (*rpart1 * (*distcc - *rpart1)))
                    + sqrt(*rpart1 * (*distcc - *rpart1) / (*rpart2 * (*distcc - *rpart2)));

  cs_real_t omega1 = pow (lphi1, 2.0) + pow (lphi2, 2.0) + alpha * lphi1 * lphi2;

  cs_real_t omega2 = pow (lphi1, 2.0) + pow (lphi2, 2.0) - alpha * lphi1 * lphi2;

  cs_real_t gamma  =  sqrt (*rpart1 * *rpart2 / (*distcc - *rpart1) / (*distcc - *rpart2))
                    * exp (1.0 / ldebye * (*rpart1 + *rpart2 - *distcc));

  *var =  2.0 * cs_math_pi * lag_pc->epseau * _free_space_permit
        * pow ((_k_boltz * tempf / charge), 2.0) * *rpart1 * *rpart2
        * (*distcc - *rpart1) * (*distcc - *rpart2)
        / (*distcc * (*distcc * (*rpart1 + *rpart2) - pow (*rpart1, 2.0) - pow (*rpart2, 2.0)))
        * (omega1 * log (1.0 + gamma) + omega2 * log (1.0 - gamma));
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the adhesion force and adhesion energy
 *
 * \param[in]  ip               particle number
 * \param[in]  tempf            thermal scalar value at current time step
 * \param[out] adhesion_energ   particle adhesion energy
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_adh(cs_lnum_t   ip,
            cs_real_t   tempf,
            cs_real_t  *adhesion_energ)
{
  /* ================================================================  */
  /* 0.    initialization */
  /* ================================================================  */

  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;
  unsigned char *part = p_set->p_buffer + p_am->extents * ip;
  cs_lagr_physico_chemical_t *lag_pc = cs_glob_lagr_physico_chemical;

  /*     step = step used to calculate the adhesion force following    */
  /*                         F = U(_d_cut_off+step)-U(_d_cut_off-step)/(2*step)     */

  cs_real_t denasp = cs_glob_lagr_reentrained_model->denasp;
  cs_real_t rayasg = cs_glob_lagr_reentrained_model->rayasg;
  cs_real_t rayasp = cs_glob_lagr_reentrained_model->rayasp;
  cs_real_t espasg = cs_glob_lagr_reentrained_model->espasg;
  cs_real_t modyeq = cs_glob_lagr_reentrained_model->modyeq;

  cs_real_t step = 1e-11;
  cs_real_t scovap = denasp * cs_math_pi * pow (rayasp, 2);
  cs_real_t scovag = cs_math_pi * pow (rayasg, 2) / pow (espasg, 2);

  cs_real_t dismin;
  /* ================================================================  */
  /* 3.    calculation of the adhesion force  */
  /* ================================================================  */

  /*     determination of the number of contacts with asperities  */
  /*     =======================================================  */

  /* *************************************/
  /* Number of large-scale asperities    */
  /* *************************************/

  cs_real_t rpart = 0.5 * cs_lagr_particle_get_real(part, p_am, CS_LAGR_DIAMETER);

  cs_real_t nmoyag = (2.0 * rpart + rayasg) / rayasg * scovag;

  if (nmoyag > 600.0) {

    cs_lnum_t tmp = 1;
    cs_real_t rtmp;

    cs_random_normal(1, &rtmp);

    tmp = (int)nmoyag + sqrt(nmoyag) * rtmp;
    tmp = CS_MAX(0, tmp);

    cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, tmp);

  }
  else {

    cs_lnum_t ntmp;

    cs_random_poisson(1, nmoyag, &ntmp);

    cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, ntmp);

  }

  cs_lnum_t nbasg = 0;
  if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) > 1) {

    nmoyag =  1.0
            +  2.0 * _d_cut_off * (2.0 * rpart + 2.0 * rayasg + 4.0 * _d_cut_off)
             / pow (rayasg, 2) * scovag;

    if (nmoyag > 600.0) {

      cs_real_t rtmp;

      cs_random_normal(1, &rtmp);

      nbasg = (int)nmoyag + sqrt (nmoyag) * rtmp;
      nbasg = CS_MAX(0, nbasg);

    }
    else {

      cs_lnum_t ntmp;

      cs_random_poisson(1, nmoyag, &ntmp);

      nbasg = ntmp;

    }

    nbasg = CS_MAX(1, nbasg);

  }
  else {

    nbasg = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES);

  }

  /* *******************************/
  /* Nb of small-scale asperities  */
  /* *******************************/

  /* 1st case: no large-scale asperities*/
  cs_lnum_t nbasp = 0;
  if (nbasg == 0) {

    cs_real_t nmoyap = (2.0 * rpart + rayasp) / rayasp * scovap;

    if (nmoyap > 600.0) {

      cs_lnum_t tmp = 1;
      cs_real_t rtmp;

      cs_random_normal(1, &rtmp);

      tmp = (int)nmoyap + sqrt(nmoyap) * rtmp;
      tmp = CS_MAX(0, tmp);

      cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, tmp);

    }
    else {

      cs_lnum_t ntmp;

      cs_random_poisson(1, nmoyap, &ntmp);

      cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, ntmp);

    }

    if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES) > 1) {

      nmoyap =  1
              +  2.0 * _d_cut_off * (2.0 * rpart + 2.0 * rayasp + 4.0 * _d_cut_off)
               / pow (rayasp, 2) * scovap;

      if (nmoyap > 600.0) {

        cs_real_t rtmp;

        cs_random_normal(1, &rtmp);

        nbasp = (int)nmoyap + sqrt (nmoyap) * rtmp;
        nbasp = CS_MAX(0, nbasp);

        cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, nbasp);

      }
      else {

        cs_lnum_t ntmp;

        cs_random_poisson(1, nmoyap, &ntmp);
        nbasp = ntmp;

      }

      nbasp = CS_MAX(1, nbasp);

    }
    else {

      nbasp = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES);

    }

    /* Determination of the minimal distance between the particle and the plate */

    dismin = rayasp * CS_MIN (1.0, cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES));

  }
 /* 2nd case: contact with large-scale asperities */
  else {

    cs_real_t paramh = 0.5 * (2.0 * rpart + rayasp) * rayasp / (rpart + rayasg);
    cs_real_t nmoyap = paramh * (2 * rayasg - paramh) / pow (rayasp, 2) * scovap;

    if (nmoyap > 600.0) {

      cs_lnum_t tmp = 1;
      cs_real_t rtmp;

      cs_random_normal(1, &rtmp);

      tmp = (int)nmoyap + sqrt(nmoyap) * rtmp;
      tmp = CS_MAX(0, tmp);

      cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, tmp);

    }
    else {

      cs_lnum_t ntmp;

      cs_random_poisson(1, nmoyap, &ntmp);

      cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, ntmp);

    }

    if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES) > 1) {

      paramh =  0.5 * (2.0 * rpart + 2 * rayasp + 4.0 * _d_cut_off)
              * 2.0 * _d_cut_off / (rpart + rayasg + rayasp + _d_cut_off);
      nmoyap = 1 + paramh * (2 * rayasg - paramh) / pow (rayasp, 2) * scovap;

      if (nmoyap > 600.0) {

        cs_real_t rtmp;

        cs_random_normal(1, &rtmp);

        nbasp = (int)nmoyap + sqrt (nmoyap) * rtmp;
        nbasp = CS_MAX(0, nbasp);

      }
      else {

        cs_lnum_t ntmp;

        cs_random_poisson(1, nmoyap, &ntmp);

        nbasp = ntmp;

      }

      nbasp = CS_MAX(1, nbasp);

    }
    else
      nbasp = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES);

    /* Mutliple contacts with large scale asperities?     */
    nbasp = nbasp * nbasg;
    cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES,
                               cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES)
                               *cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES));

    /* Determination of the minimal distance between the particle and the plate */
    dismin = rayasp * CS_MIN(1.0, nbasp * 1.0) + rayasg * CS_MIN (1.0, nbasg * 1.0);

  }

  /* Sum of {particle-plane} and {particle-asperity} energies     */

  cs_real_t udlvor[2];
  /* Interaction between the sphere and the plate  */
  for (cs_lnum_t np = 0; np < 2; np++) {

    udlvor[np] = 0.0;
    cs_real_t distp = dismin + _d_cut_off + step * (3 - 2 * (np+1));

    cs_real_t uvdwsp, uedlsp;
    _vdwsp(lag_pc->cstham, lag_pc->lambda_vdw, &distp, &rpart, &uvdwsp);
    _edlsp(&distp, &rpart, lag_pc->phi_p, lag_pc->phi_s, tempf, &uedlsp);

    udlvor[np] = (uvdwsp + uedlsp) * (1 - scovag - scovap);

  }

  cs_real_t fadhes = (udlvor[1] - udlvor[0]) / (2.0 * step);

  *adhesion_energ  = udlvor[0];

  /* Interaction between the sphere and small-scale asperities    */
  for (cs_lnum_t np = 0; np < 2; np++) {

    udlvor[np] = 0.0;

    cs_real_t distcc = _d_cut_off + step * (3 - 2 * (np + 1)) + rpart + rayasp;

    cs_real_t uvdwss, uedlss;

    _vdwss(lag_pc->cstham, lag_pc->lambda_vdw, &distcc, &rpart, &rayasp, &uvdwss);
    _edlss(&distcc, &rpart, &rayasp, lag_pc->phi_p, lag_pc->phi_s, tempf, &uedlss);

    udlvor[np] = uvdwss + uedlss;

  }

  fadhes = fadhes + (udlvor[1] - udlvor[0]) / (2.0 * step) * nbasp;

  *adhesion_energ  = *adhesion_energ + udlvor[0] * nbasp;

  /* Interaction between the sphere and large-scale asperities    */
  for (cs_lnum_t np = 0; np < 2; np++) {

    udlvor[np] = 0.0;
    cs_real_t distcc;

    if (nbasp == 0)
      distcc  = _d_cut_off + step * (3 - 2 * (np+1)) + rpart + rayasg;

    else if (nbasp > 0)
      distcc  = _d_cut_off + rayasp + step * (3 - 2 * (np+1)) + rpart + rayasg;

    cs_real_t uvdwss, uedlss;

    _vdwss(lag_pc->cstham, lag_pc->lambda_vdw, &distcc, &rpart, &rayasg, &uvdwss);
    _edlss(&distcc, &rpart, &rayasg, lag_pc->phi_p, lag_pc->phi_s, tempf, &uedlss);

    udlvor[np] = uvdwss + uedlss;

  }

  fadhes = fadhes + (udlvor[1] - udlvor[0]) / (2.0 * step) * nbasg;

  *adhesion_energ  = *adhesion_energ + udlvor[0] * nbasg;

  /* The force is negative when it is attractive   */

  if (fadhes >= 0.0)
    cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);

  else
    cs_lagr_particle_set_real(p_set, p_am, CS_LAGR_ADHESION_FORCE, -fadhes);

  /* The interaction should be negative to prevent reentrainment (attraction) */
  if (*adhesion_energ >= 0.0)
    *adhesion_energ = 0.0;

  else
    *adhesion_energ = CS_ABS(*adhesion_energ);

  /*  Calculation of adhesion torques exerted on the particle */
  cs_real_t rtmp;
  cs_random_uniform(1, &rtmp);

  cs_real_t dismom = rtmp;
  if (nbasp > 0)
    dismom =  (pow(rtmp, 1.0 / cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES)) * 2.0 - 1.0)
            * sqrt ((2.0 * rpart + rayasp) * rayasp);

  else {

    /* In the sphere-plate case, we use the deformation given by the DMT theory,
     * which is close to our approach  */

    cs_real_t omsurf = lag_pc->cstham / (24.0 * cs_math_pi * pow (_d_cut_off, 2));
    dismom = pow ((4.0 * cs_math_pi * omsurf * (pow (rpart, 2.0)) / modyeq), (1.0 / 3.0));
    /* dismom = (12.0d0 * cs_math_pi * omsurf * (rpart**2)/modyeq)**(1.0d0/3.0d0) */

  }

  dismom *= cs_lagr_particle_get_real(part, p_am, CS_LAGR_ADHESION_FORCE);

  cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, dismom);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the adhesion force and adhesion energy
 *        between two particles
 *
 * \param[in]  dpart            particle diameter
 * \param[in]  tempf            thermal scalar value at current time step
 * \param[out] adhesion_energ   particle adhesion energy
 * \param[out] adhesion_force   particle adhesion force
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_adh_pp(cs_real_t   dpart,
               cs_real_t   tempf,
               cs_real_t  *adhesion_energ,
               cs_real_t  *adhesion_force)
{
  /* ================================================================  */
  /* 0.    initialization */
  /* ================================================================  */

  cs_lagr_physico_chemical_t *lag_pc = cs_glob_lagr_physico_chemical;
  cs_lagr_clogging_model_t *lag_cm = cs_glob_lagr_clogging_model;

  cs_real_t step = 1e-11;

  /* ================================================================  */
  /* 1.    calculation of the adhesion force  */
  /* ================================================================  */

  cs_real_t rpart = 0.5 * dpart;

  /* Sum of {particle-plane} and {particle-asperity} energies     */

  cs_real_t udlvor[2];

  /* Interaction between the spheres */

  for (cs_lnum_t np = 0; np < 2; np++) {

    udlvor[np] = 0.0;

    cs_real_t distcc = _d_cut_off + step * (3 - 2 * (np + 1)) + 2.0 * rpart;

    cs_real_t uvdwss, uedlss;

    _vdwss(lag_cm->csthpp, lag_pc->lambda_vdw, &distcc, &rpart, &rpart, &uvdwss);
    _edlss(&distcc, &rpart, &rpart, lag_pc->phi_p, lag_pc->phi_p, tempf, &uedlss);

    udlvor[np] = uvdwss + uedlss;

  }

  *adhesion_force = CS_MAX( -(udlvor[1] - udlvor[0]) / (2.0 * step), 0.0);

  *adhesion_energ  = CS_MAX(-udlvor[0], 0.0);

}

END_C_DECLS
