#ifndef __CS_CONVECTION_DIFFUSION_PRIV_H__
#define __CS_CONVECTION_DIFFUSION_PRIV_H__

/*============================================================================
 * Private functions for convection-diffusion operators.
 *============================================================================*/

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_math.h"
#include "base/cs_parameters.h"  // for BC types

#include "alge/cs_convection_diffusion.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public inlined function
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Synchronize halos for scalar variables.
 *
 * parameters:
 *   m          <-- pointer to associated mesh structure
 *   halo_type  <-> halo type
 *   pvar       <-> variable
 *----------------------------------------------------------------------------*/

inline static void
cs_sync_scalar_halo(const cs_mesh_t  *m,
                    cs_halo_type_t    halo_type,
                    cs_real_t         pvar[])
{
  if (m->halo != NULL)
    cs_halo_sync_var(m->halo, halo_type, pvar);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the normalized face scalar using the specified NVD scheme.
 *
 * \param[in]   limiter     choice of the NVD scheme
 * \param[in]   nvf_p_c     normalised property of the current cell
 * \param[in]   nvf_r_f     normalised distance from the face
 * \param[in]   nvf_r_c     normalised distance from the current cell
 *
 * \return normalised face scalar value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static cs_real_t
cs_nvd_scheme_scalar(const cs_nvd_type_t  limiter,
                     const cs_real_t      nvf_p_c,
                     const cs_real_t      nvf_r_f,
                     const cs_real_t      nvf_r_c)
{
  cs_real_t nvf_p_f;

  cs_real_t beta_m, rfc, r1f, r1, r2, r3, b1, b2;

  switch (limiter) {
  case CS_NVD_GAMMA: /* Gamma scheme */
    beta_m = 0.1; /* in [0.1, 0.5] */
    rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);

    if (nvf_p_c < beta_m) {
      nvf_p_f = nvf_p_c*(1.+rfc*(1.-nvf_p_c)/beta_m);
    } else {
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = r1f*nvf_p_c+rfc;
    }

    break;

  case CS_NVD_SMART: /* SMART scheme */
    if (nvf_p_c < (nvf_r_c/3.)) {
      r1 = nvf_r_f*(1.-3.*nvf_r_c+2.*nvf_r_f);
      r2 = nvf_r_c*(1.-nvf_r_c);

      nvf_p_f = nvf_p_c*r1/r2;
    } else if (nvf_p_c <= (nvf_r_c*(1.+nvf_r_f-nvf_r_c)/nvf_r_f)) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = nvf_r_f*(r1f*nvf_p_c/nvf_r_c + rfc);
    } else {
      nvf_p_f = 1.;
    }

    break;

  case CS_NVD_CUBISTA: /* CUBISTA scheme */
    if (nvf_p_c < (3.*nvf_r_c/4.)) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);

      nvf_p_f = nvf_r_f*(1.+rfc/3.)*nvf_p_c/nvf_r_c;
    } else if (nvf_p_c <= (nvf_r_c*(1.+2.*(nvf_r_f-nvf_r_c))/(2.*nvf_r_f-nvf_r_c))) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = nvf_r_f*(r1f*nvf_p_c/nvf_r_c+rfc);
    } else {
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = 1.-.5*r1f*(1.-nvf_p_c);
    }

    break;

  case CS_NVD_SUPERBEE: /* SuperBee scheme */
    if (nvf_p_c < (nvf_r_c/(2.-nvf_r_c))) {
      nvf_p_f = (2.*nvf_r_f-nvf_r_c)*nvf_p_c/nvf_r_c;
    } else if (nvf_p_c < nvf_r_c) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = r1f*nvf_p_c+rfc;
    } else if (nvf_p_c < (nvf_r_c/nvf_r_f)) {
      nvf_p_f = nvf_r_f*nvf_p_c/nvf_r_c;
    } else {
      nvf_p_f = 1.;
    }

    break;

  case CS_NVD_MUSCL: /* MUSCL scheme */
    if (nvf_p_c < (.5*nvf_r_c)) {
      nvf_p_f = (2.*nvf_r_f-nvf_r_c)*nvf_p_c/nvf_r_c;
    } else if (nvf_p_c < (1.+nvf_r_c-nvf_r_f)) {
      nvf_p_f = nvf_p_c+nvf_r_f-nvf_r_c;
    } else {
      nvf_p_f = 1.;
    }

    break;

  case CS_NVD_MINMOD: /* MINMOD scheme */
    if (nvf_p_c < nvf_r_c) {
      nvf_p_f = nvf_r_f*nvf_p_c/nvf_r_c;
    } else {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = r1f*nvf_p_c+rfc;
    }

    break;

  case CS_NVD_CLAM: /* CLAM scheme */
    r1 = nvf_r_c*nvf_r_c-nvf_r_f;
    r2 = nvf_r_c*(nvf_r_c-1.);
    r3 = nvf_r_f-nvf_r_c;

    nvf_p_f = nvf_p_c*(r1+r3*nvf_p_c)/r2;

    break;

  case CS_NVD_STOIC: /* STOIC scheme */
    b1 = (nvf_r_c-nvf_r_f)*nvf_r_c;
    b2 = nvf_r_c+nvf_r_f+2.*nvf_r_f*nvf_r_f-4.*nvf_r_f*nvf_r_c;

    if (nvf_p_c < (b1/b2)) {
      r1 = -nvf_r_f*(1.-3.*nvf_r_c+2.*nvf_r_f);
      r2 = nvf_r_c*(nvf_r_c-1.);

      nvf_p_f = nvf_p_c*r1/r2;
    } else if (nvf_p_c < nvf_r_c) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = r1f*nvf_p_c+rfc;
    } else if (nvf_p_c < (nvf_r_c*(1.+nvf_r_f-nvf_r_c)/nvf_r_f)) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = nvf_r_f*(nvf_p_c*r1f/nvf_r_c+rfc);
    } else {
      nvf_p_f = 1.;
    }

    break;

  case CS_NVD_OSHER: /* OSHER scheme */
    if (nvf_p_c < (nvf_r_c/nvf_r_f)) {
      nvf_p_f = nvf_r_f*nvf_p_c/nvf_r_c;
    } else {
      nvf_p_f = 1.;
    }

    break;

  case CS_NVD_WASEB: /* WASEB scheme */
    r1 = nvf_r_c*nvf_r_f*(nvf_r_f-nvf_r_c);
    r2 = 2.*nvf_r_c*(1.-nvf_r_c)-nvf_r_f*(1.-nvf_r_f);

    if (nvf_p_c < (r1/r2)) {
      nvf_p_f = 2.*nvf_p_c;
    } else if (nvf_p_c <= (nvf_r_c*(1.+nvf_r_f-nvf_r_c)/nvf_r_f)) {
      rfc = (nvf_r_f-nvf_r_c)/(1.-nvf_r_c);
      r1f = (1.-nvf_r_f)/(1.-nvf_r_c);

      nvf_p_f = nvf_r_f*(nvf_p_c*r1f/nvf_r_c+rfc);
    } else {
      nvf_p_f = 1.;
    }

    break;

  default: /* Upwinding */
    nvf_p_f = nvf_p_c;
    break;
  }

  return nvf_p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the normalised face scalar using the specified NVD scheme
 *        for the case of a Volume-of-Fluid (VOF) transport equation.
 *
 * \param[in]  limiter          choice of the NVD scheme
 * \param[in]  i_face_u_normal  unit normal of face ij
 * \param[in]  nvf_p_c          normalised property of the current cell
 * \param[in]  nvf_r_f          normalised distance from the face
 * \param[in]  nvf_r_c          normalised distance from the current cell
 * \param[in]  gradv_c          gradient at central cell
 * \param[in]  c_courant        courant at central cell
 *
 * \return  normalised face scalar
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static cs_real_t
cs_nvd_vof_scheme_scalar(cs_nvd_type_t     limiter,
                         const cs_nreal_t  i_face_u_normal[3],
                         cs_real_t         nvf_p_c,
                         cs_real_t         nvf_r_f,
                         cs_real_t         nvf_r_c,
                         const cs_real_t   gradv_c[3],
                         const cs_real_t   c_courant)
{
  assert(limiter >= CS_NVD_VOF_HRIC);

  cs_real_t nvf_p_f;
  cs_real_t blend, high_order, low_order, ratio;

  /* Compute gradient angle indicator */
  cs_real_t dotp = cs::abs(cs_math_3_dot_product(gradv_c, i_face_u_normal));
  cs_real_t sgrad = cs_math_3_norm(gradv_c);
  cs_real_t denom = sgrad;

  if (limiter == CS_NVD_VOF_HRIC) {   /* M-HRIC scheme */
    /* High order scheme : Bounded Downwind */
    if (nvf_p_c <= .5) {
      high_order = 2.*nvf_p_c;
    } else {
      high_order = 1.;
    }

    /* Low order scheme : MUSCL */
    low_order = cs_nvd_scheme_scalar(CS_NVD_MUSCL,
                                     nvf_p_c,
                                     nvf_r_f,
                                     nvf_r_c);

    /* Compute the blending factor */
    if (denom <= (cs_math_epzero*dotp)) {
      blend = 1.;
    } else {
      ratio = dotp/denom;
      blend = cs_math_fmin(1., pow(ratio, .5));
    }

    /* Blending */
    nvf_p_f = blend*high_order + (1.-blend)*low_order;

    /* Extra blending due to the cell Courant number */
    if (c_courant < .7 && c_courant > .3) {
      nvf_p_f = nvf_p_f + (nvf_p_f - low_order)*(.7 - c_courant )/.4;
    } else if (c_courant >= .7) {
      nvf_p_f = low_order;
    }
  } else if (limiter == CS_NVD_VOF_CICSAM) { /* M-CICSAM scheme */
    /* High order scheme : HYPER-C + SUPERBEE */
    if (c_courant <= .3) {
      high_order = cs_math_fmin(1., nvf_p_c/(c_courant+cs_math_epzero));
    } else if (c_courant <= .6) {
      high_order = cs_math_fmin(1., nvf_p_c/.3);
    } else if (c_courant <= .7) {
      cs_real_t superbee = cs_nvd_scheme_scalar(CS_NVD_SUPERBEE,
                                                nvf_p_c,
                                                nvf_r_f,
                                                nvf_r_c);
      high_order =  10.*(  (.7-c_courant)*cs_math_fmin(1., nvf_p_c/.3)
                         + (c_courant-.6)*superbee);
    }
    else {
      high_order = cs_nvd_scheme_scalar(CS_NVD_SUPERBEE,
                                        nvf_p_c,
                                        nvf_r_f,
                                        nvf_r_c);
    }

    /* Low order scheme : MUSCL */
    low_order = cs_nvd_scheme_scalar(CS_NVD_MUSCL,
                                     nvf_p_c,
                                     nvf_r_f,
                                     nvf_r_c);

    /* Compute the blending factor */
    if (denom <= (cs_math_epzero*dotp)) {
      blend = 1.;
    } else {
      ratio = dotp/denom;
      blend = cs_math_fmin(1., pow(ratio, 2.));
    }

    /* Blending */
    nvf_p_f = blend*high_order + (1.-blend)*low_order;
  } else { /* STACS scheme */
    /* High order scheme : SUPERBEE */
    high_order = cs_nvd_scheme_scalar(CS_NVD_SUPERBEE,
                                      nvf_p_c,
                                      nvf_r_f,
                                      nvf_r_c);

    /* Low order scheme : STOIC */
    low_order = cs_nvd_scheme_scalar(CS_NVD_STOIC,
                                     nvf_p_c,
                                     nvf_r_f,
                                     nvf_r_c);

    /* Compute the blending factor */
    if (denom <= (cs_math_epzero*dotp)) {
      blend = 1.;
    } else {
      ratio = dotp/denom;
      blend = cs_math_fmin(1., pow(ratio, 4.));
    }

    /* Blending */
    nvf_p_f = blend*high_order + (1.-blend)*low_order;
  }

  return nvf_p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute slope test criteria at internal face between cell i and j.
 *
 * \param[in]     pi                value at cell i
 * \param[in]     pj                value at cell j
 * \param[in]     distf             distance IJ.Nij
 * \param[in]     i_face_u_normal   face unit normal
 * \param[in]     gradi             gradient at cell i
 * \param[in]     gradj             gradient at cell j
 * \param[in]     grdpai            upwind gradient at cell i
 * \param[in]     grdpaj            upwind gradient at cell j
 * \param[in]     i_massflux        mass flux at face (from i to j)
 * \param[out]    testij            value of slope test first criterion
 * \param[out]    tesqck            value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_slope_test(const cs_real_t   pi,
              const cs_real_t   pj,
              const cs_real_t   distf,
              const cs_nreal_t  i_face_u_normal[3],
              const cs_real_t   gradi[3],
              const cs_real_t   gradj[3],
              const cs_real_t   grdpai[3],
              const cs_real_t   grdpaj[3],
              const cs_real_t   i_massflux,
              cs_real_t        *testij,
              cs_real_t        *tesqck)
{
  cs_real_t dcc, ddi, ddj;

  /* Slope test */

  *testij =   grdpai[0]*grdpaj[0]
            + grdpai[1]*grdpaj[1]
            + grdpai[2]*grdpaj[2];

  if (i_massflux > 0.) {
    dcc =   gradi[0]*i_face_u_normal[0]
          + gradi[1]*i_face_u_normal[1]
          + gradi[2]*i_face_u_normal[2];
    ddi =   grdpai[0]*i_face_u_normal[0]
          + grdpai[1]*i_face_u_normal[1]
          + grdpai[2]*i_face_u_normal[2];
    ddj = (pj-pi)/distf;
  }
  else {
    dcc =   gradj[0]*i_face_u_normal[0]
          + gradj[1]*i_face_u_normal[1]
          + gradj[2]*i_face_u_normal[2];
    ddi = (pj-pi)/distf;
    ddj =   grdpaj[0]*i_face_u_normal[0]
          + grdpaj[1]*i_face_u_normal[1]
          + grdpaj[2]*i_face_u_normal[2];
  }

  *tesqck = cs_math_sq(dcc) - cs_math_sq(ddi-ddj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute slope test criteria at internal face between cell i and j.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     pi                value at cell i
 * \param[in]     pj                value at cell j
 * \param[in]     distf             distance IJ.Nij
 * \param[in]     i_face_u_normal   face unit normal
 * \param[in]     gradi             gradient at cell i
 * \param[in]     gradj             gradient at cell j
 * \param[in]     gradsti           upwind gradient at cell i
 * \param[in]     gradstj           upwind gradient at cell j
 * \param[in]     i_massflux        mass flux at face (from i to j)
 * \param[out]    testij            value of slope test first criterion
 * \param[out]    tesqck            value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_slope_test_strided(const cs_real_t    pi[stride],
                      const cs_real_t    pj[stride],
                      const cs_real_t    distf,
                      const cs_nreal_t   i_face_u_normal[3],
                      const cs_real_t    gradi[stride][3],
                      const cs_real_t    gradj[stride][3],
                      const cs_real_t    gradsti[stride][3],
                      const cs_real_t    gradstj[stride][3],
                      const cs_real_t    i_massflux,
                      cs_real_t         *testij,
                      cs_real_t         *tesqck)
{
  cs_real_t dcc[stride], ddi[stride], ddj[stride];

  *testij = 0.;
  *tesqck = 0.;

  /* Slope test */

  for (int i = 0; i < stride; i++) {
    *testij +=   gradsti[i][0]*gradstj[i][0]
               + gradsti[i][1]*gradstj[i][1]
               + gradsti[i][2]*gradstj[i][2];

    if (i_massflux > 0.) {
      dcc[i] =   gradi[i][0]*i_face_u_normal[0]
               + gradi[i][1]*i_face_u_normal[1]
               + gradi[i][2]*i_face_u_normal[2];
      ddi[i] =   gradsti[i][0]*i_face_u_normal[0]
               + gradsti[i][1]*i_face_u_normal[1]
               + gradsti[i][2]*i_face_u_normal[2];
      ddj[i] = (pj[i]-pi[i])/distf;
    }
    else {
      dcc[i] =   gradj[i][0]*i_face_u_normal[0]
               + gradj[i][1]*i_face_u_normal[1]
               + gradj[i][2]*i_face_u_normal[2];
      ddi[i] = (pj[i]-pi[i])/distf;
      ddj[i] =   gradstj[i][0]*i_face_u_normal[0]
               + gradstj[i][1]*i_face_u_normal[1]
               + gradstj[i][2]*i_face_u_normal[2];
    }

    *tesqck += cs_math_sq(dcc[i]) - cs_math_sq(ddi[i]-ddj[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' and J'.
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    recoi        reconstruction at cell i
 * \param[out]    recoj        reconstruction at cell j
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_compute_quantities(const cs_real_t     bldfrp,
                        const cs_rreal_3_t  diipf,
                        const cs_rreal_3_t  djjpf,
                        const cs_real_3_t   gradi,
                        const cs_real_3_t   gradj,
                        const cs_real_t     pi,
                        const cs_real_t     pj,
                        cs_real_t          *recoi,
                        cs_real_t          *recoj,
                        cs_real_t          *pip,
                        cs_real_t          *pjp)
{
  cs_real_t gradpf[3] = {0.5*(gradi[0] + gradj[0]),
                         0.5*(gradi[1] + gradj[1]),
                         0.5*(gradi[2] + gradj[2])};

  /* reconstruction only if IRCFLP = 1 */
  *recoi = bldfrp*(cs_math_3_dot_product(gradpf, diipf));
  *recoj = bldfrp*(cs_math_3_dot_product(gradpf, djjpf));
  *pip = pi + *recoi;
  *pjp = pj + *recoj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' and J'.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    recoi        reconstruction at cell i
 * \param[out]    recoj        reconstruction at cell j
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_compute_quantities_strided(const cs_real_t   bldfrp,
                                const cs_rreal_t  diipf[3],
                                const cs_rreal_t  djjpf[3],
                                const cs_real_t   gradi[stride][3],
                                const cs_real_t   gradj[stride][3],
                                const cs_real_t   pi[stride],
                                const cs_real_t   pj[stride],
                                cs_real_t         recoi[stride],
                                cs_real_t         recoj[stride],
                                cs_real_t         pip[stride],
                                cs_real_t         pjp[stride])
{
  cs_real_t dpvf[3];

  /* x-y-z components, p = u, v, w */

  for (int isou = 0; isou < stride; isou++) {

    for (int jsou = 0; jsou < 3; jsou++)
      dpvf[jsou] = 0.5*( gradi[isou][jsou]
                       + gradj[isou][jsou]);

    /* reconstruction only if IRCFLP = 1 */

    recoi[isou] = bldfrp*(cs_math_3_dot_product(dpvf, diipf));
    recoj[isou] = bldfrp*(cs_math_3_dot_product(dpvf, djjpf));

    pip[isou] = pi[isou] + recoi[isou];
    pjp[isou] = pj[isou] + recoj[isou];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute relaxed values at cell i and j.
 *
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     pia      old value at cell i
 * \param[in]     pja      old value at cell j
 * \param[in]     recoi    reconstruction at cell i
 * \param[in]     recoj    reconstruction at cell j
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pjr      relaxed value at cell j
 * \param[out]    pipr     relaxed reconstructed value at cell i
 * \param[out]    pjpr     relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_relax_c_val(const double     relaxp,
                 const cs_real_t  pia,
                 const cs_real_t  pja,
                 const cs_real_t  recoi,
                 const cs_real_t  recoj,
                 const cs_real_t  pi,
                 const cs_real_t  pj,
                 cs_real_t       *pir,
                 cs_real_t       *pjr,
                 cs_real_t       *pipr,
                 cs_real_t       *pjpr)
{
  *pir = pi/relaxp - (1.-relaxp)/relaxp * pia;
  *pjr = pj/relaxp - (1.-relaxp)/relaxp * pja;

  *pipr = *pir + recoi;
  *pjpr = *pjr + recoj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute relaxed values at cell i and j.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     pia      old value at cell i
 * \param[in]     pja      old value at cell j
 * \param[in]     recoi    reconstruction at cell i
 * \param[in]     recoj    reconstruction at cell j
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pjr      relaxed value at cell j
 * \param[out]    pipr     relaxed reconstructed value at cell i
 * \param[out]    pjpr     relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_relax_c_val_strided(cs_real_t        relaxp,
                         const cs_real_t  pia[stride],
                         const cs_real_t  pja[stride],
                         const cs_real_t  recoi[stride],
                         const cs_real_t  recoj[stride],
                         const cs_real_t  pi[stride],
                         const cs_real_t  pj[stride],
                         cs_real_t        pir[stride],
                         cs_real_t        pjr[stride],
                         cs_real_t        pipr[stride],
                         cs_real_t        pjpr[stride])
{
  for (int isou = 0; isou < stride; isou++) {
    pir[isou] = pi[isou] /relaxp - (1.-relaxp)/relaxp * pia[isou];
    pjr[isou] = pj[isou] /relaxp - (1.-relaxp)/relaxp * pja[isou];

    pipr[isou] = pir[isou] + recoi[isou];
    pjpr[isou] = pjr[isou] + recoj[isou];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using an upwind scheme.
 *
 * \param[in]     p       value at cell
 * \param[out]    pf      value at face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_upwind_f_val(const cs_real_t  p,
                cs_real_t       *pf)
{
  *pf = p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using an upwind scheme.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     p       value at cell
 * \param[out]    pf      value at face
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_upwind_f_val_strided(const cs_real_t  p[stride],
                        cs_real_t        pf[stride])
{
  for (int isou = 0; isou < stride; isou++)
    pf[isou] = p[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a centered scheme.
 *
 * \param[in]     pnd      weight
 * \param[in]     pip      (relaxed) reconstructed value at cell i
 * \param[in]     pjp      (relaxed) reconstructed value at cell j
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_centered_f_val(double      pnd,
                  cs_real_t   pip,
                  cs_real_t   pjp,
                  cs_real_t  *pf)
{
  *pf = pnd*pip + (1.-pnd)*pjp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a centered scheme.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     pnd      weight
 * \param[in]     pip      (relaxed) reconstructed value at cell i
 * \param[in]     pjp      (relaxed) reconstructed value at cell j
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_centered_f_val_strided(double           pnd,
                          const cs_real_t  pip[stride],
                          const cs_real_t  pjp[stride],
                          cs_real_t        pf[stride])
{
  for (int isou = 0; isou < stride; isou++)
    pf[isou] = pnd*pip[isou] + (1.-pnd)*pjp[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a Second Order Linear Upwind scheme.
 *
 * \param[in]     cell_cen     center of gravity coordinates of cell
 * \param[in]     i_face_cog   center of gravity coordinates of face
 * \param[in]     grad         gradient at cell
 * \param[in]     p            (relaced) value at cell
 * \param[out]    pf           face value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_solu_f_val(const cs_real_t   cell_cen[3],
              const cs_real_t   i_face_cog[3],
              const cs_real_t   grad[3],
              const cs_real_t   p,
              cs_real_t        *pf)
{
  cs_real_t df[3];

  df[0] = i_face_cog[0] - cell_cen[0];
  df[1] = i_face_cog[1] - cell_cen[1];
  df[2] = i_face_cog[2] - cell_cen[2];

  *pf = p + cs_math_3_dot_product(df, grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a Second Order Linear Upwind scheme.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     cell_cen     center of gravity coordinates of cell
 * \param[in]     i_face_cog   center of gravity coordinates of face
 * \param[in]     grad         gradient at cell
 * \param[in]     p            (relaced) value at cell
 * \param[out]    pf           face value
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_solu_f_val_strided(const cs_real_t  cell_cen[3],
                      const cs_real_t  i_face_cog[3],
                      const cs_real_t  grad[stride][3],
                      const cs_real_t  p[stride],
                      cs_real_t        pf[stride])
{
  cs_real_t df[3];

  for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
    df[jsou] = i_face_cog[jsou] - cell_cen[jsou];

  for (cs_lnum_t isou = 0; isou < stride; isou++) {
     pf[isou] = p[isou] + df[0]*grad[isou][0]
                        + df[1]*grad[isou][1]
                        + df[2]*grad[isou][2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Blend face values for a centered or SOLU scheme with face values for
 * an upwind scheme.
 *
 * \param[in]     blencp   proportion of second order scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     p        (relaxed) value at cell
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_blend_f_val(const double     blencp,
               const cs_real_t  p,
               cs_real_t       *pf)
{
  *pf = blencp * (*pf) + (1. - blencp) * p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Blend face values for a centered or SOLU scheme with face values for
 * an upwind scheme.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     blencp   proportion of second order scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     p        (relaxed) value at cell
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_blend_f_val_strided(const double     blencp,
                       const cs_real_t  p[stride],
                       cs_real_t        pf[stride])
{
  for (int isou = 0; isou < stride; isou++)
    pf[isou] = blencp*(pf[isou])+(1.-blencp)*p[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (substracting the mass accumulation from them)
 * to fluxes at face ij.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-scheme,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pifri        contribution of i to flux from i to j
 * \param[in]     pifrj        contribution of i to flux from j to i
 * \param[in]     pjfri        contribution of j to flux from i to j
 * \param[in]     pjfrj        contribution of j to flux from j to i
 * \param[in]     i_massflux   mass flux at face ij
 * \param[in]     xcppi        specific heat value if the scalar is the temperature,
 *                             1 otherwise at cell i
 * \param[in]     xcppj        specific heat value if the scalar is the temperature,
 *                             1 otherwise at cell j
 * \param[in,out] fluxij       fluxes at face ij
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_conv_flux(const int       iconvp,
               const cs_real_t thetap,
               const int       imasac,
               const cs_real_t pi,
               const cs_real_t pj,
               const cs_real_t pifri,
               const cs_real_t pifrj,
               const cs_real_t pjfri,
               const cs_real_t pjfrj,
               const cs_real_t i_massflux,
               const cs_real_t xcppi,
               const cs_real_t xcppj,
               cs_real_2_t     fluxij)
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux + fabs(i_massflux));
  fluj = 0.5*(i_massflux - fabs(i_massflux));

  fluxij[0] += iconvp*xcppi*(thetap*(flui*pifri + fluj*pjfri) - imasac*i_massflux*pi);
  fluxij[1] += iconvp*xcppj*(thetap*(flui*pifrj + fluj*pjfrj) - imasac*i_massflux*pj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (substracting the mass accumulation from them)
 * to fluxes at face ij.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 * \param[in]     i_massflux   mass flux at face ij
 * \param[in,out] fluxi       fluxes at face i
 * \param[in,out] fluxj       fluxes at face j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_conv_flux_strided(int              iconvp,
                       cs_real_t        thetap,
                       int              imasac,
                       const cs_real_t  pi[stride],
                       const cs_real_t  pj[stride],
                       const cs_real_t  pifri[stride],
                       const cs_real_t  pifrj[stride],
                       const cs_real_t  pjfri[stride],
                       const cs_real_t  pjfrj[stride],
                       cs_real_t        i_massflux,
                       cs_real_t        fluxi[stride],
                       cs_real_t        fluxj[stride])
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux + fabs(i_massflux));
  fluj = 0.5*(i_massflux - fabs(i_massflux));

  for (int isou = 0; isou < stride; isou++) {
    fluxi[isou] +=  iconvp*(  thetap*(flui*pifri[isou] + fluj*pjfri[isou])
                            - imasac*i_massflux*pi[isou]);
    fluxj[isou] +=  iconvp*(  thetap*(flui*pifrj[isou] + fluj*pjfrj[isou])
                            - imasac*i_massflux*pj[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive fluxes to fluxes at face ij.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-scheme,
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[in]     i_visc   diffusion coefficient (divided by IJ) at face ij
 * \param[in,out] fluxij   fluxes at face ij
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_diff_flux(const int        idiffp,
               const cs_real_t  thetap,
               const cs_real_t  pip,
               const cs_real_t  pjp,
               const cs_real_t  pipr,
               const cs_real_t  pjpr,
               const cs_real_t  i_visc,
               cs_real_t        fluxij[2])
{
  fluxij[0] += idiffp*thetap*i_visc*(pipr -pjp);
  fluxij[1] += idiffp*thetap*i_visc*(pip -pjpr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive fluxes to fluxes at face ij.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-scheme,
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[in]     i_visc   diffusion coefficient (divided by IJ) at face ij
 * \param[in,out] fluxi    fluxes at face i
 * \param[in,out] fluxj    fluxes at face j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_diff_flux_strided(int             idiffp,
                       cs_real_t       thetap,
                       const cs_real_t  pip[stride],
                       const cs_real_t  pjp[stride],
                       const cs_real_t  pipr[stride],
                       const cs_real_t  pjpr[stride],
                       const cs_real_t  i_visc,
                       cs_real_t        fluxi[stride],
                       cs_real_t        fluxj[stride])
{
  for (int isou = 0; isou < stride; isou++) {
    fluxi[isou] += idiffp*thetap*i_visc*(pipr[isou] -pjp[isou]);
    fluxj[isou] += idiffp*thetap*i_visc*(pip[isou] -pjpr[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and a pure upwind flux.
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pia          old value at cell i
 * \param[in]     pja          old value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 * \param[out]    pipr         relaxed reconstructed value at cell i
 * \param[out]    pjpr         relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_steady_upwind(const cs_real_t   bldfrp,
                      const cs_real_t   relaxp,
                      const cs_rreal_t  diipf[3],
                      const cs_rreal_t  djjpf[3],
                      const cs_real_t   gradi[3],
                      const cs_real_t   gradj[3],
                      const cs_real_t   pi,
                      const cs_real_t   pj,
                      const cs_real_t   pia,
                      const cs_real_t   pja,
                      cs_real_t        *pifri,
                      cs_real_t        *pifrj,
                      cs_real_t        *pjfri,
                      cs_real_t        *pjfrj,
                      cs_real_t        *pip,
                      cs_real_t        *pjp,
                      cs_real_t        *pipr,
                      cs_real_t        *pjpr)
{
  cs_real_t pir, pjr;
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);

  cs_i_relax_c_val(relaxp,
                   pia,
                   pja,
                   recoi,
                   recoj,
                   pi,
                   pj,
                   &pir,
                   &pjr,
                   pipr,
                   pjpr);

  cs_upwind_f_val(pi,
                  pifrj);
  cs_upwind_f_val(pir,
                  pifri);
  cs_upwind_f_val(pj,
                  pjfri);
  cs_upwind_f_val(pjr,
                  pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and a pure upwind flux.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pia          old value at cell i
 * \param[in]     pja          old value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 * \param[out]    pipr         relaxed reconstructed value at cell i
 * \param[out]    pjpr         relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_steady_upwind_strided(const cs_real_t   bldfrp,
                              const cs_real_t   relaxp,
                              const cs_rreal_t  diipf[3],
                              const cs_rreal_t  djjpf[3],
                              const cs_real_t   gradi[stride][3],
                              const cs_real_t   gradj[stride][3],
                              const cs_real_t   pi[stride],
                              const cs_real_t   pj[stride],
                              const cs_real_t   pia[stride],
                              const cs_real_t   pja[stride],
                              cs_real_t         pifri[stride],
                              cs_real_t         pifrj[stride],
                              cs_real_t         pjfri[stride],
                              cs_real_t         pjfrj[stride],
                              cs_real_t         pip[stride],
                              cs_real_t         pjp[stride],
                              cs_real_t         pipr[stride],
                              cs_real_t         pjpr[stride])
{
  cs_real_t pir[stride], pjr[stride];
  cs_real_t recoi[stride], recoj[stride];

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  cs_i_relax_c_val_strided<stride>(relaxp,
                                   pia,
                                   pja,
                                   recoi,
                                   recoj,
                                   pi,
                                   pj,
                                   pir,
                                   pjr,
                                   pipr,
                                   pjpr);

  cs_upwind_f_val_strided<stride>(pi, pifrj);
  cs_upwind_f_val_strided<stride>(pir, pifri);
  cs_upwind_f_val_strided<stride>(pj, pjfri);
  cs_upwind_f_val_strided<stride>(pjr, pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and a pure upwind flux.
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    pif          contribution of i to flux from i to j
 * \param[out]    pjf          contribution of j to flux from i to j
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_upwind(const cs_real_t   bldfrp,
                        const cs_rreal_t  diipf[3],
                        const cs_rreal_t  djjpf[3],
                        const cs_real_t   gradi[3],
                        const cs_real_t   gradj[3],
                        const cs_real_t   pi,
                        const cs_real_t   pj,
                        cs_real_t        *pif,
                        cs_real_t        *pjf,
                        cs_real_t        *pip,
                        cs_real_t        *pjp)
{
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);

  cs_upwind_f_val(pi, pif);
  cs_upwind_f_val(pj, pjf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and a pure upwind flux.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    pif          contribution of i to flux from i to j
 * \param[out]    pjf          contribution of j to flux from i to j
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_upwind_strided(const cs_real_t   bldfrp,
                                const cs_rreal_t  diipf[3],
                                const cs_rreal_t  djjpf[3],
                                const cs_real_t   gradi[stride][3],
                                const cs_real_t   gradj[stride][3],
                                const cs_real_t   pi[stride],
                                const cs_real_t   pj[stride],
                                cs_real_t         pif[stride],
                                cs_real_t         pjf[stride],
                                cs_real_t         pip[stride],
                                cs_real_t         pjp[stride])
{
  cs_real_t recoi[stride], recoj[stride];

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  cs_upwind_f_val_strided<stride>(pi, pif);
  cs_upwind_f_val_strided<stride>(pj, pjf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and without enabling slope tests.
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     blencp       proportion of second order scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     gradupi      gradient upwind at cell i
 * \param[in]     gradupj      gradient upwind at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pia          old value at cell i
 * \param[in]     pja          old value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 * \param[out]    pipr         relaxed reconstructed value at cell i
 * \param[out]    pjpr         relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_steady(const cs_real_t   bldfrp,
               const int         ischcp,
               const double      relaxp,
               const double      blencp,
               const cs_real_t   weight,
               const cs_real_t   cell_ceni[3],
               const cs_real_t   cell_cenj[3],
               const cs_real_t   i_face_cog[3],
               const cs_rreal_t  diipf[3],
               const cs_rreal_t  djjpf[3],
               const cs_real_t   gradi[3],
               const cs_real_t   gradj[3],
               const cs_real_t   gradupi[3],
               const cs_real_t   gradupj[3],
               const cs_real_t   pi,
               const cs_real_t   pj,
               const cs_real_t   pia,
               const cs_real_t   pja,
               cs_real_t        *pifri,
               cs_real_t        *pifrj,
               cs_real_t        *pjfri,
               cs_real_t        *pjfrj,
               cs_real_t        *pip,
               cs_real_t        *pjp,
               cs_real_t        *pipr,
               cs_real_t        *pjpr)
{
  cs_real_t pir, pjr;
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);

  cs_i_relax_c_val(relaxp,
                   pia,
                   pja,
                   recoi,
                   recoj,
                   pi,
                   pj,
                   &pir,
                   &pjr,
                   pipr,
                   pjpr);

  if (ischcp == 1) {

    /* Centered
       --------*/

    cs_centered_f_val(weight,
                      *pip,
                      *pjpr,
                      pifrj);
    cs_centered_f_val(weight,
                      *pipr,
                      *pjp,
                      pifri);
    cs_centered_f_val(weight,
                      *pipr,
                      *pjp,
                      pjfri);
    cs_centered_f_val(weight,
                      *pip,
                      *pjpr,
                      pjfrj);

  } else if (ischcp == 0) {

    /* Original SOLU
       --------------*/

    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradi,
                  pi,
                  pifrj);
    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradi,
                  pir,
                  pifri);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradj,
                  pj,
                  pjfri);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradj,
                  pjr,
                  pjfrj);

  } else {

    /* SOLU
       ----*/

    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradupi,
                  pi,
                  pifrj);
    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradupi,
                  pir,
                  pifri);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradupj,
                  pj,
                  pjfri);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradupj,
                  pjr,
                  pjfrj);

  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pifrj);
  cs_blend_f_val(blencp,
                 pir,
                 pifri);
  cs_blend_f_val(blencp,
                 pj,
                 pjfri);
  cs_blend_f_val(blencp,
                 pjr,
                 pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and without enabling slope tests.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     blencp       proportion of second order scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pia          old value at cell i
 * \param[in]     pja          old value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 * \param[out]    pipr         relaxed reconstructed value at cell i
 * \param[out]    pjpr         relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_steady_strided(cs_real_t         bldfrp,
                       int               ischcp,
                       double            relaxp,
                       double            blencp,
                       cs_real_t         weight,
                       const cs_real_t   cell_ceni[3],
                       const cs_real_t   cell_cenj[3],
                       const cs_real_t   i_face_cog[3],
                       const cs_rreal_t  diipf[3],
                       const cs_rreal_t  djjpf[3],
                       const cs_real_t   gradi[stride][3],
                       const cs_real_t   gradj[stride][3],
                       const cs_real_t   pi[stride],
                       const cs_real_t   pj[stride],
                       const cs_real_t   pia[stride],
                       const cs_real_t   pja[stride],
                       cs_real_t         pifri[stride],
                       cs_real_t         pifrj[stride],
                       cs_real_t         pjfri[stride],
                       cs_real_t         pjfrj[stride],
                       cs_real_t         pip[stride],
                       cs_real_t         pjp[stride],
                       cs_real_t         pipr[stride],
                      cs_real_t          pjpr[stride])

{
  cs_real_t pir[stride], pjr[stride];
  cs_real_t recoi[stride], recoj[stride];

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  cs_i_relax_c_val_strided<stride>(relaxp,
                                   pia,
                                   pja,
                                   recoi,
                                   recoj,
                                   pi,
                                   pj,
                                   pir,
                                   pjr,
                                   pipr,
                                   pjpr);

  if (ischcp == 1) {

    /* Centered
       --------*/

    cs_centered_f_val_strided<stride>(weight, pip, pjpr, pifrj);
    cs_centered_f_val_strided<stride>(weight, pipr, pjp, pifri);
    cs_centered_f_val_strided<stride>(weight, pipr, pjp, pjfri);
    cs_centered_f_val_strided<stride>(weight, pip, pjpr, pjfrj);

  }
  else {

    /* Second order
       ------------*/

    cs_solu_f_val_strided<stride>(cell_ceni,
                                  i_face_cog,
                                  gradi,
                                  pi,
                                  pifrj);
    cs_solu_f_val_strided<stride>(cell_ceni,
                                  i_face_cog,
                                  gradi,
                                  pir,
                                  pifri);
    cs_solu_f_val_strided<stride>(cell_cenj,
                                  i_face_cog,
                                  gradj,
                                  pj,
                                  pjfri);
    cs_solu_f_val_strided<stride>(cell_cenj,
                                  i_face_cog,
                                  gradj,
                                  pjr,
                                  pjfrj);

  }

  /* Blending
     --------*/

  cs_blend_f_val_strided<stride>(blencp, pi, pifrj);
  cs_blend_f_val_strided<stride>(blencp, pir, pifri);
  cs_blend_f_val_strided<stride>(blencp, pj, pjfri);
  cs_blend_f_val_strided<stride>(blencp, pjr, pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and without enabling slope tests.
 *
 * \param[in]     bldfrp         reconstruction blending factor
 * \param[in]     ischcp         second order convection scheme flag
 * \param[in]     blencp         proportion of second order scheme,
 *                               (1-blencp) is the proportion of upwind.
 * \param[in]     weight         geometrical weight
 * \param[in]     cell_ceni      center of gravity coordinates of cell i
 * \param[in]     cell_cenj      center of gravity coordinates of cell i
 * \param[in]     i_face_cog     center of gravity coordinates of face ij
 * \param[in]     hybrid_blend_i blending factor between SOLU and centered
 * \param[in]     hybrid_blend_j blending factor between SOLU and centered
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi          gradient at cell i
 * \param[in]     gradj          gradient at cell j
 * \param[in]     gradupi        upwind gradient at cell i
 * \param[in]     gradupj        upwind gradient at cell j
 * \param[in]     pi             value at cell i
 * \param[in]     pj             value at cell j
 * \param[out]    pif            contribution of i to flux from i to j
 * \param[out]    pjf            contribution of j to flux from i to j
 * \param[out]    pip            reconstructed value at cell i
 * \param[out]    pjp            reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady(const cs_real_t     bldfrp,
                 const int           ischcp,
                 const double        blencp,
                 const cs_real_t     weight,
                 const cs_real_3_t   cell_ceni,
                 const cs_real_3_t   cell_cenj,
                 const cs_real_3_t   i_face_cog,
                 const cs_real_t     hybrid_blend_i,
                 const cs_real_t     hybrid_blend_j,
                 const cs_rreal_3_t  diipf,
                 const cs_rreal_3_t  djjpf,
                 const cs_real_3_t   gradi,
                 const cs_real_3_t   gradj,
                 const cs_real_3_t   gradupi,
                 const cs_real_3_t   gradupj,
                 const cs_real_t     pi,
                 const cs_real_t     pj,
                 cs_real_t          *pif,
                 cs_real_t          *pjf,
                 cs_real_t          *pip,
                 cs_real_t          *pjp)
{
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);


  if (ischcp == 1) {

    /* Centered
       --------*/

    cs_centered_f_val(weight,
                      *pip,
                      *pjp,
                      pif);
    cs_centered_f_val(weight,
                      *pip,
                      *pjp,
                      pjf);

  } else if (ischcp == 0) {

    /* Legacy SOLU
       -----------*/

    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradi,
                  pi,
                  pif);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradj,
                  pj,
                  pjf);

  } else if (ischcp == 3) {

    /* Centered
       --------*/

    cs_centered_f_val(weight,
                      *pip,
                      *pjp,
                      pif);
    cs_centered_f_val(weight,
                      *pip,
                      *pjp,
                      pjf);

    /* Legacy SOLU
       -----------*/
    cs_real_t pif_up, pjf_up;
    cs_real_t hybrid_blend_interp;

    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradi,
                  pi,
                  &pif_up);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradj,
                  pj,
                  &pjf_up);

    hybrid_blend_interp = fmin(hybrid_blend_i,hybrid_blend_j);
    *pif = hybrid_blend_interp*(*pif) + (1. - hybrid_blend_interp)*pif_up;
    *pjf = hybrid_blend_interp*(*pjf) + (1. - hybrid_blend_interp)*pjf_up;

  } else {

    /* SOLU
       ----*/

    cs_solu_f_val(cell_ceni,
                  i_face_cog,
                  gradupi,
                  pi,
                  pif);
    cs_solu_f_val(cell_cenj,
                  i_face_cog,
                  gradupj,
                  pj,
                  pjf);

  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pif);
  cs_blend_f_val(blencp,
                 pj,
                 pjf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and without enabling slope tests.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp       reconstruction blending factor
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     blencp       proportion of second order scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     hybrid_blend_i blending factor between SOLU and centered
 * \param[in]     hybrid_blend_j blending factor between SOLU and centered
 * \param[in]     diipf        distance II'
 * \param[in]     djjpf        distance JJ'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    pif          contribution of i to flux from i to j
 * \param[out]    pjf          contribution of j to flux from i to j
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_strided(cs_real_t         bldfrp,
                         int               ischcp,
                         double            blencp,
                         cs_real_t         weight,
                         const cs_real_t   cell_ceni[3],
                         const cs_real_t   cell_cenj[3],
                         const cs_real_t   i_face_cog[3],
                         const cs_real_t   hybrid_blend_i,
                         const cs_real_t   hybrid_blend_j,
                         const cs_rreal_t  diipf[3],
                         const cs_rreal_t  djjpf[3],
                         const cs_real_t   gradi[stride][3],
                         const cs_real_t   gradj[stride][3],
                         const cs_real_t   pi[stride],
                         const cs_real_t   pj[stride],
                         cs_real_t         pif[stride],
                         cs_real_t         pjf[stride],
                         cs_real_t         pip[stride],
                         cs_real_t         pjp[stride])

{
  cs_real_t recoi[stride], recoj[stride];

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  if (ischcp == 1) {

    /* Centered
       --------*/

    cs_centered_f_val_strided<stride>(weight, pip, pjp, pif);
    cs_centered_f_val_strided<stride>(weight, pip, pjp, pjf);

  }
  else if (ischcp == 3) {

    /* Centered
       --------*/

    cs_centered_f_val_strided<stride>(weight, pip, pjp, pif);
    cs_centered_f_val_strided<stride>(weight, pip, pjp, pjf);

    /* SOLU
       -----*/
    cs_real_t pif_up[stride], pjf_up[stride];
    cs_real_t hybrid_blend_interp;

    cs_solu_f_val_strided<stride>(cell_ceni,
                                  i_face_cog,
                                  gradi,
                                  pi,
                                  pif_up);
    cs_solu_f_val_strided<stride>(cell_cenj,
                                  i_face_cog,
                                  gradj,
                                  pj,
                                  pjf_up);

    hybrid_blend_interp = fmin(hybrid_blend_i, hybrid_blend_j);
    for (int isou = 0; isou < stride; isou++) {
      pif[isou] =   hybrid_blend_interp      *pif[isou]
                 + (1. - hybrid_blend_interp)*pif_up[isou];
      pjf[isou] =   hybrid_blend_interp      *pjf[isou]
                 + (1. - hybrid_blend_interp)*pjf_up[isou];
    }
  }
  else {

    /* Second order
       ------------*/

    cs_solu_f_val_strided<stride>(cell_ceni,
                                  i_face_cog,
                                  gradi,
                                  pi,
                                  pif);
    cs_solu_f_val_strided<stride>(cell_cenj,
                                  i_face_cog,
                                  gradj,
                                  pj,
                                  pjf);

  }

  /* Blending
     --------*/

  cs_blend_f_val_strided<stride>(blencp, pi, pif);
  cs_blend_f_val_strided<stride>(blencp, pj, pjf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     bldfrp          reconstruction blending factor
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of second order scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell j
 * \param[in]     i_face_u_normal face unit normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     diipf           distance II'
 * \param[in]     djjpf           distance JJ'
 * \param[in]     i_massflux      mass flux at face ij
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     gradupi         upwind gradient at cell i
 * \param[in]     gradupj         upwind gradient at cell j
 * \param[in]     gradsti         slope test gradient at cell i
 * \param[in]     gradstj         slope test gradient at cell j
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     pia             old value at cell i
 * \param[in]     pja             old value at cell j
 * \param[out]    pifri           contribution of i to flux from i to j
 * \param[out]    pifrj           contribution of i to flux from j to i
 * \param[out]    pjfri           contribution of j to flux from i to j
 * \param[out]    pjfrj           contribution of j to flux from j to i
 * \param[out]    pip             reconstructed value at cell i
 * \param[out]    pjp             reconstructed value at cell j
 * \param[out]    pipr            relaxed reconstructed value at cell i
 * \param[out]    pjpr            relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_steady_slope_test(bool              *upwind_switch,
                          const int          iconvp,
                          const cs_real_t    bldfrp,
                          const int          ischcp,
                          const double       relaxp,
                          const double       blencp,
                          const double       blend_st,
                          const cs_real_t    weight,
                          const cs_real_t    i_dist,
                          const cs_real_t    cell_ceni[3],
                          const cs_real_t    cell_cenj[3],
                          const cs_nreal_t   i_face_u_normal[3],
                          const cs_real_t    i_face_cog[3],
                          const cs_rreal_t   diipf[3],
                          const cs_rreal_t   djjpf[3],
                          const cs_real_t    i_massflux,
                          const cs_real_t    gradi[3],
                          const cs_real_t    gradj[3],
                          const cs_real_t    gradupi[3],
                          const cs_real_t    gradupj[3],
                          const cs_real_t    gradsti[3],
                          const cs_real_t    gradstj[3],
                          const cs_real_t    pi,
                          const cs_real_t    pj,
                          const cs_real_t    pia,
                          const cs_real_t    pja,
                          cs_real_t         *pifri,
                          cs_real_t         *pifrj,
                          cs_real_t         *pjfri,
                          cs_real_t         *pjfrj,
                          cs_real_t         *pip,
                          cs_real_t         *pjp,
                          cs_real_t         *pipr,
                          cs_real_t         *pjpr)
{
  cs_real_t pir, pjr;
  cs_real_t recoi, recoj;
  cs_real_t testij, tesqck;

  *upwind_switch = false;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);

  cs_i_relax_c_val(relaxp,
                   pia,
                   pja,
                   recoi,
                   recoj,
                   pi,
                   pj,
                   &pir,
                   &pjr,
                   pipr,
                   pjpr);

  /* Convection slope test is needed only when iconv >0 */
  if (iconvp > 0) {
    cs_slope_test(pi,
                  pj,
                  i_dist,
                  i_face_u_normal,
                  gradi,
                  gradj,
                  gradsti,
                  gradstj,
                  i_massflux,
                  &testij,
                  &tesqck);

    if (ischcp==1) {

      /* Centered
         --------*/

      cs_centered_f_val(weight,
                        *pip,
                        *pjpr,
                        pifrj);
      cs_centered_f_val(weight,
                        *pipr,
                        *pjp,
                        pifri);
      cs_centered_f_val(weight,
                        *pipr,
                        *pjp,
                        pjfri);
      cs_centered_f_val(weight,
                        *pip,
                        *pjpr,
                        pjfrj);

    } else if (ischcp == 0) {

      /* Second order
         ------------*/

      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradi,
                    pi,
                    pifrj);
      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradi,
                    pir,
                    pifri);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradj,
                    pj,
                    pjfri);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradj,
                    pjr,
                    pjfrj);

    } else {

      /* SOLU
         -----*/

      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradupi,
                    pi,
                    pifrj);
      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradupi,
                    pir,
                    pifri);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradupj,
                    pj,
                    pjfri);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradupj,
                    pjr,
                    pjfrj);
    }


    /* Slope test: percentage of upwind
       ----------------------------------*/

    if (tesqck <= 0. || testij <= 0.) {

      cs_blend_f_val(blend_st,
                     pi,
                     pifrj);
      cs_blend_f_val(blend_st,
                     pir,
                     pifri);
      cs_blend_f_val(blend_st,
                     pj,
                     pjfri);
      cs_blend_f_val(blend_st,
                     pjr,
                     pjfrj);

      *upwind_switch = true;

    }


    /* Blending
     --------*/

    cs_blend_f_val(blencp,
                   pi,
                   pifrj);
    cs_blend_f_val(blencp,
                   pir,
                   pifri);
    cs_blend_f_val(blencp,
                   pj,
                   pjfri);
    cs_blend_f_val(blencp,
                   pjr,
                   pjfrj);

  /* If iconv=0 p*fr* are useless */
  }
  else {
    cs_upwind_f_val(pi, pifrj);
    cs_upwind_f_val(pir, pifri);
    cs_upwind_f_val(pj, pjfri);
    cs_upwind_f_val(pjr, pjfrj);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     bldfrp          reconstruction blending factor
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of second order scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_u_normal face unit normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     diipf           distance II'
 * \param[in]     djjpf           distance JJ'
 * \param[in]     i_massflux      mass flux at face ij
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     pia             old value at cell i
 * \param[in]     pja             old value at cell j
 * \param[out]    pifri           contribution of i to flux from i to j
 * \param[out]    pifrj           contribution of i to flux from j to i
 * \param[out]    pjfri           contribution of j to flux from i to j
 * \param[out]    pjfrj           contribution of j to flux from j to i
 * \param[out]    pip             reconstructed value at cell i
 * \param[out]    pjp             reconstructed value at cell j
 * \param[out]    pipr            relaxed reconstructed value at cell i
 * \param[out]    pjpr            relaxed reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_steady_slope_test_strided(bool             *upwind_switch,
                                  int               iconvp,
                                  cs_real_t         bldfrp,
                                  int               ischcp,
                                  double            relaxp,
                                  double            blencp,
                                  double            blend_st,
                                  cs_real_t         weight,
                                  cs_real_t         i_dist,
                                  const cs_real_t   cell_ceni[3],
                                  const cs_real_t   cell_cenj[3],
                                  const cs_nreal_t  i_face_u_normal[3],
                                  const cs_real_t   i_face_cog[3],
                                  const cs_rreal_t  diipf[3],
                                  const cs_rreal_t  djjpf[3],
                                  cs_real_t         i_massflux,
                                  const cs_real_t   gradi[stride][3],
                                  const cs_real_t   gradj[stride][3],
                                  const cs_real_t   grdpai[stride][3],
                                  const cs_real_t   grdpaj[stride][3],
                                  const cs_real_t   pi[stride],
                                  const cs_real_t   pj[stride],
                                  const cs_real_t   pia[stride],
                                  const cs_real_t   pja[stride],
                                  cs_real_t         pifri[stride],
                                  cs_real_t         pifrj[stride],
                                  cs_real_t         pjfri[stride],
                                  cs_real_t         pjfrj[stride],
                                  cs_real_t         pip[stride],
                                  cs_real_t         pjp[stride],
                                  cs_real_t         pipr[stride],
                                  cs_real_t         pjpr[stride])
{
  cs_real_t pir[stride], pjr[stride];
  cs_real_t recoi[stride], recoj[stride];
  cs_real_t testij, tesqck;
  int isou;

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  cs_i_relax_c_val_strided<stride>(relaxp,
                                   pia,
                                   pja,
                                   recoi,
                                   recoj,
                                   pi,
                                   pj,
                                   pir,
                                   pjr,
                                   pipr,
                                   pjpr);

  /* Convection slope test is needed only when iconv >0 */
  if (iconvp > 0) {
    cs_slope_test_strided<stride>(pi,
                                  pj,
                                  i_dist,
                                  i_face_u_normal,
                                  gradi,
                                  gradj,
                                  grdpai,
                                  grdpaj,
                                  i_massflux,
                                  &testij,
                                  &tesqck);

    for (isou = 0; isou < stride; isou++) {
      if (ischcp==1) {

        /* Centered
           --------*/

        cs_centered_f_val(weight, pip[isou], pjpr[isou], &pifrj[isou]);
        cs_centered_f_val(weight, pipr[isou], pjp[isou], &pifri[isou]);
        cs_centered_f_val(weight, pipr[isou], pjp[isou], &pjfri[isou]);
        cs_centered_f_val(weight, pip[isou], pjpr[isou], &pjfrj[isou]);

      }
      else {

        /* Second order
           ------------*/

        cs_solu_f_val(cell_ceni, i_face_cog, gradi[isou], pi[isou],
                      &pifrj[isou]);
        cs_solu_f_val(cell_ceni, i_face_cog, gradi[isou], pir[isou],
                      &pifri[isou]);
        cs_solu_f_val(cell_cenj, i_face_cog, gradj[isou], pj[isou],
                      &pjfri[isou]);
        cs_solu_f_val(cell_cenj, i_face_cog, gradj[isou], pjr[isou],
                      &pjfrj[isou]);

      }

    }

    /* Slope test: percentage of upwind
       ----------------------------------*/

    if (tesqck <= 0. || testij <= 0.) {

      cs_blend_f_val_strided<stride>(blend_st, pi, pifrj);
      cs_blend_f_val_strided<stride>(blend_st, pir, pifri);
      cs_blend_f_val_strided<stride>(blend_st, pj, pjfri);
      cs_blend_f_val_strided<stride>(blend_st, pjr, pjfrj);

      *upwind_switch = true;

    }

    /* Blending
       --------*/

    cs_blend_f_val_strided<stride>(blencp, pi, pifrj);
    cs_blend_f_val_strided<stride>(blencp, pir, pifri);
    cs_blend_f_val_strided<stride>(blencp, pj, pjfri);
    cs_blend_f_val_strided<stride>(blencp, pjr, pjfrj);

   /* If iconv=0 p*fr* are useless */
  }
  else {
    for (isou = 0; isou < stride; isou++) {
      cs_upwind_f_val(pi[isou], &pifrj[isou]);
      cs_upwind_f_val(pir[isou], &pifri[isou]);
      cs_upwind_f_val(pj[isou], &pjfri[isou]);
      cs_upwind_f_val(pjr[isou], &pjfrj[isou]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     bldfrp          reconstruction blending factor
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of second order scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_u_normal face unit normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     diipf           distance II'
 * \param[in]     djjpf           distance JJ'
 * \param[in]     i_massflux      mass flux at face ij
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     gradupi         upwind gradient at cell i
 * \param[in]     gradupj         upwind gradient at cell j
 * \param[in]     gradsti         slope test gradient at cell i
 * \param[in]     gradstj         slope test gradient at cell j
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[out]    pif             contribution of i to flux from i to j
 * \param[out]    pjf             contribution of j to flux from i to j
 * \param[out]    pip             reconstructed value at cell i
 * \param[out]    pjp             reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_slope_test(bool               *upwind_switch,
                            const int           iconvp,
                            const cs_real_t     bldfrp,
                            const int           ischcp,
                            const double        blencp,
                            const double        blend_st,
                            const cs_real_t     weight,
                            const cs_real_t     i_dist,
                            const cs_real_3_t   cell_ceni,
                            const cs_real_3_t   cell_cenj,
                            const cs_nreal_3_t  i_face_u_normal,
                            const cs_real_3_t   i_face_cog,
                            const cs_rreal_3_t  diipf,
                            const cs_rreal_3_t  djjpf,
                            const cs_real_t     i_massflux,
                            const cs_real_3_t   gradi,
                            const cs_real_3_t   gradj,
                            const cs_real_3_t   gradupi,
                            const cs_real_3_t   gradupj,
                            const cs_real_3_t   gradsti,
                            const cs_real_3_t   gradstj,
                            const cs_real_t     pi,
                            const cs_real_t     pj,
                            cs_real_t          *pif,
                            cs_real_t          *pjf,
                            cs_real_t          *pip,
                            cs_real_t          *pjp)
{
  CS_UNUSED(blend_st);

  cs_real_t recoi, recoj;
  cs_real_t testij, tesqck;

  *upwind_switch = false;

  cs_i_compute_quantities(bldfrp,
                          diipf,
                          djjpf,
                          gradi,
                          gradj,
                          pi,
                          pj,
                          &recoi,
                          &recoj,
                          pip,
                          pjp);

  /* Convection slope test is needed only when iconv >0 */
  if (iconvp > 0) {
    cs_slope_test(pi,
                  pj,
                  i_dist,
                  i_face_u_normal,
                  gradi,
                  gradj,
                  gradsti,
                  gradstj,
                  i_massflux,
                  &testij,
                  &tesqck);

    if (ischcp==1) {

      /* Centered
         --------*/

      cs_centered_f_val(weight,
                        *pip,
                        *pjp,
                        pif);
      cs_centered_f_val(weight,
                        *pip,
                        *pjp,
                        pjf);

    } else if (ischcp == 0) {

      /* Original SOLU
         --------------*/

      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradi,
                    pi,
                    pif);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradj,
                    pj,
                    pjf);

    } else {

      /* SOLU
         -----*/

      cs_solu_f_val(cell_ceni,
                    i_face_cog,
                    gradupi,
                    pi,
                    pif);
      cs_solu_f_val(cell_cenj,
                    i_face_cog,
                    gradupj,
                    pj,
                    pjf);

    }

    /* Slope test: percentage of upwind
       -------------------------------- */

    if (tesqck<=0. || testij<=0.) {

      cs_blend_f_val(blend_st,
                     pi,
                     pif);
      cs_blend_f_val(blend_st,
                     pj,
                     pjf);

      *upwind_switch = true;

    }

    /* Blending
       --------*/

    cs_blend_f_val(blencp,
                   pi,
                   pif);
    cs_blend_f_val(blencp,
                   pj,
                   pjf);

  /* If iconv=0 p*f are useless */
  } else {
    cs_upwind_f_val(pi,
                    pif);
    cs_upwind_f_val(pj,
                    pjf);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine the upwind and downwind sides of an internal face and
 *        matching cell indices.
 *
 * \param[in]     ii              index of cell (0)
 * \param[in]     jj              index of cell (1)
 * \param[in]     i_massflux      mass flux at face ij
 * \param[out]    ic              index of central cell (upwind w.r.t. the face)
 * \param[out]    id              index of downwind cell
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_central_downwind_cells(const cs_lnum_t    ii,
                          const cs_lnum_t    jj,
                          const cs_real_t    i_massflux,
                          cs_lnum_t         *ic,
                          cs_lnum_t         *id)
{
  if (i_massflux >= 0.) {
    *ic = ii;
    *id = jj;
  } else {
    *ic = jj;
    *id = ii;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the convection flux
 *        computation in case of an unsteady algorithm and using NVD schemes.
 *
 * \param[in]     limiter         choice of the NVD scheme
 * \param[in]     beta            proportion of second order scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     cell_cen_c      center of gravity coordinates of central cell
 * \param[in]     cell_cen_d      center of gravity coordinates of downwind cell
 * \param[in]     i_face_u_normal unit normal of face ij
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     gradv_c         gradient at central cell
 * \param[in]     p_c             value at central cell
 * \param[in]     p_d             value at downwind cell
 * \param[in]     local_max_c     local maximum of variable
 * \param[in]     local_min_c     local minimum of variable
 * \param[in]     courant_c       central cell Courant number
 * \param[out]    pif             contribution of i to flux from i to j
 * \param[out]    pjf             contribution of j to flux from i to j
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_nvd(const cs_nvd_type_t  limiter,
                     const double         beta,
                     const cs_real_3_t    cell_cen_c,
                     const cs_real_3_t    cell_cen_d,
                     const cs_nreal_3_t   i_face_u_normal,
                     const cs_real_3_t    i_face_cog,
                     const cs_real_3_t    gradv_c,
                     const cs_real_t      p_c,
                     const cs_real_t      p_d,
                     const cs_real_t      local_max_c,
                     const cs_real_t      local_min_c,
                     const cs_real_t      courant_c,
                     cs_real_t           *pif,
                     cs_real_t           *pjf)
{
  /* distance between face center and central cell center */
  /* Distance between face center and central cell center */
  cs_real_t dist_fc = cs_math_3_distance(cell_cen_c, i_face_cog);

  /* Unit vector and distance between central and downwind cells centers */
  cs_real_t diff[3] = {cell_cen_d[0] - cell_cen_c[0],
                       cell_cen_d[1] - cell_cen_c[1],
                       cell_cen_d[2] - cell_cen_c[2]};

  cs_real_t dist_dc = cs_math_3_norm(diff);
  cs_real_t invl = 1./dist_dc;

  cs_real_t ndc[3] = {invl*diff[0],
                      invl*diff[1],
                      invl*diff[2]};

  /* Place the upwind point on the line that joins
     the two cells on the upwind side and the same
     distance as that between the two cells */
  const cs_real_t dist_cu = dist_dc;
  const cs_real_t dist_du = dist_dc + dist_cu;

  /* Compute the property on the upwind assuming a parabolic
     variation of the property between the two cells */
  const cs_real_t gradc = cs_math_3_dot_product(gradv_c, ndc);

  const cs_real_t grad2c = ((p_d - p_c)/dist_dc - gradc)/dist_dc;

  cs_real_t p_u = p_c + (grad2c*dist_cu - gradc)*dist_cu;
  p_u = cs::max(cs::min(p_u, local_max_c), local_min_c);

  /* Compute the normalised distances */
  const cs_real_t nvf_r_f = (dist_fc+dist_cu)/dist_du;
  const cs_real_t nvf_r_c = dist_cu/dist_du;

  /* Check for the bounds of the NVD diagram and compute the face
     property according to the selected NVD scheme */
  const cs_real_t _small
    = cs_math_epzero * (cs::abs(p_u) + cs::abs(p_c) + cs::abs(p_d));

  if (cs::abs(p_d-p_u) <= _small) {
    *pif = p_c;
    *pjf = p_c;
  }
  else {
    const cs_real_t nvf_p_c = (p_c - p_u)/(p_d - p_u);

    if (nvf_p_c <= 0. || nvf_p_c >= 1.) {
      *pif = p_c;
      *pjf = p_c;
    }
    else {
      cs_real_t nvf_p_f;

      /* Highly compressive NVD scheme for VOF */
      if (limiter >= CS_NVD_VOF_HRIC) {
        nvf_p_f = cs_nvd_vof_scheme_scalar(limiter,
                                           i_face_u_normal,
                                           nvf_p_c,
                                           nvf_r_f,
                                           nvf_r_c,
                                           gradv_c,
                                           courant_c);
      } else { /* Regular NVD scheme */
        nvf_p_f = cs_nvd_scheme_scalar(limiter,
                                       nvf_p_c,
                                       nvf_r_f,
                                       nvf_r_c);
      }

      *pif = p_u + nvf_p_f*(p_d - p_u);
      *pif = cs::max(cs::min(*pif, local_max_c), local_min_c);

      cs_blend_f_val(beta,
                     p_c,
                     pif);

      *pjf = *pif;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and using slope tests.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     bldfrp          reconstruction blending factor
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of second order scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_u_normal face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     diipf           distance II'
 * \param[in]     djjpf           distance JJ'
 * \param[in]     i_massflux      mass flux at face ij
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[out]    pif             contribution of i to flux from i to j
 * \param[out]    pjf             contribution of j to flux from i to j
 * \param[out]    pip             reconstructed value at cell i
 * \param[out]    pjp             reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_i_cd_unsteady_slope_test_strided(bool             *upwind_switch,
                                    int               iconvp,
                                    cs_real_t         bldfrp,
                                    int               ischcp,
                                    double            blencp,
                                    double            blend_st,
                                    cs_real_t         weight,
                                    cs_real_t         i_dist,
                                    const cs_real_t   cell_ceni[3],
                                    const cs_real_t   cell_cenj[3],
                                    const cs_nreal_t  i_face_u_normal[3],
                                    const cs_real_t   i_face_cog[3],
                                    const cs_rreal_t  diipf[3],
                                    const cs_rreal_t  djjpf[3],
                                    cs_real_t         i_massflux,
                                    const cs_real_t   gradi[stride][3],
                                    const cs_real_t   gradj[stride][3],
                                    const cs_real_t   grdpai[stride][3],
                                    const cs_real_t   grdpaj[stride][3],
                                    const cs_real_t   pi[stride],
                                    const cs_real_t   pj[stride],
                                    cs_real_t         pif[stride],
                                   cs_real_t          pjf[stride],
                                   cs_real_t          pip[stride],
                                   cs_real_t          pjp[stride])
{
  cs_real_t recoi[stride], recoj[stride];
  cs_real_t testij, tesqck;
  int isou;

  cs_i_compute_quantities_strided<stride>(bldfrp,
                                          diipf,
                                          djjpf,
                                          gradi,
                                          gradj,
                                          pi,
                                          pj,
                                          recoi,
                                          recoj,
                                          pip,
                                          pjp);

  /* Convection slope test is needed only when iconv >0 */
  if (iconvp > 0) {
    cs_slope_test_strided<stride>(pi,
                                  pj,
                                  i_dist,
                                  i_face_u_normal,
                                  gradi,
                                  gradj,
                                  grdpai,
                                  grdpaj,
                                  i_massflux,
                                  &testij,
                                  &tesqck);

    for (isou = 0; isou < stride; isou++) {

      if (ischcp==1) {

        /* Centered
           --------*/

        cs_centered_f_val(weight,
                          pip[isou],
                          pjp[isou],
                          &pif[isou]);
        cs_centered_f_val(weight,
                          pip[isou],
                          pjp[isou],
                          &pjf[isou]);

      }
      else {

        /* Second order
           ------------*/

        cs_solu_f_val(cell_ceni,
                      i_face_cog,
                      gradi[isou],
                      pi[isou],
                      &pif[isou]);
        cs_solu_f_val(cell_cenj,
                      i_face_cog,
                      gradj[isou],
                      pj[isou],
                      &pjf[isou]);
      }

    }

    /* Slope test activated: percentage of upwind */
    if (tesqck <= 0. || testij <= 0.) {

      /* Upwind
         --------*/

      cs_blend_f_val_strided<stride>(blend_st, pi, pif);
      cs_blend_f_val_strided<stride>(blend_st, pj,  pjf);

      *upwind_switch = true;
    }

    /* Blending
       --------*/

    cs_blend_f_val_strided<stride>(blencp, pi, pif);
    cs_blend_f_val_strided<stride>(blencp, pj, pjf);

  /* If iconv=0 p*fr* are useless */
  }
  else {

    for (isou = 0; isou < stride; isou++) {
      cs_upwind_f_val(pi[isou], &pif[isou]);
      cs_upwind_f_val(pj[isou], &pjf[isou]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' at boundary cell i.
 *
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[out]    recoi    reconstruction at cell i
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_compute_quantities(const cs_rreal_t   diipb[3],
                        const cs_real_t    gradi[3],
                        const cs_real_t    bldfrp,
                        cs_real_t         *recoi)
{
  *recoi = bldfrp * (  gradi[0]*diipb[0]
                     + gradi[1]*diipb[1]
                     + gradi[2]*diipb[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' at boundary cell i.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[out]    recoi    reconstruction at cell i
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_compute_quantities_strided(const cs_rreal_t  diipb[3],
                                const cs_real_t   gradi[stride][3],
                                const cs_real_t   bldfrp,
                                cs_real_t         recoi[stride])
{
  for (int isou = 0; isou < stride; isou++) {
    recoi[isou] = bldfrp * (  gradi[isou][0]*diipb[0]
                            + gradi[isou][1]*diipb[1]
                            + gradi[isou][2]*diipb[2]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute relaxed values at boundary cell i.
 *
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[in]     recoi    reconstruction at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_relax_c_val(const double     relaxp,
                 const cs_real_t  pi,
                 const cs_real_t  pia,
                 const cs_real_t  recoi,
                 cs_real_t       *pir,
                 cs_real_t       *pipr)
{
  *pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;
  *pipr = *pir + recoi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute relaxed values at boundary cell i.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[in]     recoi    reconstruction at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_relax_c_val_strided(const double     relaxp,
                         const cs_real_t  pi[stride],
                         const cs_real_t  pia[stride],
                         const cs_real_t  recoi[stride],
                         cs_real_t        pir[stride],
                         cs_real_t        pipr[stride])
{
  for (int isou = 0; isou < stride; isou++) {
    pir[isou]  = pi[isou]/relaxp - (1.-relaxp)/relaxp*pia[isou];
    pipr[isou] = pir[isou] + recoi[isou];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux can be either an upwind flux or an
 * imposed value.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-scheme,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     inc          Not an increment flag
 * \param[in]     bc_type      type of boundary face
 * \param[in]     icvfli       imposed convective flux flag
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection operator
 * \param[in]     coefbp       implicit boundary coefficient for convection operator
 * \param[in]     coface       explicit imposed convective flux value (0 otherwise).
 * \param[in]     cofbce       implicit part of imp. conv. flux value
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in]     xcpp         specific heat value if the scalar is the temperature,
 *                             1 otherwise
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_imposed_conv_flux(int         iconvp,
                       cs_real_t   thetap,
                       int         imasac,
                       int         inc,
                       int         bc_type,
                       int         icvfli,
                       cs_real_t   pi,
                       cs_real_t   pir,
                       cs_real_t   pipr,
                       cs_real_t   coefap,
                       cs_real_t   coefbp,
                       cs_real_t   coface,
                       cs_real_t   cofbce,
                       cs_real_t   b_massflux,
                       cs_real_t   xcpp,
                       cs_real_t  *flux)
{
  cs_real_t flui, fluj, pfac;

  /* Computed convective flux */

  if (icvfli == 0) {

    /* Remove decentering for coupled faces */
    if (bc_type == CS_COUPLED_FD) {
      flui = 0.0;
      fluj = b_massflux;
    } else {
      flui = 0.5*(b_massflux +fabs(b_massflux));
      fluj = 0.5*(b_massflux -fabs(b_massflux));
    }

    pfac  = inc*coefap + coefbp*pipr;
    *flux += iconvp*xcpp*(thetap*(flui*pir + fluj*pfac) -imasac*( b_massflux*pi));

  /* Imposed convective flux */

  } else {

    pfac = inc*coface + cofbce*pipr;
    *flux += iconvp*xcpp*(-imasac*(b_massflux*pi) + thetap*(pfac));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux can be either an upwind flux or an
 * imposed value.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-scheme,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     inc          Not an increment flag
 * \param[in]     bc_type      type of boundary face
 * \param[in]     icvfli       imposed convective flux flag
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coface       explicit imposed convective flux value (0 otherwise).
 * \param[in]     cofbce       implicit part of imp. conv. flux value
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in]     pfac         value at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_imposed_conv_flux_strided(int              iconvp,
                               cs_real_t        thetap,
                               int              imasac,
                               int              inc,
                               int              bc_type,
                               int              icvfli,
                               const cs_real_t  pi[stride],
                               const cs_real_t  pir[stride],
                               const cs_real_t  pipr[stride],
                               const cs_real_t  coface[stride],
                               const cs_real_t  cofbce[stride][stride],
                               cs_real_t        b_massflux,
                               cs_real_t        pfac[stride],
                               cs_real_t        flux[stride])
{
  cs_real_t flui, fluj;

  /* Computed convective flux */

  if (icvfli == 0) {

    /* Remove decentering for coupled faces */
    if (bc_type == CS_COUPLED_FD) {
      flui = 0.0;
      fluj = b_massflux;
    }
    else {
      flui = 0.5*(b_massflux +fabs(b_massflux));
      fluj = 0.5*(b_massflux -fabs(b_massflux));
    }

    for (int isou = 0; isou < stride; isou++)
      flux[isou] += iconvp*( thetap*(flui*pir[isou] + fluj*pfac[isou])
                            - imasac*b_massflux*pi[isou]);

  /* Imposed convective flux */

  }
  else {

    for (int isou = 0; isou < stride; isou++) {
      pfac[isou]  = inc*coface[isou];
      for (int jsou = 0; jsou < stride; jsou++) {
        pfac[isou] += cofbce[isou][jsou]*pipr[jsou];
      }
      flux[isou] += iconvp*(  thetap*pfac[isou]
                            - imasac*b_massflux*pi[isou]);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux is a pure upwind flux.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-scheme,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     inc          Not an increment flag
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection operator
 * \param[in]     coefbp       implicit boundary coefficient for convection operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in]     xcpp         specific heat value if the scalar is the
 *                             temperature, 1 otherwise
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_upwind_flux(const int        iconvp,
                 const cs_real_t  thetap,
                 const int        imasac,
                 const int        inc,
                 const int        bc_type,
                 const cs_real_t  pi,
                 const cs_real_t  pir,
                 const cs_real_t  pipr,
                 const cs_real_t  coefap,
                 const cs_real_t  coefbp,
                 const cs_real_t  b_massflux,
                 const cs_real_t  xcpp,
                 cs_real_t       *flux)
{
  cs_real_t flui, fluj, pfac;

  /* Remove decentering for coupled faces */
  if (bc_type == CS_COUPLED_FD) {
    flui = 0.0;
    fluj = b_massflux;
  } else {
    flui = 0.5*(b_massflux +fabs(b_massflux));
    fluj = 0.5*(b_massflux -fabs(b_massflux));
  }

  pfac  = inc*coefap + coefbp*pipr;
  *flux += iconvp*xcpp*(thetap*(flui*pir + fluj*pfac) -imasac*( b_massflux*pi));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux is a pure upwind flux.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-scheme,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[out]    pfac         contribution from BCs
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_upwind_flux_strided(int              iconvp,
                         cs_real_t        thetap,
                         int              imasac,
                         int              bc_type,
                         const cs_real_t  pi[stride],
                         const cs_real_t  pir[stride],
                         const cs_real_t  b_massflux,
                         const cs_real_t  pfac[stride],
                         cs_real_t        flux[stride])
{
  cs_real_t flui, fluj;

  /* Remove decentering for coupled faces */
  if (bc_type == CS_COUPLED_FD) {
    flui = 0.0;
    fluj = b_massflux;
  } else {
    flui = 0.5*(b_massflux +fabs(b_massflux));
    fluj = 0.5*(b_massflux -fabs(b_massflux));
  }

  for (int isou = 0; isou < stride; isou++)
    flux[isou] += iconvp*(  thetap*(flui*pir[isou] + fluj*pfac[isou])
                          - imasac*b_massflux*pi[isou]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at boundary face.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-scheme,
 * \param[in]     inc      Not an increment flag
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     cofafp   explicit boundary coefficient for diffusion operator
 * \param[in]     cofbfp   implicit boundary coefficient for diffusion operator
 * \param[in]     b_visc   boundary face surface
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_diff_flux(const int        idiffp,
               const cs_real_t  thetap,
               const int        inc,
               const cs_real_t  pipr,
               const cs_real_t  cofafp,
               const cs_real_t  cofbfp,
               const cs_real_t  b_visc,
               cs_real_t       *flux)
{
  cs_real_t pfacd = inc*cofafp + cofbfp*pipr;
  *flux += idiffp*thetap*b_visc*pfacd;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at boundary face.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-scheme,
 * \param[in]     b_visc   boundary face surface
 * \param[in]     pfacd    value at boundary face
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_diff_flux_strided(int              idiffp,
                       cs_real_t        thetap,
                       cs_real_t        b_visc,
                       const cs_real_t  pfacd[stride],
                       cs_real_t        flux[stride])
{
  for (int isou = 0; isou < stride; isou++)
    flux[isou] += idiffp*thetap*b_visc*pfacd[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_cd_steady(const cs_real_t     bldfrp,
               const double        relaxp,
               const cs_rreal_3_t  diipb,
               const cs_real_3_t   gradi,
               const cs_real_t     pi,
               const cs_real_t     pia,
               cs_real_t          *pir,
               cs_real_t          *pipr)
{
  cs_real_t recoi;

  cs_b_compute_quantities(diipb,
                          gradi,
                          bldfrp,
                          &recoi);

  cs_b_relax_c_val(relaxp,
                   pi,
                   pia,
                   recoi,
                   pir,
                   pipr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_cd_steady_strided(cs_real_t         bldfrp,
                       double            relaxp,
                       const cs_rreal_t  diipb[3],
                       const cs_real_t   gradi[stride][3],
                       const cs_real_t   pi[stride],
                       const cs_real_t   pia[stride],
                       cs_real_t         pir[stride],
                       cs_real_t         pipr[stride])
{
  cs_real_t recoi[stride];

  cs_b_compute_quantities_strided<stride>(diipb,
                                          gradi,
                                          bldfrp,
                                          recoi);

  cs_b_relax_c_val_strided<stride>(relaxp,
                                   pi,
                                   pia,
                                   recoi,
                                   pir,
                                   pipr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of an unsteady algorithm.
 *
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[out]    pip      reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_cd_unsteady(const cs_real_t    bldfrp,
                 const cs_rreal_t   diipb[3],
                 const cs_real_t    gradi[3],
                 const cs_real_t    pi,
                 cs_real_t         *pip)
{
  cs_real_t recoi;

  cs_b_compute_quantities(diipb,
                          gradi,
                          bldfrp,
                          &recoi);

  *pip = pi + recoi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     bldfrp   reconstruction blending factor
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[out]    pip      reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_cd_unsteady_strided(cs_real_t         bldfrp,
                         const cs_rreal_t  diipb[3],
                         const cs_real_t   gradi[stride][3],
                         const cs_real_t   pi[stride],
                         cs_real_t         pip[stride])
{
  cs_real_t recoi[stride];

  cs_b_compute_quantities_strided<stride>(diipb,
                                          gradi,
                                          bldfrp,
                                          recoi);

  for (int isou = 0; isou< stride; isou++)
    pip[isou] = pi[isou] + recoi[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at an internal coupling face.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[in]     b_visc   equivalent exchange coefficient at an internal
 *                         coupling face
 * \param[in,out] fluxi    flux at internal coupling face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_b_diff_flux_coupling(int         idiffp,
                        cs_real_t   pi,
                        cs_real_t   pj,
                        cs_real_t   b_visc,
                        cs_real_t  *fluxi)
{
  *fluxi += idiffp*b_visc*(pi - pj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at an internal coupling face for a vector.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[in]     b_visc   equivalent exchange coefficient at an internal
 *                         coupling face
 * \param[in,out] fluxi    flux at internal coupling face
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
CS_F_HOST_DEVICE inline static void
cs_b_diff_flux_coupling_strided(int              idiffp,
                                const cs_real_t  pi[stride],
                                const cs_real_t  pj[stride],
                                cs_real_t        b_visc,
                                cs_real_t        fluxi[stride])
{
  for (int k = 0; k < stride; k++)
    fluxi[k] += idiffp*b_visc*(pi[k] - pj[k]);
}

/*----------------------------------------------------------------------------*/

#endif /* __CS_CONVECTION_DIFFUSION_PRIV_H__ */
