#ifndef __CS_CONVECTION_DIFFUSION_H__
#define __CS_CONVECTION_DIFFUSION_H__

/*============================================================================
 * Convection-diffusion operators.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * NVD/TVD Advection Scheme
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_NVD_GAMMA      = 0,      /* GAMMA            */
  CS_NVD_SMART      = 1,      /* SMART            */
  CS_NVD_CUBISTA    = 2,      /* CUBISTA          */
  CS_NVD_SUPERBEE   = 3,      /* SUPERBEE         */
  CS_NVD_MUSCL      = 4,      /* MUSCL            */
  CS_NVD_MINMOD     = 5,      /* MINMOD           */
  CS_NVD_CLAM       = 6,      /* CLAM             */
  CS_NVD_STOIC      = 7,      /* STOIC            */
  CS_NVD_OSHER      = 8,      /* OSHER            */
  CS_NVD_WASEB      = 9,      /* WASEB            */
  CS_NVD_VOF_HRIC   = 10,     /* M-HRIC for VOF   */
  CS_NVD_VOF_CICSAM = 11,     /* M-CICSAM for VOF */
  CS_NVD_VOF_STACS  = 12,     /* STACS for VOF    */
  CS_NVD_N_TYPES    = 13      /* number of NVD schemes */

} cs_nvd_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute slope test criteria at internal face between cell i and j.
 *
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     distf           distance IJ.Nij
 * \param[in]     srfan           face surface
 * \param[in]     i_face_normal   face normal
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     i_massflux      mass flux at face (from i to j)
 * \param[out]    testij          value of slope test first criterion
 * \param[out]    tesqck          value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_slope_test(const cs_real_t     pi,
              const cs_real_t     pj,
              const cs_real_t     distf,
              const cs_real_t     srfan,
              const cs_real_3_t   i_face_normal,
              const cs_real_3_t   gradi,
              const cs_real_3_t   gradj,
              const cs_real_3_t   grdpai,
              const cs_real_3_t   grdpaj,
              const cs_real_t     i_massflux,
              double             *testij,
              double             *tesqck)
{
  double testi, testj;
  double dcc, ddi, ddj;

  /* Slope test
     ----------*/

  testi = grdpai[0]*i_face_normal[0]
        + grdpai[1]*i_face_normal[1]
        + grdpai[2]*i_face_normal[2];
  testj = grdpaj[0]*i_face_normal[0]
        + grdpaj[1]*i_face_normal[1]
        + grdpaj[2]*i_face_normal[2];
  *testij = grdpai[0]*grdpaj[0]
          + grdpai[1]*grdpaj[1]
          + grdpai[2]*grdpaj[2];

  if (i_massflux>0.) {
    dcc = gradi[0]*i_face_normal[0]
        + gradi[1]*i_face_normal[1]
        + gradi[2]*i_face_normal[2];
    ddi = testi;
    ddj = (pj-pi)/distf *srfan;
  } else {
    dcc = gradj[0]*i_face_normal[0]
        + gradj[1]*i_face_normal[1]
        + gradj[2]*i_face_normal[2];
    ddi = (pj-pi)/distf *srfan;
    ddj = testj;
  }
  *tesqck = cs_math_sq(dcc) - cs_math_sq(ddi-ddj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute slope test criteria at internal face between cell i and j.
 *
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     distf           distance IJ.Nij
 * \param[in]     srfan           face surface
 * \param[in]     i_face_normal   face normal
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     i_massflux      mass flux at face (from i to j)
 * \param[out]    testij          value of slope test first criterion
 * \param[out]    tesqck          value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_slope_test_vector(const cs_real_3_t   pi,
                     const cs_real_3_t   pj,
                     const cs_real_t     distf,
                     const cs_real_t     srfan,
                     const cs_real_3_t   i_face_normal,
                     const cs_real_33_t  gradi,
                     const cs_real_33_t  gradj,
                     const cs_real_33_t  gradsti,
                     const cs_real_33_t  gradstj,
                     const cs_real_t     i_massflux,
                     cs_real_t          *testij,
                     cs_real_t          *tesqck)
{
  double testi[3], testj[3];
  double dcc[3], ddi[3], ddj[3];
  *testij = 0.;
  *tesqck = 0.;

  /* Slope test
     ----------*/
  for (int i = 0; i < 3; i++) {
    *testij += gradsti[i][0]*gradstj[i][0]
             + gradsti[i][1]*gradstj[i][1]
             + gradsti[i][2]*gradstj[i][2];

    testi[i] = gradsti[i][0]*i_face_normal[0]
             + gradsti[i][1]*i_face_normal[1]
             + gradsti[i][2]*i_face_normal[2];
    testj[i] = gradstj[i][0]*i_face_normal[0]
             + gradstj[i][1]*i_face_normal[1]
             + gradstj[i][2]*i_face_normal[2];

    if (i_massflux > 0.) {
      dcc[i] = gradi[i][0]*i_face_normal[0]
             + gradi[i][1]*i_face_normal[1]
             + gradi[i][2]*i_face_normal[2];
      ddi[i] = testi[i];
      ddj[i] = (pj[i]-pi[i])/distf *srfan;
    } else {
      dcc[i] = gradj[i][0]*i_face_normal[0]
             + gradj[i][1]*i_face_normal[1]
             + gradj[i][2]*i_face_normal[2];
      ddi[i] = (pj[i]-pi[i])/distf *srfan;
      ddj[i] = testj[i];
    }
  }

  *tesqck = cs_math_3_square_norm(dcc) - cs_math_3_square_distance(ddi, ddj);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief DEPRECATED Compute slope test criteria at internal face between cell i
 *        and j.
 *
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     distf           distance IJ.Nij
 * \param[in]     srfan           face surface
 * \param[in]     i_face_normal   face normal
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     i_massflux      mass flux at face (from i to j)
 * \param[out]    testij          value of slope test first criterion
 * \param[out]    tesqck          value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_slope_test_vector_old(const cs_real_3_t   pi,
                         const cs_real_3_t   pj,
                         const cs_real_t     distf,
                         const cs_real_t     srfan,
                         const cs_real_3_t   i_face_normal,
                         const cs_real_33_t  gradi,
                         const cs_real_33_t  gradj,
                         const cs_real_33_t  grdpai,
                         const cs_real_33_t  grdpaj,
                         const cs_real_t     i_massflux,
                         cs_real_t         testij[3],
                         cs_real_t         tesqck[3])
{
  double testi[3], testj[3];
  double dcc[3], ddi[3], ddj[3];

  /* Slope test
     ----------*/
  for (int isou = 0; isou < 3; isou++) {
    testi[isou] = grdpai[isou][0]*i_face_normal[0]
                + grdpai[isou][1]*i_face_normal[1]
                + grdpai[isou][2]*i_face_normal[2];
    testj[isou] = grdpaj[isou][0]*i_face_normal[0]
                + grdpaj[isou][1]*i_face_normal[1]
                + grdpaj[isou][2]*i_face_normal[2];
    testij[isou] = grdpai[isou][0]*grdpaj[isou][0]
                 + grdpai[isou][1]*grdpaj[isou][1]
                 + grdpai[isou][2]*grdpaj[isou][2];

    if (i_massflux>0.) {
      dcc[isou] = gradi[isou][0]*i_face_normal[0]
                + gradi[isou][1]*i_face_normal[1]
                + gradi[isou][2]*i_face_normal[2];
      ddi[isou] = testi[isou];
      ddj[isou] = (pj[isou]-pi[isou])/distf *srfan;
    } else {
      dcc[isou] = gradj[isou][0]*i_face_normal[0]
          + gradj[isou][1]*i_face_normal[1]
          + gradj[isou][2]*i_face_normal[2];
      ddi[isou] = (pj[isou]-pi[isou])/distf *srfan;
      ddj[isou] = testj[isou];
    }
    tesqck[isou] = cs_math_sq(dcc[isou]) - cs_math_sq(ddi[isou]-ddj[isou]);
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute slope test criteria at internal face between cell i and j.
 *
 * \param[in]     pi              value at cell i
 * \param[in]     pj              value at cell j
 * \param[in]     distf           distance IJ.Nij
 * \param[in]     srfan           face surface
 * \param[in]     i_face_normal   face normal
 * \param[in]     gradi           gradient at cell i
 * \param[in]     gradj           gradient at cell j
 * \param[in]     grdpai          upwind gradient at cell i
 * \param[in]     grdpaj          upwind gradient at cell j
 * \param[in]     i_massflux      mass flux at face (from i to j)
 * \param[out]    testij          value of slope test first criterion
 * \param[out]    tesqck          value of slope test second criterion
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_slope_test_tensor(const cs_real_6_t   pi,
                     const cs_real_6_t   pj,
                     const cs_real_t     distf,
                     const cs_real_t     srfan,
                     const cs_real_3_t   i_face_normal,
                     const cs_real_63_t  gradi,
                     const cs_real_63_t  gradj,
                     const cs_real_63_t  gradsti,
                     const cs_real_63_t  gradstj,
                     const cs_real_t     i_massflux,
                     cs_real_t          *testij,
                     cs_real_t          *tesqck)
{
  double testi[6], testj[6];
  double dcc[6], ddi[6], ddj[6];
  *testij = 0.;
  *tesqck = 0.;

  /* Slope test */

  for (int ij = 0; ij < 6; ij++) {
    *testij += gradsti[ij][0]*gradstj[ij][0]
             + gradsti[ij][1]*gradstj[ij][1]
             + gradsti[ij][2]*gradstj[ij][2];
    testi[ij] = gradsti[ij][0]*i_face_normal[0]
              + gradsti[ij][1]*i_face_normal[1]
              + gradsti[ij][2]*i_face_normal[2];
    testj[ij] = gradstj[ij][0]*i_face_normal[0]
              + gradstj[ij][1]*i_face_normal[1]
              + gradstj[ij][2]*i_face_normal[2];

    if (i_massflux > 0.) {
      dcc[ij] = gradi[ij][0]*i_face_normal[0]
              + gradi[ij][1]*i_face_normal[1]
              + gradi[ij][2]*i_face_normal[2];
      ddi[ij] = testi[ij];
      ddj[ij] = (pj[ij]-pi[ij])/distf *srfan;
    } else {
      dcc[ij] = gradj[ij][0]*i_face_normal[0]
              + gradj[ij][1]*i_face_normal[1]
              + gradj[ij][2]*i_face_normal[2];
      ddi[ij] = (pj[ij]-pi[ij])/distf *srfan;
      ddj[ij] = testj[ij];
    }

    *tesqck += cs_math_sq(dcc[ij]) - cs_math_sq(ddi[ij]-ddj[ij]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' and J'.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     pnd          geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_compute_quantities(const int          ircflp,
                        const double       pnd,
                        const cs_real_3_t  cell_ceni,
                        const cs_real_3_t  cell_cenj,
                        const cs_real_3_t  i_face_cog,
                        const cs_real_3_t  dijpf,
                        const cs_real_3_t  gradi,
                        const cs_real_3_t  gradj,
                        const cs_real_t    pi,
                        const cs_real_t    pj,
                        cs_real_t         *recoi,
                        cs_real_t         *recoj,
                        cs_real_t         *pip,
                        cs_real_t         *pjp)
{
  cs_real_t diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz;
  cs_real_t dpxf, dpyf, dpzf;//FIXME

  /* Recompute II' and JJ' at this level */

  diipfx = i_face_cog[0] - (cell_ceni[0] + (1.-pnd) * dijpf[0]);
  diipfy = i_face_cog[1] - (cell_ceni[1] + (1.-pnd) * dijpf[1]);
  diipfz = i_face_cog[2] - (cell_ceni[2] + (1.-pnd) * dijpf[2]);

  djjpfx = i_face_cog[0] -  cell_cenj[0] +     pnd  * dijpf[0];
  djjpfy = i_face_cog[1] -  cell_cenj[1] +     pnd  * dijpf[1];
  djjpfz = i_face_cog[2] -  cell_cenj[2] +     pnd  * dijpf[2];

  dpxf = 0.5*(gradi[0] + gradj[0]);
  dpyf = 0.5*(gradi[1] + gradj[1]);
  dpzf = 0.5*(gradi[2] + gradj[2]);

  /* reconstruction only if IRCFLP = 1 */
  *recoi = ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
  *recoj = ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);
  *pip = pi + *recoi;
  *pjp = pj + *recoj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' and J'.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     pnd          geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_compute_quantities_vector(const int          ircflp,
                               const double       pnd,
                               const cs_real_3_t  cell_ceni,
                               const cs_real_3_t  cell_cenj,
                               const cs_real_3_t  i_face_cog,
                               const cs_real_3_t  dijpf,
                               const cs_real_33_t gradi,
                               const cs_real_33_t gradj,
                               const cs_real_3_t  pi,
                               const cs_real_3_t  pj,
                               cs_real_t        recoi[3],
                               cs_real_t        recoj[3],
                               cs_real_t        pip[3],
                               cs_real_t        pjp[3])
{
  cs_real_3_t dijpfv, diipfv, djjpfv;
  cs_real_3_t dpvf;

  for (int jsou = 0; jsou < 3; jsou++)
    dijpfv[jsou] = dijpf[jsou];

  /* Recompute II' and JJ' at this level */
  for (int jsou = 0; jsou < 3; jsou++) {
    diipfv[jsou] =   i_face_cog[jsou]
                   - (cell_ceni[jsou] + (1.-pnd) * dijpfv[jsou]);
    djjpfv[jsou] =   i_face_cog[jsou]
                   - cell_cenj[jsou] + pnd  * dijpfv[jsou];
  }

  /* x-y-z components, p = u, v, w */

  for (int isou = 0; isou < 3; isou++) {

    for (int jsou = 0; jsou < 3; jsou++)
      dpvf[jsou] = 0.5*(  gradi[isou][jsou]
                        + gradj[isou][jsou]);

    /* reconstruction only if IRCFLP = 1 */

    recoi[isou] = ircflp*(  dpvf[0]*diipfv[0]
                          + dpvf[1]*diipfv[1]
                          + dpvf[2]*diipfv[2]);


    recoj[isou] = ircflp*(  dpvf[0]*djjpfv[0]
                          + dpvf[1]*djjpfv[1]
                          + dpvf[2]*djjpfv[2]);

    pip[isou] = pi[isou] + recoi[isou];

    pjp[isou] = pj[isou] + recoj[isou];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' and J'.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     pnd          geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_compute_quantities_tensor(const int          ircflp,
                               const double       pnd,
                               const cs_real_3_t  cell_ceni,
                               const cs_real_3_t  cell_cenj,
                               const cs_real_3_t  i_face_cog,
                               const cs_real_3_t  dijpf,
                               const cs_real_63_t gradi,
                               const cs_real_63_t gradj,
                               const cs_real_6_t  pi,
                               const cs_real_6_t  pj,
                               cs_real_t        recoi[6],
                               cs_real_t        recoj[6],
                               cs_real_t        pip[6],
                               cs_real_t        pjp[6])
{
  cs_real_3_t dijpfv, diipfv, djjpfv;
  cs_real_3_t dpvf;

  for (int jsou = 0; jsou < 3; jsou++)
    dijpfv[jsou] = dijpf[jsou];

  /* Recompute II' and JJ' at this level */

  for (int jsou = 0; jsou < 3; jsou++) {
    diipfv[jsou] =   i_face_cog[jsou]
                   - (cell_ceni[jsou] + (1.-pnd) * dijpfv[jsou]);
    djjpfv[jsou] =   i_face_cog[jsou]
                   - cell_cenj[jsou] + pnd  * dijpfv[jsou];
  }

  /* x-y-z components, p = u, v, w */

  for (int isou = 0; isou < 6; isou++) {

    for (int jsou = 0; jsou < 3; jsou++)
      dpvf[jsou] = 0.5*( gradi[isou][jsou]
                       + gradj[isou][jsou]);

    /* reconstruction only if IRCFLP = 1 */

    recoi[isou] = ircflp*( dpvf[0]*diipfv[0]
                         + dpvf[1]*diipfv[1]
                         + dpvf[2]*diipfv[2]);


    recoj[isou] = ircflp*( dpvf[0]*djjpfv[0]
                         + dpvf[1]*djjpfv[1]
                         + dpvf[2]*djjpfv[2]);

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

inline static void
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

inline static void
cs_i_relax_c_val_vector(const double       relaxp,
                        const cs_real_3_t  pia,
                        const cs_real_3_t  pja,
                        const cs_real_3_t  recoi,
                        const cs_real_3_t  recoj,
                        const cs_real_3_t  pi,
                        const cs_real_3_t  pj,
                        cs_real_t       pir[3],
                        cs_real_t       pjr[3],
                        cs_real_t       pipr[3],
                        cs_real_t       pjpr[3])
{
  for (int isou = 0; isou < 3; isou++) {
    pir[isou] = pi[isou] /relaxp - (1.-relaxp)/relaxp * pia[isou];
    pjr[isou] = pj[isou] /relaxp - (1.-relaxp)/relaxp * pja[isou];

    pipr[isou] = pir[isou] + recoi[isou];
    pjpr[isou] = pjr[isou] + recoj[isou];
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

inline static void
cs_i_relax_c_val_tensor(const double       relaxp,
                        const cs_real_6_t  pia,
                        const cs_real_6_t  pja,
                        const cs_real_6_t  recoi,
                        const cs_real_6_t  recoj,
                        const cs_real_6_t  pi,
                        const cs_real_6_t  pj,
                        cs_real_t          pir[6],
                        cs_real_t          pjr[6],
                        cs_real_t          pipr[6],
                        cs_real_t          pjpr[6])
{
  for (int isou = 0; isou < 6; isou++) {
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

inline static void
cs_upwind_f_val(const cs_real_t  p,
                cs_real_t       *pf)
{
  *pf = p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using an upwind scheme.
 *
 * \param[in]     p       value at cell
 * \param[out]    pf      value at face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_upwind_f_val_vector(const cs_real_3_t  p,
                       cs_real_t          pf[3])
{
  for (int isou = 0; isou < 3; isou++)
    pf[isou] = p[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using an upwind scheme.
 *
 * \param[in]     p       value at cell
 * \param[out]    pf      value at face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_upwind_f_val_tensor(const cs_real_6_t  p,
                       cs_real_t          pf[6])
{
  for (int isou = 0; isou < 6; isou++)
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

inline static void
cs_centered_f_val(const double     pnd,
                  const cs_real_t  pip,
                  const cs_real_t  pjp,
                  cs_real_t       *pf)
{
  *pf = pnd*pip + (1.-pnd)*pjp;
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

inline static void
cs_centered_f_val_vector(const double       pnd,
                         const cs_real_3_t  pip,
                         const cs_real_3_t  pjp,
                         cs_real_t          pf[3])
{
  for (int isou = 0; isou < 3; isou++)
    pf[isou] = pnd*pip[isou] + (1.-pnd)*pjp[isou];
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

inline static void
cs_centered_f_val_tensor(const double       pnd,
                         const cs_real_6_t  pip,
                         const cs_real_6_t  pjp,
                         cs_real_t          pf[6])
{
  for (int isou = 0; isou < 6; isou++)
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

inline static void
cs_solu_f_val(const cs_real_3_t  cell_cen,
              const cs_real_3_t  i_face_cog,
              const cs_real_3_t  grad,
              const cs_real_t    p,
              cs_real_t         *pf)
{
  cs_real_3_t df;

  df[0] = i_face_cog[0] - cell_cen[0];
  df[1] = i_face_cog[1] - cell_cen[1];
  df[2] = i_face_cog[2] - cell_cen[2];

  *pf = p + cs_math_3_dot_product(df, grad);
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

inline static void
cs_solu_f_val_vector(const cs_real_3_t   cell_cen,
                     const cs_real_3_t   i_face_cog,
                     const cs_real_33_t  grad,
                     const cs_real_3_t   p,
                     cs_real_t           pf[3])
{
  cs_real_3_t df;

  for (int jsou = 0; jsou < 3; jsou++)
    df[jsou] = i_face_cog[jsou] - cell_cen[jsou];

  for (int isou = 0; isou < 3; isou++) {
     pf[isou] = p[isou] + df[0]*grad[isou][0]
                        + df[1]*grad[isou][1]
                        + df[2]*grad[isou][2];

  }
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

inline static void
cs_solu_f_val_tensor(const cs_real_3_t   cell_cen,
                     const cs_real_3_t   i_face_cog,
                     const cs_real_63_t  grad,
                     const cs_real_6_t   p,
                     cs_real_t           pf[6])
{
  cs_real_3_t df;

  for (int jsou = 0; jsou < 3; jsou++)
    df[jsou] = i_face_cog[jsou] - cell_cen[jsou];

  for (int isou = 0; isou < 6; isou++) {
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
 * \param[in]     blencp   proportion of centered or SOLU scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     p        (relaxed) value at cell
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

inline static void
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
 * \param[in]     blencp   proportion of centered or SOLU scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     p        (relaxed) value at cell
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_blend_f_val_vector(const double       blencp,
                      const cs_real_3_t  p,
                      cs_real_t          pf[3])
{
  for (int isou = 0; isou < 3; isou++)
    pf[isou] = blencp*(pf[isou])+(1.-blencp)*p[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Blend face values for a centered or SOLU scheme with face values for
 * an upwind scheme.
 *
 * \param[in]     blencp   proportion of centered or SOLU scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     p        (relaxed) value at cell
 * \param[out]    pf       face value
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_blend_f_val_tensor(const double       blencp,
                      const cs_real_6_t  p,
                      cs_real_t          pf[6])
{
  for (int isou = 0; isou < 6; isou++)
    pf[isou] = blencp*(pf[isou])+(1.-blencp)*p[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (substracting the mass accumulation from them)
 * to fluxes at face ij.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-schema,
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

inline static void
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
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap        weighting coefficient for the theta-schema,
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

inline static void
cs_i_conv_flux_vector(const int         iconvp,
                      const cs_real_t   thetap,
                      const int         imasac,
                      const cs_real_3_t pi,
                      const cs_real_3_t pj,
                      const cs_real_3_t pifri,
                      const cs_real_3_t pifrj,
                      const cs_real_3_t pjfri,
                      const cs_real_3_t pjfrj,
                      const cs_real_t   i_massflux,
                      cs_real_t         fluxi[3],
                      cs_real_t         fluxj[3])
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux + fabs(i_massflux));
  fluj = 0.5*(i_massflux - fabs(i_massflux));

  for (int isou = 0; isou < 3; isou++) {

    fluxi[isou] +=  iconvp*(  thetap*(flui*pifri[isou] + fluj*pjfri[isou])
                            - imasac*i_massflux*pi[isou]);
    fluxj[isou] +=  iconvp*(  thetap*(flui*pifrj[isou] + fluj*pjfrj[isou])
                            - imasac*i_massflux*pj[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (substracting the mass accumulation from them)
 * to fluxes at face ij.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap        weighting coefficient for the theta-schema,
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

inline static void
cs_i_conv_flux_tensor(const int         iconvp,
                      const cs_real_t   thetap,
                      const int         imasac,
                      const cs_real_6_t pi,
                      const cs_real_6_t pj,
                      const cs_real_6_t pifri,
                      const cs_real_6_t pifrj,
                      const cs_real_6_t pjfri,
                      const cs_real_6_t pjfrj,
                      const cs_real_t   i_massflux,
                      cs_real_t         fluxi[6],
                      cs_real_t         fluxj[6])
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux + fabs(i_massflux));
  fluj = 0.5*(i_massflux - fabs(i_massflux));

  for (int isou = 0; isou < 6; isou++) {
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
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[in]     i_visc   diffusion coefficient (divided by IJ) at face ij
 * \param[in,out] fluxij   fluxes at face ij
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_diff_flux(const int       idiffp,
               const cs_real_t thetap,
               const cs_real_t pip,
               const cs_real_t pjp,
               const cs_real_t pipr,
               const cs_real_t pjpr,
               const cs_real_t i_visc,
               cs_real_2_t     fluxij)
{
  fluxij[0] += idiffp*thetap*i_visc*(pipr -pjp);
  fluxij[1] += idiffp*thetap*i_visc*(pip -pjpr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive fluxes to fluxes at face ij.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[in]     i_visc   diffusion coefficient (divided by IJ) at face ij
 * \param[in,out] fluxi    fluxes at face i
 * \param[in,out] fluxj    fluxes at face j
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_diff_flux_vector(const int         idiffp,
                      const cs_real_t   thetap,
                      const cs_real_3_t pip,
                      const cs_real_3_t pjp,
                      const cs_real_3_t pipr,
                      const cs_real_3_t pjpr,
                      const cs_real_t   i_visc,
                      cs_real_t         fluxi[3],
                      cs_real_t         fluxj[3])
{
  for (int isou = 0; isou < 3; isou++) {
    fluxi[isou] += idiffp*thetap*i_visc*(pipr[isou] -pjp[isou]);
    fluxj[isou] += idiffp*thetap*i_visc*(pip[isou] -pjpr[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive fluxes to fluxes at face ij.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[in]     i_visc   diffusion coefficient (divided by IJ) at face ij
 * \param[in,out] fluxi    fluxes at face i
 * \param[in,out] fluxj    fluxes at face j
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_diff_flux_tensor(const int         idiffp,
                      const cs_real_t   thetap,
                      const cs_real_6_t pip,
                      const cs_real_6_t pjp,
                      const cs_real_6_t pipr,
                      const cs_real_6_t pjpr,
                      const cs_real_t   i_visc,
                      cs_real_t         fluxi[6],
                      cs_real_t         fluxj[6])
{
  for (int isou = 0; isou < 6; isou++) {
    fluxi[isou] += idiffp*thetap*i_visc*(pipr[isou] -pjp[isou]);
    fluxj[isou] += idiffp*thetap*i_visc*(pip[isou] -pjpr[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and a pure upwind flux.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady_upwind(const int          ircflp,
                      const cs_real_t    relaxp,
                      const cs_real_t    weight,
                      const cs_real_3_t  cell_ceni,
                      const cs_real_3_t  cell_cenj,
                      const cs_real_3_t  i_face_cog,
                      const cs_real_3_t  dijpf,
                      const cs_real_3_t  gradi,
                      const cs_real_3_t  gradj,
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

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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
 * \param[in]     ircflp       recontruction flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady_upwind_vector(const int          ircflp,
                             const cs_real_t    relaxp,
                             const cs_real_t    weight,
                             const cs_real_3_t  cell_ceni,
                             const cs_real_3_t  cell_cenj,
                             const cs_real_3_t  i_face_cog,
                             const cs_real_3_t  dijpf,
                             const cs_real_33_t gradi,
                             const cs_real_33_t gradj,
                             const cs_real_3_t  pi,
                             const cs_real_3_t  pj,
                             const cs_real_3_t  pia,
                             const cs_real_3_t  pja,
                             cs_real_t          pifri[3],
                             cs_real_t          pifrj[3],
                             cs_real_t          pjfri[3],
                             cs_real_t          pjfrj[3],
                             cs_real_t          pip[3],
                             cs_real_t          pjp[3],
                             cs_real_t          pipr[3],
                             cs_real_t          pjpr[3])
{
  cs_real_3_t pir, pjr;
  cs_real_3_t recoi, recoj;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_vector(relaxp,
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

  cs_upwind_f_val_vector(pi,
                         pifrj);
  cs_upwind_f_val_vector(pir,
                         pifri);
  cs_upwind_f_val_vector(pj,
                         pjfri);
  cs_upwind_f_val_vector(pjr,
                         pjfrj);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and a pure upwind flux.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady_upwind_tensor(const int          ircflp,
                             const cs_real_t    relaxp,
                             const cs_real_t    weight,
                             const cs_real_3_t  cell_ceni,
                             const cs_real_3_t  cell_cenj,
                             const cs_real_3_t  i_face_cog,
                             const cs_real_3_t  dijpf,
                             const cs_real_63_t gradi,
                             const cs_real_63_t gradj,
                             const cs_real_6_t  pi,
                             const cs_real_6_t  pj,
                             const cs_real_6_t  pia,
                             const cs_real_6_t  pja,
                             cs_real_t          pifri[6],
                             cs_real_t          pifrj[6],
                             cs_real_t          pjfri[6],
                             cs_real_t          pjfrj[6],
                             cs_real_t          pip[6],
                             cs_real_t          pjp[6],
                             cs_real_t          pipr[6],
                             cs_real_t          pjpr[6])
{
  cs_real_6_t pir, pjr;
  cs_real_6_t recoi, recoj;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_tensor(relaxp,
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

  cs_upwind_f_val_tensor(pi,
                         pifrj);
  cs_upwind_f_val_tensor(pir,
                         pifri);
  cs_upwind_f_val_tensor(pj,
                         pjfri);
  cs_upwind_f_val_tensor(pjr,
                         pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and a pure upwind flux.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_unsteady_upwind(const int          ircflp,
                        const cs_real_t    weight,
                        const cs_real_3_t  cell_ceni,
                        const cs_real_3_t  cell_cenj,
                        const cs_real_3_t  i_face_cog,
                        const cs_real_3_t  dijpf,
                        const cs_real_3_t  gradi,
                        const cs_real_3_t  gradj,
                        const cs_real_t    pi,
                        const cs_real_t    pj,
                        cs_real_t         *pif,
                        cs_real_t         *pjf,
                        cs_real_t         *pip,
                        cs_real_t         *pjp)
{
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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
 * \param[in]     ircflp       recontruction flag
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_unsteady_upwind_vector(const int           ircflp,
                               const cs_real_t     weight,
                               const cs_real_3_t   cell_ceni,
                               const cs_real_3_t   cell_cenj,
                               const cs_real_3_t   i_face_cog,
                               const cs_real_3_t   dijpf,
                               const cs_real_33_t  gradi,
                               const cs_real_33_t  gradj,
                               const cs_real_3_t   pi,
                               const cs_real_3_t   pj,
                               cs_real_t           pif[3],
                               cs_real_t           pjf[3],
                               cs_real_t           pip[3],
                               cs_real_t           pjp[3])
{
  cs_real_3_t recoi, recoj;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_upwind_f_val_vector(pi, pif);
  cs_upwind_f_val_vector(pj, pjf);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and a pure upwind flux.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_unsteady_upwind_tensor(const int           ircflp,
                               const cs_real_t     weight,
                               const cs_real_3_t   cell_ceni,
                               const cs_real_3_t   cell_cenj,
                               const cs_real_3_t   i_face_cog,
                               const cs_real_3_t   dijpf,
                               const cs_real_63_t  gradi,
                               const cs_real_63_t  gradj,
                               const cs_real_6_t   pi,
                               const cs_real_6_t   pj,
                               cs_real_t           pif[6],
                               cs_real_t           pjf[6],
                               cs_real_t           pip[6],
                               cs_real_t           pjp[6])
{
  cs_real_6_t recoi, recoj;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_upwind_f_val_tensor(pi, pif);
  cs_upwind_f_val_tensor(pj, pjf);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and without enabling slope tests.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady(const int          ircflp,
               const int          ischcp,
               const double       relaxp,
               const double       blencp,
               const cs_real_t    weight,
               const cs_real_3_t  cell_ceni,
               const cs_real_3_t  cell_cenj,
               const cs_real_3_t  i_face_cog,
               const cs_real_3_t  dijpf,
               const cs_real_3_t  gradi,
               const cs_real_3_t  gradj,
               const cs_real_3_t  gradupi,
               const cs_real_3_t  gradupj,
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

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady_vector(const int           ircflp,
                      const int           ischcp,
                      const double        relaxp,
                      const double        blencp,
                      const cs_real_t     weight,
                      const cs_real_3_t   cell_ceni,
                      const cs_real_3_t   cell_cenj,
                      const cs_real_3_t   i_face_cog,
                      const cs_real_3_t   dijpf,
                      const cs_real_33_t  gradi,
                      const cs_real_33_t  gradj,
                      const cs_real_3_t   pi,
                      const cs_real_3_t   pj,
                      const cs_real_3_t   pia,
                      const cs_real_3_t   pja,
                      cs_real_t           pifri[3],
                      cs_real_t           pifrj[3],
                      cs_real_t           pjfri[3],
                      cs_real_t           pjfrj[3],
                      cs_real_t           pip[3],
                      cs_real_t           pjp[3],
                      cs_real_t           pipr[3],
                      cs_real_t           pjpr[3])
{
  cs_real_3_t pir, pjr;
  cs_real_3_t recoi, recoj;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_vector(relaxp,
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

    cs_centered_f_val_vector(weight,
                             pip,
                             pjpr,
                             pifrj);
    cs_centered_f_val_vector(weight,
                             pipr,
                             pjp,
                             pifri);
    cs_centered_f_val_vector(weight,
                             pipr,
                             pjp,
                             pjfri);
    cs_centered_f_val_vector(weight,
                             pip,
                             pjpr,
                             pjfrj);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val_vector(cell_ceni,
                         i_face_cog,
                         gradi,
                         pi,
                         pifrj);
    cs_solu_f_val_vector(cell_ceni,
                         i_face_cog,
                         gradi,
                         pir,
                         pifri);
    cs_solu_f_val_vector(cell_cenj,
                         i_face_cog,
                         gradj,
                         pj,
                         pjfri);
    cs_solu_f_val_vector(cell_cenj,
                         i_face_cog,
                         gradj,
                         pjr,
                         pjfrj);

  }

  /* Blending
     --------*/
  cs_blend_f_val_vector(blencp,
                        pi,
                        pifrj);
  cs_blend_f_val_vector(blencp,
                        pir,
                        pifri);
  cs_blend_f_val_vector(blencp,
                        pj,
                        pjfri);
  cs_blend_f_val_vector(blencp,
                        pjr,
                        pjfrj);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and without enabling slope tests.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     relaxp       relaxation coefficient
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_steady_tensor(const int           ircflp,
                      const int           ischcp,
                      const double        relaxp,
                      const double        blencp,
                      const cs_real_t     weight,
                      const cs_real_3_t   cell_ceni,
                      const cs_real_3_t   cell_cenj,
                      const cs_real_3_t   i_face_cog,
                      const cs_real_3_t   dijpf,
                      const cs_real_63_t  gradi,
                      const cs_real_63_t  gradj,
                      const cs_real_6_t   pi,
                      const cs_real_6_t   pj,
                      const cs_real_6_t   pia,
                      const cs_real_6_t   pja,
                      cs_real_t           pifri[6],
                      cs_real_t           pifrj[6],
                      cs_real_t           pjfri[6],
                      cs_real_t           pjfrj[6],
                      cs_real_t           pip[6],
                      cs_real_t           pjp[6],
                      cs_real_t           pipr[6],
                      cs_real_t           pjpr[6])

{
  cs_real_6_t pir, pjr;
  cs_real_6_t recoi, recoj;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_tensor(relaxp,
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

    cs_centered_f_val_tensor(weight,
                             pip,
                             pjpr,
                             pifrj);
    cs_centered_f_val_tensor(weight,
                             pipr,
                             pjp,
                             pifri);
    cs_centered_f_val_tensor(weight,
                             pipr,
                             pjp,
                             pjfri);
    cs_centered_f_val_tensor(weight,
                             pip,
                             pjpr,
                             pjfrj);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val_tensor(cell_ceni,
                         i_face_cog,
                         gradi,
                         pi,
                         pifrj);
    cs_solu_f_val_tensor(cell_ceni,
                         i_face_cog,
                         gradi,
                         pir,
                         pifri);
    cs_solu_f_val_tensor(cell_cenj,
                         i_face_cog,
                         gradj,
                         pj,
                         pjfri);
    cs_solu_f_val_tensor(cell_cenj,
                         i_face_cog,
                         gradj,
                         pjr,
                         pjfrj);

  }

  /* Blending
     --------*/

  cs_blend_f_val_tensor(blencp,
                        pi,
                        pifrj);
  cs_blend_f_val_tensor(blencp,
                        pir,
                        pifri);
  cs_blend_f_val_tensor(blencp,
                        pj,
                        pjfri);
  cs_blend_f_val_tensor(blencp,
                        pjr,
                        pjfrj);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and without enabling slope tests.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     gradupi      upwind gradient at cell i
 * \param[in]     gradupj      upwind gradient at cell j
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

inline static void
cs_i_cd_unsteady(const int          ircflp,
                 const int          ischcp,
                 const double       blencp,
                 const cs_real_t    weight,
                 const cs_real_3_t  cell_ceni,
                 const cs_real_3_t  cell_cenj,
                 const cs_real_3_t  i_face_cog,
                 const cs_real_3_t  dijpf,
                 const cs_real_3_t  gradi,
                 const cs_real_3_t  gradj,
                 const cs_real_3_t  gradupi,
                 const cs_real_3_t  gradupj,
                 const cs_real_t    pi,
                 const cs_real_t    pj,
                 cs_real_t         *pif,
                 cs_real_t         *pjf,
                 cs_real_t         *pip,
                 cs_real_t         *pjp)
{
  cs_real_t recoi, recoj;

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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

    /* Original SOLU
       ------------*/

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
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[out]    pif          contribution of i to flux from i to j
 * \param[out]    pjf          contribution of j to flux from j to i
 * \param[out]    pip          reconstructed value at cell i
 * \param[out]    pjp          reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_cd_unsteady_vector(const int           ircflp,
                        const int           ischcp,
                        const double        blencp,
                        const cs_real_t     weight,
                        const cs_real_3_t   cell_ceni,
                        const cs_real_3_t   cell_cenj,
                        const cs_real_3_t   i_face_cog,
                        const cs_real_3_t   dijpf,
                        const cs_real_33_t  gradi,
                        const cs_real_33_t  gradj,
                        const cs_real_3_t   pi,
                        const cs_real_3_t   pj,
                        cs_real_t           pif[3],
                        cs_real_t           pjf[3],
                        cs_real_t           pip[3],
                        cs_real_t           pjp[3])

{
  cs_real_3_t recoi, recoj;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
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

    cs_centered_f_val_vector(weight,
                             pip,
                             pjp,
                             pif);
    cs_centered_f_val_vector(weight,
                             pip,
                             pjp,
                             pjf);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val_vector(cell_ceni,
                         i_face_cog,
                         gradi,
                         pi,
                         pif);
    cs_solu_f_val_vector(cell_cenj,
                         i_face_cog,
                         gradj,
                         pj,
                         pjf);

  }

  /* Blending
     --------*/

  cs_blend_f_val_vector(blencp,
                        pi,
                        pif);
  cs_blend_f_val_vector(blencp,
                        pj,
                        pjf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of an unsteady algorithm and without enabling slope tests.
 *
 * \param[in]     ircflp       recontruction flag
 * \param[in]     ischcp       second order convection scheme flag
 * \param[in]     blencp       proportion of centered or SOLU scheme,
 *                             (1-blencp) is the proportion of upwind.
 * \param[in]     weight       geometrical weight
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     dijpf        distance I'J'
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

inline static void
cs_i_cd_unsteady_tensor(const int           ircflp,
                        const int           ischcp,
                        const double        blencp,
                        const cs_real_t     weight,
                        const cs_real_3_t   cell_ceni,
                        const cs_real_3_t   cell_cenj,
                        const cs_real_3_t   i_face_cog,
                        const cs_real_3_t   dijpf,
                        const cs_real_63_t  gradi,
                        const cs_real_63_t  gradj,
                        const cs_real_6_t   pi,
                        const cs_real_6_t   pj,
                        cs_real_t           pif[6],
                        cs_real_t           pjf[6],
                        cs_real_t           pip[6],
                        cs_real_t           pjp[6])

{
  cs_real_6_t recoi, recoj;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
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

    cs_centered_f_val_tensor(weight,
                             pip,
                             pjp,
                             pif);
    cs_centered_f_val_tensor(weight,
                             pip,
                             pjp,
                             pjf);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val_tensor(cell_ceni,
                         i_face_cog,
                         gradi,
                         pi,
                         pif);
    cs_solu_f_val_tensor(cell_cenj,
                         i_face_cog,
                         gradj,
                         pj,
                         pjf);

  }

  /* Blending
     --------*/

  cs_blend_f_val_tensor(blencp,
                        pi,
                        pif);
  cs_blend_f_val_tensor(blencp,
                        pj,
                        pjf);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell j
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_steady_slope_test(bool              *upwind_switch,
                          const int          iconvp,
                          const int          ircflp,
                          const int          ischcp,
                          const double       relaxp,
                          const double       blencp,
                          const double       blend_st,
                          const cs_real_t    weight,
                          const cs_real_t    i_dist,
                          const cs_real_t    i_face_surf,
                          const cs_real_3_t  cell_ceni,
                          const cs_real_3_t  cell_cenj,
                          const cs_real_3_t  i_face_normal,
                          const cs_real_3_t  i_face_cog,
                          const cs_real_3_t  dijpf,
                          const cs_real_t    i_massflux,
                          const cs_real_3_t  gradi,
                          const cs_real_3_t  gradj,
                          const cs_real_3_t  gradupi,
                          const cs_real_3_t  gradupj,
                          const cs_real_3_t  gradsti,
                          const cs_real_3_t  gradstj,
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
  cs_real_t distf, srfan, testij, tesqck;

  distf = i_dist;
  srfan = i_face_surf;

  *upwind_switch = false;

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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
                  distf,
                  srfan,
                  i_face_normal,
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


    /* Slope test: Pourcentage of upwind
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
  } else {
    cs_upwind_f_val(pi,
                    pifrj);
    cs_upwind_f_val(pir,
                    pifri);
    cs_upwind_f_val(pj,
                    pjfri);
    cs_upwind_f_val(pjr,
                    pjfrj);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_steady_slope_test_vector_old(bool               *upwind_switch,
                                     const int           iconvp,
                                     const int           ircflp,
                                     const int           ischcp,
                                     const double        relaxp,
                                     const double        blencp,
                                     const cs_real_t     weight,
                                     const cs_real_t     i_dist,
                                     const cs_real_t     i_face_surf,
                                     const cs_real_3_t   cell_ceni,
                                     const cs_real_3_t   cell_cenj,
                                     const cs_real_3_t   i_face_normal,
                                     const cs_real_3_t   i_face_cog,
                                     const cs_real_3_t   dijpf,
                                     const cs_real_t     i_massflux,
                                     const cs_real_33_t  gradi,
                                     const cs_real_33_t  gradj,
                                     const cs_real_33_t  grdpai,
                                     const cs_real_33_t  grdpaj,
                                     const cs_real_3_t   pi,
                                     const cs_real_3_t   pj,
                                     const cs_real_3_t   pia,
                                     const cs_real_3_t   pja,
                                     cs_real_t           pifri[3],
                                     cs_real_t           pifrj[3],
                                     cs_real_t           pjfri[3],
                                     cs_real_t           pjfrj[3],
                                     cs_real_t           pip[3],
                                     cs_real_t           pjp[3],
                                     cs_real_t           pipr[3],
                                     cs_real_t           pjpr[3])
{
  cs_real_3_t pir, pjr;
  cs_real_3_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_3_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_vector(relaxp,
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
    cs_slope_test_vector_old(pi,
                             pj,
                             distf,
                             srfan,
                             i_face_normal,
                             gradi,
                             gradj,
                             grdpai,
                             grdpaj,
                             i_massflux,
                             testij,
                             tesqck);

    for (isou = 0; isou < 3; isou++) {
      if (tesqck[isou]<=0. || testij[isou]<=0.) {

        /* Upwind
           --------*/

        cs_upwind_f_val(pi[isou],
                        &pifrj[isou]);
        cs_upwind_f_val(pir[isou],
                        &pifri[isou]);
        cs_upwind_f_val(pj[isou],
                        &pjfri[isou]);
        cs_upwind_f_val(pjr[isou],
                        &pjfrj[isou]);

        *upwind_switch = true;

      } else {

        if (ischcp==1) {

          /* Centered
             --------*/

          cs_centered_f_val(weight,
                            pip[isou],
                            pjpr[isou],
                            &pifrj[isou]);
          cs_centered_f_val(weight,
                            pipr[isou],
                            pjp[isou],
                            &pifri[isou]);
          cs_centered_f_val(weight,
                            pipr[isou],
                            pjp[isou],
                            &pjfri[isou]);
          cs_centered_f_val(weight,
                            pip[isou],
                            pjpr[isou],
                            &pjfrj[isou]);

        } else {

          /* Second order
             ------------*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog,
                        gradi[isou],
                        pi[isou],
                        &pifrj[isou]);
          cs_solu_f_val(cell_ceni,
                        i_face_cog,
                        gradi[isou],
                        pir[isou],
                        &pifri[isou]);
          cs_solu_f_val(cell_cenj,
                        i_face_cog,
                        gradj[isou],
                        pj[isou],
                        &pjfri[isou]);
          cs_solu_f_val(cell_cenj,
                        i_face_cog,
                        gradj[isou],
                        pjr[isou],
                        &pjfrj[isou]);

        }
      }
    }

    /* Blending
       --------*/
    cs_blend_f_val_vector(blencp,
                          pi,
                          pifrj);
    cs_blend_f_val_vector(blencp,
                          pir,
                          pifri);
    cs_blend_f_val_vector(blencp,
                          pj,
                          pjfri);
    cs_blend_f_val_vector(blencp,
                          pjr,
                          pjfrj);

  /* If iconv=0 p*fr* are useless */
  } else {
    for (isou = 0; isou < 3; isou++) {
        cs_upwind_f_val(pi[isou],
                        &pifrj[isou]);
        cs_upwind_f_val(pir[isou],
                        &pifri[isou]);
        cs_upwind_f_val(pj[isou],
                        &pjfri[isou]);
        cs_upwind_f_val(pjr[isou],
                        &pjfrj[isou]);
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
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_steady_slope_test_vector(bool               *upwind_switch,
                                 const int           iconvp,
                                 const int           ircflp,
                                 const int           ischcp,
                                 const double        relaxp,
                                 const double        blencp,
                                 const double        blend_st,
                                 const cs_real_t     weight,
                                 const cs_real_t     i_dist,
                                 const cs_real_t     i_face_surf,
                                 const cs_real_3_t   cell_ceni,
                                 const cs_real_3_t   cell_cenj,
                                 const cs_real_3_t   i_face_normal,
                                 const cs_real_3_t   i_face_cog,
                                 const cs_real_3_t   dijpf,
                                 const cs_real_t     i_massflux,
                                 const cs_real_33_t  gradi,
                                 const cs_real_33_t  gradj,
                                 const cs_real_33_t  grdpai,
                                 const cs_real_33_t  grdpaj,
                                 const cs_real_3_t   pi,
                                 const cs_real_3_t   pj,
                                 const cs_real_3_t   pia,
                                 const cs_real_3_t   pja,
                                 cs_real_t           pifri[3],
                                 cs_real_t           pifrj[3],
                                 cs_real_t           pjfri[3],
                                 cs_real_t           pjfrj[3],
                                 cs_real_t           pip[3],
                                 cs_real_t           pjp[3],
                                 cs_real_t           pipr[3],
                                 cs_real_t           pjpr[3])
{
  cs_real_3_t pir, pjr;
  cs_real_3_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_vector(relaxp,
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
    cs_slope_test_vector(pi,
                         pj,
                         distf,
                         srfan,
                         i_face_normal,
                         gradi,
                         gradj,
                         grdpai,
                         grdpaj,
                         i_massflux,
                         &testij,
                         &tesqck);

    for (isou = 0; isou < 3; isou++) {
      if (ischcp==1) {

        /* Centered
           --------*/

        cs_centered_f_val(weight,
                          pip[isou],
                          pjpr[isou],
                          &pifrj[isou]);
        cs_centered_f_val(weight,
                          pipr[isou],
                          pjp[isou],
                          &pifri[isou]);
        cs_centered_f_val(weight,
                          pipr[isou],
                          pjp[isou],
                          &pjfri[isou]);
        cs_centered_f_val(weight,
                          pip[isou],
                          pjpr[isou],
                          &pjfrj[isou]);

      } else {

        /* Second order
           ------------*/

        cs_solu_f_val(cell_ceni,
                      i_face_cog,
                      gradi[isou],
                      pi[isou],
                      &pifrj[isou]);
        cs_solu_f_val(cell_ceni,
                      i_face_cog,
                      gradi[isou],
                      pir[isou],
                      &pifri[isou]);
        cs_solu_f_val(cell_cenj,
                      i_face_cog,
                      gradj[isou],
                      pj[isou],
                      &pjfri[isou]);
        cs_solu_f_val(cell_cenj,
                      i_face_cog,
                      gradj[isou],
                      pjr[isou],
                      &pjfrj[isou]);

      }

    }

    /* Slope test: Pourcentage of upwind
       ----------------------------------*/

    if (tesqck <= 0. || testij <= 0.) {
      cs_blend_f_val_vector(blend_st,
                            pi,
                            pifrj);
      cs_blend_f_val_vector(blend_st,
                            pir,
                            pifri);
      cs_blend_f_val_vector(blend_st,
                            pj,
                            pjfri);
      cs_blend_f_val_vector(blend_st,
                            pjr,
                            pjfrj);

      *upwind_switch = true;
    }


    /* Blending
       --------*/
    cs_blend_f_val_vector(blencp,
                          pi,
                          pifrj);
    cs_blend_f_val_vector(blencp,
                          pir,
                          pifri);
    cs_blend_f_val_vector(blencp,
                          pj,
                          pjfri);
    cs_blend_f_val_vector(blencp,
                          pjr,
                          pjfrj);

  /* If iconv=0 p*fr* are useless */
  } else {
    for (isou = 0; isou < 3; isou++) {
        cs_upwind_f_val(pi[isou],
                        &pifrj[isou]);
        cs_upwind_f_val(pir[isou],
                        &pifri[isou]);
        cs_upwind_f_val(pj[isou],
                        &pjfri[isou]);
        cs_upwind_f_val(pjr[isou],
                        &pjfrj[isou]);
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
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     relaxp          relaxation coefficient
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_steady_slope_test_tensor(bool               *upwind_switch,
                                 const int           iconvp,
                                 const int           ircflp,
                                 const int           ischcp,
                                 const double        relaxp,
                                 const double        blencp,
                                 const double        blend_st,
                                 const cs_real_t     weight,
                                 const cs_real_t     i_dist,
                                 const cs_real_t     i_face_surf,
                                 const cs_real_3_t   cell_ceni,
                                 const cs_real_3_t   cell_cenj,
                                 const cs_real_3_t   i_face_normal,
                                 const cs_real_3_t   i_face_cog,
                                 const cs_real_3_t   dijpf,
                                 const cs_real_t     i_massflux,
                                 const cs_real_63_t  gradi,
                                 const cs_real_63_t  gradj,
                                 const cs_real_63_t  grdpai,
                                 const cs_real_63_t  grdpaj,
                                 const cs_real_6_t   pi,
                                 const cs_real_6_t   pj,
                                 const cs_real_6_t   pia,
                                 const cs_real_6_t   pja,
                                 cs_real_t           pifri[6],
                                 cs_real_t           pifrj[6],
                                 cs_real_t           pjfri[6],
                                 cs_real_t           pjfrj[6],
                                 cs_real_t           pip[6],
                                 cs_real_t           pjp[6],
                                 cs_real_t           pipr[6],
                                 cs_real_t           pjpr[6])
{
  cs_real_6_t pir, pjr;
  cs_real_6_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
                                 gradi,
                                 gradj,
                                 pi,
                                 pj,
                                 recoi,
                                 recoj,
                                 pip,
                                 pjp);

  cs_i_relax_c_val_tensor(relaxp,
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
    cs_slope_test_tensor(pi,
                         pj,
                         distf,
                         srfan,
                         i_face_normal,
                         gradi,
                         gradj,
                         grdpai,
                         grdpaj,
                         i_massflux,
                         &testij,
                         &tesqck);

    for (isou = 0; isou < 6; isou++) {
      if (ischcp==1) {

        /* Centered
           --------*/

        cs_centered_f_val(weight,
                          pip[isou],
                          pjpr[isou],
                          &pifrj[isou]);
        cs_centered_f_val(weight,
                          pipr[isou],
                          pjp[isou],
                          &pifri[isou]);
        cs_centered_f_val(weight,
                          pipr[isou],
                          pjp[isou],
                          &pjfri[isou]);
        cs_centered_f_val(weight,
                          pip[isou],
                          pjpr[isou],
                          &pjfrj[isou]);

      } else {

        /* Second order
           ------------*/

        cs_solu_f_val(cell_ceni,
                      i_face_cog,
                      gradi[isou],
                      pi[isou],
                      &pifrj[isou]);
        cs_solu_f_val(cell_ceni,
                      i_face_cog,
                      gradi[isou],
                      pir[isou],
                      &pifri[isou]);
        cs_solu_f_val(cell_cenj,
                      i_face_cog,
                      gradj[isou],
                      pj[isou],
                      &pjfri[isou]);
        cs_solu_f_val(cell_cenj,
                      i_face_cog,
                      gradj[isou],
                      pjr[isou],
                      &pjfrj[isou]);

      }

    }

    /* Slope test: Pourcentage of upwind
       ----------------------------------*/

    if (tesqck <= 0. || testij <= 0.) {

      cs_blend_f_val_tensor(blend_st,
                            pi,
                            pifrj);
      cs_blend_f_val_tensor(blend_st,
                            pir,
                            pifri);
      cs_blend_f_val_tensor(blend_st,
                            pj,
                            pjfri);
      cs_blend_f_val_tensor(blend_st,
                            pjr,
                            pjfrj);

      *upwind_switch = true;

    }


    /* Blending
       --------*/

    cs_blend_f_val_tensor(blencp,
                          pi,
                          pifrj);
    cs_blend_f_val_tensor(blencp,
                          pir,
                          pifri);
    cs_blend_f_val_tensor(blencp,
                          pj,
                          pjfri);
    cs_blend_f_val_tensor(blencp,
                          pjr,
                          pjfrj);

   /* If iconv=0 p*fr* are useless */
  } else {
    for (isou = 0; isou < 6; isou++) {
        cs_upwind_f_val(pi[isou],
                        &pifrj[isou]);
        cs_upwind_f_val(pir[isou],
                        &pifri[isou]);
        cs_upwind_f_val(pj[isou],
                        &pjfri[isou]);
        cs_upwind_f_val(pjr[isou],
                        &pjfrj[isou]);
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
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_unsteady_slope_test(bool              *upwind_switch,
                            const int          iconvp,
                            const int          ircflp,
                            const int          ischcp,
                            const double       blencp,
                            const double       blend_st,
                            const cs_real_t    weight,
                            const cs_real_t    i_dist,
                            const cs_real_t    i_face_surf,
                            const cs_real_3_t  cell_ceni,
                            const cs_real_3_t  cell_cenj,
                            const cs_real_3_t  i_face_normal,
                            const cs_real_3_t  i_face_cog,
                            const cs_real_3_t  dijpf,
                            const cs_real_t    i_massflux,
                            const cs_real_3_t  gradi,
                            const cs_real_3_t  gradj,
                            const cs_real_3_t  gradupi,
                            const cs_real_3_t  gradupj,
                            const cs_real_3_t  gradsti,
                            const cs_real_3_t  gradstj,
                            const cs_real_t    pi,
                            const cs_real_t    pj,
                            cs_real_t         *pif,
                            cs_real_t         *pjf,
                            cs_real_t         *pip,
                            cs_real_t         *pjp)
{
  CS_UNUSED(blend_st);

  cs_real_t recoi, recoj;
  cs_real_t distf, srfan, testij, tesqck;

  distf = i_dist;
  srfan = i_face_surf;

  *upwind_switch = false;

  cs_i_compute_quantities(ircflp,
                          weight,
                          cell_ceni,
                          cell_cenj,
                          i_face_cog,
                          dijpf,
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
                  distf,
                  srfan,
                  i_face_normal,
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

    /* Slope test: Pourcentage of upwind
       ----------------------------------*/

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
 * \brief DEPRECATED Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_unsteady_slope_test_vector_old(bool               *upwind_switch,
                                       const int           iconvp,
                                       const int           ircflp,
                                       const int           ischcp,
                                       const double        blencp,
                                       const cs_real_t     weight,
                                       const cs_real_t     i_dist,
                                       const cs_real_t     i_face_surf,
                                       const cs_real_3_t   cell_ceni,
                                       const cs_real_3_t   cell_cenj,
                                       const cs_real_3_t   i_face_normal,
                                       const cs_real_3_t   i_face_cog,
                                       const cs_real_3_t   dijpf,
                                       const cs_real_t     i_massflux,
                                       const cs_real_33_t  gradi,
                                       const cs_real_33_t  gradj,
                                       const cs_real_33_t  grdpai,
                                       const cs_real_33_t  grdpaj,
                                       const cs_real_3_t   pi,
                                       const cs_real_3_t   pj,
                                       cs_real_t           pif[3],
                                       cs_real_t           pjf[3],
                                       cs_real_t           pip[3],
                                       cs_real_t           pjp[3])
{
  cs_real_3_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_3_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
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
    cs_slope_test_vector_old(pi,
                             pj,
                             distf,
                             srfan,
                             i_face_normal,
                             gradi,
                             gradj,
                             grdpai,
                             grdpaj,
                             i_massflux,
                             testij,
                             tesqck);

    /* FIXME: slope test should be done for the vector and not component by
     * component. This is conserved for compatibility only. */
    for (isou = 0; isou < 3; isou++) {
      if (tesqck[isou]<=0. || testij[isou]<=0.) {

        /* Upwind
           --------*/

        cs_upwind_f_val(pi[isou],
                        &pif[isou]);
        cs_upwind_f_val(pj[isou],
                        &pjf[isou]);

        *upwind_switch = true;

      } else {

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

        } else {

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
    }

    /* Blending
       --------*/
    cs_blend_f_val_vector(blencp,
                          pi,
                          pif);
    cs_blend_f_val_vector(blencp,
                          pj,
                          pjf);

  /* If iconv=0 p*f are useless */
  } else {

    for (isou = 0; isou < 3; isou++) {
      cs_upwind_f_val(pi[isou],
                      &pif[isou]);
      cs_upwind_f_val(pj[isou],
                      &pjf[isou]);

    }
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_unsteady_slope_test_vector(bool               *upwind_switch,
                                   const int           iconvp,
                                   const int           ircflp,
                                   const int           ischcp,
                                   const double        blencp,
                                   const double        blend_st,
                                   const cs_real_t     weight,
                                   const cs_real_t     i_dist,
                                   const cs_real_t     i_face_surf,
                                   const cs_real_3_t   cell_ceni,
                                   const cs_real_3_t   cell_cenj,
                                   const cs_real_3_t   i_face_normal,
                                   const cs_real_3_t   i_face_cog,
                                   const cs_real_3_t   dijpf,
                                   const cs_real_t     i_massflux,
                                   const cs_real_33_t  gradi,
                                   const cs_real_33_t  gradj,
                                   const cs_real_33_t  grdpai,
                                   const cs_real_33_t  grdpaj,
                                   const cs_real_3_t   pi,
                                   const cs_real_3_t   pj,
                                   cs_real_t           pif[3],
                                   cs_real_t           pjf[3],
                                   cs_real_t           pip[3],
                                   cs_real_t           pjp[3])
{
  cs_real_3_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_vector(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
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
    cs_slope_test_vector(pi,
                         pj,
                         distf,
                         srfan,
                         i_face_normal,
                         gradi,
                         gradj,
                         grdpai,
                         grdpaj,
                         i_massflux,
                         &testij,
                         &tesqck);

    for (isou = 0; isou < 3; isou++) {
      if (ischcp == 1) {

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

      } else {

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

    /* Slope test: Pourcentage of upwind
       ----------------------------------*/

    if (tesqck <= 0. || testij <= 0.) {

      cs_blend_f_val_vector(blend_st,
                            pi,
                            pif);
      cs_blend_f_val_vector(blend_st,
                            pj,
                            pjf);

      *upwind_switch = true;

    }


    /* Blending
       --------*/
    cs_blend_f_val_vector(blencp,
                          pi,
                          pif);
    cs_blend_f_val_vector(blencp,
                          pj,
                          pjf);

  /* If iconv=0 p*f are useless */
  } else {

    for (isou = 0; isou < 3; isou++) {
      cs_upwind_f_val(pi[isou],
                      &pif[isou]);
      cs_upwind_f_val(pj[isou],
                      &pjf[isou]);

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
 * \param[in]     iconvp          convection flag
 * \param[in]     ircflp          recontruction flag
 * \param[in]     ischcp          second order convection scheme flag
 * \param[in]     blencp          proportion of centered or SOLU scheme,
 *                                (1-blencp) is the proportion of upwind.
 * \param[in]     blend_st        proportion of centered or SOLU scheme,
 *                                when the slope test is activated
 *                                (1-blend_st) is the proportion of upwind.
 * \param[in]     weight          geometrical weight
 * \param[in]     i_dist          distance IJ.Nij
 * \param[in]     i_face_surf     face surface
 * \param[in]     cell_ceni       center of gravity coordinates of cell i
 * \param[in]     cell_cenj       center of gravity coordinates of cell i
 * \param[in]     i_face_normal   face normal
 * \param[in]     i_face_cog      center of gravity coordinates of face ij
 * \param[in]     dijpf           distance I'J'
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

inline static void
cs_i_cd_unsteady_slope_test_tensor(bool               *upwind_switch,
                                   const int           iconvp,
                                   const int           ircflp,
                                   const int           ischcp,
                                   const double        blencp,
                                   const double        blend_st,
                                   const cs_real_t     weight,
                                   const cs_real_t     i_dist,
                                   const cs_real_t     i_face_surf,
                                   const cs_real_3_t   cell_ceni,
                                   const cs_real_3_t   cell_cenj,
                                   const cs_real_3_t   i_face_normal,
                                   const cs_real_3_t   i_face_cog,
                                   const cs_real_3_t   dijpf,
                                   const cs_real_t     i_massflux,
                                   const cs_real_63_t  gradi,
                                   const cs_real_63_t  gradj,
                                   const cs_real_63_t  grdpai,
                                   const cs_real_63_t  grdpaj,
                                   const cs_real_6_t   pi,
                                   const cs_real_6_t   pj,
                                   cs_real_t           pif[6],
                                   cs_real_t           pjf[6],
                                   cs_real_t           pip[6],
                                   cs_real_t           pjp[6])
{
  cs_real_6_t recoi, recoj;
  cs_real_t distf, srfan;
  cs_real_t testij, tesqck;
  int isou;

  distf = i_dist;
  srfan = i_face_surf;

  cs_i_compute_quantities_tensor(ircflp,
                                 weight,
                                 cell_ceni,
                                 cell_cenj,
                                 i_face_cog,
                                 dijpf,
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
    cs_slope_test_tensor(pi,
                         pj,
                         distf,
                         srfan,
                         i_face_normal,
                         gradi,
                         gradj,
                         grdpai,
                         grdpaj,
                         i_massflux,
                         &testij,
                         &tesqck);

    for (isou = 0; isou < 6; isou++) {

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

      } else {

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

    /* Slope test activated: poucentage of upwind */
    if (tesqck <= 0. || testij <= 0.) {

        /* Upwind
           --------*/

      cs_blend_f_val_tensor(blend_st,
                            pi,
                            pif);
      cs_blend_f_val_tensor(blend_st,
                            pj,
                            pjf);

      *upwind_switch = true;
    }


    /* Blending
       --------*/

    cs_blend_f_val_tensor(blencp,
                          pi,
                          pif);
    cs_blend_f_val_tensor(blencp,
                          pj,
                          pjf);

  /* If iconv=0 p*fr* are useless */
  } else {

    for (isou = 0; isou < 6; isou++) {
      cs_upwind_f_val(pi[isou],
                      &pif[isou]);
      cs_upwind_f_val(pj[isou],
                      &pjf[isou]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' at boundary cell i.
 *
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     ircflp   recontruction flag
 * \param[out]    recoi    reconstruction at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_compute_quantities(const cs_real_3_t  diipb,
                        const cs_real_3_t  gradi,
                        const int          ircflp,
                        cs_real_t         *recoi)
{
  *recoi = ircflp * (  gradi[0]*diipb[0]
                     + gradi[1]*diipb[1]
                     + gradi[2]*diipb[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' at boundary cell i.
 *
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     ircflp   recontruction flag
 * \param[out]    recoi    reconstruction at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_compute_quantities_vector(const cs_real_3_t   diipb,
                               const cs_real_33_t  gradi,
                               const int           ircflp,
                               cs_real_t           recoi[3])
{
  for (int isou = 0; isou < 3; isou++) {
    recoi[isou] = ircflp * (gradi[isou][0]*diipb[0]
        + gradi[isou][1]*diipb[1]
        + gradi[isou][2]*diipb[2]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct values in I' at boundary cell i.
 *
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     ircflp   recontruction flag
 * \param[out]    recoi    reconstruction at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_compute_quantities_tensor(const cs_real_3_t   diipb,
                               const cs_real_63_t  gradi,
                               const int           ircflp,
                               cs_real_t           recoi[6])
{
  for (int isou = 0; isou < 6; isou++) {
    recoi[isou] = ircflp * (gradi[isou][0]*diipb[0]
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

inline static void
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
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[in]     recoi    reconstruction at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_relax_c_val_vector(const double       relaxp,
                        const cs_real_3_t  pi,
                        const cs_real_3_t  pia,
                        const cs_real_3_t  recoi,
                        cs_real_t          pir[3],
                        cs_real_t          pipr[3])
{
  for (int isou = 0; isou < 3; isou++) {
    pir[isou]  = pi[isou]/relaxp - (1.-relaxp)/relaxp*pia[isou];
    pipr[isou] = pir[isou] + recoi[isou];
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

inline static void
cs_b_relax_c_val_tensor(const double       relaxp,
                        const cs_real_6_t  pi,
                        const cs_real_6_t  pia,
                        const cs_real_6_t  recoi,
                        cs_real_t          pir[6],
                        cs_real_t          pipr[6])
{
  for (int isou = 0; isou < 6; isou++) {
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
 * \param[in]     thetap       weighting coefficient for the theta-schema,
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

inline static void
cs_b_imposed_conv_flux(int         iconvp,
                       cs_real_t   thetap,
                       int         imasac,
                       int         inc,
                       cs_int_t    bc_type,
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
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-schema,
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
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_imposed_conv_flux_vector(int              iconvp,
                              cs_real_t        thetap,
                              int              imasac,
                              int              inc,
                              cs_int_t         bc_type,
                              int              icvfli,
                              const cs_real_t  pi[restrict 3],
                              const cs_real_t  pir[restrict 3],
                              const cs_real_t  pipr[restrict 3],
                              const cs_real_t  coefap[restrict 3],
                              const cs_real_t  coefbp[restrict 3][3],
                              const cs_real_t  coface[restrict 3],
                              const cs_real_t  cofbce[restrict 3][3],
                              cs_real_t        b_massflux,
                              cs_real_t        flux[restrict 3])
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
    for (int isou = 0; isou < 3; isou++) {
      pfac  = inc*coefap[isou];
      for (int jsou = 0; jsou < 3; jsou++) {
        pfac += coefbp[isou][jsou]*pipr[jsou];
      }
      flux[isou] += iconvp*( thetap*(flui*pir[isou] + fluj*pfac)
                           - imasac*b_massflux*pi[isou]);
    }

  /* Imposed convective flux */

  } else {

    for (int isou = 0; isou < 3; isou++) {
      pfac  = inc*coface[isou];
      for (int jsou = 0; jsou < 3; jsou++) {
        pfac += cofbce[isou][jsou]*pipr[jsou];
      }
      flux[isou] += iconvp*( thetap*pfac
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
 * \param[in]     thetap       weighting coefficient for the theta-schema,
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

inline static void
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
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-schema,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     inc          Not an increment flag
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefa        explicit boundary coefficient for convection
 *                             operator
 * \param[in]     coefb        implicit boundary coefficient for convection
 *                             operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_upwind_flux_vector(const int          iconvp,
                        const cs_real_t    thetap,
                        const int          imasac,
                        const int          inc,
                        const int          bc_type,
                        const cs_real_3_t  pi,
                        const cs_real_3_t  pir,
                        const cs_real_3_t  pipr,
                        const cs_real_3_t  coefa,
                        const cs_real_33_t coefb,
                        const cs_real_t    b_massflux,
                        cs_real_t          flux[3])
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
  for (int isou = 0; isou < 3; isou++) {
    pfac  = inc*coefa[isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      pfac += coefb[isou][jsou]*pipr[jsou];
    }
    flux[isou] += iconvp*( thetap*(flui*pir[isou] + fluj*pfac)
                         - imasac*b_massflux*pi[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux is a pure upwind flux.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     thetap       weighting coefficient for the theta-schema,
 * \param[in]     imasac       take mass accumulation into account?
 * \param[in]     inc          Not an increment flag
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefa        explicit boundary coefficient for convection
 *                             operator
 * \param[in]     coefb        implicit boundary coefficient for convection
 *                             operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_upwind_flux_tensor(const int          iconvp,
                        const cs_real_t    thetap,
                        const int          imasac,
                        const int          inc,
                        const int          bc_type,
                        const cs_real_6_t  pi,
                        const cs_real_6_t  pir,
                        const cs_real_6_t  pipr,
                        const cs_real_6_t  coefa,
                        const cs_real_66_t coefb,
                        const cs_real_t    b_massflux,
                        cs_real_t          flux[6])
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
  for (int isou = 0; isou < 6; isou++) {
    pfac  = inc*coefa[isou];
    for (int jsou = 0; jsou < 6; jsou++) {
      pfac += coefb[isou][jsou]*pipr[jsou];
    }
    flux[isou] += iconvp*( thetap*(flui*pir[isou] + fluj*pfac)
                         - imasac*b_massflux*pi[isou]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at boundary face.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     inc      Not an increment flag
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     cofafp   explicit boundary coefficient for diffusion operator
 * \param[in]     cofbfp   implicit boundary coefficient for diffusion operator
 * \param[in]     b_visc   boundary face surface
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
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
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     inc      Not an increment flag
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     cofaf    explicit boundary coefficient for diffusion operator
 * \param[in]     cofbf    implicit boundary coefficient for diffusion operator
 * \param[in]     b_visc   boundary face surface
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_diff_flux_vector(const int          idiffp,
                      const cs_real_t    thetap,
                      const int          inc,
                      const cs_real_3_t  pipr,
                      const cs_real_3_t  cofaf,
                      const cs_real_33_t cofbf,
                      const cs_real_t    b_visc,
                      cs_real_t          flux[3])
{
  cs_real_t pfacd ;
  for (int isou = 0; isou < 3; isou++) {
    pfacd  = inc*cofaf[isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      pfacd += cofbf[isou][jsou]*pipr[jsou];
    }
    flux[isou] += idiffp*thetap*b_visc*pfacd;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at boundary face.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     thetap   weighting coefficient for the theta-schema,
 * \param[in]     inc      Not an increment flag
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     cofaf    explicit boundary coefficient for diffusion operator
 * \param[in]     cofbf    implicit boundary coefficient for diffusion operator
 * \param[in]     b_visc   boundary face surface
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_diff_flux_tensor(const int          idiffp,
                      const cs_real_t    thetap,
                      const int          inc,
                      const cs_real_6_t  pipr,
                      const cs_real_6_t  cofaf,
                      const cs_real_66_t cofbf,
                      const cs_real_t    b_visc,
                      cs_real_t          flux[6])
{
  cs_real_t pfacd ;
  for (int isou = 0; isou < 6; isou++) {
    pfacd  = inc*cofaf[isou];
    for (int jsou = 0; jsou < 6; jsou++) {
      pfacd += cofbf[isou][jsou]*pipr[jsou];
    }
    flux[isou] += idiffp*thetap*b_visc*pfacd;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * \param[in]     ircflp   recontruction flag
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_steady(const int          ircflp,
               const double       relaxp,
               const cs_real_3_t  diipb,
               const cs_real_3_t  gradi,
               const cs_real_t    pi,
               const cs_real_t    pia,
               cs_real_t         *pir,
               cs_real_t         *pipr)
{
  cs_real_t recoi;

  cs_b_compute_quantities(diipb,
                          gradi,
                          ircflp,
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
 * \param[in]     ircflp   recontruction flag
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_steady_vector(const int          ircflp,
                      const double       relaxp,
                      const cs_real_3_t  diipb,
                      const cs_real_33_t gradi,
                      const cs_real_3_t  pi,
                      const cs_real_3_t  pia,
                      cs_real_t          pir[3],
                      cs_real_t          pipr[3])
{
  cs_real_3_t recoi;

  cs_b_compute_quantities_vector(diipb,
                                 gradi,
                                 ircflp,
                                 recoi);

  cs_b_relax_c_val_vector(relaxp,
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
 * \param[in]     ircflp   recontruction flag
 * \param[in]     relaxp   relaxation coefficient
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pia      old value at cell i
 * \param[out]    pir      relaxed value at cell i
 * \param[out]    pipr     relaxed reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_steady_tensor(const int          ircflp,
                      const double       relaxp,
                      const cs_real_3_t  diipb,
                      const cs_real_63_t gradi,
                      const cs_real_6_t  pi,
                      const cs_real_6_t  pia,
                      cs_real_t          pir[6],
                      cs_real_t          pipr[6])
{
  cs_real_6_t recoi;

  cs_b_compute_quantities_tensor(diipb,
                                 gradi,
                                 ircflp,
                                 recoi);

  cs_b_relax_c_val_tensor(relaxp,
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
 * \param[in]     ircflp   recontruction flag
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[out]    pip      reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_unsteady(const int          ircflp,
                 const cs_real_3_t  diipb,
                 const cs_real_3_t  gradi,
                 const cs_real_t    pi,
                 cs_real_t         *pip)
{
  cs_real_t recoi;

  cs_b_compute_quantities(diipb,
                          gradi,
                          ircflp,
                          &recoi);

  *pip = pi + recoi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * \param[in]     ircflp   recontruction flag
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[out]    pip      reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_unsteady_vector(const int          ircflp,
                        const cs_real_3_t  diipb,
                        const cs_real_33_t gradi,
                        const cs_real_3_t  pi,
                        cs_real_t          pip[3])
{
  cs_real_3_t recoi;

  cs_b_compute_quantities_vector(diipb,
                                 gradi,
                                 ircflp,
                                 recoi);

  for (int isou = 0; isou < 3; isou++)
    pip[isou] = pi[isou] + recoi[isou];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of boundary face values for the flux computation in
 * case of a steady algorithm.
 *
 * \param[in]     ircflp   recontruction flag
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[out]    pip      reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_unsteady_tensor(const int          ircflp,
                        const cs_real_3_t  diipb,
                        const cs_real_63_t gradi,
                        const cs_real_6_t  pi,
                        cs_real_t          pip[6])
{
  cs_real_6_t recoi;

  cs_b_compute_quantities_tensor(diipb,
                                 gradi,
                                 ircflp,
                                 recoi);

  for(int isou = 0; isou< 6; isou++)
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

inline static void
cs_b_diff_flux_coupling(const int        idiffp,
                        const cs_real_t  pi,
                        const cs_real_t  pj,
                        const cs_real_t  b_visc,
                        cs_real_t       *fluxi)
{
  *fluxi += idiffp*b_visc*(pi - pj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at an internal coupling face for a vector.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[in]     b_visc   equivalent exchange coefficient at an internal
 *                         coupling face
 * \param[in,out] fluxi    flux at internal coupling face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_diff_flux_coupling_vector(const int          idiffp,
                               const cs_real_3_t  pi,
                               const cs_real_3_t  pj,
                               const cs_real_t    b_visc,
                               cs_real_t          fluxi[3])
{
  for (int k = 0; k < 3; k++)
    fluxi[k] += idiffp*b_visc*(pi[k] - pj[k]);
}


/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmas, ITRMAS)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwgrp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_t                visel[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmav, ITRMAV)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwgrp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrp, ITRGRP)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_t                visel[],
 cs_real_t                diverg[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrv, ITRGRV)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                diverg[]
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     f_id         field id
 * \param[in]     inc          Not an increment flag
 * \param[in]     halo_type    halo type
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     coefap       boundary condition array for the variable
 *                             (explicit part)
 * \param[in]     coefbp       boundary condition array for the variable
 *                             (implicit part)
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

void
cs_slope_test_gradient(int                     f_id,
                       int                     inc,
                       cs_halo_type_t          halo_type,
                       const cs_real_3_t      *grad,
                       cs_real_3_t            *grdpa,
                       const cs_real_t        *pvar,
                       const cs_real_t        *coefap,
                       const cs_real_t        *coefbp,
                       const cs_real_t        *i_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient in order to cope with SOLU schemes
 *        observed in the litterature.
 *
 * \param[in]     f_id         field index
 * \param[in]     inc          Not an increment flag
 * \param[in]     halo_type    halo type
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     coefap       boundary condition array for the variable
 *                             (explicit part)
 * \param[in]     coefbp       boundary condition array for the variable
 *                             (implicit part)
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

void
cs_upwind_gradient(const int                     f_id,
                   const int                     inc,
                   const cs_halo_type_t          halo_type,
                   const cs_real_t               coefap[],
                   const cs_real_t               coefbp[],
                   const cs_real_t               i_massflux[],
                   const cs_real_t               b_massflux[],
                   const cs_real_t     *restrict pvar,
                   cs_real_3_t         *restrict grdpa);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     inc          Not an increment flag
 * \param[in]     halo_type    halo type
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     coefa        boundary condition array for the variable
 *                             (explicit part)
 * \param[in]     coefb        boundary condition array for the variable
 *                             (implicit part)
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

void
cs_slope_test_gradient_vector(const int              inc,
                              const cs_halo_type_t   halo_type,
                              const cs_real_33_t    *grad,
                              cs_real_33_t          *grdpa,
                              const cs_real_3_t     *pvar,
                              const cs_real_3_t     *coefa,
                              const cs_real_33_t    *coefb,
                              const cs_real_t       *i_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     inc          Not an increment flag
 * \param[in]     halo_type    halo type
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     coefa        boundary condition array for the variable
 *                             (explicit part)
 * \param[in]     coefb        boundary condition array for the variable
 *                             (implicit part)
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

void
cs_slope_test_gradient_tensor(const int               inc,
                              const cs_halo_type_t    halo_type,
                              const cs_real_63_t     *grad,
                              cs_real_63_t           *grdpa,
                              const cs_real_6_t      *pvar,
                              const cs_real_6_t      *coefa,
                              const cs_real_66_t     *coefb,
                              const cs_real_t        *i_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a coefficient for blending that ensures the positivity
 *  of the scalar.
 *
 * \param[in]     f_id         field id
 * \param[in]     inc          "not an increment" flag
 * \param[in]     rovsdt       rho * volume / dt
 */
/*----------------------------------------------------------------------------*/

void
cs_max_limiter_building(int              f_id,
                        int              inc,
                        const cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                   (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar(int                       idtvar,
                               int                       f_id,
                               const cs_var_cal_opt_t    var_cal_opt,
                               int                       icvflb,
                               int                       inc,
                               int                       iccocg,
                               int                       imasac,
                               cs_real_t       *restrict pvar,
                               const cs_real_t *restrict pvara,
                               const cs_int_t            icvfli[],
                               const cs_real_t           coefap[],
                               const cs_real_t           coefbp[],
                               const cs_real_t           cofafp[],
                               const cs_real_t           cofbfp[],
                               const cs_real_t           i_massflux[],
                               const cs_real_t           b_massflux[],
                               const cs_real_t           i_visc[],
                               const cs_real_t           b_visc[],
                               cs_real_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_vector(int                         idtvar,
                               int                         f_id,
                               const cs_var_cal_opt_t      var_cal_opt,
                               int                         icvflb,
                               int                         inc,
                               int                         ivisep,
                               int                         imasac,
                               cs_real_3_t       *restrict pvar,
                               const cs_real_3_t *restrict pvara,
                               const cs_int_t              icvfli[],
                               const cs_real_3_t           coefav[],
                               const cs_real_33_t          coefbv[],
                               const cs_real_3_t           cofafv[],
                               const cs_real_33_t          cofbfv[],
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t             secvif[],
                               cs_real_3_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     coefa         boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefb         boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofaf         boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbf         boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_tensor(int                         idtvar,
                               int                         f_id,
                               const cs_var_cal_opt_t      var_cal_opt,
                               int                         icvflb,
                               int                         inc,
                               int                         imasac,
                               cs_real_6_t       *restrict pvar,
                               const cs_real_6_t *restrict pvara,
                               const cs_real_6_t           coefa[],
                               const cs_real_66_t          coefb[],
                               const cs_real_6_t           cofaf[],
                               const cs_real_66_t          cofbf[],
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               cs_real_6_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 * equation of a scalar field \f$ \varia \f$ such as the temperature.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs + \sum_{\fij \in \Facei{\celli}}      \left(
 *        C_p\dot{m}_\ij \varia_\fij
 *      - \lambda_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * \f$ Rhs \f$ has already been initialized before calling bilsct!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_thermal(int                       idtvar,
                                int                       f_id,
                                const cs_var_cal_opt_t    var_cal_opt,
                                int                       inc,
                                int                       iccocg,
                                int                       imasac,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_real_t           coefap[],
                                const cs_real_t           coefbp[],
                                const cs_real_t           cofafp[],
                                const cs_real_t           cofbfp[],
                                const cs_real_t           i_massflux[],
                                const cs_real_t           b_massflux[],
                                const cs_real_t           i_visc[],
                                const cs_real_t           b_visc[],
                                const cs_real_t           xcpp[],
                                cs_real_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before
 *   calling cs_anisotropic_diffusion_scalar!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
                                (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_scalar(int                       idtvar,
                                int                       f_id,
                                const cs_var_cal_opt_t    var_cal_opt,
                                int                       inc,
                                int                       iccocg,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_real_t           coefap[],
                                const cs_real_t           coefbp[],
                                const cs_real_t           cofafp[],
                                const cs_real_t           cofbfp[],
                                const cs_real_t           i_visc[],
                                const cs_real_t           b_visc[],
                                cs_real_6_t     *restrict viscel,
                                const cs_real_2_t         weighf[],
                                const cs_real_t           weighb[],
                                cs_real_t       *restrict rhs);

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensorial
 * diffusivity for a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling diftnv!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_generalized_diffusion_vector(int                         idtvar,
                                int                         f_id,
                                const cs_var_cal_opt_t      var_cal_opt,
                                int                         inc,
                                int                         ivisep,
                                cs_real_3_t       *restrict pvar,
                                const cs_real_3_t *restrict pvara,
                                const cs_real_3_t           coefav[],
                                const cs_real_33_t          coefbv[],
                                const cs_real_3_t           cofafv[],
                                const cs_real_33_t          cofbfv[],
                                const cs_real_33_t          i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             secvif[],
                                cs_real_3_t       *restrict rhs);

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensorial
 * diffusivity for a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling diftnv!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_vector(int                         idtvar,
                                int                         f_id,
                                const cs_var_cal_opt_t      var_cal_opt,
                                int                         inc,
                                int                         ivisep,
                                cs_real_3_t       *restrict pvar,
                                const cs_real_3_t *restrict pvara,
                                const cs_real_3_t           coefav[],
                                const cs_real_33_t          coefbv[],
                                const cs_real_3_t           cofafv[],
                                const cs_real_33_t          cofbfv[],
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             secvif[],
                                cs_real_6_t       *restrict viscel,
                                const cs_real_2_t           weighf[],
                                const cs_real_t             weighb[],
                                cs_real_3_t       *restrict rhs);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before
 *   calling cs_anisotropic_diffusion_scalar!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefa         boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefb         boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofaf         boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbf         boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_tensor(int                         idtvar,
                                int                         f_id,
                                const cs_var_cal_opt_t      var_cal_opt,
                                int                         inc,
                                cs_real_6_t       *restrict pvar,
                                const cs_real_6_t *restrict pvara,
                                const cs_real_6_t           coefa[],
                                const cs_real_66_t          coefb[],
                                const cs_real_6_t           cofaf[],
                                const cs_real_66_t          cofbf[],
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                cs_real_6_t       *restrict viscel,
                                const cs_real_2_t           weighf[],
                                const cs_real_t             weighb[],
                                cs_real_6_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the face mass flux with the face pressure (or pressure
 * increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     visel         viscosity by cell
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_diffusion_potential(const int                 f_id,
                            const cs_mesh_t          *m,
                            cs_mesh_quantities_t     *fvq,
                            int                       init,
                            int                       inc,
                            int                       imrgra,
                            int                       iccocg,
                            int                       nswrgp,
                            int                       imligp,
                            int                       iphydp,
                            int                       iwgrp,
                            int                       iwarnp,
                            double                    epsrgp,
                            double                    climgp,
                            double                    extrap,
                            cs_real_3_t     *restrict frcxt,
                            cs_real_t       *restrict pvar,
                            const cs_real_t           coefap[],
                            const cs_real_t           coefbp[],
                            const cs_real_t           cofafp[],
                            const cs_real_t           cofbfp[],
                            const cs_real_t           i_visc[],
                            const cs_real_t           b_visc[],
                            cs_real_t       *restrict visel,
                            cs_real_t       *restrict i_massflux,
                            cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the pressure gradient term to the mass flux
 * in case of anisotropic diffusion of the pressure field \f$ P \f$.
 *
 * More precisely, the mass flux side \f$ \dot{m}_\fij \f$ is updated as
 * follows:
 * \f[
 * \dot{m}_\fij = \dot{m}_\fij -
 *              \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
                                    (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_potential(const int                 f_id,
                                        const cs_mesh_t          *m,
                                        cs_mesh_quantities_t     *fvq,
                                        int                       init,
                                        int                       inc,
                                        int                       imrgra,
                                        int                       iccocg,
                                        int                       nswrgp,
                                        int                       imligp,
                                        int                       ircflp,
                                        int                       iphydp,
                                        int                       iwgrp,
                                        int                       iwarnp,
                                        double                    epsrgp,
                                        double                    climgp,
                                        double                    extrap,
                                        cs_real_3_t     *restrict frcxt,
                                        cs_real_t       *restrict pvar,
                                        const cs_real_t           coefap[],
                                        const cs_real_t           coefbp[],
                                        const cs_real_t           cofafp[],
                                        const cs_real_t           cofbfp[],
                                        const cs_real_t           i_visc[],
                                        const cs_real_t           b_visc[],
                                        cs_real_6_t     *restrict viscel,
                                        const cs_real_2_t         weighf[],
                                        const cs_real_t           weighb[],
                                        cs_real_t       *restrict i_massflux,
                                        cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the cell mass flux divergence with the face pressure (or
 * pressure increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     visel         viscosity by cell
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_diffusion_potential(const int                 f_id,
                       const cs_mesh_t          *m,
                       cs_mesh_quantities_t     *fvq,
                       int                       init,
                       int                       inc,
                       int                       imrgra,
                       int                       iccocg,
                       int                       nswrgp,
                       int                       imligp,
                       int                       iphydp,
                       int                       iwarnp,
                       double                    epsrgp,
                       double                    climgp,
                       double                    extrap,
                       cs_real_3_t     *restrict frcxt,
                       cs_real_t       *restrict pvar,
                       const cs_real_t           coefap[],
                       const cs_real_t           coefbp[],
                       const cs_real_t           cofafp[],
                       const cs_real_t           cofbfp[],
                       const cs_real_t           i_visc[],
                       const cs_real_t           b_visc[],
                       cs_real_t                 visel[],
                       cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the divergence of the mass flux due to the
 * pressure gradient (routine analog to cs_anisotropic_diffusion_scalar).
 *
 * More precisely, the divergence of the mass flux side
 * \f$ \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij \f$ is updated as follows:
 * \f[
 * \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  = \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  - \sum_{\fij \in \Facei{\celli}}
 *    \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
                                     (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] diverg        divergence of the mass flux
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_potential(const int                 f_id,
                                   const cs_mesh_t          *m,
                                   cs_mesh_quantities_t     *fvq,
                                   int                       init,
                                   int                       inc,
                                   int                       imrgra,
                                   int                       iccocg,
                                   int                       nswrgp,
                                   int                       imligp,
                                   int                       ircflp,
                                   int                       iphydp,
                                   int                       iwarnp,
                                   double                    epsrgp,
                                   double                    climgp,
                                   double                    extrap,
                                   cs_real_3_t     *restrict frcxt,
                                   cs_real_t       *restrict pvar,
                                   const cs_real_t           coefap[],
                                   const cs_real_t           coefbp[],
                                   const cs_real_t           cofafp[],
                                   const cs_real_t           cofbfp[],
                                   const cs_real_t           i_visc[],
                                   const cs_real_t           b_visc[],
                                   cs_real_6_t     *restrict viscel,
                                   const cs_real_2_t         weighf[],
                                   const cs_real_t           weighb[],
                                   cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CONVECTION_DIFFUSION_H__ */
