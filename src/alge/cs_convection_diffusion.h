#ifndef __CS_CONVECTION_DIFFUSION_H__
#define __CS_CONVECTION_DIFFUSION_H__

/*============================================================================
 * Convection-diffusion operators.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

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
 * \param[out]    testij          value of slope test first criterium
 * \param[out]    tesqck          value of slope test second criterium
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
  *tesqck = pow(dcc, 2.) - pow(ddi-ddj, 2.);
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
  cs_real_t dpxf, dpyf, dpzf;

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
 * \brief Copy reconstructed or not cell values for consistency with relaxation
 * case.
 *
 * \param[in]     pi     value at cell i
 * \param[in]     pj     value at cell j
 * \param[in]     pip    reconstructed value at cell i
 * \param[in]     pjp    reconstructed value at cell j
 * \param[out]    pir    value at cell i
 * \param[out]    pjr    value at cell j
 * \param[out]    pipr   reconstructed value at cell i
 * \param[out]    pjpr   reconstructed value at cell j
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_no_relax_c_val(const cs_real_t   pi,
                    const cs_real_t   pj,
                    const cs_real_t   pip,
                    const cs_real_t   pjp,
                    cs_real_t        *pir,
                    cs_real_t        *pjr,
                    cs_real_t        *pipr,
                    cs_real_t        *pjpr)
{
  *pir = pi;
  *pjr = pj;

  *pipr = pip;
  *pjpr = pjp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using an upwind scheme.
 *
 * \param[in]     pi      value at cell i
 * \param[in]     pj      value at cell j
 * \param[in]     pir     relaxed value at cell i
 * \param[in]     pjr     relaxed value at cell j
 * \param[out]    pifri   contribution of i to flux from i to j
 * \param[out]    pifrj   contribution of i to flux from j to i
 * \param[out]    pjfri   contribution of j to flux from i to j
 * \param[out]    pjfrj   contribution of j to flux from j to i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_upwind_f_val(const cs_real_t  pi,
                const cs_real_t  pj,
                const cs_real_t  pir,
                const cs_real_t  pjr,
                cs_real_t       *pifri,
                cs_real_t       *pifrj,
                cs_real_t       *pjfri,
                cs_real_t       *pjfrj)
{
  *pifri = pir;
  *pifrj = pi;
  *pjfri = pj;
  *pjfrj = pjr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a centered scheme.
 *
 * \param[in]     pnd      weight
 * \param[in]     pip      reconstructed value at cell i
 * \param[in]     pjp      reconstructed value at cell j
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     pjpr     relaxed reconstructed value at cell j
 * \param[out]    pifri    contribution of i to flux from i to j
 * \param[out]    pifrj    contribution of i to flux from j to i
 * \param[out]    pjfri    contribution of j to flux from i to j
 * \param[out]    pjfrj    contribution of j to flux from j to i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_centered_f_val(const double     pnd,
                  const cs_real_t  pip,
                  const cs_real_t  pjp,
                  const cs_real_t  pipr,
                  const cs_real_t  pjpr,
                  cs_real_t       *pifri,
                  cs_real_t       *pifrj,
                  cs_real_t       *pjfri,
                  cs_real_t       *pjfrj)
{
  *pifri = pnd*pipr + (1.-pnd)*pjp;
  *pjfri = *pifri;
  *pifrj = pnd*pip  + (1.-pnd)*pjpr;
  *pjfrj = *pifrj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare value at face ij by using a Second Order Linear Upwind scheme.
 *
 * \param[in]     cell_ceni    center of gravity coordinates of cell i
 * \param[in]     cell_cenj    center of gravity coordinates of cell i
 * \param[in]     i_face_cog   center of gravity coordinates of face ij
 * \param[in]     gradi        gradient at cell i
 * \param[in]     gradj        gradient at cell j
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pjr          relaxed value at cell j
 * \param[out]    pifri        contribution of i to flux from i to j
 * \param[out]    pifrj        contribution of i to flux from j to i
 * \param[out]    pjfri        contribution of j to flux from i to j
 * \param[out]    pjfrj        contribution of j to flux from j to i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_solu_f_val(const cs_real_3_t  cell_ceni,
              const cs_real_3_t  cell_cenj,
              const cs_real_3_t  i_face_cog,
              const cs_real_3_t  gradi,
              const cs_real_3_t  gradj,
              const cs_real_t    pi,
              const cs_real_t    pj,
              const cs_real_t    pir,
              const cs_real_t    pjr,
              cs_real_t         *pifri,
              cs_real_t         *pifrj,
              cs_real_t         *pjfri,
              cs_real_t         *pjfrj)
{
  cs_real_t difx, dify, difz, djfx, djfy, djfz;

  difx = i_face_cog[0] - cell_ceni[0];
  dify = i_face_cog[1] - cell_ceni[1];
  difz = i_face_cog[2] - cell_ceni[2];
  djfx = i_face_cog[0] - cell_cenj[0];
  djfy = i_face_cog[1] - cell_cenj[1];
  djfz = i_face_cog[2] - cell_cenj[2];

  /* leave reconstruction of PIF and PJF even if IRCFLP=0
     otherwise, it is the same as using upwind */
  *pifri = pir + difx*gradi[0]+dify*gradi[1]+difz*gradi[2];
  *pifrj = pi  + difx*gradi[0]+dify*gradi[1]+difz*gradi[2];
  *pjfrj = pjr + djfx*gradj[0]+djfy*gradj[1]+djfz*gradj[2];
  *pjfri = pj  + djfx*gradj[0]+djfy*gradj[1]+djfz*gradj[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Blend face values for a centered or SOLU scheme with face values for
 * an upwind scheme.
 *
 * \param[in]     blencp   proportion of centered or SOLU scheme,
 *                         (1-blencp) is the proportion of upwind.
 * \param[in]     pi       value at cell i
 * \param[in]     pj       value at cell j
 * \param[in]     pir      relaxed value at cell i
 * \param[in]     pjr      relaxed value at cell j
 * \param[out]    pifri    contribution of i to flux from i to j
 * \param[out]    pifrj    contribution of i to flux from j to i
 * \param[out]    pjfri    contribution of j to flux from i to j
 * \param[out]    pjfrj    contribution of j to flux from j to i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_blend_f_val(const double     blencp,
               const cs_real_t  pi,
               const cs_real_t  pj,
               const cs_real_t  pir,
               const cs_real_t  pjr,
               cs_real_t       *pifri,
               cs_real_t       *pifrj,
               cs_real_t       *pjfri,
               cs_real_t       *pjfrj)
{
  *pifri = blencp*(*pifri)+(1.-blencp)*pir;
  *pifrj = blencp*(*pifrj)+(1.-blencp)*pi;
  *pjfri = blencp*(*pjfri)+(1.-blencp)*pj;
  *pjfrj = blencp*(*pjfrj)+(1.-blencp)*pjr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (substracting the mass accumulation from them)
 * to fluxes at face ij.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     pifri        contribution of i to flux from i to j
 * \param[in]     pifrj        contribution of i to flux from j to i
 * \param[in]     pjfri        contribution of j to flux from i to j
 * \param[in]     pjfrj        contribution of j to flux from j to i
 * \param[in]     i_massflux   mass flux at face ij
 * \param[in,out] fluxij       fluxes at face ij
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_i_conv_flux(const int       iconvp,
               const cs_real_t pi,
               const cs_real_t pj,
               const cs_real_t pifri,
               const cs_real_t pifrj,
               const cs_real_t pjfri,
               const cs_real_t pjfrj,
               const cs_real_t i_massflux,
               cs_real_2_t     fluxij)
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux + fabs(i_massflux));
  fluj = 0.5*(i_massflux - fabs(i_massflux));

  fluxij[0] +=  iconvp*(flui*pifri + fluj*pjfri - i_massflux*pi);
  fluxij[1] +=  iconvp*(flui*pifrj + fluj*pjfrj - i_massflux*pj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes (cons. formulation) to fluxes at face ij.
 *
 * \param[in]     iconvp       convection flag
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
cs_i_conv_flux_cons(const int       iconvp,
                    const cs_real_t pifri,
                    const cs_real_t pifrj,
                    const cs_real_t pjfri,
                    const cs_real_t pjfrj,
                    const cs_real_t i_massflux,
                    const cs_real_t  xcppi,
                    const cs_real_t  xcppj,
                    cs_real_2_t     fluxij)
{
  cs_real_t flui, fluj;

  flui = 0.5*(i_massflux +fabs(i_massflux));
  fluj = 0.5*(i_massflux -fabs(i_massflux));

  fluxij[0] +=  iconvp*xcppi*(flui*pifri + fluj*pjfri);
  fluxij[1] +=  iconvp*xcppj*(flui*pifrj + fluj*pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive fluxes to fluxes at face ij.
 *
 * \param[in]     idiffp   diffusion flag
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
             const cs_real_t pip,
             const cs_real_t pjp,
             const cs_real_t pipr,
             const cs_real_t pjpr,
             const cs_real_t i_visc,
             cs_real_2_t     fluxij)
{
  fluxij[0] += idiffp*i_visc*(pipr -pjp);
  fluxij[1] += idiffp*i_visc*(pip -pjpr);
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
                  pj,
                  pir,
                  pjr,
                  pifri,
                  pifrj,
                  pjfri,
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

  cs_i_no_relax_c_val(pi,
                      pj,
                      *pip,
                      *pjp,
                      &pir,
                      &pjr,
                      pipr,
                      pjpr);

  cs_upwind_f_val(pi,
                  pj,
                  pir,
                  pjr,
                  pifri,
                  pifrj,
                  pjfri,
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
                      *pjp,
                      *pipr,
                      *pjpr,
                      pifri,
                      pifrj,
                      pjfri,
                      pjfrj);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val(cell_ceni,
                  cell_cenj,
                  i_face_cog,
                  gradi,
                  gradj,
                  pi,
                  pj,
                  pir,
                  pjr,
                  pifri,
                  pifrj,
                  pjfri,
                  pjfrj);

  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pj,
                 pir,
                 pjr,
                 pifri,
                 pifrj,
                 pjfri,
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
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
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
                 const cs_real_t    pi,
                 const cs_real_t    pj,
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

  cs_i_no_relax_c_val(pi,
                      pj,
                      *pip,
                      *pjp,
                      &pir,
                      &pjr,
                      pipr,
                      pjpr);

  if (ischcp == 1) {

    /* Centered
       --------*/

    cs_centered_f_val(weight,
                      *pip,
                      *pjp,
                      *pipr,
                      *pjpr,
                      pifri,
                      pifrj,
                      pjfri,
                      pjfrj);

  } else {

    /* Second order
       ------------*/

    cs_solu_f_val(cell_ceni,
                  cell_cenj,
                  i_face_cog,
                  gradi,
                  gradj,
                  pi,
                  pj,
                  pir,
                  pjr,
                  pifri,
                  pifrj,
                  pjfri,
                  pjfrj);

  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pj,
                 pir,
                 pjr,
                 pifri,
                 pifrj,
                 pjfri,
                 pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a steady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
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
cs_i_cd_steady_slope_test(bool              *upwind_switch,
                          const int          ircflp,
                          const int          ischcp,
                          const double       relaxp,
                          const double       blencp,
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
                          const cs_real_3_t  grdpai,
                          const cs_real_3_t  grdpaj,
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

  cs_slope_test(pi,
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

  if (tesqck<=0. || testij<=0.) {

    /* Upwind
       --------*/

    cs_upwind_f_val(pi,
                    pj,
                    pir,
                    pjr,
                    pifri,
                    pifrj,
                    pjfri,
                    pjfrj);

    *upwind_switch = true;

  } else {

    if (ischcp==1) {

      /* Centered
         --------*/

      cs_centered_f_val(weight,
                        *pip,
                        *pjp,
                        *pipr,
                        *pjpr,
                        pifri,
                        pifrj,
                        pjfri,
                        pjfrj);

    } else {

      /* Second order
         ------------*/

      cs_solu_f_val(cell_ceni,
                    cell_cenj,
                    i_face_cog,
                    gradi,
                    gradj,
                    pi,
                    pj,
                    pir,
                    pjr,
                    pifri,
                    pifrj,
                    pjfri,
                    pjfrj);

    }
  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pj,
                 pir,
                 pjr,
                 pifri,
                 pifrj,
                 pjfri,
                 pjfrj);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle preparation of internal face values for the fluxes computation
 * in case of a unsteady algorithm and using slope tests.
 *
 * \param[out]    upwind_switch   slope test result
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
cs_i_cd_unsteady_slope_test(bool              *upwind_switch,
                            const int          ircflp,
                            const int          ischcp,
                            const double       blencp,
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
                            const cs_real_3_t  grdpai,
                            const cs_real_3_t  grdpaj,
                            const cs_real_t    pi,
                            const cs_real_t    pj,
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

  cs_i_no_relax_c_val(pi,
                      pj,
                      *pip,
                      *pjp,
                      &pir,
                      &pjr,
                      pipr,
                      pjpr);

  cs_slope_test(pi,
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

  if (tesqck<=0. || testij<=0.) {

    /* Upwind
       --------*/

    cs_upwind_f_val(pi,
                    pj,
                    pir,
                    pjr,
                    pifri,
                    pifrj,
                    pjfri,
                    pjfrj);

    *upwind_switch = true;

  } else {

    if (ischcp==1) {

      /* Centered
         --------*/

      cs_centered_f_val(weight,
                        *pip,
                        *pjp,
                        *pipr,
                        *pjpr,
                        pifri,
                        pifrj,
                        pjfri,
                        pjfrj);

    } else {

      /* Second order
         ------------*/

      cs_solu_f_val(cell_ceni,
                    cell_cenj,
                    i_face_cog,
                    gradi,
                    gradj,
                    pi,
                    pj,
                    pir,
                    pjr,
                    pifri,
                    pifrj,
                    pjfri,
                    pjfrj);

    }
  }

  /* Blending
     --------*/

  cs_blend_f_val(blencp,
                 pi,
                 pj,
                 pir,
                 pjr,
                 pifri,
                 pifrj,
                 pjfri,
                 pjfrj);
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
 * \brief Copy values at bounadry cell i for consistency with relaxation case.
 *
 * \param[in]     pi       value at cell i
 * \param[in]     recoi    reconstruction at cell i
 * \param[out]    pir      value at cell i
 * \param[out]    pipr     reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_no_relax_c_val(const cs_real_t  pi,
                    const cs_real_t  recoi,
                    cs_real_t       *pir,
                    cs_real_t       *pipr)
{
  *pir  = pi;
  *pipr = pi + recoi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux can be either an upwind flux or an
 * imposed value.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     inc          Not an increment flag
 * \param[in]     ifaccp       indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
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
cs_b_imposed_conv_flux(const int        iconvp,
                       const int        inc,
                       const int        ifaccp,
                       const cs_int_t   bc_type,
                       const int        icvfli,
                       const cs_real_t  pi,
                       const cs_real_t  pir,
                       const cs_real_t  pipr,
                       const cs_real_t  coefap,
                       const cs_real_t  coefbp,
                       const cs_real_t  coface,
                       const cs_real_t  cofbce,
                       const cs_real_t  b_massflux,
                       cs_real_t       *flux)
{
  cs_real_t flui, fluj, pfac;

  /* Computed convective flux */

  if (icvfli == 0) {

    /* Remove decentering for coupled faces */
    if (ifaccp==1 && bc_type==CS_COUPLED) {
      flui = 0.0;
      fluj = b_massflux;
    } else {
      flui = 0.5*(b_massflux +fabs(b_massflux));
      fluj = 0.5*(b_massflux -fabs(b_massflux));
    }

    pfac  = inc*coefap + coefbp*pipr;
    *flux += iconvp*(flui*pir + fluj*pfac - b_massflux*pi);

  /* Imposed convective flux */

  } else {

    pfac = inc*coface + cofbce*pipr;
    *flux += iconvp*(-b_massflux*pi + pfac);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (conservative formulation) to flux at boundary
 * face. The convective flux can be either an upwind flux or an imposed value.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     inc          Not an increment flag
 * \param[in]     ifaccp       indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     bc_type      type of boundary face
 * \param[in]     icvfli       imposed convective flux flag
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection operator
 * \param[in]     coefbp       implicit boundary coefficient for convection operator
 * \param[in]     coface       explicit imposed convective flux value (0 otherwise).
 * \param[in]     cofbce       implicit part of imp. conv. flux value
 * \param[in]     xcpp         specific heat value if the scalar is the temperature,
 *                             1 otherwise
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_imposed_conv_flux_cons(const int        iconvp,
                            const int        inc,
                            const int        ifaccp,
                            const cs_int_t   bc_type,
                            const int        icvfli,
                            const cs_real_t  pir,
                            const cs_real_t  pipr,
                            const cs_real_t  coefap,
                            const cs_real_t  coefbp,
                            const cs_real_t  coface,
                            const cs_real_t  cofbce,
                            const cs_real_t  xcpp,
                            const cs_real_t  b_massflux,
                            cs_real_t       *flux)
{
  cs_real_t flui, fluj, pfac;

  /* Computed convective flux */

  if (icvfli == 0) {

    /* Remove decentering for coupled faces */
    if (ifaccp==1 && bc_type==CS_COUPLED) {
      flui = 0.0;
      fluj = b_massflux;
    } else {
      flui = 0.5*(b_massflux +fabs(b_massflux));
      fluj = 0.5*(b_massflux -fabs(b_massflux));
    }

    pfac  = inc*coefap + coefbp*pipr;
    *flux += iconvp*xcpp*(flui*pir + fluj*pfac);

  /* Imposed convective flux */

  } else {

    pfac = inc*coface + cofbce*pipr;
    *flux += iconvp*pfac;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (substracting the mass accumulation from it) to
 * flux at boundary face. The convective flux is a pure upwind flux.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     inc          Not an increment flag
 * \param[in]     ifaccp       indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection operator
 * \param[in]     coefbp       implicit boundary coefficient for convection operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_upwind_flux(const int        iconvp,
                 const int        inc,
                 const int        ifaccp,
                 const int        bc_type,
                 const cs_real_t  pi,
                 const cs_real_t  pir,
                 const cs_real_t  pipr,
                 const cs_real_t  coefap,
                 const cs_real_t  coefbp,
                 const cs_real_t  b_massflux,
                 cs_real_t       *flux)
{
  cs_real_t flui, fluj, pfac;

  /* Remove decentering for coupled faces */
  if (ifaccp==1 && bc_type==CS_COUPLED) {
    flui = 0.0;
    fluj = b_massflux;
  } else {
    flui = 0.5*(b_massflux +fabs(b_massflux));
    fluj = 0.5*(b_massflux -fabs(b_massflux));
  }

  pfac  = inc*coefap + coefbp*pipr;
  *flux += iconvp*(flui*pir + fluj*pfac - b_massflux*pi);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux (conservative formulation) to flux at boundary
 * face. The convective flux is a pure upwind flux.
 *
 * \param[in]     iconvp       convection flag
 * \param[in]     inc          Not an increment flag
 * \param[in]     ifaccp       indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pir          relaxed value at cell i
 * \param[in]     pipr         relaxed reconstructed value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection
 *                             operator
 * \param[in]     coefbp       implicit boundary coefficient for convection
 *                             operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in]     xcpp         specific heat value if the scalar is the
 *                             temperature, 1 otherwise
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_upwind_flux_cons(const int        iconvp,
                      const int        inc,
                      const int        ifaccp,
                      const cs_int_t   bc_type,
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
  if (ifaccp==1 && bc_type==CS_COUPLED) {
    flui = 0.0;
    fluj = b_massflux;
  } else {
    flui = 0.5*(b_massflux +fabs(b_massflux));
    fluj = 0.5*(b_massflux -fabs(b_massflux));
  }

  pfac  = inc*coefap + coefbp*pipr;
  *flux += iconvp*xcpp*(flui*pir + fluj*pfac);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add diffusive flux to flux at boundary face.
 *
 * \param[in]     idiffp   diffusion flag
 * \param[in]     inc      Not an increment flag
 * \param[in]     pipr     relaxed reconstructed value at cell i
 * \param[in]     cofafp   explicit boundary coefficient for diffusion operator
 * \param[in]     cofbfp   implicit boundary coefficient for diffusion operator
 * \param[in]     b_visc   mass flux at boundary face
 * \param[in,out] flux     flux at boundary face
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_diff_flux(const int        idiffp,
               const int        inc,
               const cs_real_t  pipr,
               const cs_real_t  cofafp,
               const cs_real_t  cofbfp,
               const cs_real_t  b_visc,
               cs_real_t       *flux)
{
  cs_real_t pfacd = inc*cofafp + cofbfp*pipr;
  *flux += idiffp*b_visc*pfacd;
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
 * case of an unsteady algorithm.
 *
 * \param[in]     ircflp   recontruction flag
 * \param[in]     diipb    distance I'I'
 * \param[in]     gradi    gradient at cell i
 * \param[in]     pi       value at cell i
 * \param[in]     pir      value at cell i (for consistancy)
 * \param[out]    pipr     reconstructed value at cell i
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_b_cd_unsteady(const int          ircflp,
                 const cs_real_3_t  diipb,
                 const cs_real_3_t  gradi,
                 const cs_real_t    pi,
                 cs_real_t         *pir,
                 cs_real_t         *pipr)
{
  cs_real_t recoi;

  cs_b_compute_quantities(diipb,
                          gradi,
                          ircflp,
                          &recoi);

  cs_b_no_relax_c_val(pi,
                      recoi,
                      pir,
                      pipr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * \param[in]     f_id         field index
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

inline static void
cs_slope_test_gradient(const int                     f_id,
                       const int                     inc,
                       const cs_halo_type_t          halo_type,
                       cs_real_3_t                  *grad,
                       cs_real_3_t                  *grdpa,
                       cs_real_t                    *pvar,
                       const cs_real_t              *coefap,
                       const cs_real_t              *coefbp,
                       const cs_real_t              *i_massflux)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  cs_lnum_t cell_id, face_id, ii, jj;
  int g_id, t_id;

  cs_real_t difx,dify,difz,djfx,djfy,djfz;
  cs_real_t diipbx, diipby, diipbz;
  cs_real_t pfac, pfac1, pfac2, pfac3, unsvol;
  cs_real_t pif,pjf;

  for (g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for private(face_id, ii, jj, difx, dify, difz, \
                                    djfx, djfy, djfz, pif, pjf, pfac,  \
                                    pfac1, pfac2, pfac3)
    for (t_id = 0; t_id < n_i_threads; t_id++) {
      for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        ii = i_face_cells[face_id][0];
        jj = i_face_cells[face_id][1];

        difx = i_face_cog[face_id][0] - cell_cen[ii][0];
        dify = i_face_cog[face_id][1] - cell_cen[ii][1];
        difz = i_face_cog[face_id][2] - cell_cen[ii][2];
        djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
        djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
        djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

        pif =  pvar[ii]
             + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
        pjf =  pvar[jj]
             + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

        pfac = pjf;
        if (i_massflux[face_id] > 0.) pfac = pif;

        pfac1 = pfac*i_face_normal[face_id][0];
        pfac2 = pfac*i_face_normal[face_id][1];
        pfac3 = pfac*i_face_normal[face_id][2];

        grdpa[ii][0] = grdpa[ii][0] + pfac1;
        grdpa[ii][1] = grdpa[ii][1] + pfac2;
        grdpa[ii][2] = grdpa[ii][2] + pfac3;

        grdpa[jj][0] = grdpa[jj][0] - pfac1;
        grdpa[jj][1] = grdpa[jj][1] - pfac2;
        grdpa[jj][2] = grdpa[jj][2] - pfac3;

      }
    }
  }

  for (g_id = 0; g_id < n_b_groups; g_id++) {
# pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz, pfac) \
                          if(m->n_b_faces > CS_THR_MIN)
    for (t_id = 0; t_id < n_b_threads; t_id++) {
      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        ii = b_face_cells[face_id];

        diipbx = diipb[face_id][0];
        diipby = diipb[face_id][1];
        diipbz = diipb[face_id][2];
        pfac =   inc*coefap[face_id]
               + coefbp[face_id] * (pvar[ii] + diipbx*grad[ii][0]
                                             + diipby*grad[ii][1]
                                             + diipbz*grad[ii][2]);
        grdpa[ii][0] = grdpa[ii][0] + pfac*b_face_normal[face_id][0];
        grdpa[ii][1] = grdpa[ii][1] + pfac*b_face_normal[face_id][1];
        grdpa[ii][2] = grdpa[ii][2] + pfac*b_face_normal[face_id][2];

      }
    }
  }

# pragma omp parallel for private(unsvol)
  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    unsvol = 1./cell_vol[cell_id];

    grdpa[cell_id][0] = grdpa[cell_id][0]*unsvol;
    grdpa[cell_id][1] = grdpa[cell_id][1]*unsvol;
    grdpa[cell_id][2] = grdpa[cell_id][2]*unsvol;

  }

  /* Synchronization for parallelism or periodicity */

  if (halo != NULL) {
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)grdpa, 3);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)grdpa, 3);

    /* Gradient periodicity of rotation for Reynolds stress components */
    if (cs_glob_mesh->have_rotation_perio > 0 && f_id != -1)
      cs_gradient_perio_process_rij(&f_id, grdpa);
  }

}

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc2, BILSC2)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   icvflb,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 const cs_int_t          *const   ifaccp,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_int_t                   bc_type[],
 const cs_int_t                   icvfli[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 cs_real_t                        rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc4, BILSC4)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   icvflb,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   ifaccp,
 const cs_int_t          *const   ivisep,
 cs_real_3_t                      pvar[],
 const cs_real_3_t                pvara[],
 const cs_int_t                   bc_type[],
 const cs_int_t                   icvfli[],
 const cs_real_3_t                coefav[],
 const cs_real_33_t               coefbv[],
 const cs_real_3_t                cofafv[],
 const cs_real_33_t               cofbfv[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  secvif[],
 cs_real_3_t                      rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_thermal
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsct, BILSCT)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 const cs_int_t          *const   ifaccp,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_int_t                   bc_type[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  xcpp[],
 cs_real_t                        rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (diften, DIFTEN)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 cs_real_6_t                      viscel[],
 const cs_real_2_t                weighf[],
 const cs_real_t                  weighb[],
 cs_real_t                        rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (diftnv, DIFTNV)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   ifaccp,
 const cs_int_t          *const   ivisep,
 cs_real_3_t                      pvar[],
 const cs_real_3_t                pvara[],
 const cs_int_t                   bc_type[],
 const cs_real_3_t                coefav[],
 const cs_real_33_t               coefbv[],
 const cs_real_3_t                cofafv[],
 const cs_real_33_t               cofbfv[],
 const cs_real_33_t               i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  secvif[],
 cs_real_3_t                      rhs[]
);

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
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmav, ITRMAV)
(
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
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                diverg[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrv, ITRGRV)
(
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
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
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
                               int                       ifaccp,
                               cs_real_t       *restrict pvar,
                               const cs_real_t *restrict pvara,
                               const cs_int_t            bc_type[],
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
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_type       boundary condition type
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
                               int                         ifaccp,
                               int                         ivisep,
                               cs_real_3_t       *restrict pvar,
                               const cs_real_3_t *restrict pvara,
                               const cs_int_t              bc_type[],
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
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
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
                                int                       ifaccp,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_int_t            bc_type[],
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
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
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
cs_anisotropic_diffusion_vector(int                         idtvar,
                                int                         f_id,
                                const cs_var_cal_opt_t      var_cal_opt,
                                int                         inc,
                                int                         ifaccp,
                                int                         ivisep,
                                cs_real_3_t       *restrict pvar,
                                const cs_real_3_t *restrict pvara,
                                const cs_int_t              bc_type[],
                                const cs_real_3_t           coefav[],
                                const cs_real_33_t          coefbv[],
                                const cs_real_3_t           cofafv[],
                                const cs_real_33_t          cofbfv[],
                                const cs_real_33_t          i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             secvif[],
                                cs_real_3_t       *restrict rhs);

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
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
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
                            const cs_real_t           viselx[],
                            const cs_real_t           visely[],
                            const cs_real_t           viselz[],
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
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_potential(const cs_mesh_t          *m,
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
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
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
                       const cs_real_t           viselx[],
                       const cs_real_t           visely[],
                       const cs_real_t           viselz[],
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
cs_anisotropic_diffusion_potential(const cs_mesh_t          *m,
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
