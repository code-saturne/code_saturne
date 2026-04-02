/*============================================================================
 * Convection-diffusion operators for stead algorithm (deprecated).
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

/*-----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <cmath>
#include <chrono>

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_halo.h"
#include "base/cs_log.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_profiling.h"
#include "base/cs_timer.h"
#include "base/cs_velocity_pressure.h"

#include "alge/cs_blas.h"
#include "alge/cs_gradient.h"
#include "alge/cs_gradient_boundary.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_convection_diffusion.h"
#include "alge/cs_convection_diffusion_priv.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_convection_diffusion_steady.cpp
 *
 * \brief Convection-diffusion operators.
 *
 * Please refer to the
 * <a href="../../theory.pdf#conv-diff"><b>convection-diffusion</b></a> section
 *  of the theory guide for more information.
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
                   pifri,
                   pjfrj,
                   pipr,
                   pjpr);

  *pifrj = pi;
  *pjfri = pj;
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
                                   pifri,
                                   pjfrj,
                                   pipr,
                                   pjpr);

  for (int l = 0; l < stride; l++) {
    pifrj[l] = pi[l];
    pjfri[l] = pj[l];
  }
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
    *pifrj = pi;
    *pifri = pir;
    *pjfri = pj;
    *pjfrj = pjr;
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

    if (ischcp==1) {

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

      cs_solu_f_val_strided<stride>(cell_ceni, i_face_cog, gradi, pi, pifrj);
      cs_solu_f_val_strided<stride>(cell_ceni, i_face_cog, gradi, pir, pifri);
      cs_solu_f_val_strided<stride>(cell_cenj, i_face_cog, gradj, pj, pjfri);
      cs_solu_f_val_strided<stride>(cell_cenj, i_face_cog, gradj, pjr, pjfrj);

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
    for (int l = 0; l < stride; l++) {
      pifrj[l] = pi[l];
      pifri[l] = pir[l];
      pjfri[l] = pj[l];
      pjfrj[l] = pjr[l];
    }
  }
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

/*----------------------------------------------------------------------------
 * Return the equivalent heat transfer coefficient. If both terms are
 * below a given tolerance, 0. is returned.
 *
 * parameters:
 *   h1     <-- first exchange coefficient
 *   h2     <-- second exchange coefficient
 *
 * return:
 *   value of equivalent exchange coefficient
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_calc_heq(cs_real_t h1,
          cs_real_t h2)
{
  const cs_real_t h_eps = 1.e-12;

  cs_real_t heq = 0.;
  if (h1 + h2 > h_eps)
    heq = h1 * h2 / (h1 + h2);

  return heq;
}

/*----------------------------------------------------------------------------
 * Synchronize strided gradient ghost cell values.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   on_device,     <-- is data on device (GPU) ?
 *   halo_type      <-- halo type (extended or not)
 *   grad           --> gradient of a variable
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sync_strided_gradient_halo(const cs_mesh_t         *m,
                            cs_halo_type_t           halo_type,
                            [[maybe_unused]] bool    on_device,
                            cs_real_t (*restrict grad)[stride][3])
{
#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_d(m->halo, halo_type, CS_REAL_TYPE, stride*3,
                   (cs_real_t *)grad);
  else
#endif
    cs_halo_sync(m->halo, halo_type, CS_REAL_TYPE, stride*3,
                 (cs_real_t *)grad);

  if (m->have_rotation_perio) {
#if defined(HAVE_ACCEL)
    if (on_device)
      cs_sync_d2h((void  *)grad);
#endif
    if (stride == 1)
      cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)grad, 3);
    else if (stride == 3)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
    else if (stride == 6)
      cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                           halo_type,
                                           (cs_real_t *)grad);
#if defined(HAVE_ACCEL)
    if (on_device)
      cs_sync_h2d((void  *)grad);
#endif
  }
}

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
 * Please refer to the
 * <a href="../../theory.pdf#bilsc2"><b>bilsc2</b></a> section of the
 * theory guide for more information.
 *
 * \param[in]     f             pointer to field, or null
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$), or nullptr
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_steady_scalar
(
  const cs_field_t           *f,
  const cs_equation_param_t  &eqp,
  bool                        icvflb,
  int                         inc,
  cs_real_t         *restrict pvar,
  const cs_real_t   *restrict pvara,
  const int                   icvfli[],
  const cs_field_bc_coeffs_t *bc_coeffs,
  const cs_real_t             i_massflux[],
  const cs_real_t             b_massflux[],
  const cs_real_t             i_visc[],
  const cs_real_t             b_visc[],
  const cs_real_t             xcpp[],
  cs_real_t         *restrict rhs
)
{
  CS_PROFILE_FUNC_RANGE();

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofafp = bc_coeffs->af;
  const cs_real_t *cofbfp = bc_coeffs->bf;
  /* Precomputed boundary face value */
  const cs_real_t *val_f = bc_coeffs->val_f;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* steady case not ported to GPU */
  cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

  /* Local variables */

  const int f_id = (f != nullptr) ? f->id : -1;
  char var_name[64];

  int w_stride = 1;

  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *df_limiter = nullptr;

  cs_real_t *gweight = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Internal coupling variables */
  cs_real_t *pvar_local = nullptr;
  cs_real_t *pvar_distant = nullptr;
  cs_real_t *df_limiter_local = nullptr;
  int coupling_id = -1;
  cs_lnum_t n_local = 0, n_distant = 0;
  const cs_lnum_t *faces_local = nullptr, *faces_distant = nullptr;
  cs_internal_coupling_t *cpl = nullptr;

  /* Initialization */

  /* Allocate work arrays */

  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  /* Limiters */
  int df_limiter_id = eqp.diffusion_limiter_id;
  if (df_limiter_id > -1)
    df_limiter = cs_field_by_id(df_limiter_id)->val;

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar convection-diffusion]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  bool pure_upwind = (blencp > 0.) ? false : true;

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  bool is_ischcp
    = (   (xcpp != nullptr && ischcp == 0)
       || (xcpp == nullptr && (ischcp == 0 || ischcp == 3 || ischcp == 4)));

  if (   (idiffp != 0 && ircflp == 1)
      || (   iconvp != 0 && pure_upwind == false
          && (ircflp == 1 || isstpp == 0 || is_ischcp))) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)eqp.imligr,
                                    eqp.epsrgr,
                                    eqp.climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    pvar,
                                    gweight, /* Weighted gradient */
                                    grad,
                                    nullptr);

  }
  else {

    h_ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      grad[c_id][0] = 0.;
      grad[c_id][1] = 0.;
      grad[c_id][2] = 0.;
    });
  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && pure_upwind == 0) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC(gradst, n_cells_ext, cs_real_3_t);

      cs_slope_test_gradient(f_id,
                             ctx,
                             (const cs_real_3_t *)grad,
                             gradst,
                             pvar,
                             val_f,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2 || (xcpp != nullptr && ischcp == 4)) {

      CS_MALLOC(gradup, n_cells_ext, cs_real_3_t);

      h_ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        gradup[c_id][0] = 0.;
        gradup[c_id][1] = 0.;
        gradup[c_id][2] = 0.;
      });

      cs_upwind_gradient(ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         pvar,
                         gradup);

    }

  }

  /* ======================================================================
     Contribution from interior faces
     ======================================================================*/

  short *i_upwind = nullptr;
  if (eqp.verbosity >= 2 && iconvp == 1 && pure_upwind == false) {
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* Pure upwind flux
     ================ */

  if (pure_upwind) {

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (xcpp != nullptr) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_2_t fluxij = {0.,0.};

      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

      cs_i_cd_steady_upwind(bldfrp,
                            relaxp,
                            diipf[face_id],
                            djjpf[face_id],
                            grad[ii],
                            grad[jj],
                            pvar[ii],
                            pvar[jj],
                            pvara[ii],
                            pvara[jj],
                            &pifri,
                            &pifrj,
                            &pjfri,
                            &pjfrj,
                            &pip,
                            &pjp,
                            &pipr,
                            &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     cpi,
                     cpj,
                     fluxij);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc[face_id],
                     fluxij);

      if (ii < n_cells)
        rhs[ii] -= fluxij[0];
      if (jj < n_cells)
        rhs[jj] += fluxij[1];

    });

    /* Flux with no slope test or Min/Max Beta limiter
       =============================================== */

  }

  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (xcpp != nullptr) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_2_t fluxij = {0.,0.};

      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii],
                                 df_limiter[jj]), 0.);

      cs_i_cd_steady(bldfrp,
                     ischcp,
                     relaxp,
                     blencp,
                     weight[face_id],
                     cell_cen[ii],
                     cell_cen[jj],
                     i_face_cog[face_id],
                     diipf[face_id],
                     djjpf[face_id],
                     grad[ii],
                     grad[jj],
                     gradup[ii],
                     gradup[jj],
                     pvar[ii],
                     pvar[jj],
                     pvara[ii],
                     pvara[jj],
                     &pifri,
                     &pifrj,
                     &pjfri,
                     &pjfrj,
                     &pip,
                     &pjp,
                     &pipr,
                     &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     cpi,
                     cpj,
                     fluxij);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc[face_id],
                     fluxij);

      if (ii < n_cells)
        rhs[ii] -= fluxij[0];
      if (jj < n_cells)
        rhs[jj] += fluxij[1];

    });

    /* Flux with slope test
       ==================== */

  }

  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (xcpp != nullptr) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_2_t fluxij = {0., 0.};

      bool upwind_switch = false;
      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady_slope_test(&upwind_switch,
                                iconvp,
                                bldfrp,
                                ischcp,
                                relaxp,
                                blencp,
                                blend_st,
                                weight[face_id],
                                i_dist[face_id],
                                cell_cen[ii],
                                cell_cen[jj],
                                i_face_u_normal[face_id],
                                i_face_cog[face_id],
                                diipf[face_id],
                                djjpf[face_id],
                                i_massflux[face_id],
                                grad[ii],
                                grad[jj],
                                gradup[ii],
                                gradup[jj],
                                gradst[ii],
                                gradst[jj],
                                pvar[ii],
                                pvar[jj],
                                pvara[ii],
                                pvara[jj],
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     cpi,
                     cpj,
                     fluxij);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc[face_id],
                     fluxij);

      if (upwind_switch) {

        /* in parallel, face will be counted by one and only one rank */
        if (i_upwind != nullptr && ii < n_cells)
          i_upwind[face_id] = 1;

        if (v_slope_test != nullptr) {
          cs_real_t q_d_vol_ii
            =   std::abs(i_massflux[face_id])
            * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
          cs_real_t q_d_vol_jj
            =   std::abs(i_massflux[face_id])
            * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

          v_slope_test[ii] += q_d_vol_ii;
          v_slope_test[jj] += q_d_vol_jj;
        }

      }

      if (ii < n_cells)
        rhs[ii] -= fluxij[0];
      if (jj < n_cells)
        rhs[jj] += fluxij[1];

    });

  } /* End test on pure upwind */

  if (iwarnp >= 2 && iconvp == 1) {

    if (eqp.verbosity >= 2) {
      cs_gnum_t n_upwind = 0;

      if (i_upwind != nullptr) {
        ctx.parallel_for_reduce_sum
          (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
           (cs_lnum_t i,
            CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
            sum += (cs_gnum_t)i_upwind[i];
          });

        ctx.wait();
        cs_parall_counter(&n_upwind, 1);
      }
      else if (pure_upwind)
        n_upwind = m->n_g_i_faces;

      bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);
    }

  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == false || xcpp != nullptr) {

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t cpi = 1.0;
      if (xcpp != nullptr)
        cpi = xcpp[ii];

      cs_real_t fluxi = 0.;
      cs_real_t pir, pipr;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady(bldfrp,
                     relaxp,
                     diipb[face_id],
                     grad[ii],
                     pvar[ii],
                     pvara[ii],
                     &pir,
                     &pipr);

      /* Compute face value for gradient and diffusion for the
         steady case (relaxation value in iprime) */
      cs_real_t val_f_steady = inc * coefap[face_id] + coefbp[face_id] * pipr;
      cs_real_t flux_steady = inc * cofafp[face_id] + cofbfp[face_id] * pipr;

      cs_b_upwind_flux(iconvp,
                       1.,
                       1,
                       bc_type[face_id],
                       pvar[ii],
                       pir,
                       val_f_steady,
                       b_massflux[face_id],
                       cpi,
                       &fluxi);

      cs_b_diff_flux(idiffp,
                     1., /* thetap */
                     flux_steady,
                     b_visc[face_id],
                     &fluxi);

      rhs[ii] -= fluxi;

    });

    /* The scalar is internal_coupled and an implicit contribution
     * is required */
    if (icoupl > 0) {
      //TODO update BCs
      /* Prepare data for sending */
      CS_MALLOC(pvar_distant, n_distant, cs_real_t);

      for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
        cs_lnum_t face_id = faces_distant[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t pip, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[jj], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[jj],
                       pvar[jj],
                       pvara[jj],
                       &pip,
                       &pipr);
        pvar_distant[ii] = pipr;
      }

      /* Receive data */
      CS_MALLOC(pvar_local, n_local, cs_real_t);
      cs_internal_coupling_exchange_var(cpl,
                                        1, /* Dimension */
                                        pvar_distant,
                                        pvar_local);

      /* Exchange diffusion limiter */
      if (df_limiter != nullptr) {
        CS_MALLOC(df_limiter_local, n_local, cs_real_t);
        cs_internal_coupling_exchange_var(cpl,
                                          1, /* Dimension */
                                          df_limiter,
                                          df_limiter_local);
      }

      /* Flux contribution */
      assert(f != nullptr);
      cs_real_t *hintp = f->bc_coeffs->h_int_tot;
      cs_real_t *hextp = f->bc_coeffs->rcodcl2;
      for (cs_lnum_t ii = 0; ii < n_local; ii++) {
        cs_lnum_t face_id = faces_local[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t surf = b_face_surf[face_id];
        cs_real_t pip, pipr, pjpr;
        cs_real_t fluxi = 0.;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(cs::min(df_limiter_local[ii],
                                   df_limiter[jj]),
                           0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[jj],
                       pvar[jj],
                       pvara[jj],
                       &pip,
                       &pipr);

        pjpr = pvar_local[ii];

        cs_real_t hint = hintp[face_id];
        cs_real_t hext = hextp[face_id];
        cs_real_t heq = _calc_heq(hint, hext)*surf;

        cs_b_diff_flux_coupling(idiffp,
                                pipr,
                                pjpr,
                                heq,
                                &fluxi);

        rhs[jj] -= thetap * fluxi;
      };

      CS_FREE(pvar_local);
      /* Sending structures are no longer needed */
      CS_FREE(pvar_distant);
      if (df_limiter != nullptr)
        CS_FREE(df_limiter_local);
    }

    /* Boundary convective flux is imposed at some faces
       (tagged in icvfli array) */
  }
  else if (icvflb && xcpp == nullptr) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t fluxi = 0.;
      cs_real_t pir, pipr;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady(bldfrp,
                     relaxp,
                     diipb[face_id],
                     grad[ii],
                     pvar[ii],
                     pvara[ii],
                     &pir,
                     &pipr);

      /* Compute face value for gradient and diffusion for the
         steady case (relaxation value in iprime) */
        cs_real_t val_f_steady = inc * coefap[face_id] + coefbp[face_id] * pipr;
        cs_real_t flux_steady = inc * cofafp[face_id] + cofbfp[face_id] * pipr;

        cs_b_imposed_conv_flux(iconvp,
                               1.,
                               1,
                               inc,
                               bc_type[face_id],
                               icvfli[face_id],
                               pvar[ii],
                               pir,
                               pipr,
                               coface[face_id],
                               cofbce[face_id],
                               b_massflux[face_id],
                               1., /* xcpp */
                               val_f_steady,
                               &fluxi);

        cs_b_diff_flux(idiffp,
                       1., /* thetap */
                       flux_steady,
                       b_visc[face_id],
                       &fluxi);

      rhs[ii] -= fluxi;

    });
  }

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(gradup);
  CS_FREE(gradst);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update face flux with convection contribution of a standard transport
 * equation of a scalar field \f$ \varia \f$.
 *
 * <a name="cs_face_convection_scalar"></a>
 *
 * \f[
 * C_\ij = \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 * \f]
 *
 * \param[in]     f_id          pointer to field id, or nullptr
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] i_conv_flux   scalar convection flux at interior faces
 * \param[in,out] b_conv_flux   scalar convection flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_convection_steady_scalar
(
  const cs_field_t           *f,
  const cs_equation_param_t   eqp,
  int                         icvflb,
  int                         inc,
  cs_real_t         *restrict pvar,
  const cs_real_t   *restrict pvara,
  const int                   icvfli[],
  const cs_field_bc_coeffs_t *bc_coeffs,
  const cs_real_t             i_massflux[],
  const cs_real_t             b_massflux[],
  cs_real_t                   i_conv_flux[][2],
  cs_real_t                   b_conv_flux[]
)
{
  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;
  /* Precomputed boundary face value */
  const cs_real_t *val_f = bc_coeffs->val_f;

  const int iconvp = eqp.iconv;
  const int nswrgp = eqp.nswrgr;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal  = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* steady case not ported to GPU */
  cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

  /* Local variables */

  char var_name[64];
  const int f_id = (f != nullptr) ? f->id : -1;

  int w_stride = 1;

  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;

  cs_real_t *df_limiter = nullptr;
  cs_real_t *gweight = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Allocate work arrays */

  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  /* Slope limiters */
  int df_limiter_id = eqp.diffusion_limiter_id;
  if (df_limiter_id > -1)
    df_limiter = cs_field_by_id(df_limiter_id)->val;

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      CS_MALLOC(local_max, n_cells_ext, cs_real_t);
      CS_MALLOC(local_min, n_cells_ext, cs_real_t);
      cs_field_local_extrema_scalar(f_id,
                                    halo_type,
                                    local_max,
                                    local_min);
    }

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar face flux from convection]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  bool pure_upwind = (blencp > 0.) ? false : true;

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  if (   iconvp != 0 && pure_upwind == false
      && (ischcp == 0 || ircflp == 1 || isstpp == 0 || ischcp == 4)) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    pvar,
                                    gweight, /* Weighted gradient */
                                    grad,
                                    nullptr);

  }
  else {

    h_ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      grad[c_id][0] = 0.;
      grad[c_id][1] = 0.;
      grad[c_id][2] = 0.;
    });

  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && pure_upwind == false) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC(gradst, n_cells_ext, cs_real_3_t);

      cs_slope_test_gradient(f_id,
                             ctx,
                             (const cs_real_3_t *)grad,
                             gradst,
                             pvar,
                             val_f,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2) {

      CS_MALLOC(gradup, n_cells_ext, cs_real_3_t);

      h_ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        gradup[c_id][0] = 0.;
        gradup[c_id][1] = 0.;
        gradup[c_id][2] = 0.;
      });

      cs_upwind_gradient(ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         pvar,
                         gradup);

    }

  }

  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  short *i_upwind = nullptr;
  if (eqp.verbosity >= 2 && iconvp == 1 && pure_upwind == false) {
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* Pure upwind flux
     ================ */

  if (pure_upwind) {

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

      cs_i_cd_steady_upwind(bldfrp,
                            relaxp,
                            diipf[face_id],
                            djjpf[face_id],
                            grad[ii],
                            grad[jj],
                            pvar[ii],
                            pvar[jj],
                            pvara[ii],
                            pvara[jj],
                            &pifri,
                            &pifrj,
                            &pjfri,
                            &pjfrj,
                            &pip,
                            &pjp,
                            &pipr,
                            &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     1., /* xcpp */
                     1., /* xcpp */
                     i_conv_flux[face_id]);
    });
  }

  /* Flux with no slope test or Min/Max Beta limiter
     =============================================== */

  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady(bldfrp,
                     ischcp,
                     relaxp,
                     blencp,
                     weight[face_id],
                     cell_cen[ii],
                     cell_cen[jj],
                     i_face_cog[face_id],
                     diipf[face_id],
                     djjpf[face_id],
                     grad[ii],
                     grad[jj],
                     gradup[ii],
                     gradup[jj],
                     pvar[ii],
                     pvar[jj],
                     pvara[ii],
                     pvara[jj],
                     &pifri,
                     &pifrj,
                     &pjfri,
                     &pjfrj,
                     &pip,
                     &pjp,
                     &pipr,
                     &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     1., /* xcpp */
                     1., /* xcpp */
                     i_conv_flux[face_id]);

    });

  }

  /* Flux with slope test or NVD/TVD limiter
     ======================================= */

  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      bool upwind_switch = false;
      cs_real_t pifri, pjfri, pifrj, pjfrj;
      cs_real_t pip, pjp, pipr, pjpr;

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady_slope_test(&upwind_switch,
                                iconvp,
                                bldfrp,
                                ischcp,
                                relaxp,
                                blencp,
                                blend_st,
                                weight[face_id],
                                i_dist[face_id],
                                cell_cen[ii],
                                cell_cen[jj],
                                i_face_u_normal[face_id],
                                i_face_cog[face_id],
                                diipf[face_id],
                                djjpf[face_id],
                                i_massflux[face_id],
                                grad[ii],
                                grad[jj],
                                gradup[ii],
                                gradup[jj],
                                gradst[ii],
                                gradst[jj],
                                pvar[ii],
                                pvar[jj],
                                pvara[ii],
                                pvara[jj],
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     1,
                     pvar[ii],
                     pvar[jj],
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_massflux[face_id],
                     1., /* xcpp */
                     1., /* xcpp */
                     i_conv_flux[face_id]);

      if (upwind_switch) {

        /* in parallel, face will be counted by one and only one rank */
        if (i_upwind != nullptr && ii < n_cells)
          i_upwind[face_id] = 1;

        if (v_slope_test != nullptr) {
          cs_real_t q_d_vol_ii
            =   std::abs(i_massflux[face_id])
            * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
          cs_real_t q_d_vol_jj
            =   std::abs(i_massflux[face_id])
            * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

          v_slope_test[ii] += q_d_vol_ii;
          v_slope_test[jj] += q_d_vol_jj;
        }

      }

    });

  } /* End pure_upwind */

  if (iwarnp >= 2 && iconvp == 1) {

    if (eqp.verbosity >= 2) {
      cs_gnum_t n_upwind = 0;

      if (i_upwind != nullptr) {
        ctx.parallel_for_reduce_sum
          (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
           (cs_lnum_t i,
            CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
            sum += (cs_gnum_t)i_upwind[i];
          });

        ctx.wait();
        cs_parall_counter(&n_upwind, 1);
      }
      else if (pure_upwind)
        n_upwind = m->n_g_i_faces;

      bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);
    }

  }

  /* ======================================================================
    ---> Contribution from boundary faces
    ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pir, pipr;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady(bldfrp,
                     relaxp,
                     diipb[face_id],
                     grad[ii],
                     pvar[ii],
                     pvara[ii],
                     &pir,
                     &pipr);

      cs_b_upwind_flux(iconvp,
                       1.,
                       1,
                       inc,
                       bc_type[face_id],
                       pvar[ii],
                       pir,
                       pipr,
                       coefap[face_id],
                       coefbp[face_id],
                       b_massflux[face_id],
                       1., /* xcpp */
                       &(b_conv_flux[face_id]));

    });

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */
  }
  else if (icvflb == 1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pir, pipr;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady(bldfrp,
                     relaxp,
                     diipb[face_id],
                     grad[ii],
                     pvar[ii],
                     pvara[ii],
                     &pir,
                     &pipr);

      /* Compute face value for gradient and diffusion for the
         steady case (relaxation value in iprime) */
      cs_real_t val_f_steady = inc * coefap[face_id] + coefbp[face_id] * pipr;

      cs_b_imposed_conv_flux(iconvp,
                             1.,
                             1,
                             inc,
                             bc_type[face_id],
                             icvfli[face_id],
                             pvar[ii],
                             pir,
                             pipr,
                             coface[face_id],
                             cofbce[face_id],
                             b_massflux[face_id],
                             1., /* xcpp */
                             val_f_steady,
                             &(b_conv_flux[face_id]));

    });
  }

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(gradup);
  CS_FREE(gradst);
  CS_FREE(local_max);
  CS_FREE(local_min);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector or tensor field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ or \f$ \tens{Rhs} \f$
 * is updated as * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Or in the tensor case:
 *
 * \f[
 *  \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
 * \f]
 *
 * \deprecated
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - rhs has already been initialized!
 * - mind the minus sign
 *
 * \param[in]      f             pointer to field, or nullptr
 * \param[in]      name          pointer to associated field or array name
 * \param[in]      eqp           equation parameters
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                                - 0 upwind scheme at all boundary faces
 *                                - 1 imposed flux at some boundary faces
 * \param[in]     inc            indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]      pvar          solved velocity (current time step)
 * \param[in]      pvara         solved velocity (previous time step)
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                                - 0 upwind scheme
 *                                - 1 imposed flux
 * \param[in]      bc_coeffs     boundary conditions structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                                at border faces for the r.h.s.
 * \param[in]      grad          associated gradient
 * \param[in, out] rhs           right hand side \f$ \vect{Rhs} \f$
 *                               or \f$ \tens{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_convection_diffusion_steady_strided
(
  cs_field_t                   *f,
  const char                   *var_name,
  const cs_equation_param_t    &eqp,
  int                           icvflb,
  int                           inc,
  cs_real_t                   (*pvar)[stride],
  const cs_real_t             (*pvara)[stride],
  const int                     icvfli[],
  const cs_field_bc_coeffs_t   *bc_coeffs,
  const cs_real_t               i_massflux[],
  const cs_real_t               b_massflux[],
  const cs_real_t               i_visc[],
  const cs_real_t               b_visc[],
  cs_real_t          (*restrict grad)[stride][3],
  cs_real_t          (*restrict rhs)[stride]
)
{
  using grad_t = cs_real_t[stride][3];
  using var_t = cs_real_t[stride];
  using b_t = cs_real_t[stride][stride];

  const var_t *coefa = (const var_t  *)bc_coeffs->a;
  const b_t   *coefb = (const b_t *)bc_coeffs->b;
  const var_t *cofaf = (const var_t *)bc_coeffs->af;
  const b_t   *cofbf = (const b_t *)bc_coeffs->bf;

  /* precomputed boundary face value */
  const var_t *val_f = (const var_t  *)bc_coeffs->val_f;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  cs_real_2_t *i_f_face_factor = nullptr;
  cs_real_t *b_f_face_factor = nullptr;

  /* Flux limiter */
  cs_real_t *df_limiter = nullptr;
  int df_limiter_id = eqp.diffusion_limiter_id;
  if (df_limiter_id > -1)
    df_limiter = cs_field_by_id(df_limiter_id)->val;

  /* Discontinuous porous treatment */
  if (cs_glob_porous_model == 3 && f == CS_F_(vel)) {
    i_f_face_factor = fvq->i_f_face_factor;
    b_f_face_factor = fvq->b_f_face_factor;
  }

  const var_t *coface = nullptr;
  const b_t *cofbce = nullptr;

  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* steady case not ported to GPU */
  cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

  /*==========================================================================*/

  /* 1. Initialization */

  /* Allocate work arrays */

  grad_t *grdpa = nullptr;

  /* Slope limiters */

  if (iwarnp >= 2 && iconvp == 1) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  int pure_upwind = (blencp > 0.) ? 0 : 1;

  /* ======================================================================
     ---> Compute uncentered gradient grdpa for the slope test
     ======================================================================*/

  if (iconvp > 0 && pure_upwind == 0 && isstpp == 0) {
    CS_MALLOC_HD(grdpa, n_cells_ext, grad_t, cs_alloc_mode);

    cs_slope_test_gradient_strided<stride>(ctx,
                                           (const grad_t *)grad,
                                           grdpa,
                                           pvar,
                                           val_f,
                                           i_massflux);
    ctx.wait();
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  short *i_upwind = nullptr;
  if (iwarnp >= 2 && iconvp == 1) {
    cs_lnum_t n_i_faces = m->n_i_faces;
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    h_ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* Pure upwind flux
     ================ */

  if (pure_upwind == 1) {

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      /* in parallel, face will be counted by one and only one rank */
      if (i_upwind != nullptr && ii < n_cells) {
        i_upwind[face_id] = 1;
      }

      cs_real_t fluxi[stride], fluxj[stride] ;
      for (cs_lnum_t isou =  0; isou < stride; isou++) {
        fluxi[isou] = 0;
        fluxj[isou] = 0;
      }

      var_t pip, pjp, pipr, pjpr;
      var_t pifri, pifrj, pjfri, pjfrj;
      var_t _pi, _pj, _pia, _pja;

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = pvar[ii][i];
        _pj[i]  = pvar[jj][i];
        _pia[i] = pvara[ii][i];
        _pja[i] = pvara[jj][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (i_f_face_factor != nullptr) {
        const cs_nreal_t *n = i_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady_upwind_strided<stride>(bldfrp,
                                            relaxp,
                                            diipf[face_id],
                                            djjpf[face_id],
                                            grad[ii],
                                            grad[jj],
                                            _pi,
                                            _pj,
                                            _pia,
                                            _pja,
                                            pifri,
                                            pifrj,
                                            pjfri,
                                            pjfrj,
                                            pip,
                                            pjp,
                                            pipr,
                                            pjpr);

      cs_i_conv_flux_strided<stride>(iconvp,
                                     1.,
                                     1,
                                     pvar[ii],
                                     pvar[jj],
                                     pifri,
                                     pifrj,
                                     pjfri,
                                     pjfrj,
                                     i_massflux[face_id],
                                     fluxi,
                                     fluxj);

      cs_i_diff_flux_strided<stride>(idiffp,
                                     1.,
                                     pip,
                                     pjp,
                                     pipr,
                                     pjpr,
                                     i_visc[face_id],
                                     fluxi,
                                     fluxj);

      if (ii < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[ii][isou] -= fluxi[isou];
      }
      if (jj < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[jj][isou] += fluxj[isou];
      }

    });

  }

  /* --> Flux with no slope test
     ============================*/

  else if (isstpp == 1) {

    if (ischcp < 0 || ischcp == 2 || ischcp > 3) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      var_t fluxi, fluxj;
      for (cs_lnum_t isou =  0; isou < stride; isou++) {
        fluxi[isou] = 0;
        fluxj[isou] = 0;
      }
      var_t pip, pjp, pipr, pjpr;
      var_t pifri, pifrj, pjfri, pjfrj;
      var_t _pi, _pj, _pia, _pja;

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = pvar[ii][i];
        _pj[i]  = pvar[jj][i];
        _pia[i] = pvara[ii][i];
        _pja[i] = pvara[jj][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (i_f_face_factor != nullptr) {
        const cs_nreal_t *n = i_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady_strided<stride>(bldfrp,
                                     ischcp,
                                     relaxp,
                                     blencp,
                                     weight[face_id],
                                     cell_cen[ii],
                                     cell_cen[jj],
                                     i_face_cog[face_id],
                                     diipf[face_id],
                                     djjpf[face_id],
                                     grad[ii],
                                     grad[jj],
                                     _pi,
                                     _pj,
                                     _pia,
                                     _pja,
                                     pifri,
                                     pifrj,
                                     pjfri,
                                     pjfrj,
                                     pip,
                                     pjp,
                                     pipr,
                                     pjpr);

      cs_i_conv_flux_strided<stride>(iconvp,
                                     1.,
                                     1,
                                     pvar[ii],
                                     pvar[jj],
                                     pifri,
                                     pifrj,
                                     pjfri,
                                     pjfrj,
                                     i_massflux[face_id],
                                     fluxi,
                                     fluxj);

      cs_i_diff_flux_strided<stride>(idiffp,
                                     1.,
                                     pip,
                                     pjp,
                                     pipr,
                                     pjpr,
                                     i_visc[face_id],
                                     fluxi,
                                     fluxj);

      if (ii < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[ii][isou] -= fluxi[isou];
      }
      if (jj < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[jj][isou] += fluxj[isou];
      }

    });

  }

  /* Flux with slope test
     ==================== */

  else {

    if (ischcp < 0 || ischcp > 1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      var_t fluxi, fluxj;
      for (cs_lnum_t isou =  0; isou < stride; isou++) {
        fluxi[isou] = 0;
        fluxj[isou] = 0;
      }
      var_t pip, pjp, pipr, pjpr;
      var_t pifri, pifrj, pjfri, pjfrj;
      bool upwind_switch = false;
      var_t _pi, _pj, _pia, _pja;

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = pvar[ii][i];
        _pj[i]  = pvar[jj][i];
        _pia[i] = pvara[ii][i];
        _pja[i] = pvara[jj][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (i_f_face_factor != nullptr) {
        const cs_nreal_t *n = i_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      cs_i_cd_steady_slope_test_strided<stride>(&upwind_switch,
                                                iconvp,
                                                bldfrp,
                                                ischcp,
                                                relaxp,
                                                blencp,
                                                blend_st,
                                                weight[face_id],
                                                i_dist[face_id],
                                                cell_cen[ii],
                                                cell_cen[jj],
                                                i_face_u_normal[face_id],
                                                i_face_cog[face_id],
                                                diipf[face_id],
                                                djjpf[face_id],
                                                i_massflux[face_id],
                                                grad[ii],
                                                grad[jj],
                                                grdpa[ii],
                                                grdpa[jj],
                                                _pi,
                                                _pj,
                                                _pia,
                                                _pja,
                                                pifri,
                                                pifrj,
                                                pjfri,
                                                pjfrj,
                                                pip,
                                                pjp,
                                                pipr,
                                                pjpr);

      cs_i_conv_flux_strided<stride>(iconvp,
                                     1.,
                                     1,
                                     pvar[ii],
                                     pvar[jj],
                                     pifri,
                                     pifrj,
                                     pjfri,
                                     pjfrj,
                                     i_massflux[face_id],
                                     fluxi,
                                     fluxj);

      cs_i_diff_flux_strided<stride>(idiffp,
                                     1.,
                                     pip,
                                     pjp,
                                     pipr,
                                     pjpr,
                                     i_visc[face_id],
                                     fluxi,
                                     fluxj);

      if (ii < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[ii][isou] -= fluxi[isou];
      }
      if (jj < n_cells) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          rhs[jj][isou] += fluxj[isou];
      }

    });

  } /* iupwin */

  if (iwarnp >= 2 && i_upwind != nullptr) {
    cs_gnum_t n_upwind = 0;
    const cs_lnum_t n_i_faces = m->n_i_faces;

    if (i_upwind != nullptr) {
      h_ctx.parallel_for_reduce_sum
        (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
         (cs_lnum_t i,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
        sum += (cs_gnum_t)i_upwind[i];
      });

      h_ctx.wait();
      cs_parall_counter(&n_upwind, 1);
    }

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE(i_upwind);
  }

  /* ======================================================================
     Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      var_t fluxi;
      for (cs_lnum_t isou =  0; isou < stride; isou++) {
        fluxi[isou] = 0;
      }
      var_t pir, pipr, _val_f, _flux_d;
      var_t _pi, _pia;
      for (int i = 0; i < stride; i++) {
        _pi[i]  = pvar[ii][i];
        _pia[i] = pvara[ii][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (b_f_face_factor != nullptr) {
        const cs_nreal_t *n = i_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pia);
      }

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady_strided<stride>(bldfrp,
                                     relaxp,
                                     diipb[face_id],
                                     (const cs_real_3_t *)grad[ii],
                                     _pi,
                                     _pia,
                                     pir,
                                     pipr);

      /* Compute face value for gradient and diffusion for the
         steady case (relaxation value in iprime) */
      for (cs_lnum_t i = 0; i < stride; i++) {
        _val_f[i] = inc * coefa[face_id][i];
        _flux_d[i] = inc * cofaf[face_id][i];

        for (cs_lnum_t j = 0; j < stride; j++) {
          _val_f[i] += coefb[face_id][j][i] * pipr[j];
          _flux_d[i] += cofbf[face_id][j][i] * pipr[j];
        }
      }

      cs_b_upwind_flux_strided<stride>(iconvp,
                                       1., /* thetap */
                                       1, /* imasac */
                                       bc_type[face_id],
                                       _pi,
                                       pir,
                                       b_massflux[face_id],
                                       _val_f,
                                       fluxi);

      cs_b_diff_flux_strided<stride>(idiffp,
                                     1., /* thetap */
                                     b_visc[face_id],
                                     _flux_d,
                                     fluxi);

      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        rhs[ii][isou] -= fluxi[isou];
      } /* isou */

    });

  }

  /* Boundary convective flux imposed at some faces (tags in icvfli array) */

  else if (icvflb == 1 && icvfli != nullptr) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f != nullptr) {
      coface = (const var_t *)(f->bc_coeffs->ac);
      cofbce = (const b_t *)(f->bc_coeffs->bc);
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      var_t fluxi, pir, pipr, _val_f, flux;

      for (cs_lnum_t isou =  0; isou < stride; isou++) {
        fluxi[isou] = 0;
      }
      var_t _pi, _pia;

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = pvar[ii][i];
        _pia[i] = pvara[ii][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (b_f_face_factor != nullptr) {
        const cs_nreal_t *n = b_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pia);
      }

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_steady_strided<stride>(bldfrp,
                                     relaxp,
                                     diipb[face_id],
                                     grad[ii],
                                     _pi,
                                     _pia,
                                     pir,
                                     pipr);

      /* Compute face value for gradient and diffusion for the
         steady case (relaxation value in iprime) */
      for (cs_lnum_t i = 0; i < stride; i++) {
        _val_f[i] = inc * coefa[face_id][i];
        flux[i] = inc * cofaf[face_id][i];

        for (cs_lnum_t j = 0; j < stride; j++) {
          _val_f[i] += coefb[face_id][j][i] * pipr[j];
          flux[i] += cofbf[face_id][j][i] * pipr[j];
        }
      }

      cs_b_imposed_conv_flux_strided<stride>(iconvp,
                                             1., /* thetap */
                                             1., /* imasac */
                                             inc,
                                             bc_type[face_id],
                                             icvfli[face_id],
                                             pvar[ii],
                                             pir,
                                             pipr,
                                             coface[face_id],
                                             cofbce[face_id],
                                             b_massflux[face_id],
                                             _val_f,
                                             fluxi);

      cs_b_diff_flux_strided<stride>(idiffp,
                                     1., /* thetap */
                                     b_visc[face_id],
                                     flux,
                                     fluxi);

      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        rhs[ii][isou] -= fluxi[isou];
      }

    });

  }

  /* Free memory */
  CS_FREE(grdpa);
}

// Force instanciation

template void
cs_convection_diffusion_steady_strided
(
  cs_field_t                   *f,
  const char                   *var_name,
  const cs_equation_param_t    &eqp,
  int                           icvflb,
  int                           inc,
  cs_real_t                   (*pvar)[3],
  const cs_real_t             (*pvara)[3],
  const int                     icvfli[],
  const cs_field_bc_coeffs_t   *bc_coeffs,
  const cs_real_t               i_massflux[],
  const cs_real_t               b_massflux[],
  const cs_real_t               i_visc[],
  const cs_real_t               b_visc[],
  cs_real_t          (*restrict grad)[3][3],
  cs_real_t          (*restrict rhs)[3]
);

template void
cs_convection_diffusion_steady_strided
(
  cs_field_t                   *f,
  const char                   *var_name,
  const cs_equation_param_t    &eqp,
  int                           icvflb,
  int                           inc,
  cs_real_t                   (*pvar)[6],
  const cs_real_t             (*pvara)[6],
  const int                     icvfli[],
  const cs_field_bc_coeffs_t   *bc_coeffs,
  const cs_real_t               i_massflux[],
  const cs_real_t               b_massflux[],
  const cs_real_t               i_visc[],
  const cs_real_t               b_visc[],
  cs_real_t          (*restrict grad)[6][3],
  cs_real_t          (*restrict rhs)[6]
);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
