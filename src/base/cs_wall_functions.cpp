/*============================================================================
 * Wall functions
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "turb/cs_turbulence_model.h"
#include "cdo/cs_domain.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_wall_functions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_functions.cpp
        Wall functions descriptor and computation.
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_wall_functions_t

  \brief wall functions descriptor.

  Members of this wall functions descriptor are publicly accessible, to allow
  for concise syntax, as it is expected to be used in many places.

  \var  cs_wall_functions_t::iwallf
        Indicates the type of wall function used for the velocity
        boundary conditions on a frictional wall.\n
        - CS_WALL_F_DISABLED: no wall functions
        - CS_WALL_F_1SCALE_POWER: one scale of friction velocities (power law)
        - CS_WALL_F_1SCALE_LOG: one scale of friction velocities (log law)
        - CS_WALL_F_2SCALES_LOG: two scales of friction velocities (log law)
        - CS_WALL_F_SCALABLE_2SCALES_LOG: two scales of friction velocities
          (log law) (scalable wall functions)
        - CS_WALL_F_2SCALES_VDRIEST: two scales of friction velocities
          (mixing length based on V. Driest analysis)
        - CS_WALL_F_2SCALES_SMOOTH_ROUGH: wall function unifying rough and smooth
          friction regimes
        - CS_WALL_F_2SCALES_CONTINUOUS: All \f$ y^+ \f$  for low Reynolds models\n
        \ref iwallf is initialised to CS_WALL_F_1SCALE_LOG for \ref model = 10,
          40, 41 or 70
        (mixing length, LES and Spalart Allmaras).\n
        \ref iwallf is initialised to CS_WALL_F_DISABLED for \ref model = 0, 32,
          50 or 51\n
        \ref iwallf is initialised to CS_WALL_F_2SCALES_LOG for \ref model = 20,
          21, 30, 31 or 60
        (\f$k-\epsilon\f$, \f$R_{ij}-\epsilon\f$ LRR, \f$R_{ij}-\epsilon\f$ SSG
          and \f$k-\omega\f$ SST models).\n
        The v2f model (\ref model=50) is not designed to use wall functions
        (the mesh must be low Reynolds).\n
        The value \ref iwallf = CS_WALL_F_2SCALES_LOG is not compatible with
          \ref model=0, 10, 40
        or 41 (laminar, mixing length and LES).\n
        Concerning the \f$k-\epsilon\f$ and \f$R_{ij}-\epsilon\f$ models, the
        two-scales model is usually at least as satisfactory as the one-scale
        model.\n
        The scalable wall function allows to virtually shift the wall when
        necessary in order to be always in a logarithmic layer. It is used to make up for
        the problems related to the use of High-Reynolds models on very refined
        meshes.\n
        Useful if \ref model is different from 50.
  \var  cs_wall_functions_t::iwalfs
        wall functions for scalar
        - CS_WALL_F_S_ARPACI_LARSEN: three layers (Arpaci and Larsen) or two layers (Prandtl-Taylor) for
             Prandtl number smaller than 0.1
        - CS_WALL_F_S_VDRIEST: consistant with the 2 scales wall function for velocity using Van
             Driest mixing length
        - CS_WALL_F_S_LOUIS: default wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.
        - CS_WALL_F_S_MONIN_OBUKHOV: Monin Obukhov wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.

  \var  cs_wall_functions_t::ypluli
        limit value of \f$y^+\f$ for the viscous sublayer

        \ref ypluli depends on the chosen wall function: it is initialized to
        10.88 for the scalable wall function (\ref iwallf=CS_WALL_F_SCALABLE_2SCALES_LOG), otherwise it is
        initialized to \f$1/\kappa\approx 2,38\f$. In LES, \ref ypluli is taken
        by default to be 10.88. Always useful.

*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Wall functions structure and associated pointers

   ypluli is set to -grand*10. If the user has not changed this value,
   its value is modified in cs_parameters_*_complete (10.88 with invariant wall
   laws, 1/kappa otherwise). */

static cs_wall_functions_t  _wall_functions =
{
  .iwallf = CS_WALL_F_UNSET,
  .iwalfs = CS_WALL_F_S_UNSET,
  .ypluli = -1e13
};

const cs_wall_functions_t  * cs_glob_wall_functions = &_wall_functions;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *! \brief Provide access to cs_glob_wall_functions
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_wall_functions_t *
cs_get_glob_wall_functions(void)
{
  return &_wall_functions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the friction velocity and \f$y^+\f$ / \f$u^+\f$.
 *
 * \param[in]     iwallf        wall function type
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     t_visc        turbulent kinematic viscosity
 * \param[in]     vel           wall projected cell center velocity
 * \param[in]     y             wall distance
 * \param[in]     rough_d       roughness length scale
 *                              (not sand grain roughness)
 * \param[in]     rnnb          \f$\vec{n}.(\tens{R}\vec{n})\f$
 * \param[in]     kinetic_en    turbulent kinetic energy (cell center)
 * \param[in,out] iuntur        indicator: 0 in the viscous sublayer
 * \param[in,out] nsubla        counter of cell in the viscous sublayer
 * \param[in,out] nlogla        counter of cell in the log-layer
 * \param[out]    ustar         friction velocity
 * \param[out]    uk            friction velocity
 * \param[out]    yplus         dimensionless distance to the wall
 * \param[out]    ypup          yplus projected vel ratio
 * \param[out]    cofimp        \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                              turbulence production
 * \param[out]    dplus         dimensionless shift to the wall for scalable
 *                              wall functions
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_velocity(cs_wall_f_type_t  iwallf,
                           cs_real_t         l_visc,
                           cs_real_t         t_visc,
                           cs_real_t         vel,
                           cs_real_t         y,
                           cs_real_t         rough_d,
                           cs_real_t         rnnb,
                           cs_real_t         kinetic_en,
                           int              *iuntur,
                           cs_gnum_t        *nsubla,
                           cs_gnum_t        *nlogla,
                           cs_real_t        *ustar,
                           cs_real_t        *uk,
                           cs_real_t        *yplus,
                           cs_real_t        *ypup,
                           cs_real_t        *cofimp,
                           cs_real_t        *dplus)
{
  cs_real_t lmk;
  bool wf = true;

  /* Pseudo shift of the wall, 0 by default */
  *dplus = 0.;

  /* Activation of wall function by default */
  *iuntur = 1;

  /* Sand Grain roughness */
  cs_real_t sg_rough = rough_d * exp(cs_turb_xkappa*cs_turb_cstlog_rough);

  switch (iwallf) {
  case CS_WALL_F_DISABLED:
    cs_wall_functions_disabled(l_visc,
                               t_visc,
                               vel,
                               y,
                               iuntur,
                               nsubla,
                               nlogla,
                               ustar,
                               uk,
                               yplus,
                               dplus,
                               ypup,
                               cofimp);
    break;
  case CS_WALL_F_1SCALE_POWER:
    cs_wall_functions_1scale_power(l_visc,
                                   vel,
                                   y,
                                   iuntur,
                                   nsubla,
                                   nlogla,
                                   ustar,
                                   uk,
                                   yplus,
                                   ypup,
                                   cofimp);
    break;
  case CS_WALL_F_1SCALE_LOG:
    cs_wall_functions_1scale_log(l_visc,
                                 vel,
                                 y,
                                 iuntur,
                                 nsubla,
                                 nlogla,
                                 ustar,
                                 uk,
                                 yplus,
                                 ypup,
                                 cofimp);
    break;
  case CS_WALL_F_2SCALES_LOG:
    cs_wall_functions_2scales_log(l_visc,
                                  t_visc,
                                  vel,
                                  y,
                                  kinetic_en,
                                  iuntur,
                                  nsubla,
                                  nlogla,
                                  ustar,
                                  uk,
                                  yplus,
                                  ypup,
                                  cofimp);
    break;
  case CS_WALL_F_SCALABLE_2SCALES_LOG:
    cs_wall_functions_2scales_scalable(l_visc,
                                       t_visc,
                                       vel,
                                       y,
                                       kinetic_en,
                                       iuntur,
                                       nsubla,
                                       nlogla,
                                       ustar,
                                       uk,
                                       yplus,
                                       dplus,
                                       ypup,
                                       cofimp);
    break;
  case CS_WALL_F_2SCALES_VDRIEST:
   cs_wall_functions_2scales_vdriest(rnnb,
                                      l_visc,
                                      vel,
                                      y,
                                      kinetic_en,
                                      iuntur,
                                      nsubla,
                                      nlogla,
                                      ustar,
                                      uk,
                                      yplus,
                                      ypup,
                                      cofimp,
                                      &lmk,
                                      sg_rough,
                                      wf);
    break;
  case CS_WALL_F_2SCALES_SMOOTH_ROUGH:
    cs_wall_functions_2scales_smooth_rough(l_visc,
                                           t_visc,
                                           vel,
                                           y,
                                           rough_d,
                                           kinetic_en,
                                           iuntur,
                                           nsubla,
                                           nlogla,
                                           ustar,
                                           uk,
                                           yplus,
                                           dplus,
                                           ypup,
                                           cofimp);
    break;
  case CS_WALL_F_2SCALES_CONTINUOUS:
    cs_wall_functions_2scales_continuous(rnnb,
                                         l_visc,
                                         t_visc,
                                         vel,
                                         y,
                                         kinetic_en,
                                         iuntur,
                                         nsubla,
                                         nlogla,
                                         ustar,
                                         uk,
                                         yplus,
                                         ypup,
                                         cofimp);
    break;
  default:
    break;
  }

  /* To be coherent with a wall function, clip it to 0 */
  *cofimp = cs::max(*cofimp, 0.);
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Compute the correction of the exchange coefficient between the
 *         fluid and the wall for a turbulent flow.
 *
 *  This is function of the dimensionless
 *  distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.
 *
 *  Then the return coefficient reads:
 *  \f[
 *  h_{tur} = Pr \dfrac{y^+}{T^+}
 *  \f]
 *
 * \param[in]     iwalfs        type of wall functions for scalar
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     prl           laminar Prandtl number
 * \param[in]     prt           turbulent Prandtl number
 * \param[in]     rough_t       scalar roughness lenghth scale
 * \param[in]     uk            velocity scale based on TKE
 * \param[in]     yplus         dimensionless distance to the wall
 * \param[in]     dplus         dimensionless distance for scalable
 *                              wall functions
 * \param[out]    htur          corrected exchange coefficient
 * \param[out]    yplim         value of the limit for \f$ y^+ \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_scalar(cs_wall_f_s_type_t  iwalfs,
                         cs_real_t           l_visc,
                         cs_real_t           prl,
                         cs_real_t           prt,
                         cs_real_t           rough_t,
                         cs_real_t           uk,
                         cs_real_t           yplus,
                         cs_real_t           dplus,
                         cs_real_t          *htur,
                         cs_real_t          *yplim)
{
  switch (iwalfs) {
  case CS_WALL_F_S_ARPACI_LARSEN:
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  case CS_WALL_F_S_VDRIEST:
    cs_wall_functions_s_vdriest(prl,
                                prt,
                                yplus,
                                htur);
    break;
  case CS_WALL_F_S_SMOOTH_ROUGH:
    cs_wall_functions_s_smooth_rough(l_visc,
                                     prl,
                                     prt,
                                     rough_t,
                                     uk,
                                     yplus,
                                     dplus,
                                     htur);
    break;
  default:
    /* TODO Monin Obukhov or Louis atmospheric wall function
     * must be adapted to smooth wall functions.
     * Arpaci and Larsen wall functions are put as in previous versions of
     * code_Saturne.
     * */
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
