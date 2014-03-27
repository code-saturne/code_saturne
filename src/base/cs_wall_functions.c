/*============================================================================
 * Wall functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

#include "cs_config.h"
#include "cs_defs.h"

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_log.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_functions.c
        Wall functions descriptor and computation.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_wall_functions_t

  \brief wall functions descriptor.

  Members of this wall functions descriptor are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_wall_functions_t::ideuch
        wall functions
        - 0: one scale of friction velocities
        - 1: two scale of friction velocities
        - 2: scalable wall functions
  \var  cs_wall_functions_t::iwallt
        exchange coefficient correlation
        - 0: not use by default
        - 1: exchange coefficient computed with a correlation
  \var  cs_wall_functions_t::ilogpo
        wall functions with
        - 0: a power lay (deprecated)
        - 1: a log lay
*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* wall functions structure and associated pointer */

static cs_wall_functions_t  _wall_functions =
  {-999, 0, 1};

const cs_wall_functions_t  *cs_glob_wall_functions = &_wall_functions;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_wall_functions_get_pointers(int     **ideuch,
                                  int     **iwallt,
                                  int     **ilogpo);

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Power law: Werner & Wengle
 *
 * parameters:
 *   ypluli   <-- y+ limit
 *   apow     <-- Coef of the Wener law
 *   bpow     <-- Coef of the Wener law
 *   dpow     <-- Coef of the Wener law
 *   l_visc   <-- kinematic viscosity
 *   vel      <-- wall projected cell center velocity
 *   y        <-- wall distance
 *   iuntur   <-- indicator: 0 in the viscous sublayer
 *   nsubla   <-- counter of cell in the viscous sublayer
 *   nlogla   <-- counter of cell in the log-layer
 *   ustar    --> friction velocity
 *   uk       --> friction velocity
 *   yplus    --> non-dimension wall distance
 *   ypup     --> yplus projected vel ratio
 *   cofimp   --> |U_F|/|U_I^p| to ensure a good turbulence production
 *----------------------------------------------------------------------------*/

static void
_1scale_power_law(cs_real_t   ypluli,
                  cs_real_t   apow,
                  cs_real_t   bpow,
                  cs_real_t   dpow,
                  cs_real_t   l_visc,
                  cs_real_t   vel,
                  cs_real_t   y,
                  int        *iuntur,
                  cs_lnum_t  *nsubla,
                  cs_lnum_t  *nlogla,
                  cs_real_t  *ustar,
                  cs_real_t  *uk,
                  cs_real_t  *yplus,
                  cs_real_t  *ypup,
                  cs_real_t  *cofimp)
{
  const double ydvisc =  y / l_visc;

  /* Compute the friction velocity ustar */

  *ustar = pow((vel/(apow * pow(ydvisc, bpow))), dpow);
  *uk = *ustar;
  *yplus = *ustar * ydvisc;

  /* In the viscous sub-layer: U+ = y+ */
  if (*yplus <= ypluli) {

    *ustar = sqrt(vel / ydvisc);
    *yplus = *ustar * ydvisc;
    *uk = *ustar;
    *ypup = 1.;
    *cofimp = 0.;

    /* Disable the wall funcion count the cell in the viscous sub-layer */
    *iuntur = 0;
    *nsubla += 1;

  /* In the log layer */
  } else {
    *ypup = pow(vel, 2. * dpow-1.) / pow(apow, 2. * dpow);
    *cofimp = 1. + bpow * pow(*ustar, bpow + 1. - 1./dpow)
                 * (pow(2., bpow - 1.) - 2.);

    /* Count the cell in the log layer */
    *nlogla += 1;

  }
}

/*----------------------------------------------------------------------------
 * Log law: piecewise linear and log, with one velocity scale base on the
 * friction.
 *
 * parameters:
 *   ifac   <-- face number
 *   xkappa <-- Von Karman constant
 *   cstlog <-- Log law constant
 *   ypluli <-- y+ limit
 *   apow   <-- Coef of the Wener law
 *   bpow   <-- Coef of the Wener law
 *   dpow   <-- Coef of the Wener law
 *   l_visc <-- kinematic viscosity
 *   vel    <-- wall projected cell center velocity
 *   y      <-- wall distance
 *   iuntur <-- indicator: 0 in the viscous sublayer
 *   nsubla <-- counter of cell in the viscous sublayer
 *   nlogla <-- counter of cell in the log-layer
 *   ustar  --> friction velocity
 *   uk     --> friction velocity
 *   yplus  --> non-dimension wall distance
 *   ypup   --> yplus projected vel ratio
 *   cofimp --> |U_F|/|U_I^p| to ensure a good turbulence production
 *----------------------------------------------------------------------------*/

static void
_1scale_log_law(cs_lnum_t    ifac,
                cs_real_t    xkappa,
                cs_real_t    cstlog,
                cs_real_t    ypluli,
                cs_real_t    apow,
                cs_real_t    bpow,
                cs_real_t    dpow,
                cs_real_t    l_visc,
                cs_real_t    vel,
                cs_real_t    y,
                int         *iuntur,
                cs_lnum_t   *nsubla,
                cs_lnum_t   *nlogla,
                cs_real_t   *ustar,
                cs_real_t   *uk,
                cs_real_t   *yplus,
                cs_real_t   *ypup,
                cs_real_t   *cofimp)
{
  double ustarwer, ustarmin, ustaro, ydvisc;
  double eps = 0.001;
  int niter_max = 100;
  int iter = 0;
  double reynolds;

  /* Compute the local Reynolds number */

  ydvisc = y / l_visc;
  reynolds = vel * ydvisc;

  /*
   * Compute the friction velocity ustar
   */

  /* In the viscous sub-layer: U+ = y+ */
  if (reynolds <= ypluli * ypluli) {

    *ustar = sqrt(vel / ydvisc);
    *yplus = *ustar * ydvisc;
    *uk = *ustar;
    *ypup = 1.;
    *cofimp = 0.;

    /* Disable the wall funcion count the cell in the viscous sub-layer */
    *iuntur = 0;
    *nsubla += 1;

  /* In the log layer */
  } else {

    /* The initial value is Wener or the minimun ustar to ensure convergence */
    ustarwer = pow(fabs(vel) / apow / pow(ydvisc, bpow), dpow);
    ustarmin = exp(-cstlog * xkappa)/ydvisc;
    ustaro = CS_MAX(ustarwer, ustarmin);
    *ustar = (xkappa * vel + ustaro)
           / (log(ydvisc * ustaro) + xkappa * cstlog + 1.);

    /* Iterative solving */
    for (iter = 0;   iter < niter_max
                  && fabs(*ustar - ustaro) >= eps * ustaro; iter++) {
      ustaro = *ustar;
      *ustar = (xkappa * vel + ustaro)
             / (log(ydvisc * ustaro) + xkappa * cstlog + 1.);
    }

    if (iter >= niter_max) {
      bft_printf(_("WARNING: non-convergence in the computation\n"
                   "******** of the friction velocity\n\n"
                   "face number: %d \n"
                   "friction vel: %f \n" ), ifac, *ustar);
    }

    *uk = *ustar;
    *yplus = *ustar * ydvisc;
    *ypup = *yplus / (log(*yplus) / xkappa + cstlog);
    *cofimp = 1. - *ypup / xkappa * 1.5 / *yplus;

    /* Count the cell in the log layer */
    *nlogla += 1;

  }

}

/*----------------------------------------------------------------------------
 * Log law: picewise linear and log, with two velocity scales base on the
 * friction and the TKE.
 *
 * parameters:
 *   xkappa     <-- Von Karman constant
 *   cstlog     <-- Log law constant
 *   cmu025     <-- C_mu^1/4
 *   ypluli     <-- y+ limit
 *   l_visc     <-- kinematic viscosity
 *   t_visc     <-- turbulent kinematic viscosity
 *   vel        <-- wall projected cell center velocity
 *   y          <-- wall distance
 *   kinetic_en <-- turbulent kinetic energy
 *   iuntur     <-- indicator: 0 in the viscous sublayer
 *   nsubla     <-- counter of cell in the viscous sublayer
 *   nlogla     <-- counter of cell in the log-layer
 *   ustar      --> friction velocity
 *   uk         --> friction velocity
 *   yplus      --> non-dimension wall distance
 *   ypup       --> yplus projected vel ratio
 *   cofimp     --> |U_F|/|U_I^p| to ensure a good turbulence production
 *----------------------------------------------------------------------------*/

static void
_2scales_log_law(cs_real_t   xkappa,
                 cs_real_t   cstlog,
                 cs_real_t   cmu025,
                 cs_real_t   ypluli,
                 cs_real_t   l_visc,
                 cs_real_t   t_visc,
                 cs_real_t   vel,
                 cs_real_t   y,
                 cs_real_t   kinetic_en,
                 int        *iuntur,
                 cs_lnum_t  *nsubla,
                 cs_lnum_t  *nlogla,
                 cs_real_t  *ustar,
                 cs_real_t  *uk,
                 cs_real_t  *yplus,
                 cs_real_t  *ypup,
                 cs_real_t  *cofimp)
{
  double rcprod, ml_visc, Re, g;

  /* Compute the friction velocity ustar */

  /* Blending for very low values of k */
  Re = sqrt(kinetic_en) * y / l_visc;
  g = exp(-Re/11.);

  *uk = sqrt( (1.-g) * cmu025 * cmu025 * kinetic_en
            + g * l_visc * vel / y);

  *yplus = *uk * y / l_visc;

  /* log layer */
  if (*yplus > ypluli) {

    *ustar = vel / (log(*yplus) / xkappa + cstlog);
    *ypup = *yplus / (log(*yplus) / xkappa + cstlog);
    /* Mixing length viscosity */
    ml_visc = xkappa * l_visc * *yplus;
    rcprod = CS_MIN(xkappa, CS_MAX(1., sqrt(ml_visc / t_visc)) / *yplus);
    *cofimp = 1. - *ypup / xkappa * ( 2. * rcprod - 1. / (2. * *yplus));

    *nlogla += 1;

  /* viscous sub-layer */
  } else {

    if (*yplus > 1.e-12) {
      *ustar = fabs(vel / *yplus); /* FIXME remove that: its is here only to
                                      be fully equivalent to the former code. */
    } else {
      *ustar = 0.;
    }
    *ypup = 1.;
    *cofimp = 0.;

    *iuntur = 0;
    *nsubla += 1;

  }
}

/*----------------------------------------------------------------------------
 * Scalable wall function: shift the wall if "y+ < y+lim".
 *
 * parameters:
 *   xkappa     <-- Von Karman constant
 *   cstlog     <-- Log law constant
 *   cmu025     <-- C_mu^1/4
 *   ypluli     <-- y+ limit
 *   l_visc     <-- kinematic viscosity
 *   t_visc     <-- turbulent kinematic viscosity
 *   vel        <-- wall projected cell center velocity
 *   y          <-- wall distance
 *   kinetic_en <-- turbulent kinetic energy
 *   iuntur     <-- indicator: 0 in the viscous sublayer
 *   nsubla     <-- counter of cell in the viscous sublayer
 *   nlogla     <-- counter of cell in the log-layer
 *   ustar      --> friction velocity
 *   uk         --> friction velocity
 *   yplus      --> non-dimension wall distance
 *   dplus      --> non-dimension shift
 *   ypup       --> yplus projected vel ratio
 *   cofimp     --> |U_F|/|U_I^p| to ensure a good turbulence production
 *----------------------------------------------------------------------------*/

static void
_2scales_scalable_wallfunction(cs_real_t   xkappa,
                               cs_real_t   cstlog,
                               cs_real_t   cmu025,
                               cs_real_t   ypluli,
                               cs_real_t   l_visc,
                               cs_real_t   t_visc,
                               cs_real_t   vel,
                               cs_real_t   y,
                               cs_real_t   kinetic_en,
                               int        *iuntur,
                               cs_lnum_t  *nsubla,
                               cs_lnum_t  *nlogla,
                               cs_real_t  *ustar,
                               cs_real_t  *uk,
                               cs_real_t  *yplus,
                               cs_real_t  *dplus,
                               cs_real_t  *ypup,
                               cs_real_t  *cofimp
)
{
  double rcprod, ml_visc, Re, g;
  /* Compute the friction velocity ustar */

  /* Blending for very low values of k */
  Re = sqrt(kinetic_en) * y / l_visc;
  g = exp(-Re/11.);

  *uk = sqrt( (1.-g) * cmu025 * cmu025 * kinetic_en
            + g * l_visc * vel / y);

  *yplus = *uk * y / l_visc;

  /* Compute the friction velocity ustar */
  *uk = cmu025 * sqrt(kinetic_en);
  *yplus = *uk * y / l_visc;

  /* Log layer */
  if (*yplus > ypluli) {

    *dplus = 0.;

    *nlogla += 1;

  /* Viscous sub-layer and therefore shift */
  } else {

    *dplus = ypluli - *yplus;
    *yplus = ypluli;

    /* Count the cell as if it was in the viscous sub-layer */
    *nsubla += 1;

  }

  /* Mixing length viscosity */
  ml_visc = xkappa * l_visc * *yplus;
  rcprod = CS_MIN(xkappa, CS_MAX(1., sqrt(ml_visc / t_visc)) / *yplus);

  *ustar = vel / (log(*yplus) / xkappa + cstlog);
  *ypup = (*yplus - *dplus) / (log(*yplus) / xkappa + cstlog);
  *cofimp = 1. - *ypup / xkappa * (2. * rcprod - 1. / (2. * *yplus - *dplus));
}

/*----------------------------------------------------------------------------
 * No wall fucntions
 *
 * parameters:
 *   ypluli     <-- y+ limit
 *   l_visc     <-- kinematic viscosity
 *   t_visc     <-- turbulent kinematic viscosity
 *   vel        <-- wall projected cell center velocity
 *   y          <-- wall distance
 *   iuntur     <-- indicator: 0 in the viscous sublayer
 *   nsubla     <-- counter of cell in the viscous sublayer
 *   nlogla     <-- counter of cell in the log-layer
 *   ustar      --> friction velocity
 *   uk         --> friction velocity
 *   yplus      --> non-dimension wall distance
 *   dplus      --> non-dimension shift
 *   ypup       --> yplus projected vel ratio
 *   cofimp     --> |U_F|/|U_I^p| to ensure a good turbulence production
 *----------------------------------------------------------------------------*/

static void
_no_wallfunction(cs_real_t   ypluli,
                 cs_real_t   l_visc,
                 cs_real_t   t_visc,
                 cs_real_t   vel,
                 cs_real_t   y,
                 int        *iuntur,
                 cs_lnum_t  *nsubla,
                 cs_lnum_t  *nlogla,
                 cs_real_t  *ustar,
                 cs_real_t  *uk,
                 cs_real_t  *yplus,
                 cs_real_t  *dplus,
                 cs_real_t  *ypup,
                 cs_real_t  *cofimp)
{
  /* Compute the friction velocity ustar */

  *ustar = sqrt(vel * l_visc / y);
  *yplus = *ustar * y / l_visc;
  *uk = *ustar;
  *ypup = l_visc / (l_visc + t_visc);
  *cofimp = 0.;
  *iuntur = 0;

  if (*yplus <= ypluli) {

    /* Disable the wall funcion count the cell in the viscous sub-layer */
    *nsubla += 1;

  } else {

    /* Count the cell as if it was in the viscous sub-layer */
    *nsubla += 1;

  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to members of the wall functions structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ideuch --> pointer to cs_glob_wall_functions->ideuch
 *   iwallt --> pointer to cs_glob_wall_functions->iwallt
 *   ilogpo --> pointer to cs_glob_wall_functions->ilogpo
 *----------------------------------------------------------------------------*/

void
cs_f_wall_functions_get_pointers(int     **ideuch,
                                 int     **iwallt,
                                 int     **ilogpo)
{
  *ideuch = &(_wall_functions.ideuch);
  *iwallt = &(_wall_functions.iwallt);
  *ilogpo = &(_wall_functions.ilogpo);
}

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (wallfunctions, WALLFUNCTIONS)
(
 const cs_int_t   *const iwallf,
 const cs_lnum_t  *const ifac,
 const cs_real_t  *const xkappa,
 const cs_real_t  *const cstlog,
 const cs_real_t  *const cmu025,
 const cs_real_t  *const ypluli,
 const cs_real_t  *const apow,
 const cs_real_t  *const bpow,
 const cs_real_t  *const dpow,
 const cs_real_t  *const l_visc,
 const cs_real_t  *const t_visc,
 const cs_real_t  *const vel,
 const cs_real_t  *const y,
 const cs_real_t  *const kinetic_en,
       cs_int_t         *iuntur,
       cs_lnum_t        *nsubla,
       cs_lnum_t        *nlogla,
       cs_real_t        *ustar,
       cs_real_t        *uk,
       cs_real_t        *yplus,
       cs_real_t        *ypup,
       cs_real_t        *cofimp,
       cs_real_t        *dplus
)
{
  cs_wall_functions_velocity(*iwallf,
                             *ifac,
                             *xkappa,
                             *cstlog,
                             *cmu025,
                             *ypluli,
                             *apow,
                             *bpow,
                             *dpow,
                             *l_visc,
                             *t_visc,
                             *vel,
                             *y,
                             *kinetic_en,
                             iuntur,
                             nsubla,
                             nlogla,
                             ustar,
                             uk,
                             yplus,
                             ypup,
                             cofimp,
                             dplus);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (hturbp, HTURBP)
(
 const cs_real_t  *const prl,
 const cs_real_t  *const prt,
 const cs_real_t  *const ckarm,
 const cs_real_t  *const yplus,
 const cs_real_t  *const dplus,
       cs_real_t        *htur,
       cs_real_t        *yplim
)
{
  cs_wall_functions_scalar(*prl,
                           *prt,
                           *ckarm,
                           *yplus,
                           *dplus,
                           htur,
                           yplim);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*! \brief  Compute the friction velocity and \f$y^+\f$ / \f$u^+\f$.

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     iwallf        wall function type
 * \param[in]     ifac          face number
 * \param[in]     xkappa        Von Karman constant
 * \param[in]     cstlog        Log law constant
 * \param[in]     cmu025        \f$ C_{\mu}^{1/4} \f$
 * \param[in]     ypluli        \f$y^+\f$ limit
 * \param[in]     apow          Coef of the Wener law
 * \param[in]     bpow          Coef of the Wener law
 * \param[in]     dpow          Coef of the Wener law
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     t_visc        turbulent kinematic viscosity
 * \param[in]     vel           wall projected cell center velocity
 * \param[in]     y             wall distance
 * \param[in]     kinetic_en    turbulente kinetic energy
 * \param[in]     iuntur        indicator: 0 in the viscous sublayer
 * \param[in]     nsubla        counter of cell in the viscous sublayer
 * \param[in]     nlogla        counter of cell in the log-layer
 * \param[out]    ustar         friction velocity
 * \param[out]    uk            friction velocity
 * \param[out]    yplus         non-dimension wall distance
 * \param[out]    ypup          yplus projected vel ratio
 * \param[out]    cofimp        \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                              turbulence production
 * \param[out]    dplus         dimensionless shift to the wall for scalable
 *                              wall functions
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_velocity(int          iwallf,
                           cs_lnum_t    ifac,
                           cs_real_t    xkappa,
                           cs_real_t    cstlog,
                           cs_real_t    cmu025,
                           cs_real_t    ypluli,
                           cs_real_t    apow,
                           cs_real_t    bpow,
                           cs_real_t    dpow,
                           cs_real_t    l_visc,
                           cs_real_t    t_visc,
                           cs_real_t    vel,
                           cs_real_t    y,
                           cs_real_t    kinetic_en,
                           int         *iuntur,
                           cs_lnum_t   *nsubla,
                           cs_lnum_t   *nlogla,
                           cs_real_t   *ustar,
                           cs_real_t   *uk,
                           cs_real_t   *yplus,
                           cs_real_t   *ypup,
                           cs_real_t   *cofimp,
                           cs_real_t   *dplus)
{
  /* Pseudo shift of the wall, 0 by default */
  *dplus = 0.;

  /* Activation of wall function by default */
  *iuntur = 1;

  if (iwallf == 3) {

    _1scale_power_law(ypluli,
                      apow,
                      bpow,
                      dpow,
                      l_visc,
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

  } else if (iwallf == 0) {

    _1scale_log_law(ifac,
                    xkappa,
                    cstlog,
                    ypluli,
                    apow,
                    bpow,
                    dpow,
                    l_visc,
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

  } else if (iwallf == 1) {

    _2scales_log_law(xkappa,
                     cstlog,
                     cmu025,
                     ypluli,
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

  } else if (iwallf == 2) {

    _2scales_scalable_wallfunction(xkappa,
                                   cstlog,
                                   cmu025,
                                   ypluli,
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
                                   dplus,
                                   ypup,
                                   cofimp);

  } else if (iwallf == 4) {

    _no_wallfunction(ypluli,
                     l_visc,
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
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief Compute the correction of the exchange coefficient between the fluid
  and the wall for a turbulent flow.

  This is function of the dimensionless
  distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.

  Then the return coefficient reads:
  \f[
  h_{tur} = Pr \dfrac{y^+}{T^+}
  \f]

  This coefficient is computed thanks to a similarity model between
  dynamic viscous sub-layer and themal sub-layer.

  \f$ T^+ \f$ is computed as follows:

  - For a laminar Prandtl number smaller than 0.1 (such as liquid metals),
    the standard model with two sub-layers (Prandtl-Taylor) is used.

  - For a laminar Prandtl number larger than 0.1 (such as liquids and gaz),
    a model with three sub-layers (Arpaci-Larsen) is used.

  The final exchange coefficient is:
  \f[
  h = \dfrac{K}{\centip \centf} h_{tur}
  \f]

*/
/*-----------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     prl           laminar Prandtl number
 * \param[in]     prt           turbulent Prandtl number
 * \param[in]     ckarm         Von Karman constant
 * \param[in]     yplus         dimensionless distance to the wall
 * \param[in]     dplus         dimensionless distance for scalable
 *                              wall functions
 * \param[out]    htur          corrected exchange coefficient
 * \param[out]    yplim         value of the limit for \f$ y^+ \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_scalar(double  prl,
                         double  prt,
                         double  ckarm,
                         double  yplus,
                         double  dplus,
                         double  *htur,
                         double  *yplim)
{

  /* Local variables */

  double tplus;
  double beta2,a2;
  double yp2;
  double prlm1;

  const double epzero = 1.e-12;

  /*==========================================================================*/

  /*==========================================================================
    1. Initializations
    ==========================================================================*/

  /*==========================================================================*/

  (*htur) = CS_MAX(yplus-dplus,epzero)/CS_MAX(yplus,epzero);

  prlm1 = 0.1;

  /*==========================================================================
    2. Compute htur for small Prandtl numbers
    ==========================================================================*/

  if (prl <= prlm1) {
    (*yplim)   = prt/(prl*ckarm);
    if (yplus > (*yplim)) {
      tplus = prl*(*yplim) + prt/ckarm * log(yplus/(*yplim));
      (*htur) = prl*(yplus-dplus)/tplus;
    }

    /*========================================================================
      3. Compute htur for the model with three sub-layers
      ========================================================================*/

  } else {
    yp2   = ckarm*1000./prt;
    yp2   = sqrt(yp2);
    (*yplim)   = pow(1000./prl,1./3.);

    a2 = 15.*pow(prl,2./3.);
    beta2 = a2 - 500./ pow(yp2,2);

    if (yplus >= (*yplim) && yplus < yp2) {
      tplus = a2 - 500./(yplus*yplus);
      (*htur) = prl*(yplus-dplus)/tplus;
    }

    if (yplus >= yp2) {
      tplus = beta2 + prt/ckarm*log(yplus/yp2);
      (*htur) = prl*(yplus-dplus)/tplus;
    }

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
