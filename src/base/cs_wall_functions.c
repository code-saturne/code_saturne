/*============================================================================
 * Wall functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
                  cs_int_t   *iuntur,
                  cs_int_t   *nsubla,
                  cs_int_t   *nlogla,
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
_1scale_log_law(cs_int_t     ifac,
                cs_real_t    xkappa,
                cs_real_t    cstlog,
                cs_real_t    ypluli,
                cs_real_t    apow,
                cs_real_t    bpow,
                cs_real_t    dpow,
                cs_real_t    l_visc,
                cs_real_t    vel,
                cs_real_t    y,
                cs_int_t    *iuntur,
                cs_int_t    *nsubla,
                cs_int_t    *nlogla,
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
                 cs_int_t   *iuntur,
                 cs_int_t   *nsubla,
                 cs_int_t   *nlogla,
                 cs_real_t  *ustar,
                 cs_real_t  *uk,
                 cs_real_t  *yplus,
                 cs_real_t  *ypup,
                 cs_real_t  *cofimp)
{
  double rcprod, ml_visc;

  /* Compute the friction velocity ustar */
  *uk = cmu025 * sqrt(kinetic_en);
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
                               cs_int_t   *iuntur,
                               cs_int_t   *nsubla,
                               cs_int_t   *nlogla,
                               cs_real_t  *ustar,
                               cs_real_t  *uk,
                               cs_real_t  *yplus,
                               cs_real_t  *dplus,
                               cs_real_t  *ypup,
                               cs_real_t  *cofimp
)
{
  double rcprod, ml_visc;

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
                 cs_int_t   *iuntur,
                 cs_int_t   *nsubla,
                 cs_int_t   *nlogla,
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
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the friction velocity and y+/u+
 *----------------------------------------------------------------------------*/

void CS_PROCF (wallfunctions, WALLFUNCTIONS)
(
 const cs_int_t   *const iwallf,      /* <-- wall function type               */
 const cs_int_t   *const ifac,        /* <-- face number                      */
 const cs_real_t  *const xkappa,      /* <-- Von Karman constant              */
 const cs_real_t  *const cstlog,      /* <-- Log law constant                 */
 const cs_real_t  *const cmu025,      /* <-- C_mu^1/4                         */
 const cs_real_t  *const ypluli,      /* <-- y+ limit                         */
 const cs_real_t  *const apow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const bpow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const dpow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const l_visc,      /* <-- kinematic viscosity              */
 const cs_real_t  *const t_visc,      /* <-- turbulent kinematic viscosity    */
 const cs_real_t  *const vel,         /* <-- wall projected
                                             cell center velocity             */
 const cs_real_t  *const y,           /* <-- wall distance                    */
 const cs_real_t  *const kinetic_en,  /* <-- turbulente kinetic energy        */
       cs_int_t         *iuntur,      /* <-- indicator:
                                              0 in the viscous sublayer       */
       cs_int_t         *nsubla,      /* <-- counter of cell in the viscous
                                             sublayer                         */
       cs_int_t         *nlogla,      /* <-- counter of cell in the log-layer */
       cs_real_t        *ustar,       /* --> friction velocity                */
       cs_real_t        *uk,          /* --> friction velocity                */
       cs_real_t        *yplus,       /* --> non-dimension wall distance      */
       cs_real_t        *ypup,        /* --> yplus projected vel ratio        */
       cs_real_t        *cofimp,      /* --> |U_F|/|U_I^p| to ensure a good
                                             turbulence production            */
       cs_real_t        *dplus        /* --> dimensionless shift to the wall
                                             for scalable wall functions      */
)
{
  /* Pseudo shift of the wall, 0 by default */
  *dplus = 0.;

  /* Activation of wall function by default */
  *iuntur = 1;

  if (*iwallf == 3) {

    _1scale_power_law(*ypluli,
                      *apow,
                      *bpow,
                      *dpow,
                      *l_visc,
                      *vel,
                      *y,
                      iuntur,
                      nsubla,
                      nlogla,
                      ustar,
                      uk,
                      yplus,
                      ypup,
                      cofimp);

  } else if (*iwallf == 0) {

    _1scale_log_law(*ifac,
                    *xkappa,
                    *cstlog,
                    *ypluli,
                    *apow,
                    *bpow,
                    *dpow,
                    *l_visc,
                    *vel,
                    *y,
                     iuntur,
                     nsubla,
                     nlogla,
                     ustar,
                     uk,
                     yplus,
                     ypup,
                     cofimp);

  } else if (*iwallf == 1) {

    _2scales_log_law(*xkappa,
                     *cstlog,
                     *cmu025,
                     *ypluli,
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
                      cofimp);

  } else if (*iwallf == 2) {

    _2scales_scalable_wallfunction(*xkappa,
                                   *cstlog,
                                   *cmu025,
                                   *ypluli,
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
                                    dplus,
                                    ypup,
                                    cofimp);

  } else if (*iwallf == 4) {

    _no_wallfunction(*ypluli,
                     *l_visc,
                     *t_visc,
                     *vel,
                     *y,
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
